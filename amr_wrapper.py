#!/usr/bin/env python

__version__ = "0.0.1"

import os
import sys
import re
import glob
import argparse
import textwrap
import operator
import subprocess
import datetime
import pandas as pd
from Bio import SeqIO

from vsnp_fastq_quality import FASTQ_Quality
from spades_assembly import Spades_Assembly
from spades_stats_parse import Spades_Stats
from kraken2_run import Kraken2_Identification
from seqsero2_wrapper import SeqSero2
from mlst_wrapper import MLST
from abricate_wrapper import Abricate
from amrfinder_wrapper import AMR_Finder
from quality_scaling import Quality_Scaling
from latex_reporter import AMR_Latex_Report

class Run_AMR_Wrapper:
    ''' 
    '''
    def __init__(self, **kwargs):
        self.cwd = os.getcwd()
        self.FASTQ_R1 = kwargs.get('FASTQ_R1', None)
        self.FASTQ_R2 = kwargs.get('FASTQ_R2', None)
        self.FASTA = kwargs.get('FASTA', None)
        self.abricate_report = kwargs.get('abricate_report', None)
        self.abricate_depth = kwargs.get('abricate_depth', 0)
        self.abricate_coverage = kwargs.get('abricate_coverage', 75)
        self.debug = kwargs.get('debug', False)
        if self.FASTA:
            self.assembly = self.FASTA
            self.sample_name = re.sub('[_.].*', '', self.FASTA)
        else:
            self.sample_name = re.sub('[_.].*', '', self.FASTQ_R1)

    def run_kraken(self,):
        if self.FASTA:
            kraken2 = Kraken2_Identification(FASTA=self.FASTA, db_contents=None, directory="kraken2")
            report, output = kraken2.kraken2_run()
        else:
            kraken2 = Kraken2_Identification(FASTQ_R1=self.FASTQ_R1, FASTQ_R2=self.FASTQ_R2, db_contents=None, directory="kraken2")
            report, output = kraken2.kraken2_run()
        self.krona_html = kraken2.krona_make_graph(report, output)
        kraken2.bracken(report, output)

    def quality(self,):
        self.fq = FASTQ_Quality(self.FASTQ_R1, self.FASTQ_R2)
        self.fq.get_quality()

    def spades_assembly(self,):
        if self.FASTQ_R1:
            assemble = Spades_Assembly(self.FASTQ_R1, self.FASTQ_R2)
            self.assembly = assemble.assembly_file
            self.stats = Spades_Stats(self.assembly)
            self.stats.write_stats(self.stats, self.fq)
        else:
            self.stats = Spades_Stats(self.assembly)
            self.stats.write_stats(self.stats)
        
    def seqsero2_wrapper(self,):
        self.seqsero2 = SeqSero2(self.FASTQ_R1, self.FASTQ_R2)
        self.seqsero2.run()

    def mlst_wrapper(self,):
        self.mlst = MLST(self.assembly)
        self.mlst.run()

    def abricate_wrapper(self,):
        self.abricate = Abricate(self.assembly, self.abricate_depth, self.abricate_coverage)
        self.abricate.run()

    def amrfinder_wrapper(self,):
        self.amr_finder = AMR_Finder(self.assembly)
        self.amr_finder.run()
    
    def calculate_quality_scaling(self,):
        self.quality_scaling = Quality_Scaling()
        if self.FASTQ_R1 and self.FASTQ_R2:
            if self.mlst.mlst_size_lookup:
                genome_size = int(self.mlst.mlst_size_lookup) * 1000000 #mlst_size_lookup is in Mb
                size_method = "Based on MLST identification"
            else:
                genome_size=self.stats.total_contig_lengths
                size_method = "SPAdes total contig length"
            fastq_scaling_variable, genome_coverage_depth = self.quality_scaling.fastq_scaling(
                genome_size=genome_size,
                read1_reads_gt_q30=self.fq.read1.reads_gt_q30,
                read2_reads_gt_q30=self.fq.read2.reads_gt_q30,
                sampling_size=self.fq.read1.sampling_size,
                read1_total_read_count=self.fq.read1.total_read_count,
                read2_total_read_count=self.fq.read2.total_read_count, 
                read1_read_average=self.fq.read1.read_average,
                read2_read_average=self.fq.read2.read_average,)

            self.fastq_scaling_variable = fastq_scaling_variable
            self.genome_size = genome_size
            self.genome_coverage_depth = genome_coverage_depth
            self.coverage_method = "Calculated by number of reads x 240/Genome Length"
            self.size_method = size_method
        else: # when just FASTA supplied to script
            if self.mlst.mlst_size_lookup:
                genome_size = int(self.mlst.mlst_size_lookup) * 1000000  #mlst_size_lookup is in Mb
                size_method = "Based on MLST identification"
            else:
                genome_size = self.stats.total_contig_lengths
                size_method = "SPAdes total contig length"
            self.genome_size = genome_size
            self.fastq_scaling_variable = None
            self.genome_coverage_depth = self.stats.mean_coverage
            self.coverage_method = "SPAdes estimated average coverage"
            self.size_method = size_method

        assembly_scaling_variable = self.quality_scaling.assembly_scaling(
            longest_contig=self.stats.longest_contig,
            greater_one_kb_count = self.stats.greater_one_kb_count,
            contig_count = self.stats.contig_count,
            n50=self.stats.n50,
            l50 = self.stats.l50,
            rgl = self.rgl,
            stat_total_contig_lengths=self.stats.total_contig_lengths,
            stat_contig_count=self.stats.contig_count,
            genome_coverage_depth=self.genome_coverage_depth)  #this is the only connection to quality_scaling.fastq_scaling()
            
        self.assembly_scaling_variable = assembly_scaling_variable

    def genome_gt_1kb(self,):
        assembled_fasta = SeqIO.parse(self.assembly, "fasta")
        measure_length = 0
        for contig in assembled_fasta:
            if len(contig) > 1000:
                measure_length = measure_length + len(contig)
        gl = (measure_length/self.stats.total_contig_lengths)*100
        self.rgl = round(gl, 2)

    def run(self,):
        self.run_kraken()
        if self.FASTQ_R1:
            self.quality()
            self.seqsero2_wrapper()
        self.spades_assembly()
        self.mlst_wrapper()
        self.abricate_wrapper()
        self.amrfinder_wrapper()
        self.genome_gt_1kb()
        self.calculate_quality_scaling()


if __name__ == "__main__": # execute if directly access by the interpreter
    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\
    ---------------------------------------------------------
    Usage:
        amr_wrapper.py -r1 *_R1*fastq.gz -r2 *_R2*fastq.gz
        Default setting will include FASTQ metrics, Assembly metrics, SeqSero2 results (if salmonella), MLST, and AMRFinder results.
            -a will add Abricate
            -m will just minimize the report to just FASTQ metrics, Assembly metrics, SeqSero2 results (if salmonella) and MLST.

    Additional:  nahln_amr_updates.sh
    If NAHLM samples... separate after running amr_wrapper.sh on all samples.  Separate as EC (ecoli), MH, SIG (staph intermedius group) then run nahln_amr_updates.sh on each group.  Before running make sure the assemblies (.fasta) are in the sample folders, which they should be.
    mkdir EC MH SIG; mv EC-* EC; mv MH-* MH; mv SIG-* SIG
    each directory... nahln_amr_updates.sh

    '''), epilog='''---------------------------------------------------------''')
    
    parser.add_argument('-r1', '--FASTQ_R1', action='store', dest='FASTQ_R1', required=False, default=None, help='R1 FASTQ gz file')
    parser.add_argument('-r2', '--FASTQ_R2', action='store', dest='FASTQ_R2', required=False, default=None, help='R2 FASTQ gz file')
    parser.add_argument('-f', '--FASTA', action='store', dest='FASTA', required=False, default=None, help='Assembly FASTA file')
    parser.add_argument('-m', '--mininum_report', action='store_true', dest='mininum_report', default=False, help='OPTIONAL: Only include FASTQ, Assembly SeqSero2 and MLST on report')
    parser.add_argument('-a', '--abricate_report', action='store_true', dest='abricate_report', default=False, help='OPTIONAL: include abricate results in report')
    parser.add_argument('-b', '--abricate_depth', action='store', dest='abricate_depth', default=0, help='OPTIONAL: percent average depth cutoff for abricate, aka: --mincov, cvb use -b 50')
    parser.add_argument('-c', '--abricate_coverage', action='store', dest='abricate_coverage', default=75, help='OPTIONAL: percent genome coverage cutoff for abricate, aka: --minid, cvb use -c 90')
    parser.add_argument('-d', '--debug', action='store_true', dest='debug', default=False, help='keep files for debugging.')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.basename(__file__)}: version {__version__}')
    args = parser.parse_args()
    
    print(f'\n{os.path.basename(__file__)} SET ARGUMENTS:')
    print(args)
        
    run_amr_wrapper = Run_AMR_Wrapper(FASTQ_R1=args.FASTQ_R1, FASTQ_R2=args.FASTQ_R2, FASTA=args.FASTA, abricate_report=args.abricate_report, abricate_depth=args.abricate_depth, abricate_coverage=args.abricate_coverage, debug=args.debug)
    run_amr_wrapper.run()

    if args.FASTQ_R1 and args.FASTQ_R2:
        read1_fastq = run_amr_wrapper.fq.read1.fastq
        read2_fastq = run_amr_wrapper.fq.read2.fastq
        read1_file_size = run_amr_wrapper.fq.read1.file_size
        read2_file_size = run_amr_wrapper.fq.read2.file_size
        read1_read_average = run_amr_wrapper.fq.read1.read_average
        read2_read_average = run_amr_wrapper.fq.read2.read_average
        read1_reads_gt_q30 = run_amr_wrapper.fq.read1.reads_gt_q30
        read2_reads_gt_q30 = run_amr_wrapper.fq.read2.reads_gt_q30
        sampling_size = run_amr_wrapper.fq.read1.sampling_size
        seqsero2_serotype = run_amr_wrapper.seqsero2.serotype
        seqsero2_antigenic = run_amr_wrapper.seqsero2.antigenic
        seqsero2_subspecies=run_amr_wrapper.seqsero2.subspecies
        seqserocomment=run_amr_wrapper.seqsero2.seqserocomment
        seqsero_file=run_amr_wrapper.seqsero2.seqsero_file
    elif args.FASTQ_R1 and not args.FASTQ_R2:
        read1_fastq = run_amr_wrapper.fq.read1.fastq
        read2_fastq = None
        read1_file_size = run_amr_wrapper.fq.read1.file_size
        read2_file_size = None
        read1_read_average = run_amr_wrapper.fq.read1.read_average
        read2_read_average = None
        read1_reads_gt_q30 = run_amr_wrapper.fq.read1.reads_gt_q30
        read2_reads_gt_q30 = None
        sampling_size = run_amr_wrapper.fq.read1.sampling_size
        seqsero2_serotype = run_amr_wrapper.seqsero2.serotype
        seqsero2_antigenic = run_amr_wrapper.seqsero2.antigenic
        seqsero2_subspecies=run_amr_wrapper.seqsero2.subspecies
        seqserocomment=run_amr_wrapper.seqsero2.seqserocomment
        seqsero_file=run_amr_wrapper.seqsero2.seqsero_file
    elif args.FASTA: #not available when FASTA input, set to None before calling function
        read1_fastq = None
        read2_fastq = None
        read1_file_size = None
        read2_file_size = None
        read1_read_average = None
        read2_read_average = None
        read1_reads_gt_q30 = None
        read2_reads_gt_q30 = None
        sampling_size = None
        seqval = None
        seqsero2_serotype = None
        seqsero2_antigenic = None
        seqsero2_subspecies = None
        seqserocomment = None
        seqsero_file = None
        
    print(f'Building reports...')
    amr_latex_report = AMR_Latex_Report(
        fastq_scaling_variable=run_amr_wrapper.fastq_scaling_variable,
        assembly_scaling_variable=run_amr_wrapper.assembly_scaling_variable,
        genome_size=run_amr_wrapper.genome_size,
        genome_coverage_depth=run_amr_wrapper.genome_coverage_depth,
        coverage_method=run_amr_wrapper.coverage_method,
        size_method=run_amr_wrapper.size_method,
        abricate_ab_version = run_amr_wrapper.abricate.abricate_ab_version,
        amr_version=run_amr_wrapper.amr_finder.version,
        read1_fastq = read1_fastq,
        read2_fastq = read2_fastq)

    amr_latex_report.latex_document(
        sample_name = run_amr_wrapper.stats.sample_name,
        read1_fastq = read1_fastq,
        read2_fastq = read2_fastq,
        read1_file_size = read1_file_size,
        read2_file_size = read2_file_size,
        read1_read_average = read1_read_average,
        read2_read_average = read2_read_average,
        read1_reads_gt_q30 = read1_reads_gt_q30,
        read2_reads_gt_q30=read2_reads_gt_q30,
        sampling_size = sampling_size,
        stat_contig_count=run_amr_wrapper.stats.contig_count,
        stat_mean_coverage = run_amr_wrapper.stats.mean_coverage,
        stat_total_contig_lengths = run_amr_wrapper.stats.total_contig_lengths,
        stat_longest_contig = run_amr_wrapper.stats.longest_contig,
        stat_greater_one_kb_count = run_amr_wrapper.stats.greater_one_kb_count,
        stat_n50 = run_amr_wrapper.stats.n50,
        stat_l50=run_amr_wrapper.stats.l50,
        spades_version = run_amr_wrapper.stats.spades_version,
        mlst_file = run_amr_wrapper.mlst.mlst_file,
        mlst_scheme = run_amr_wrapper.mlst.mlst_scheme,
        mlst_st=run_amr_wrapper.mlst.mlst_type,
        mlst_detail = run_amr_wrapper.mlst.mlst_detail,
        mlst_version = run_amr_wrapper.mlst.version,
        mlst_species_lookup = run_amr_wrapper.mlst.mlst_species_lookup,
        seqsero2_serotype = seqsero2_serotype,
        seqsero2_antigenic = seqsero2_antigenic,
        seqsero2_subspecies = seqsero2_subspecies,
        seqserocomment = seqserocomment,
        seqsero_file = seqsero_file,
        amrfinder_file=run_amr_wrapper.amr_finder.amrfinder_file,
        mininum_report = args.mininum_report,
        abricate_report = args.abricate_report,
        abricate_depth = run_amr_wrapper.abricate.depth,
        abricate_coverage=run_amr_wrapper.abricate.coverage,
        ab_ncbi_file = run_amr_wrapper.abricate.ab_ncbi_file,
        ab_resfinder_file = run_amr_wrapper.abricate.ab_resfinder_file,
        abricate_ncbi_version_date = run_amr_wrapper.abricate.abricate_ncbi_version_date,
        abricate_res_version_date=run_amr_wrapper.abricate.abricate_res_version_date,
        abricate_ncbi_seq_number=run_amr_wrapper.abricate.abricate_ncbi_seq_number,
        abricate_res_seq_number=run_amr_wrapper.abricate.abricate_res_seq_number,
        rgl=run_amr_wrapper.rgl)
    
    if not args.debug:
        for each_file in ('*.png', '*.svg', '*.tex', '*out', '*.aux', '*.log'):
            for each in glob.glob(each_file):
                os.remove(each) 

    print("done")

# Created March 2020 by Tod Stuber
