#!/usr/bin/env python

import os
import sys
import re
from operator import itemgetter
import numpy as np
import pandas as pd
import argparse
import textwrap
import time
from datetime import datetime
from collections import Counter
from Bio import Seq
from Bio import SeqIO

from vsnp_fastq_quality import FASTQ_Quality

class bcolors:
    PURPLE = '\033[95m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    ENDC = '\033[0m'

class Spades_Stats:

    def __init__(self, fasta_in):
        base_name = os.path.basename(fasta_in)
        self.sample_name = re.sub('[._].*', '', base_name)
        records = list(SeqIO.parse(fasta_in, "fasta"))
        cov_length_list=[]
        contig_count = 0
        coverage_list=[]
        length_list=[]
        small_contigs=[]
        greater_one_kb=[]
        mid_size = []
        for rec in records:
            header = rec.description
            try:
                coverage_value = header.split('_')[5]
            except IndexError:
                coverage_value = 1
            coverage_value = int(float(coverage_value))
            cov_length_list.append({'name': rec.description, 'cov': coverage_value, 'length': len(rec)})
            coverage_list.append(coverage_value)
            length_list.append(len(rec))
            contig_count += 1
            if len(rec) <= 300:
                small_contigs.append(len(rec))
            elif len(rec) >= 1000:
                greater_one_kb.append(len(rec))
            else:
                mid_size.append(len(rec))
        total_contig_lengths = int(sum(length_list))

        #Calculate mean coverage
        normalized_list = []
        for rec in records:
            header = rec.description
            try:
                coverage_value = header.split('_')[5]
            except IndexError:
                coverage_value = 1
            coverage_value = int(float(coverage_value))
            normalized_list.append((len(rec) / total_contig_lengths) * coverage_value)
        mean_coverage = sum(normalized_list)

        #N50 calculation
        all_len = sorted(length_list, reverse=True)
        csum = np.cumsum(all_len)
        n2 = int(sum(length_list)/2)
        csumn2 = min(csum[csum >= n2])
        ind = np.where(csum == csumn2)
        self.n50 = all_len[int(ind[0])] # n50 smallest size contig which, along with the larger contigs, contain half of sequence of a particular genome
        self.l50 = int(ind[0][0]) + 1 # l50 smallest number of contigs whose length sum makes up half of genome

        self.cov_length_list = cov_length_list
        self.contig_count = int(contig_count)
        self.longest_contig = int(max(length_list))
        self.total_contig_lengths = total_contig_lengths
        self.mean_coverage = mean_coverage
        self.small_contigs_count = len(small_contigs)
        self.greater_one_kb_count = len(greater_one_kb)
        self.mid_size = len(mid_size)
        self.spades_version = os.popen("spades.py --version").readlines()[0]
        
    def print_by_coverage(self,):
        for each_dict in sorted(self.cov_length_list, key=itemgetter('cov')):
            print(f'{each_dict["cov"]:,}X  {each_dict["length"]:,}  {each_dict["name"]}')

    def print_by_length(self,):
        for each_dict in sorted(self.cov_length_list, key=itemgetter('length')):
            print(f'{each_dict["length"]:,}  {each_dict["cov"]:,}X  {each_dict["name"]}')

    def write_stats(self, stats, fq=None):
        #stats to excel
        ts = time.time()
        st = datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H-%M-%S')
        sample_name = self.sample_name
        if fq:
            coverage_title = 'Mean Coverage: read count * read size / total assembly length'
        else:
            coverage_title = 'Mean Coverage: SPAdes reporting'
        df = pd.DataFrame(index=[sample_name], \
            columns=[ \
                'Read 1', 'R1 File Size', 'R1 Total Reads', 'R1 Mean Length', 'R1 Mean Quality', 'R1 Passing Q30', \
                'Read 2', 'R2 File Size', 'R2 Total Reads', 'R2 Mean Length', 'R2 Mean Quality', 'R2 Passing Q30', \
                'Assembly Contig Count', '<300bp Count', '301-999bp Count', '>1kb Count', 'Total Length', 'Longest Contig', 'N50', coverage_title, ])
        try:
            df.at[sample_name, 'Read 1'] = f'{fq.read1.fastq}'
            df.at[sample_name, 'R1 File Size'] = f'{fq.read1.file_size}'
            df.at[sample_name, 'R1 Total Reads'] = f'{fq.read1.total_read_count:,}'
            df.at[sample_name, 'R1 Mean Length'] = f'{fq.read1.length_mean:.1f}'
            df.at[sample_name, 'R1 Mean Quality'] = f'{fq.read1.read_average:.1f}'
            df.at[sample_name, 'R1 Passing Q30'] = f'{fq.read1.reads_gt_q30/fq.read1.sampling_size:0.1%}'
            all_reads = fq.read1.total_read_count
        except AttributeError:
            df.at[sample_name, 'Read 1'] = 'NA'
            df.at[sample_name, 'R1 File Size'] = 'NA'
            df.at[sample_name, 'R1 Total Reads'] = 'NA'
            df.at[sample_name, 'R1 Mean Length'] = 'NA'
            df.at[sample_name, 'R1 Mean Quality'] = 'NA'
            df.at[sample_name, 'R1 Passing Q30'] = 'NA'
        try: 
            df.at[sample_name, 'Read 2'] = f'{fq.read2.fastq}'
            df.at[sample_name, 'R2 File Size'] = f'{fq.read2.file_size}'
            df.at[sample_name, 'R2 Total Reads'] = f'{fq.read2.total_read_count:,}'
            df.at[sample_name, 'R2 Mean Length'] = f'{fq.read2.length_mean:.1f}'
            df.at[sample_name, 'R2 Mean Quality'] = f'{fq.read2.read_average:.1f}'
            df.at[sample_name, 'R2 Passing Q30'] = f'{fq.read2.reads_gt_q30/fq.read2.sampling_size:0.1%}'
            all_reads = all_reads + fq.read2.total_read_count
        except AttributeError:
            df.at[sample_name, 'Read 2'] = 'NA'
            df.at[sample_name, 'R2 File Size'] = 'NA'
            df.at[sample_name, 'R2 Total Reads'] = 'NA'
            df.at[sample_name, 'R2 Mean Length'] = 'NA'
            df.at[sample_name, 'R2 Mean Quality'] = 'NA'
            df.at[sample_name, 'R2 Passing Q30'] = 'NA'
        df.at[sample_name, 'Assembly Contig Count'] = f'{stats.contig_count:,}'
        df.at[sample_name, '<300bp Count'] = f'{stats.small_contigs_count:,}'
        df.at[sample_name, '301-999bp Count'] = f'{stats.mid_size:,}'
        df.at[sample_name, '>1kb Count'] = f'{stats.greater_one_kb_count:,}'
        df.at[sample_name, 'Total Length'] = f'{stats.total_contig_lengths:,}'
        df.at[sample_name, 'Longest Contig'] = f'{stats.longest_contig:,}'
        df.at[sample_name, 'N50'] = f'{stats.n50:,}'
        if fq: #calculating from FASTQ reads more accurate than spades reporting
            mean_coverage = ((fq.read1.total_read_count * fq.read1.length_mean)*2)/stats.total_contig_lengths
            df.at[sample_name, coverage_title] = f'{mean_coverage:,.1f}X'
        else: #when no FASTQ info just use spades reporting
            df.at[sample_name, coverage_title] = f'{stats.mean_coverage:,.1f}X'
            mean_coverage = stats.mean_coverage
        df.index.name = 'sample'
        df.to_excel(f'{sample_name}_{st}_stats.xlsx')
        self.stats_excel = f'{os.getcwd()}/{sample_name}_{st}_stats.xlsx'

        print(f'Contig count: {bcolors.YELLOW}{stats.contig_count:,}{bcolors.ENDC}, Contig length counts <|301-999bp|>: {bcolors.RED}{stats.small_contigs_count:,}{bcolors.ENDC}|{bcolors.BLUE}{stats.mid_size:,}{bcolors.ENDC}|{bcolors.GREEN}{stats.greater_one_kb_count:,}{bcolors.ENDC}, Longest contig: {bcolors.GREEN}{stats.longest_contig:,}{bcolors.ENDC}, Total length: {bcolors.BLUE}{stats.total_contig_lengths:,}{bcolors.ENDC}, N50: {bcolors.UNDERLINE}{stats.n50:,}{bcolors.ENDC}, {coverage_title}: {bcolors.YELLOW}{mean_coverage:,.1f}X{bcolors.ENDC}')

class Spades_Parse:

    def __init__(self, fasta_in, coverage_threshold):
        self.records = list(SeqIO.parse(fasta_in, "fasta"))
        file_name = fasta_in.replace(' ', '_')
        self.file_name = re.sub('[._].*', '', file_name)
        self.coverage_threshold = coverage_threshold

    def write_out(self, in_list):
        final_collection = []
        final_list = list(in_list.elements())
        for rec in self.records:
            for out in final_list:
                if out == rec.id:
                    final_collection.append(rec)
        final_collection_count = len(final_collection)
        print(f'The number of genomes in the final collection {final_collection_count}')
        SeqIO.write(final_collection, f'{self.file_name}_selected.fasta', "fasta")
        

    def greater_than_coverage(self, ):
        gtlist = []
        for rec in self.records:
            header = rec.description
            coverage_value = header.split('_')[5]
            coverage_value = int(float(coverage_value))
            if coverage_value > int(self.coverage_threshold):
                gtlist.append(rec.id)
            else:
                print(f'REMOVED: -g sequence length: {len(rec.seq)}, ID: {str(rec.id)}')
        gtlist = Counter(gtlist)
        self.write_out(gtlist)

    def less_than_coverage(self, ):
        ltlist = []
        for rec in self.records:
            header = rec.description
            coverage_value = header.split('_')[5]
            coverage_value = int(float(coverage_value))
            if coverage_value < int(self.coverage_threshold):
                ltlist.append(rec.id)
            else:
                print(f'REMOVED: -l sequence length: {len(rec.seq)}, ID: {str(rec.id)}')
        if args.greater_than_coverage:
            ltlist = Counter(ltlist)
            all_pass = (all_pass & ltlist)
        else:
            ltlist = Counter(ltlist)
        self.write_out(ltlist)

if __name__ == "__main__": # execute if directly access by the interpreter
    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

    ---------------------------------------------------------
    check_downloaded_genomes.py --> $ check_downloaded_genomes.py -s "bos tarus"


    '''), epilog='''---------------------------------------------------------''')

    parser.add_argument('-f', '--fasta_in', action='store', dest='fasta_in', required=True, help='REQUIRED: FASTA file to be processed')
    parser.add_argument('-r1', '--FASTQ_R1', action='store', dest='FASTQ_R1', required=False, default=None, help='R1 FASTQ gz file')
    parser.add_argument('-r2', '--FASTQ_R2', action='store', dest='FASTQ_R2', required=False, default=None, help='R2 FASTQ gz file')
    parser.add_argument('-x', '--write_stats', action='store_true', dest='write_stats', help='OPTIONAL: write stats to Excel')
    parser.add_argument('-c', '--coverage_sizes', action='store_true', dest='coverage_sizes', help='OPTIONAL: sort on coverage and print')
    parser.add_argument('-s', '--length_sizes', action='store_true', dest='length_sizes', help='OPTIONAL: sort on lenth and print')
    parser.add_argument('-g', '--greater_than_coverage', action='store', dest='greater_than_coverage', help='OPTIONAL: KEEP those GREATER THAN this value')
    parser.add_argument('-l', '--less_than_coverage', action='store', dest='less_than_coverage', help='OPTIONAL: KEEP those LESS THAN this value')

    args = parser.parse_args()
    print ("\nSET ARGUMENTS: ")
    print (args)

    fasta_in = args.fasta_in
    stats = Spades_Stats(fasta_in)
    print(f'Contig count: {bcolors.YELLOW}{stats.contig_count:,}{bcolors.ENDC}, Contig length counts <|301-999bp|>: {bcolors.RED}{stats.small_contigs_count:,}{bcolors.ENDC}|{bcolors.BLUE}{stats.mid_size:,}{bcolors.ENDC}|{bcolors.GREEN}{stats.greater_one_kb_count:,}{bcolors.ENDC}, Longest contig: {bcolors.GREEN}{stats.longest_contig:,}{bcolors.ENDC}, Total length: {bcolors.BLUE}{stats.total_contig_lengths:,}{bcolors.ENDC}, N50: {bcolors.UNDERLINE}{stats.n50:,}{bcolors.ENDC}, Mean coverage: {bcolors.YELLOW}{stats.mean_coverage:,.1f}X{bcolors.ENDC}')

    if args.coverage_sizes:
        stats.print_by_coverage()
    if args.length_sizes:
        stats.print_by_length()

    if args.greater_than_coverage:
        parse = Spades_Parse(fasta_in, args.greater_than_coverage)
        parse.greater_than_coverage()
    if args.less_than_coverage:
        parse = Spades_Parse(fasta_in, args.less_than_coverage)
        parse.less_than_coverage()

    # Use FASTQ to calculate coverage based on (read counts * read sized/total assembly length)
    if args.write_stats:
        print(f'{bcolors.RED}Stats writing to Excel{bcolors.ENDC}')
        if args.FASTQ_R1:
            fq = FASTQ_Quality(args.FASTQ_R1, args.FASTQ_R2)
            fq.get_quality()
            stats.write_stats(stats, fq)
        else:
            stats.write_stats(stats, fq=None)

# Created February 2020 by Tod Stuber
