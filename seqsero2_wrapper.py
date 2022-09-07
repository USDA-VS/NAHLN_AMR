#!/usr/bin/env python

__version__ = "0.0.1"

import os
import sys
import re
import argparse
import textwrap
import pandas as pd

class SeqSero2:
    ''' 
    '''

    def __init__(self, FASTQ_R1, FASTQ_R2):
        self.cwd = os.getcwd()
        self.FASTQ_R1 = FASTQ_R1
        self.FASTQ_R2 = FASTQ_R2
        self.sample_name = re.sub('[_.].*', '', self.FASTQ_R1)
        # Set defaults to None
        self.results = None
        self.serotype = None
        self.antigenic = None
        self.subspecies = None
        self.version = None
        self.seqserocomment = None
        self.seqsero_file = None

    def run(self,):
        print(f'\n*** Running SeqSero2...')
        if self.FASTQ_R2 == 'None':
            os.system(f'SeqSero2_package.py -d seqsero -t 3 -i {s.FASTQ_R1}')
        else:
            os.system(f'SeqSero2_package.py -d seqsero -t 2 -i {self.FASTQ_R1} {self.FASTQ_R2}')
        if os.path.isfile(f'{self.cwd}/seqsero/SeqSero_result.txt'):
            self.seqsero_file = f'{self.cwd}/seqsero/SeqSero_result.txt'
            seqsero_results=[] 
            with open(f'{self.cwd}/seqsero/SeqSero_result.txt', 'r') as infile: 
                for line in infile: 
                    line = line.strip('\n') 
                    line_split = line.split('\t') 
                    seqsero_results.append(line_split)
            self.results = seqsero_results
            self.serotype = seqsero_results[7][1]
            self.antigenic = seqsero_results[6][1]
            self.subspecies = seqsero_results[5][1]
            try:
                self.seqserocomment = seqsero_results[8][0] + " " + seqsero_results[8][1]
            except IndexError:
                self.seqserocomment = seqsero_results[8][0]
            self.version = os.popen("SeqSero2_package.py --version").readlines()[0].split()[1]
        else:
            print(f'\t####### {self.sample_name} #########\n\tSEQSERO2 FAILED\n\t####### {self.sample_name} #########')
            self.seqsero_results = None


if __name__ == "__main__": # execute if directly access by the interpreter
    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

    ---------------------------------------------------------
    SEQSERO2 OVERVIEW:
        - Salmonella serotypes are defined by two surface structures, O antigen and H antigen. 2,600 Salmonella serotypes have been described in the White-Kauffmann-Le Minor scheme
        - Molecular methods for serotype determination are based on genes responsible for serotype antigens, including the rfb  gene cluster (O antigen), fliC, and fljB (H antigen).
    Workflow:
        - Summary: Align -> Assemble -> BLAST -> cross-reference top hits to schemes
        - Workflow assembles serotype determinants, i.e., wzx and wzy for O antigen identification and fliC and fljB for H antigen identification.  All sequencing reads from a query genome were mapped to the serotype determinant database using BWA-MEM.  Reads that had been mapped to alleles in the database were extracted and assembled by SPAdes.  The resulting contigs were then aligned back to the serotype determinant database by BLAST.  The type of O or H antigen of the query genome was determined by the type of the corresponding allele that yielded the highest BLAST similarity. The final serotype of the query genome was called by consulting the White-Kauffmann-Le Minor scheme.
    Options:
        SeqSero2_package.py -h
    Modes:
        1) Allele micro-assembly (default) using FASTQs, best results
        2) kmer based using FASTQs
        3) kmer based using FASTA
    Contamination:
        - Genome contaminated with a second Salmonella serotype might be detected by the presence of 2 or more O antigen or 3 or more H antigen calls.
        - Micro-assembly workflow provides additional information regarding sequence matches to the SeqSero serotype determinant databases, which can assist in detecting atypical results due to interserotype contamination or to divergent flagellin alleles.  Besides serotype prediction, the microassembly workflow reports all assembled O and H antigen alleles, along with their sequence similarity scores to alleles in the serotype determinant databases (number of base matches/length of allele) in order to judge how good the matches are.  A genome contaminated with a second Salmonella serotype might be detected by the presence of 2 or more O antigen or 3 or more H antigen calls.  Sample contamination detection using the microassembly workflow is a good quality control tool for both in-house sequencing and public data sourcing.  It should be noted that SeqSero2 will not detect contamination from the same serotype or from a non-Salmonella organism. Also, contamination detection requires the presence of serotype determinant sequences, e.g., rfb, fliC, or fljB, from the contaminant genome. At low levels of contamination, such reads may be absent.
    Datbase:
        /software/opt/spack/linux-centos7-cascadelake/gcc-9.2.0/seqsero2-1.0.2-icaywckz4zwcwpnwjuawqhkn32n24wnw/seqsero2_db

    '''), epilog='''---------------------------------------------------------''')
    
    parser.add_argument('-r1', '--FASTQ_R1', action='store', dest='FASTQ_R1', required=True, help='R1 FASTQ gz file')
    parser.add_argument('-r2', '--FASTQ_R2', action='store', dest='FASTQ_R2', required=True, default=None, help='R2 FASTQ gz file')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.basename(__file__)}: version {__version__}')
    args = parser.parse_args()
    
    print(f'\n{os.path.basename(__file__)} SET ARGUMENTS:')
    print(args)
    print("\n")

    seqsero2 = SeqSero2(args.FASTQ_R1, args.FASTQ_R2)
    seqsero2.run()

    # print(f'done')

# Created March 2020 by Tod Stuber
