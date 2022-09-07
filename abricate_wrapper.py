#!/usr/bin/env python

__version__ = "0.0.1"

import os
import sys
import re
import argparse
import textwrap
import pandas as pd

class Abricate:
    ''' 
    '''

    def __init__(self, assembly, depth, coverage):
        self.cwd = os.getcwd()
        self.assembly = assembly
        self.depth = depth
        self.coverage = coverage
        basename = os.path.basename(self.assembly)
        self.sample_name = re.sub('[_.].*', '', basename)

    def tab_to_excel(self, tab_file):
        excel_file = tab_file.replace('tab', 'xlsx')
        df = pd.read_csv(tab_file, sep='\t')
        df.to_excel(excel_file)

    def run(self,):
        os.mkdir("abricate")
        os.chdir("abricate")
        sample_name = self.sample_name
        if os.path.isabs(self.assembly):
            assembly = f'{self.assembly}'
        else:
            assembly = f'{self.cwd}/{self.assembly}'
        print("\n*** Running abricate resfinder (default)...")
        os.system(f'abricate --db resfinder --mincov {self.depth} --minid {self.coverage} {assembly} > {sample_name}-resfinder.tab')
        self.tab_to_excel(f'{sample_name}-resfinder.tab')
        print("\n*** Running abricate NCBI...")
        os.system(f'abricate --db ncbi --mincov {self.depth} --minid {self.coverage} {assembly} > {sample_name}-ncbi.tab')
        self.tab_to_excel(f'{sample_name}-ncbi.tab')
        print("\n*** Running abricate plasmidfinder...")
        os.system(f'abricate --db plasmidfinder --mincov {self.depth} --minid {self.coverage} {assembly} > {sample_name}-plasmidfinder.tab')
        self.tab_to_excel(f'{sample_name}-plasmidfinder.tab')
        self.ab_resfinder_file = f'{self.cwd}/abricate/{sample_name}-resfinder.tab'
        self.ab_ncbi_file = f'{self.cwd}/abricate/{sample_name}-ncbi.tab'
        self.ab_plasmid_file = f'{self.cwd}/abricate/{sample_name}-plasmidfinder.tab'
        self.abricate_ab_version = os.popen("abricate --version").readlines()[0].split()[1]
        self.abricate_ncbi_version_date = os.popen("abricate --list").readlines()[4].split()[3] 
        self.abricate_res_version_date = os.popen("abricate --list").readlines()[5].split()[3]
        self.abricate_ncbi_seq_number = os.popen("abricate --list").readlines()[4].split()[1] 
        self.abricate_res_seq_number = os.popen("abricate --list").readlines()[5].split()[1]
        os.chdir(self.cwd)

if __name__ == "__main__": # execute if directly access by the interpreter
    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

    ---------------------------------------------------------
    Place description

    '''), epilog='''---------------------------------------------------------''')
    
    parser.add_argument('-a', '--assembly', action='store', dest='assembly', required=True, help='FASTA assembly file')
    parser.add_argument('-d', '--abricate_depth', action='store', dest='abricate_depth', default=0, help='OPTIONAL: percent average depth cutoff for abricate, aka: --mincov, cvb use -a 50')
    parser.add_argument('-c', '--abricate_coverage', action='store', dest='abricate_coverage', default=75, help='OPTIONAL: percent genome coverage cutoff for abricate, aka: --minid, cvb use -b 90')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.basename(__file__)}: version {__version__}')
    args = parser.parse_args()
    
    print(f'\n{os.path.basename(__file__)} SET ARGUMENTS:')
    print(args)
    print("\n")

    abricate = Abricate(args.assembly, args.abricate_depth, args.abricate_coverage)
    abricate.run()

# Created March 2020 by Tod Stuber