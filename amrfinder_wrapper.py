#!/usr/bin/env python

__version__ = "0.0.1"

import os
import sys
import re
import argparse
import textwrap
import pandas as pd

class AMR_Finder:
    ''' 
    '''

    def __init__(self, assembly):
        self.cwd = os.getcwd()
        self.assembly = assembly
        basename = os.path.basename(self.assembly)
        self.sample_name = re.sub('[_.].*', '', basename)

    def tab_to_excel(self, tab_file):
        excel_file = tab_file.replace('tab', 'xlsx')
        df = pd.read_csv(tab_file, sep='\t')
        df.to_excel(excel_file)

    def run(self,):
        sample_name = self.sample_name
        os.mkdir("amrfinder")
        os.chdir("amrfinder")
        print("\n*** Running amrfinder...")
        if os.path.isabs(self.assembly):
            assembly = f'{self.assembly}'
        else:
            assembly = f'{self.cwd}/{self.assembly}'
        # os.system(f'amrfinder.pl -n {assembly} >> {sample_name}-amrfinder.tab')
        os.system(f'amrfinder --nucleotide {assembly} --output {sample_name}-amrfinder.tab')
        self.tab_to_excel(f'{sample_name}-amrfinder.tab')
        self.amrfinder_file = f'{self.cwd}/amrfinder/{sample_name}-amrfinder.tab'
        os.chdir(self.cwd)
        self.version = os.popen("amrfinder --version").readlines()[0].split()[0]
        

if __name__ == "__main__": # execute if directly access by the interpreter
    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

    ---------------------------------------------------------
    Place description

    '''), epilog='''---------------------------------------------------------''')
    
    parser.add_argument('-a', '--assembly', action='store', dest='assembly', required=True, help='FASTA assembly file')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.basename(__file__)}: version {__version__}')
    args = parser.parse_args()
    
    print(f'\n{os.path.basename(__file__)} SET ARGUMENTS:')
    print(args)
    print("\n")

    amr_finder = AMR_Finder(args.assembly)
    amr_finder.run()

# Created March 2020 by Tod Stuber