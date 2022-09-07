#!/usr/bin/env python

__version__ = "0.0.1"

import os
import sys
import shutil
import re
import glob
import subprocess
import argparse
import textwrap
import pandas as pd

from amr_file_setup import Setup
from amr_file_setup import bcolors
from amr_file_setup import Banner
from amr_file_setup import Latex_Report
from amr_file_setup import Excel_Stats

class Sourmash(Setup):
    ''' 
    '''

    def __init__(self, FASTA=None, debug=False):
        '''
        sourmash sketch dna *fasta
        sourmash search Bcanis_B-REF-BC-RM666.fasta.sig /project/bioinformatic_databases/sourmash/genbank-k31.sbt.json --save-match search.findings.sig -o search.csv
        sourmash compare -k 31 search.findings.sig -o compare --csv compare_full_distance_matrix.csv
        sourmash plot --subsample 50 --labels --csv plot_downsample_distance_matrix.csv compare

        fastANI #one to many
        https://github.com/ParBLiSS/FastANI
        realpath *fasta > sample_list.txt
        fastANI -q LR9-5583-USVI-K-Ms-Survey.fasta --rl sample_list.txt -o fastani.out
        '''
        Setup.__init__(self, FASTA=FASTA, debug=debug)
        self.print_run_time('Sourmash')
        subprocess.run(["sourmash", "sketch", "dna", FASTA, ], capture_output=True)
        fasta_basename = os.path.basename(FASTA)
        subprocess.run(["sourmash", "search", f'{fasta_basename}.sig', "/project/bioinformatic_databases/sourmash/genbank-k31.sbt.json", "--save-match", "search.findings.sig", "-o", f'{self.sample_name}_search.csv', '--containment'], )
        self.sourmash_df = pd.read_csv(f'{self.sample_name}_search.csv')
        dir = 'sourmash'
        if not os.path.exists(dir):
            os.makedirs(dir)
        files_list = []
        for files in ('*sig', '*search.csv'):
            files_list.extend(glob.glob(files))
        for each in files_list:
            shutil.move(each, dir)
    
    def latex(self, tex):
        blast_banner = Banner("Sourmash Sequence Similarity")
        print(r'\begin{table}[ht!]', file=tex)
        print(r'\begin{adjustbox}{width=1\textwidth}', file=tex)
        print(r'\begin{center}', file=tex)
        print('\includegraphics[scale=1]{' + blast_banner.banner + '}', file=tex)
        print(r'\end{center}', file=tex)
        print(r'\end{adjustbox}', file=tex)
        print(r'\begin{adjustbox}{width=1\textwidth}', file=tex)
        print(r'\begin{tabular}{l|l}', file=tex)
        print(r'Similarity & ID \\', file=tex)
        print(r'\hline', file=tex)
        count=0
        if self.sourmash_df.size == 0:
            print('Sourmash - No Data Output & Sourmash - No Data Output \\\\', file=tex)
        else:
            for row in self.sourmash_df.itertuples():
                count+=1
                if count <= 10:
                    percetage = f'{row[1]:.1%}'
                    print(percetage.replace("%", "\%") + ' & ' + row[2].replace("_", "\_") + ' \\\\', file=tex)
                    print(r'\hline', file=tex)
        print(r'\end{tabular}', file=tex)
        print(r'\end{adjustbox}', file=tex)
        print(r'\\', file=tex)
        print(r'\end{table}', file=tex)

    def excel(self, excel_dict):
        count=0
        for row in self.sourmash_df.itertuples():
            count+=1
            if count <= 1:
                excel_dict[f'Sourmash Sequence Similarity'] = f'{row[1]:.1%}:{row[2]}'

if __name__ == "__main__": # execute if directly access by the interpreter
    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

    ---------------------------------------------------------
    usda_sourmash.py -f *fasta

    '''), epilog='''---------------------------------------------------------''')
    parser.add_argument('-f', '--fasta', action='store', dest='FASTA', default=None, help='provide assembly if just stats are needed.  Assumed SPAdes assembly')
    parser.add_argument('-d', '--debug', action='store_true', dest='debug', default=False, help='keep temp file')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.basename(__file__)}: version {__version__}')
    args = parser.parse_args()
    
    print(f'\n{os.path.basename(__file__)} SET ARGUMENTS:')
    print(args)
    print("\n")

    sourmash = Sourmash(FASTA=args.FASTA,)

    #Latex report
    latex_report = Latex_Report(sourmash.sample_name)
    sourmash.latex(latex_report.tex)
    latex_report.latex_ending()

    #Excel Stats
    excel_stats = Excel_Stats(sourmash.sample_name)
    sourmash.excel(excel_stats.excel_dict)
    excel_stats.post_excel()

    temp_dir = './temp'
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    files_grab = []
    for files in ('*.aux', '*.log', '*tex', '*png', '*out'):
        files_grab.extend(glob.glob(files))
    for each in files_grab:
        shutil.move(each, temp_dir)

    if args.debug is False:
        shutil.rmtree(temp_dir)

# Created 2021 by Tod Stuber
