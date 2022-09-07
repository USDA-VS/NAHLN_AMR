#!/usr/bin/env python

__version__ = "0.0.1"

import os
import sys
import re
import json
import argparse
import textwrap
import pandas as pd
import time
from datetime import datetime

class MLST:
    ''' 
    '''

    def __init__(self, assembly):
        self.cwd = os.getcwd()
        self.assembly = assembly
        base_name = os.path.basename(self.assembly)
        self.sample_name = re.sub('[_.].*', '', base_name)
    
    def run(self,):
        os.mkdir("mlst")
        os.chdir("mlst")
        sample_name = self.sample_name
        print("\n*** Running MLST...")
        if os.path.isabs(self.assembly):
            assembly = f'{self.assembly}'
        else:
            assembly = f'{self.cwd}/{self.assembly}'
        os.system(f'mlst {assembly} > {self.sample_name}_mlst.txt')
        mlst_detail=[]
        mlst_file = f'{self.cwd}/mlst/{self.sample_name}_mlst.txt'
        if os.path.isfile(mlst_file) and not os.stat(mlst_file).st_size == 0:
            self.mlst_file = mlst_file
            mlst_file = self.mlst_file
            with open(mlst_file, 'r') as mlst_out:
                whole_file = mlst_out.readlines()
                mlst_out_list = whole_file[0].split()
                self.mlst_scheme = mlst_out_list[1]
                self.mlst_type = mlst_out_list[2]
                for mlst_call in mlst_out_list[3:]:
                    mlst_detail.append(mlst_call)
                self.mlst_detail = mlst_detail

        df = pd.DataFrame(index=[sample_name], columns=[ 'MLST scheme', 'MLST type', 'MLST detail', ])
        df.at[sample_name, 'MLST scheme'] = f'{self.mlst_scheme}'
        df.at[sample_name, 'MLST type'] = f'{self.mlst_type}'
        df.at[sample_name, 'MLST detail'] = f'{self.mlst_detail}'
        df.index.name = 'sample'
        ts = time.time()
        st = datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H-%M-%S')
        df.to_excel(f'{sample_name}_{st}_mlst.xlsx')
        self.mlst_excel = f'{os.getcwd()}/{sample_name}_{st}_mlst.xlsx'
        real_path = os.path.dirname(os.path.realpath(__file__))
        lookup_genome_size = real_path + "/lookup_genome_size.json"
        with open(lookup_genome_size) as infile:
            lookup_dict = json.load(infile)
        try:
            self.mlst_species_lookup = lookup_dict.get(self.mlst_scheme, {}).get("Scientific Name", None)
            self.mlst_size_lookup = lookup_dict.get(self.mlst_scheme, {}).get("Approx Genome Size (Mb)", None)
        except AttributeError:
            self.mlst_species_lookup = None
            self.mlst_size_lookup = None

        self.version = os.popen("mlst --version").readlines()[0].split()[1] 
        os.chdir(self.cwd)


if __name__ == "__main__": # execute if directly access by the interpreter
    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

    ---------------------------------------------------------
    MLST OVERVIEW:
    Summary:
        Allele IDs (as FASTAs) are used to make a BLAST database.  Top hits used to find pubMLST scheme and allele IDs are reported.
    Usage: 
        mlst contigs.fasta
    Output: 
        contigs.fa  neisseria  11149  abcZ(672) adk(3) aroE(4) fumC(3) gdh(8) pdhC(4) pgm(6)
            1) filename
            2) matching PubMLST scheme name
            3) ST (sequence type)
            4) allele IDs
    Show available schemes:
        mlst --longlist
        Database of schemes: ~/.conda/envs/tod/db/pubmlst
            each scheme folder has a text file from making the final ST call based on allele ID findings.  For each gene type there is a .tfa file listing all all allele IDs for that gene, ie abcZ(671, 2, 3)... etc.
        Schemes can be mapped to genus/specie identification via scheme using ~/.conda/envs/tod/db/scheme_species_map.tab
    Allele ID coding:
        n = 100% match
        ~n = novel full length allele similar to n
        n? = partial match to known allele
        - = allele missing
        n, m = multiple alleles
    Update MLST Database:
        Follow instructions here to update: https://github.com/tseemann/mlst
        Script are here: ~/.conda/envs/tod/scripts/mlst-download_pub_mlst

    '''), epilog='''---------------------------------------------------------''')
    
    parser.add_argument('-a', '--assembly', action='store', dest='assembly', required=True, help='FASTA assembly file')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.basename(__file__)}: version {__version__}')
    args = parser.parse_args()
    
    print(f'\n{os.path.basename(__file__)} SET ARGUMENTS:')
    print(args)

    mlst = MLST(args.assembly)
    mlst.run()

    print('done')

# Created March 2020 by Tod Stuber