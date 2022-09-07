#!/usr/bin/env python

import os
import subprocess
import re
import time
import sys
import re
import shutil
import gzip
import glob
import argparse
import textwrap
import numpy as np
import pandas as pd
import humanize
import json
from numpy import mean
from collections import defaultdict
from collections import Counter
from collections import OrderedDict
import time
from datetime import datetime

from vsnp_fastq_quality import FASTQ_Quality
from spades_stats_parse import Spades_Stats

class Spades_Assembly:
    ''' 
    Assemble reads using Spades assembler.
    Paired or single reads
    '''

    def __init__(self, FASTQ_R1, FASTQ_R2, debug=False):
        '''
        See -h
        '''
        FASTQ_list = [FASTQ_R1, FASTQ_R2]
        FASTQ_list = [x for x in FASTQ_list if x is not None]  # remove None when single read
        sample_name = re.sub('[._].*', '', FASTQ_list[0])
        cwd = os.getcwd()

        print(f'SPAdes Running...')
        self.assembly_file = f'{sample_name}_spades_assembly.fasta'
        if len(FASTQ_list) == 2:
            subprocess.run(["spades.py", "-1", FASTQ_list[0], "-2", FASTQ_list[1], "-o", "spades_assembly"], capture_output=True) # default -t 16
        elif len(FASTQ_list) == 1: # assume iontorrent
            subprocess.run(["spades.py", "--iontorrent", "-s", FASTQ_list[0], "-o", "spades_assembly"], capture_output=True)
        else:
            print(f'\n### Must have either single or paired read set.\n')
            sys.exit(0)
        if os.path.exists(f'{cwd}/spades_assembly/scaffolds.fasta'):
            shutil.copy2(f'{cwd}/spades_assembly/scaffolds.fasta', f'{cwd}/{sample_name}.fasta')
            self.assembly_file = f'{cwd}/{sample_name}.fasta'
        else:
            print(f'\n### SPAdes did not complete, see log\n')
            sys.exit(0)
        
        if not debug:
            shutil.rmtree(f'{cwd}/spades_assembly')
        

if __name__ == "__main__": # execute if directly access by the interpreter

    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

        ---------------------------------------------------------
        from spades_assembly import Spades_Assembly
        assemble = Spades_Assembly("01-3941_kp_q25_R1.fastq.gz", "01-3941_kp_q25_R2.fastq.gz")
        assemble.assembly_file

        '''), epilog='''---------------------------------------------------------''')

    parser.add_argument('-r1', '--FASTQ_R1', action='store', dest='FASTQ_R1', required=True, help='Required: provide R1 FASTQ gz file')
    parser.add_argument('-r2', '--FASTQ_R2', action='store', dest='FASTQ_R2', required=False, default=None, help='Optional: provide R2 FASTQ gz file')
    parser.add_argument('-d', '--debug', action='store_true', dest='debug', help='keep spades output directory')

    args = parser.parse_args()
    print ("\nSET ARGUMENTS: ")
    print (args)
    print("")

    fq = FASTQ_Quality(args.FASTQ_R1, args.FASTQ_R2)
    fq.get_quality()
    assemble = Spades_Assembly(args.FASTQ_R1, args.FASTQ_R2, args.debug)
    stats = Spades_Stats(assemble.assembly_file)
    stats.write_stats(stats, fq)

# Created November 2020 by Tod Stuber