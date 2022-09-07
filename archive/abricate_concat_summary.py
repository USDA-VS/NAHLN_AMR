#!/usr/bin/env python

import re
import glob
import pandas as pd
from sys import argv
import time
from datetime import datetime
import argparse
import textwrap

parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

---------------------------------------------------------

Usage:
$ mkdir ~/temp; rm -r ~/temp/*
$ find /bioinfo11/TStuber/Results/bi/data/ -name "abricate-ncbi*.xlsx" -exec cp {} ~/temp \;
$ cd ~/temp
~/temp $ abricate_ncbi_concat.py
~/temp $ cp abricate-ncbi-concat_summary*xlsx /bioinfo11/TStuber/Results/bi/stats

'''), epilog='''---------------------------------------------------------''')

args = parser.parse_args()
print ("\nSET ARGUMENTS: ")
print (args)
print("")

ts = time.time()
st = datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H-%M-%S')
databases = ["resfinder", "ncbi", "plasmidfinder"]


for database in databases:

    excel_files = glob.glob("abricate-" + database + "*.xlsx")
    file_out = 'abricate-' + database + '-concat_summary' + "_" + st + '.xlsx'

    frames = []
    sname = database + "_summary"
    for excel_file in excel_files:
        df = pd.read_excel(excel_file, sheet_name=0, index_col='#FILE')
        frames.append(df)

    results = pd.concat(frames)
    results.sort_index(inplace=True)
    cols = list(results.columns.values)
    cols.pop(cols.index('NUM_FOUND'))
    results = results[['NUM_FOUND'] + cols]
    results.to_excel(file_out, sheet_name='summary_' + st)

