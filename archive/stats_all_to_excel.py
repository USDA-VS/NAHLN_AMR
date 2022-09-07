#!/usr/bin/env python

import re
import glob
import pandas as pd
from sys import argv
import time
from datetime import datetime

assembly_files = glob.glob("*_merged.txt")
all_assembly_stats = "batch_allmerged.txt"
ts = time.time()
st = datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H-%M-%S')

for assembly_file in assembly_files:
    name=re.sub('_merged.txt', '', assembly_file)[:31]
    file_out = 'stat-' + name + '_' + st + '.xlsx'
    writer = pd.ExcelWriter(file_out)
    df = pd.read_csv(assembly_file, sep='\t', index_col='Sample')
    df.to_excel(writer, name)
    df = pd.read_csv(all_assembly_stats, sep='\t', index_col='Sample')
    df.to_excel(writer, 'batch_assembly_stats')
    print ("Text file, {}, has been placed into Excel file {}" .format(assembly_file, file_out))
    writer.save()

