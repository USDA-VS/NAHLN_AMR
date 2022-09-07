#!/usr/bin/env python

import re
import glob
import pandas as pd
from sys import argv
import time
from datetime import datetime

#Working directory must be ISEScan prediction folder

isescan_files = glob.glob("*")
#remove orf files from list
isescan_files = [file for file in isescan_files if not re.search(r'orf', file)]
ts = time.time()
st = datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H-%M-%S')
file_name=re.sub('[._].*', '', isescan_files[0])

file_out = 'isescan-' + file_name + '-' + st + '.xlsx'
writer = pd.ExcelWriter(file_out)

#enter list items individually to order Excel tabs
tab_in = [file for file in isescan_files if re.search(r'sum', file)][0]
tab_name=tab_in[-31:] #keep the last 31 characters in string
df = pd.read_csv(tab_in, delim_whitespace=True, index_col='seqid')
df.to_excel(writer, tab_name)
print("Done: {}\n\n" .format(tab_in))

tab_in = [file for file in isescan_files if re.search(r'raw', file)][0]
tab_name=tab_in[-31:] #keep the last 31 characters in string
df = pd.read_csv(tab_in, delim_whitespace=True, index_col='seqID')
df.to_excel(writer, tab_name)
print("Done: {}\n\n" .format(tab_in))

tab_in = [file for file in isescan_files if re.search(r'out', file)][0]
tab_name=tab_in[-31:] #keep the last 31 characters in string
df = pd.read_csv(tab_in, delim_whitespace=True, index_col='seqID')
df.to_excel(writer, tab_name)
print("Done: {}\n\n" .format(tab_in))

tab_in = [file for file in isescan_files if re.search(r'gff', file)][0]
tab_name=tab_in[-31:] #keep the last 31 characters in string
df = pd.read_csv(tab_in, delim_whitespace=True)
df.to_excel(writer, tab_name)
print("Done: {}\n\n" .format(tab_in))

tab_in = [file for file in isescan_files if re.search(r'is.fna', file)][0]
tab_name=tab_in[-31:] #keep the last 31 characters in string
df = pd.read_csv(tab_in, delim_whitespace=True, header=None)
df.to_excel(writer, tab_name)
print("Done: {}\n\n" .format(tab_in))

writer.save()



