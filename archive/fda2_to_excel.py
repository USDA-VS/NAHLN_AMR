#!/usr/bin/env python

import re
import pandas as pd
from sys import argv
import time
from datetime import datetime

ts = time.time()
st = datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H-%M-%S')


file = ['class.txt','all_sort.out','blast5050','output_file_ed21.txt','study_genelist.txt']

file_out = 'fda-' + st + '.xlsx'

writer = pd.ExcelWriter(file_out)
for index, f in enumerate(file, start = 0):
    name=re.sub('[.].*', '', file[index])
    #df = pd.read_csv(f, sep='\t')
    if f == 'blast5050' or f == 'all_sort.out':
        df = pd.read_csv(f, sep='\t', header=None, index_col='qseqid', names=['qseqid','sseqid','pident','length','qlen','slen','qstart','qend','sstart','send','gap','bitscore'])
        print ("Column headers added to {}" .format(f))
    elif f == 'study_genelist.txt':
        df = pd.read_csv(f, sep='\t', header=None, index_col='isolate', names=['isolate','gene list'])
        print ("Column headers added to {}" .format(f)) 
    else:
        df = pd.read_csv(f, sep='\t', index_col=0)
    df.to_excel(writer, name)
    print ("Text file, {}, has been placed into Excel file {}" .format(f, file_out))

#Argument must be fda class file
#usage: fda_to_excel.py class.txt

#fda_class_file = argv[1]
#name=re.sub('[_.].*', '', fda_class_file)
#file_out = 'fda-' + name + '_' + st + '.xlsx'
#df = pd.read_csv(fda_class_file, sep='\t', index_col=0)
#df.to_excel(file_out)