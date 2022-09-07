#!/usr/bin/env python

import re
import glob
import pandas as pd
from sys import argv
import time
from datetime import datetime

databases = ["resfinder", "ncbi", "plasmidfinder", "amrfinder"]
ts = time.time()
st = datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H-%M-%S')

for database in databases:

    text_files = glob.glob("*" + database + "*txt")
    text_files.sort()
    tab_files = glob.glob("*" + database + "*tab")
    tab_files.sort()

    #df['BrandName'].replace(['ABC', 'AB'], 'A')

    file_out = 'abricate-' + database + "_" + st + '.xlsx'
    #Card file is the .txt file which is tab delimited
    writer = pd.ExcelWriter(file_out)
    for text_file in text_files:
        name=re.sub('[.].*', '', text_file)[:31]
        df = pd.read_csv(text_file, sep='\t', index_col='#FILE')
        df.to_excel(writer, name)
        print ("Text file, {}, has been placed into Excel file {}" .format(text_file, file_out))

    for tab_file in tab_files:
        in_file = tab_file
        name=re.sub('[.].*', '', tab_file)[:31]
        df = pd.read_csv(tab_file, sep='\t', index_col='SEQUENCE')
        df.replace(to_replace="^=", value="'=", regex=True, inplace=True)
        df.to_excel(writer, name)
        print ("Text file, {}, has been placed into Excel file {}" .format(tab_file, file_out))

    writer.save()
