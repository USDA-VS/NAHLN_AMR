#!/usr/bin/env python

import re
import pandas as pd
from sys import argv
import time
from datetime import datetime

ts = time.time()
st = datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H-%M-%S')

#Card file is the .txt file which is tab delimited
card_file = argv[1]
name=re.sub('[_.].*', '', card_file)
df = pd.read_csv(card_file, sep='\t', index_col='Contig')
df.to_excel('card-' + name + '_' + st + '.xlsx')

print ("Text file, {}, has been placed into Excel" .format(card_file))
