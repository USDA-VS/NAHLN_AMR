#!/usr/bin/env python

import os
import shutil
import re
import gzip
import time
import glob
from datetime import datetime
import argparse
import textwrap
from numpy import mean
from collections import Counter
from concurrent import futures
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator

quality_key = {'!':'0', '"':'1', '#':'2', '$':'3', '%':'4', '&':'5', "'":'6', '(':'7', ')':'8', '*':'9', '+':'10', ',':'11', '-':'12', '.':'13', '/':'14', '0':'15', '1':'16', '2':'17', '3':'18', '4':'19', '5':'20', '6':'21', '7':'22', '8':'23', '9':'24', ':':'25', ';':'26', '<':'27', '=':'28', '>':'29', '?':'30', '@':'31', 'A':'32', 'B':'33', 'C':'34', 'D':'35', 'E':'36', 'F':'37', 'G':'38', 'H':'39', 'I':'40'}


parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

    ---------------------------------------------------------

    Usage: fastq_quality.py -r1 *_R1*gz -r2 *_R2*gz
    Usage: fastq_quality.py -r1 *_R1*gz -r2 *_R2*gz -w # to write passing reads to file

    '''), epilog='''---------------------------------------------------------''')

parser.add_argument('-r1', '--FASTQ_R1', action='store', dest='FASTQ_R1',  required=True, help='Required: provide R1 FASTQ gz file')
parser.add_argument('-r2', '--FASTQ_R2', action='store', dest='FASTQ_R2',  required=True, help='Required: provide R2 FASTQ gz file')
parser.add_argument('-q', '--quality', action='store', dest='qcutoff',  required=False, default=30, help='Quality cutoff')
parser.add_argument('-w', '--write_reads', action='store_true', dest='write_reads', help='Optional: If given reads with an average phred quality >= 30 are written to a file')

args = parser.parse_args()
print ("\nSET ARGUMENTS: ")
print (args)
print("")

qcutoff = int(args.qcutoff)

ts = time.time()
st = datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H-%M-%S')
startTime = datetime.now()

# def efqp_sub_def(base_q):
#     q_list = []
#     mean_quality_list=[]
#     mean_high = 0
#     mean_low = 0 
#     count = 0
#     q_list.append(int(quality_key[base_q]))
#     mean_q = int(mean(q_list))
#     if mean_q >= 30:
#         mean_high+=1
#     else:
#         mean_low+=1
#     mean_quality_list.append(mean_q)
#     count+=1
#     return mean_quality_list, mean_high, mean_low, count

# def even_faster_qual_parser(fastq_file):
#     size_list=[]
#     handle = gzip.open(fastq_file, "rt")
#     for title, seq, qual in FastqGeneralIterator(handle):
#         qual = list(qual)
#         size_list.append(len(seq))
#         with futures.ProcessPoolExecutor() as pool:
#             for mean_quality_list, mean_high, mean_low, count in pool.map(efqp_sub_def, qual):
#                 pass
#     return mean_quality_list, size_list, mean_high, mean_low, count

# def faster_qual_parser(fastq_file):
#     mean_quality_list=[]
#     size_list=[]
#     mean_high = 0
#     mean_low = 0 
#     count = 0
#     handle = gzip.open(fastq_file, "rt")
#     for title, seq, qual in FastqGeneralIterator(handle):
#         q_list = []
#         for base_q in qual:
#             q_list.append(int(quality_key[base_q]))
#         mean_q = int(mean(q_list))
#         if mean_q >= 30:
#             mean_high+=1
#         else:
#                 mean_low+=1
#         mean_quality_list.append(mean_q)
#         size_list.append(len(seq))
#         count+=1
#     return mean_quality_list, size_list, mean_high, mean_low, count

# good_reads = (rec for rec in SeqIO.parse(select_on_header_fastq, "fastq") if mean(rec.letter_annotations["phred_quality"]) > 0)

# mean_quality_list = (int(mean(rec.letter_annotations["phred_quality"])) for rec in SeqIO.parse(handle, "fastq"))
# mean_quality_list = (int(mean(rec.letter_annotations["phred_quality"])) for rec in SeqIO.parse(handle, "fastq"))

# def efqp_sub_def(base_q):
#     q_list = []
#     mean_quality_list=[]
#     mean_high = 0
#     mean_low = 0 
#     count = 0
#     q_list.append(int(quality_key[base_q]))
#     mean_q = int(mean(q_list))
#     if mean_q >= 30:
#         mean_high+=1
#     else:
#         mean_low+=1
#     mean_quality_list.append(mean_q)
#     count+=1
#     return mean_quality_list, mean_high, mean_low, count

# def try_again_mean_quality_size(fastq_file):
#     mean_quality_list=[]
#     size_list=[]
#     mean_high = 0
#     mean_low = 0 
#     # handle = gzip.open(fastq_file, "rt")
#     count = 0
#     fq_dict = SeqIO.index(fastq_file, "fastq")
#     for seq_record in fq_dict.values():
#         read_mean = int(mean(seq_record.letter_annotations["phred_quality"]))
#         if read_mean >= 30:
#             mean_high+=1
#         else:
#             mean_low+=1
#         mean_quality_list.append(read_mean)
#         count+=1
#     return mean_quality_list, mean_high, mean_low, count

# def another_mean_quality_size(fastq_file):
#     mean_quality_list=[]
#     size_list=[]
#     mean_high = 0
#     mean_low = 0 
#     handle = gzip.open(fastq_file, "rt")
#     count = 0
#     read_means = (int(mean(rec.letter_annotations["phred_quality"])) for rec in SeqIO.parse(handle, "fastq"))
#     for read_mean in read_means:
#         if read_mean >= 30:
#             mean_high+=1
#         else:
#             mean_low+=1
#         mean_quality_list.append(read_mean)
#         count+=1
#     return mean_quality_list, mean_high, mean_low, count

def mean_quality_size(fastq_file):
    mean_quality_list=[]
    size_list=[]
    mean_high = 0
    mean_low = 0 
    handle = gzip.open(fastq_file, "rt")
    samplename = re.sub('[._].*', '', fastq_file)
    count = 0
    seq_record_list = []
    for rec in SeqIO.parse(handle, "fastq"):
        # if count < 100:
        mean_q = int(mean(rec.letter_annotations["phred_quality"]))
        if mean_q >= qcutoff:
            mean_high+=1
            seq_record_list.append(rec)
        else:
            mean_low+=1
        mean_quality_list.append(mean_q)
        size_list.append(len(rec))
        count+=1
    if args.write_reads:
        pair_type = re.search('_R[1,2]', fastq_file)
        SeqIO.write(seq_record_list, samplename + "_q" + str(qcutoff) + pair_type.group(0) + ".fastq", "fastq")
    return mean_quality_list, size_list, mean_high, mean_low, count
        # if count < 100:

def sizeof_fmt(num, suffix='B'):
    for unit in ['','K','M','G','T','P','E','Z']:
        if abs(num) < 1024.0:
            return "%3.1f%s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f%s%s" % (num, 'Yi', suffix)

def pairing_intergrety(fastq1, fastq2):
    fastq1_header_list = []
    handle = open(fastq1, "r")
    for title, seq, qual in FastqGeneralIterator(handle):
        fastq1_header_list.append(title.rsplit()[0])
    handle.close()

    fastq2_header_list = []
    handle = open(fastq2, "r")
    for title, seq, qual in FastqGeneralIterator(handle):
        fastq2_header_list.append(title.rsplit()[0])
    handle.close()

    print("\nFinding unbalanced headers")
    unbalanced_headers = list(set(fastq1_header_list).symmetric_difference(set(fastq2_header_list)))

    print("Balancing R1")
    handle = open(fastq1, "r")
    write_out_fastq_p = open("out1", "w")
    for title, seq, qual in FastqGeneralIterator(handle):
        if title.rsplit()[0] not in unbalanced_headers:
            print("@{}\n{}\n+\n{}" .format(title, seq, qual), file=write_out_fastq_p)
    handle.close()
    write_out_fastq_p.close()
    shutil.move("out1", fastq1)

    print("Balancing R2")
    handle = open(fastq2, "r")
    write_out_fastq_p = open("out2", "w")
    for title, seq, qual in FastqGeneralIterator(handle):
        if title.rsplit()[0] not in unbalanced_headers:
            print("@{}\n{}\n+\n{}" .format(title, seq, qual), file=write_out_fastq_p)
    handle.close()
    write_out_fastq_p.close()
    shutil.move("out2", fastq2)

FASTQ_R1 = args.FASTQ_R1
FASTQ_R2 = args.FASTQ_R2

sample_name = re.sub('_.*', '', os.path.basename(FASTQ_R1))
R1size = sizeof_fmt(os.path.getsize(FASTQ_R1))
R2size = sizeof_fmt(os.path.getsize(FASTQ_R2))

stat_summary = {}
stat_summary["time_stamp"] = st
stat_summary["sample_name"] = sample_name
stat_summary["R1size"] = R1size
stat_summary["R2size"] = R2size

print("Getting mean for {}" .format(FASTQ_R1))
mean_quality_list1, size_list1, mean_high1, mean_low1, count1 = mean_quality_size(FASTQ_R1)
# mean_quality_list1, size_list1, mean_high1, mean_low1, count1 = faster_qual_parser(FASTQ_R1) #not faster
# mean_quality_list1, size_list1, mean_high1, mean_low1, count1 = even_faster_qual_parser(FASTQ_R1) #not faster
# mean_quality_list1, mean_high1, mean_low1, count1 = try_again_mean_quality_size(FASTQ_R1)
print("Getting mean for {}" .format(FASTQ_R2))
mean_quality_list2, size_list2, mean_high2, mean_low2, count2 = mean_quality_size(FASTQ_R2)
# mean_quality_list2, size_list2, mean_high2, mean_low2, count2 = faster_qual_parser(FASTQ_R2) #not faster
# mean_quality_list2, size_list2, mean_high2, mean_low2, count2 = even_faster_qual_parser(FASTQ_R2) #not faster
# mean_quality_list2, mean_high2, mean_low2, count2 = try_again_mean_quality_size(FASTQ_R2)

if args.write_reads:
    q_fastq1 = glob.glob("*_q" + str(qcutoff) + "*_R1*")
    q_fastq2 = glob.glob("*_q" + str(qcutoff) + "*_R2*")
    pairing_intergrety(q_fastq1[0], q_fastq2[0])

mean_quality_list = mean_quality_list1 + mean_quality_list2
size_list = size_list1 + size_list2
mean_quality_bin = Counter(mean_quality_list)
size_bin = Counter(size_list)

mean_quality1 = mean(mean_quality_list1)
mean_quality2 = mean(mean_quality_list2)
combined_mean_quality = mean(mean_quality1 + mean_quality2)/2

size_list1 = mean(size_list1)
size_list2 = mean(size_list2)
combined_size_list = mean(size_list1 + size_list2)/2

mean_high = mean_high1 + mean_high2
mean_low = mean_low1 + mean_low2

count = count1 + count2

print("\n\nsize_list1: {:.1f}" .format(size_list1))
print("size_list2: {:.1f}" .format(size_list2))
print("combined_size_list: {:.1f}" .format(combined_size_list))

print("\nmean_quality1: {:.1f}" .format(mean_quality1))
print("mean_quality2: {:.1f}" .format(mean_quality2))
print("combined_mean_quality: {:.1f}" .format(combined_mean_quality))

freq_greater_than_qcutoff = mean_high/count
freq_less_than_qcutoff = mean_low/count

print("\n--Q" + str(qcutoff) + " cutoff--")
print("Greater or equal:")
print("\tRead1: {:.1%}" .format(mean_high1/count1))
print("\tRead2: {:.1%}" .format(mean_high2/count2))
print("\tAverage: {:.1%}" .format(freq_greater_than_qcutoff))
print("\nAverage less than: {:.1%}" .format(freq_less_than_qcutoff))

#print("\tReads selected on TB complex taxon number {:,} -- {:.1%}" .format(selected_on_tax, perc))
runtime = (datetime.now() - startTime)
print ("\n\nruntime: %s:  \n" % runtime)

with open(sample_name + "_stat_quality" + str(qcutoff) + ".txt", "w") as write_out:
    print("{}\t{}\t{}\t{}\t{}\t{}" .format("mean_quality1", "mean_quality2", "q" + str(qcutoff) + "pass1", "q" + str(qcutoff) + "pass2", "count1", "count2", "readlen1", "readlen2", "mean_readlen"), file=write_out)
    print("{:.1f}\t{:.1f}\t{:.1%}\t{:.1%}\t{}\t{}" .format(mean_quality1, mean_quality2, mean_high1/count1, mean_high2/count2, count1, count2, size_list1, size_list2, combined_size_list), file=write_out)

print("Done")
