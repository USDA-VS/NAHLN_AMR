#!/usr/bin/env python
import os
import re
import glob
import pandas as pd
from sys import argv
import time
from datetime import datetime
import argparse
import textwrap


def merge_files(sample_name, assembly_stats, quality_file):
    ts = time.time()
    st = datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H-%M-%S')

    merged_df = {}
    quality_file_df = pd.read_csv(quality_file, sep='\t')
    quality_file_df_T = quality_file_df.T
    for idx, row in quality_file_df_T.iterrows():
        merged_df[idx] = row[0]

    assembly_df = pd.read_csv(assembly_stats, sep='\t')
    assembly_df_T = assembly_df.T
    for idx, row in assembly_df_T.iterrows():
        merged_df[idx] = row[0]

    file_out = sample_name + '_merged.txt'
    with open(file_out, "w") as write_out:
        print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}" .format("Sample", "File_R1", "Size_R1", "mean_quality1", "q30pass1", "Count_R1", "File_R2", "Size_R2", "mean_quality2", "q30pass1", "Count_R2", "Number of scaffolds", "Total size of scaffolds", "Longest scaffold", "Number of scaffolds > 1K nt", "N50 scaffold length", "L50 scaffold count"), file=write_out)
        print(merged_df.get("Sample", "n/a"), end="\t", file=write_out)
        print(merged_df.get("File_R1", "n/a"), end="\t", file=write_out)
        print(merged_df.get("Size_R1", "n/a"), end="\t", file=write_out)
        print(merged_df.get("mean_quality1", "n/a"), end="\t", file=write_out)
        print(merged_df.get("q30pass1", "n/a"), end="\t", file=write_out)
        print(merged_df.get("count1", "n/a"), end="\t", file=write_out)
        print(merged_df.get("File_R2", "n/a"), end="\t", file=write_out)
        print(merged_df.get("Size_R2", "n/a"), end="\t", file=write_out)
        print(merged_df.get("mean_quality2", "n/a"), end="\t", file=write_out)
        print(merged_df.get("q30pass2", "n/a"), end="\t", file=write_out)
        print(merged_df.get("count2", "n/a"), end="\t", file=write_out)
        print(merged_df.get("Number of scaffolds", "n/a"), end="\t", file=write_out)
        print(merged_df.get("Total size of scaffolds", "n/a"), end="\t", file=write_out)
        print(merged_df.get("Longest scaffold", "n/a"), end="\t", file=write_out)
        print(merged_df.get("Number of scaffolds > 1K nt", "n/a"), end="\t", file=write_out)
        print(merged_df.get("N50 scaffold length", "n/a"), end="\t", file=write_out)
        print(merged_df.get("L50 scaffold count", "n/a"), file=write_out)

    print("\nDone")

if __name__ == "__main__": 
    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

        ---------------------------------------------------------

        Usage: merge_quality_files.py -r1 *_R1*gz -r2 *_R2*gz

        '''), epilog='''---------------------------------------------------------''')

    parser.add_argument('-a', '--assemble_stats', action='store', dest='assemble_stats',  required=True, help='Required: SPAdes assembly output from simple_spades.sh')
    parser.add_argument('-q', '--quality_stats', action='store', dest='quality_stats',  required=True, help='Required: simple text file output from fastq_quality.py')
    args = parser.parse_args()
    print ("\nSET ARGUMENTS: ")
    print (args)
    print("")

    assembly_stats = args.assemble_stats
    quality_file = args.quality_stats
    sample_name = re.sub('_.*', '', os.path.basename(quality_file))

    merge_files(sample_name, assembly_stats, quality_file)
