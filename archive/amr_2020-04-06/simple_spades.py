#!/usr/bin/env python

import os
import sys
import argparse
import textwrap
from Bio import SeqIO

def trim_reads(sample_name, read1, read2):
    real_path = os.path.dirname(os.path.realpath(__file__))
    nextera_file = real_path + "/nextera.fa.gz" #file must be in same diretory as this script
    R1out = sample_name + "_R1_trimmed.fastq.gz"
    R2out = sample_name + "_R2_trimmed.fastq.gz"
    os.system(f'bbduk.sh -Xmx80g in1={read1} in2={read2} ref={nextera_file} ktrim=r k=23 mink=11 hdist=1 qtrim=lr trimq=5 minlen=36 out1={R1out} out2={R2out} stats=trim_stats.txt qchist=qc_by_base.txt threads=auto showspeed=f')
    return(R1out, R2out)

def simple_spades(sample_name, read1, read2):
    cwd = os.getcwd()
    kmer = "21,33,55,77,99,127"
    print(f'SPAdes assembling {sample_name}...')
    os.system(f'spades.py -k {kmer} --careful -1 {read1} -2 {read2} -o spades_output &> /dev/null')
    assembled_fasta = SeqIO.parse(f'{cwd}/spades_output/scaffolds.fasta', "fasta")
    return assembled_fasta

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

    ---------------------------------------------------------

    '''), epilog='''---------------------------------------------------------''')

    parser.add_argument('-n', '--name', action='store', dest='name', required=True, help='sample name')
    parser.add_argument('-r1', '--read1', action='store', dest='read1', required=True, help='R1 read')
    parser.add_argument('-r2', '--read2', action='store', dest='read2', required=True, help='R2 read')
    parser.add_argument('-t', '--trim', action='store_true', dest='trim', help='trim reads with bbduk')

    args = parser.parse_args()
    read1 = args.read1
    read2 = args.read2

    if args.trim:
        read1, read2 = trim_reads(args.name, read1, read2)

    assembled_fasta = simple_spades(args.name, read1, read2)
    SeqIO.write(assembled_fasta, "spades_scaffolds.fasta", "fasta")
