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

parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

    ---------------------------------------------------------

    Usage: merge_quality_files.py -r1 *_R1*gz -r2 *_R2*gz

    '''), epilog='''---------------------------------------------------------''')

parser.add_argument('-c', '--coverage', action='store', dest='coverage', type=int, required=True, help='Required: Minimum percent coverage set for ABRicate.')
parser.add_argument('-f', '--identity', action='store', dest='identity', type=int, required=True, help='Required: Minimum percent identity set for ABRicate.')
parser.add_argument('-i', '--isolate', action='store', dest='sample_name', required=True, help='Required: Isolate/sample name to appear on the report.')
args = parser.parse_args()
print ("\nSET ARGUMENTS: ")
print (args)
print("")


sample_name = args.sample_name

ts = time.time()
st = datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H-%M-%S')

per = " \%"
cov = '%s %s' % (args.coverage, str(per))
ident = '%s %s' % (args.identity, str(per))

####Abricate
ab_version = os.popen("abricate --version").readlines()[0]
ncbi = os.popen("abricate --list|sed '5!d'|tr '\t' '&'").readlines()[0]
res = os.popen("abricate --list|sed '7!d'|tr '\t' '&'").readlines()[0]


tex_document = sample_name + "-run_info_" + st + ".tex"

write_out = open(tex_document, 'w')

 #################################
            # LATEX

print(r'\documentclass[a4paper,12pt]{article}', file=write_out)
print(r'\usepackage[margin=0.5in]{geometry}', file=write_out)
print(r'\usepackage{graphicx}', file=write_out)
print(r'\usepackage[table]{xcolor}', file=write_out)
print(r'\usepackage{hyperref}', file=write_out)
print(r'\usepackage{tabularx}', file=write_out)
print(r'\usepackage{float}', file=write_out)
print(r'\usepackage{multirow}', file=write_out)
print(r'\usepackage{charter}', file=write_out)
print(r'\usepackage{mdwlist}', file=write_out)
print(r'\usepackage{fancyhdr}', file=write_out)
print(r'\usepackage[lastpage,user]{zref}', file=write_out)
print(r'\setlength{\headheight}{30pt}', file=write_out)
print(r'\setlength{\footskip}{10pt}', file=write_out)
print(r'\fancyhead[L]{\includegraphics[scale=0.15]{/home/shared/usdalogo.png}} \fancyhead[R]{\textbf{Isolate ID:}{%s}}' % (sample_name), file=write_out)
print(r'\lfoot{\today}', file=write_out)
print(r'\cfoot{\thepage\ of \zpageref{LastPage}}', file=write_out)
print(r'\pagestyle{fancy}', file=write_out)
print(r'\begin{document}', file=write_out)
print(r'\vspace{15mm}', file=write_out)
print(r'', file=write_out)
print(r'\centerline{\Huge{\bf{Antimicrobial Resistance Analysis}}}', file=write_out)
print(r'', file=write_out)
print(r'\vspace{10mm}', file=write_out)
print(r'\noindent\makebox[\linewidth]{\rule{14cm}{0.4pt}}', file=write_out)
print(r'\vspace{-3mm}', file=write_out)
print(r'\begin{center}', file=write_out)
print(r'\begin{minipage}{6in}', file=write_out)
print(r'Results were obtained using \href{https://github.com/tseemann/abricate}{ABRicate}.  ABRicate locates genes based on nucleotide BLAST searching a database, which is specified in the filename. Searches return results when the minimum percent coverage and percent identity are met.\\', file=write_out)
print(r'\end{minipage}', file=write_out)
print(r'\end{center}', file=write_out)
print(r'', file=write_out)
print(r'\vspace{-3mm}', file=write_out)
print(r'', file=write_out)
print(r'\noindent\makebox[\linewidth]{\rule{14cm}{0.4pt}}', file=write_out)
print(r'\vspace{5mm}', file=write_out)
print(r'', file=write_out)
print(r'\begin{center}', file=write_out)
print(r'{\Large Software Version}\\', file=write_out)
print(r'\vspace{3mm}', file=write_out)
print(r'\begin{tabular}{c}', file=write_out)
print(r'\hline', file=write_out)
print(r'%s \\' % ab_version, file=write_out)
print(r'\hline', file=write_out)
print(r'\end{tabular}', file=write_out)
print(r'\vspace{10mm}', file=write_out)
print(r'', file=write_out)
print(r'{\Large Database Versions} \\', file=write_out)
print(r'\vspace{5mm}', file=write_out)
print(r'', file=write_out)
print(r'\begin{tabular}{l|c|c}', file=write_out)
print(r'\hline', file=write_out)
print(r'\textbf{Database} & \textbf{Sequences} & \textbf{Date Updated} \\', file=write_out)
print(r'\hline', file=write_out)
print(r' %s \\' % ncbi, file=write_out)
print(r' %s \\' % res, file=write_out)
print(r'\hline', file=write_out)
print(r'\end{tabular}', file=write_out)
print(r'\end{center}', file=write_out)
print(r'', file=write_out)
print(r'\vspace{3mm}', file=write_out)
print(r'\begin{center}', file=write_out)
print(r'{\Large Paramaters}\\', file=write_out)
print(r'\vspace{5mm}', file=write_out)
print(r'\begin{tabular}{l|c}', file=write_out)
print(r'\hline', file=write_out)
print(r'Minimum coverage & %s \\ ' % (cov), file=write_out)
print(r'Minimum identity & %s \\' % (ident), file=write_out)
print(r'\hline', file=write_out)
print(r'\end{tabular}', file=write_out)
print(r'\vspace{5mm}', file=write_out)
print(r'\end{center}', file=write_out)
print(r'', file=write_out)
print(r'\vspace{5mm}', file=write_out)
print(r'\noindent', file=write_out)
print(r'\textbf{\Large Summary Output}\\', file=write_out)
print(r'\noindent\makebox[\linewidth]{\rule{19cm}{0.4pt}}', file=write_out)
print(r'This abbreviated section appears in the first tab of each Excel workbook allowing you to see basic information on all samples in the analysis.\\', file=write_out)
print(r'', file=write_out)
print(r'\begin{center}', file=write_out)
print(r'\begin{minipage}{6in}', file=write_out)
print(r'\begin{basedescript}{\desclabelstyle{\nextlinelabel}}', file=write_out)
print(r'\item[NUM\_FOUND] The number of distinct genes identified in the analysis of the sample.  This is NOT the total number of hits if gene duplicates are identified, it is the number of distinct genes identified.', file=write_out)
print(r'\item[Gene List and Percent Coverage] A gene list with percent coverage is also given by isolate.  If multiple identifications of the same gene were made in a single isolate then the percent coverage information is given in a colon separated list in the order the identifications are listed in the full analysis.', file=write_out)
print(r'\end{basedescript}', file=write_out)
print(r'\end{minipage}', file=write_out)
print(r'\end{center}', file=write_out)
print(r'', file=write_out)
print(r'\newpage', file=write_out)
print(r'', file=write_out)
print(r'\noindent', file=write_out)
print(r'\textbf{\Large Full Analysis Output}\\', file=write_out)
print(r'\noindent\makebox[\linewidth]{\rule{19cm}{0.4pt}}', file=write_out)
print(r'Details of the full analysis are located on subsequent tabs of the workbook.  Each tab represents a single sample.\\', file=write_out)
print(r'', file=write_out)
print(r'\begin{center}', file=write_out)
print(r'\begin{minipage}{6in}', file=write_out)
print(r'\begin{basedescript}{\desclabelstyle{\nextlinelabel}}', file=write_out)
print(r'\item[SEQUENCE] The analysis is ran on scaffolds(also known as nodes or contigs) output by the assembler. Node length and coverage are included with the name.', file=write_out)
print(r'\item[FILE] This field refers to the original file output by the analysis prior to the workbook creation.', file=write_out)
print(r'\item[START] This is the position in the scaffold where the alignment with the given reference gene starts.', file=write_out)
print(r'\item[END] This is the position of the scaffold where the alignment with the given reference gene ends.', file=write_out)
print(r'\item[GENE] The reference gene to which the alignment with the scaffold is performed.', file=write_out)
print(r'\item[COVERAGE] The positions (given as a range) of the reference gene that align with the scaffold divided by the length of the gene sequence. (aligned length/gene length)', file=write_out)
print(r'\item[COVERAGE\_MAP] This gives an overview of the alignment relative to the reference gene. ....-no alignment ====-alignment /-gap in alignment.', file=write_out)
print(r'\item[GAPS] The number of gaps in the alignment of the reference gene to the scaffold.', file=write_out)
print(r'\item[\%COVERAGE] The percent of the reference gene covered by the alignment with the scaffold.', file=write_out)
print(r'\item[\%IDENTITY] The percent identity of the reference gene with the scaffold.', file=write_out)
print(r'\item[DATABASE] The gene database that was used for comparison to the scaffolds.', file=write_out)
print(r'\item[ACCESSION] The accession number of the reference gene in the database used for the comparison.', file=write_out)
print(r'\item[PRODUCT] A description of the gene product provided from the selected database.', file=write_out)
print(r'\end{basedescript}', file=write_out)
print(r'\end{minipage}', file=write_out)
print(r'\end{center}', file=write_out)
print(r'', file=write_out)

print(r'\end{document}', file=write_out)

write_out.close()

os.system("pdflatex -interaction=nonstopmode {}" .format(tex_document))
os.system("pdflatex -interaction=nonstopmode {}" .format(tex_document))

