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

parser.add_argument('-a', '--assemble_stats', action='store', dest='assemble_stats',  required=True, help='Required: SPAdes assembly output from simple_spades.sh')
parser.add_argument('-q', '--quality_stats', action='store', dest='quality_stats',  required=True, help='Required: simple text file output from fastq_quality.py')
parser.add_argument('-s', '--seqsero', action='store', dest='seqsero',  required=False, help='Optional: Seqsero_result.txt text file output by SeqSero')
parser.add_argument('-m', '--mlst', action='store', dest='mlst', required=False, help='Optional: MLST text file output by mlst against pubMLST')
parser.add_argument('-l', '--lcontigs', action='store', dest='lcontigs', required=True, help='Required: Sum of the length of the contigs greater than 1000 nucleotides- used for calculatig percent of the genome on those contigs.')
parser.add_argument('-c', '--coverage', action='store', dest='coverage', type=int, required=True, help='Required: Minimum percent coverage set for ABRicate.')
parser.add_argument('-f', '--identity', action='store', dest='identity', type=int, required=True, help='Required: Minimum percent identity set for ABRicate.')
args = parser.parse_args()
print ("\nSET ARGUMENTS: ")
print (args)
print("")

assembly_stats = args.assemble_stats
quality_file = args.quality_stats
seqsero_file = args.seqsero
mlst_file = args.mlst
#mlstnames_file = os.path.dirname(os.path.abspath(mlst_names.txt))
lcontigs_file = args.lcontigs

sample_name = re.sub('_.*', '', quality_file)

ts = time.time()
st = datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H-%M-%S')

####Sequencing and Assembly Stats
merged_df = {}
quality_file_df = pd.read_csv(quality_file, sep='\t')
quality_file_df_T = quality_file_df.T
for idx, row in quality_file_df_T.iterrows():
    merged_df[idx] = row[0]

assembly_df = pd.read_csv(assembly_stats, sep='\t')
assembly_df_T = assembly_df.T
for idx, row in assembly_df_T.iterrows():
    merged_df[idx] = row[0]

filer1 = merged_df.get("File_R1", "n/a")
f1end=filer1.find('_R1_')
r1name=filer1.replace("_", "-")
r1name=r1name.replace("L001\_", "")
filer2 = merged_df.get("File_R2", "n/a")
f2end=filer2.find('_R2_')
r2name=filer2.replace("_", "-")
r2name=r2name.replace("L001", "")
sizer1 = merged_df.get("Size_R1", "n/a")
sizer2 = merged_df.get("Size_R2", "n/a")
quality1 = merged_df.get("mean_quality1", "n/a")
quality2 = merged_df.get("mean_quality2", "n/a")
q30 = merged_df.get("q30pass1", "n/a")
q30 = q30[0:4]
per = " \%"
q30p1 = '%s %s' % (q30, str(per))
q30 = merged_df.get("q30pass2", "n/a")
q30 = q30[0:4]
q30p2 = '%s %s' % (q30, str(per))
numscaff = merged_df.get("Number of scaffolds", "n/a")
sizescaff = "{:,}".format(merged_df.get("Total size of scaffolds", "n/a"))
sscaff = merged_df.get("Total size of scaffolds", "n/a")
longscaff = "{:,}".format(merged_df.get("Longest scaffold", "n/a"))
gt1k = merged_df.get("Number of scaffolds > 1K nt", "n/a")
n50len = merged_df.get("N50 scaffold length", "n/a")
l50cnt = merged_df.get("L50 scaffold count", "n/a")


cov = '%s %s' % (args.coverage, str(per))
ident = '%s %s' % (args.identity, str(per))

####SeqSero Data
if seqsero_file is not None:
    seqsero_indx = pd.read_csv(seqsero_file, sep='\t', skiprows=5, index_col=0, header=None)
    predant = seqsero_indx.iloc[0][1]
    predsubsp = seqsero_indx.iloc[1][1]
    predst = seqsero_indx.iloc[2][1]
    seqserocomment = seqsero_indx.tail(1).index.values[0]

####MLST Data
if mlst_file is not None:
    mlst_version = os.popen("mlst --version").readlines()[0]
    try:
        mlst_data = pd.read_csv(mlst_file, sep='\t', header=None)

        mlst_id = mlst_data.iloc[0,1]
        mlst_id = mlst_id.partition('_')

        if mlst_id[0] == '-':
            mlst_sp = "Not Identified"
        elif mlst_id[0] == 'senterica':
            mlst_sp = "Salmonella enterica"
        elif mlst_id[0] == "-":
            mlst_sp = "Not Identified"
            mlst_st = "Not Identified"
            mlst_sch = "Not Identified"
        elif mlst_id[0] == "ecoli":
            mlst_sp = "Escherichia coli"
        elif mlst_id[0] == 'yersinia':
            mlst_sp = 'Yersinia sp.'
        elif mlst_id[0] == "pmultocida":
            mlst_sp = "Pasturella multocida"
        elif mlst_id[0] == "bordetella":
            mlst_sp = "Bordetella spp."
        elif mlst_id[0] == "brachyspira":
            mlst_sp = "Brachyspira spp."
        elif mlst_id[0] == "brucella":
            mlst_sp = "Brucella spp."
        elif mlst_id[0] == "campylobacter":
            mlst_sp = "Campylobacter spp."
        elif mlst_id[0] == "cbotulinum":
            mlst_sp = "Clostridium botulinum"
        elif mlst_id[0] == "cdifficile":
            mlst_sp = "Clostridium difficile"
        elif mlst_id[0] == "cdiphtheriae":
            mlst_sp = "Corynebacterium diphtheriae"
        elif mlst_id[0] == "cfreundii":
            mlst_sp = "Citrobacter freundii"
        elif mlst_id[0] == "cronobacter":
            mlst_sp = "Cronobacter spp."
        elif mlst_id[0] == "csepticum":
            mlst_sp = "Clostridium septicum"
        elif mlst_id[0] == "cdifficile":
            mlst_sp = "Clostridium difficile"
        elif mlst_id[0] == "ecloacae":
            mlst_sp = "Enterobacter cloacae"
        elif mlst_id[0] == "edwardsiella":
            mlst_sp = "Edwardsiella spp."
        elif mlst_id[0] == "ecoli_2":
            mlst_sp = "Escherichia coli"
        elif mlst_id[0] == "hinfluenzae":
            mlst_sp = "Haemophilus influenzae"
        elif mlst_id[0] == "hparasuis":
            mlst_sp = "Haemophilus parasuis"
        elif mlst_id[0] == "hpylori":
            mlst_sp = "Helicobacter pylori"
        elif mlst_id[0] == "koxytoca":
            mlst_sp = "Klebsiella oxytoca"
        elif mlst_id[0] == "kpneumoniae":
            mlst_sp = "Klebsiella pneumoniae"
        elif mlst_id[0] == "leptospira":
            mlst_sp = "Leptospira spp."
        elif mlst_id[0] == "lmonocytogenes":
            mlst_sp = "Lysteria monocytogenes"
        elif mlst_id[0] == "mabscessus":
            mlst_sp = "Mycobacterium abscessus complex"
        elif mlst_id[0] == "magalactiae":
            mlst_sp = "Mycoplasma agalactiae"
        elif mlst_id[0] == "mbovis":
            mlst_sp = "Mycoplasma bovis"
        elif mlst_id[0] == "mhaemolytica":
            mlst_sp = "Mannheimia haemolytica"
        elif mlst_id[0] == "mhyopneumoniae":
            mlst_sp = "Mycoplasma hyopneumoniae"
        elif mlst_id[0] == "shaemolyticus":
            mlst_sp = "Staphylococcus haemolyticus"
        elif mlst_id[0] == "spneumoniae":
            mlst_sp = "Streptococcus pneumoniae"
        elif mlst_id[0] == "spseudintermedius":
            mlst_sp = "Staphylococcus pseudintermedius"
        elif mlst_id[0] == "taylorella":
            mlst_sp = "Taylorella spp."
        elif mlst_id[0] == "vcholerae":
            mlst_sp = "Vibrio cholerae"            
        elif mlst_id[0] == "vibrio":
            mlst_sp = "Vibrio spp."
        elif mlst_id[0] == "yersinia":
            mlst_sp = "Yersinia spp."
        elif mlst_id[0] == "ypseudotuberculosis":
            mlst_sp = "pseudotuberculosis"   
        else:
            mlst_id = mlst_id[0]
            mlst_genus = mlst_id[:1].upper()
            mlst_species = mlst_id[1:]
            mlst_sp = mlst_genus + ". " + mlst_species

        mlst_sch = mlst_id[2]
        mlst_st = mlst_data.iloc[0,2]
        mlst_g1 = mlst_data.iloc[0,3]
        mlst_g2 = mlst_data.iloc[0,4]
        mlst_g3 = mlst_data.iloc[0,5]
        mlst_g4 = mlst_data.iloc[0,6]
        mlst_g5 = mlst_data.iloc[0,7]
        mlst_g6 = mlst_data.iloc[0,8]
        mlst_g7 = mlst_data.iloc[0,9]
    except IndexError:
        mlst_sp = "Not Identified"
        mlst_st = "Not Identified"
 
   #mlst_version = os.popen("mlst --version").readlines()[0]
####Percent of Genome on contigs > 1k nt
lcontig_data = pd.read_csv(lcontigs_file, sep='\t', header=None)

gl = (lcontig_data[0][0]/sscaff)*100
rgl = round(gl, 2)
glc = '%s %s' % (str(rgl), str(per))

####Abricate
ab_version = os.popen("abricate --version").readlines()[0]
ncbi = os.popen("abricate --list|sed '5!d'|tr '\t' '&'").readlines()[0]
res = os.popen("abricate --list|sed '7!d'|tr '\t' '&'").readlines()[0]

##NCBI
ab_ncbi_file = glob.glob("*-ncbi-report.tab")


##Resfinder
ab_resfinder_file = glob.glob("*-resfinder-report.tab")


####AMRFinder
cwl_version = os.popen("/usr/local/bin/amr_finder/amrfinder -v| sed '1!d' |sed 's/:/\&/g'").readlines()[0]
dock_version = os.popen("/usr/local/bin/amr_finder/amrfinder -v| sed '3!d' |sed 's/container:/container\&/g'").readlines()[0]
amrdb_version = os.popen("/usr/local/bin/amr_finder/amrfinder -v| sed '4!d' | sed 's/: /\&/g'").readlines()[0]

amrfinder_file = glob.glob("*-amrfinder-report.tab")



tex_document = sample_name + ".tex"

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
print(r'\usepackage[lastpage,abspage,user]{zref}', file=write_out)
print(r'\usepackage{pdflscape}', file=write_out)
print(r'\usepackage{rotating}', file=write_out)
print(r'\setlength{\headheight}{30pt}', file=write_out)
print(r'\setlength{\footskip}{10pt}', file=write_out)
print(r'\fancyhead[L]{\includegraphics[scale=0.15]{/home/shared/usdalogo.png}} \fancyhead[R]{\textbf{Isolate ID:}{%s}}' % (sample_name), file=write_out)
print(r'\lfoot{\today}', file=write_out)
print(r'\cfoot{\thepage\ of \zref[abspage]{LastPage}}', file=write_out) 
print(r'\pagestyle{fancy}', file=write_out)
print(r'\thispagestyle{plain}', file=write_out)
print(r'\renewcommand{\thepage}{Page \arabic{page}}', file=write_out)
print(r'\includegraphics[scale=0.2]{/home/shared/usdalogo.png}', file=write_out)
print(r'\begin{document}', file=write_out)
print(r'', file=write_out)
print(r'\vspace{15mm}', file=write_out)
print(r'\centerline{\LARGE Whole Genome Sequencing Report}', file=write_out)
print(r'', file=write_out)
print(r'\vspace{10mm}', file=write_out)
print(r'', file=write_out)
print(r'{\large \today}', file=write_out)
print(r'', file=write_out)
print(r'\vspace{5mm}', file=write_out)
print(r'', file=write_out)
print(r'\textbf{Sample ID:} {\large %s}' % (sample_name), file=write_out) 
print(r'\vspace{5mm}', file=write_out)
print(r'', file=write_out)
print(r'\textbf{Sequencing Technology}', file=write_out)
print(r'\vspace{2mm}', file=write_out)
print(r'', file=write_out)
print(r'\begin{tabular}{ p{11cm} }', file=write_out)
print(r'Nextera XT DNA Library Preparation \\', file=write_out)
print(r'MiSeq 2 x 250 Read Generation \\', file=write_out)
print(r'\end{tabular}', file=write_out)
print(r'', file=write_out)
print(r'\vspace{5mm}', file=write_out)
print(r'', file=write_out)
print(r'\textbf{File Stats}', file=write_out)
print(r'\vspace{2mm}', file=write_out)
print(r'', file=write_out)
print(r'\begin{tabular}{ l | p{7cm} | p{7cm} }', file=write_out)
print(r'\hline', file=write_out)
print(r'Filename &  %s & %s \\' % (r1name, r2name), file=write_out )
print(r'\hline', file=write_out)
print(r'File size & %s & %s \\' % (sizer1, sizer2), file=write_out)
print(r'Mean Read Score & %s & %s \\' % (quality1, quality2), file=write_out)
print(r'Q30 Passing & %s & %s \\' % (q30p1, q30p2), file=write_out)
print(r'\hline', file=write_out)
print(r'\end{tabular}', file=write_out)
print(r'', file=write_out)
print(r'\vspace{5mm}', file=write_out)
print(r'', file=write_out)
print(r'\textbf{Assembly}', file=write_out)
print(r'\vspace{2mm}', file=write_out)
print(r'', file=write_out)
print(r'\begin{tabularx}{\columnwidth}{X|X|X|X|X|X|X}', file=write_out)
print(r'\hline', file=write_out)
print(r'Scaffolds & Total length & Longest scaffold & Scaffolds \textgreater 1K nt & Genome \textgreater 1K nt & N50 & L50 \\', file=write_out)
print(r'\hline', file=write_out)
print(r'%s & %s & %s & %s & %s & %s & %s \\' % (numscaff, sizescaff, longscaff, gt1k, glc, n50len, l50cnt), file=write_out)
print(r'\hline', file=write_out)
print(r'\end{tabularx}', file=write_out)
print(r'', file=write_out)
print(r'{\footnotesize\ De novo assembly performed using \href{http://cab.spbu.ru/software/spades/} {SPAdes.} }\\', file=write_out)

if mlst_file is not None:
    print(r'\vspace{5mm}', file=write_out)
    print(r'', file=write_out)
    print(r'\textbf{Multi Locus Sequence Typing (MLST)}', file=write_out)
    print(r'\vspace{2mm}', file=write_out)
    print(r'', file=write_out)
    print(r'\begin{tabular}{ p{8cm} | p{8cm} }', file=write_out)
    print(r'\hline', file=write_out)

    if mlst_sch == '':
        print(r'Organism ID: \textbf{\large %s} & Sequence type: \textbf{\large %s} \\' % (mlst_sp, mlst_st), file=write_out)
    else:
        print(r'Organism ID: \textbf{\large %s} & Schema \& Sequence type: \textbf{\large %s-%s} \\' % (mlst_sp, mlst_sch, mlst_st), file=write_out)

    if mlst_st == 'Not Identified':
        print(r'\hline', file=write_out)
        print(r'\end{tabular}', file=write_out)
        print(r'', file=write_out)
        print(r'\begin{tabular}{ p{18cm} }', file=write_out)
        print(r'{\footnotesize\ Data obatined using %s.  Software website: \href{https://github.com/tseemann/mlst} {https://github.com/tseemann/mlst}}\\' % mlst_version, file=write_out)
        print(r'{\footnotesize\ Information on MLST schemes and allelic profiles can be found at \href{https://pubmlst.org} {pubMLST.org} }\\', file=write_out) 
        print(r'\end{tabular}', file=write_out)
        print(r'', file=write_out)
    else:
        print(r'\hline', file=write_out)
        print(r'\end{tabular}', file=write_out)
        print(r'', file=write_out)
        print(r'\begin{tabular}{ p{2cm} | p{2cm} | p{2cm} | p{2cm} | p{2cm} | p{2cm} | p{2cm} }', file=write_out)
        print(r'%s & %s & %s & %s & %s & %s & %s \\' % (mlst_g1, mlst_g2, mlst_g3, mlst_g4, mlst_g5, mlst_g6, mlst_g7), file=write_out)
        print(r'\hline', file=write_out)
        print(r'\end{tabular}', file=write_out)
        print(r'', file=write_out)
        print(r'\begin{tabular}{ p{18cm} }', file=write_out)
        print(r'{\footnotesize\ Data obatined using %s.  Software website: \href{https://github.com/tseemann/mlst} {https://github.com/tseemann/mlst}}\\' % mlst_version, file=write_out)
        print(r'{\footnotesize\ Information on MLST schemes and allelic profiles can be found at \href{https://pubmlst.org} {pubMLST.org} }\\', file=write_out) 
        print(r'\end{tabular}', file=write_out)
        print(r'', file=write_out)

        if mlst_id[0] == "senterica":  
            print(r'\vspace{7mm}', file=write_out)
            print(r'\textbf{Serotyping for Salmonella Isolates}', file=write_out)
            print(r'\vspace{2mm}', file=write_out)
            print(r'', file=write_out)
            print(r'\begin{tabular}{ p{5cm} | p{12cm} }', file=write_out)
            print(r'\hline', file=write_out)
            print(r'Predicted serotype(s) & \textbf{\large %s} \\' % predst, file=write_out)
            print(r'Predicted antigenic profile & %s \\' % predant, file=write_out)             
            print(r'Predicted subspecies & %s \\' % predsubsp, file=write_out)
            print(r'\hline', file=write_out)
            print(r'\end{tabular}', file=write_out)
            print(r'', file=write_out)
            print(r'\begin{tabular}{ p{18cm} }', file=write_out)
            print(r'Comments: %s \\' % seqserocomment, file=write_out)
            print(r'\hline', file=write_out)
            print(r'{\footnotesize\ Data obtained using SeqSero.  Software website: \href{https://github.com/denglab/SeqSero2} {https://github.com/denglab/SeqSero2} }\\', file=write_out)
            print(r'\end{tabular}', file=write_out)
            print(r'', file=write_out)

if mlst_file is None and seqsero_file is not None:
    print(r'\vspace{7mm}', file=write_out)
    print(r'\textbf{Serotyping for Salmonella Isolates}', file=write_out)
    print(r'\vspace{2mm}', file=write_out)
    print(r'', file=write_out)
    print(r'\begin{tabular}{ p{5cm} | p{12cm} }', file=write_out)
    print(r'\hline', file=write_out)
    print(r'Predicted serotype(s) & \textbf{\large %s} \\' % predst, file=write_out)
    print(r'Predicted antigenic profile & %s \\' % predant, file=write_out)
    print(r'Predicted subspecies & %s \\' % predsubsp, file=write_out)
    print(r'\hline', file=write_out)
    print(r'\end{tabular}', file=write_out)
    print(r'', file=write_out)
    print(r'\begin{tabular}{ p{18cm} }', file=write_out)
    print(r'Comments: %s \\' % seqserocomment, file=write_out)
    print(r'\hline', file=write_out)
    print(r'{\footnotesize\ Data obtained using SeqSero.  Software website: \href{https://github.com/denglab/SeqSero2} {https://github.com/denglab/SeqSero2} }\\', file=write_out)
    print(r'\end{tabular}', file=write_out)
    print(r'', file=write_out)

print(r'', file=write_out)

#AMR START

print(r'\newpage', file=write_out)
print(r'', file=write_out)
print(r'\vspace{15mm}', file=write_out)
print(r'', file=write_out)
print(r'\centerline{\Huge{\bf{Antimicrobial Resistance Analysis}}}', file=write_out)
print(r'', file=write_out)
print(r'\vspace{5mm}', file=write_out)
print(r'\noindent\makebox[\linewidth]{\rule{14cm}{0.4pt}}', file=write_out)
print(r'\vspace{-3mm}', file=write_out)
print(r'\begin{center}', file=write_out)
print(r'\begin{minipage}{6in}', file=write_out)
print(r'Results were obtained using \href{https://github.com/tseemann/abricate}{ABRicate} and \href{https://github.com/ncbi/pipelines/tree/master/amr\_finder}{AMRFinder}.  ABRicate locates genes based on nucleotide BLAST searching a database, which is specified in the filename. Searches return results when the minimum percent coverage and percent identity are met.  AMRFinder uses BLASTX to search a hierarchy of gene families with predetermined cutoffs.\\', file=write_out)
print(r'\end{minipage}', file=write_out)
print(r'\end{center}', file=write_out)
print(r'', file=write_out)
print(r'\vspace{-3mm}', file=write_out)
print(r'', file=write_out)
print(r'\noindent\makebox[\linewidth]{\rule{14cm}{0.4pt}}', file=write_out)
print(r'\vspace{3mm}', file=write_out)
print(r'', file=write_out)
print(r'\begin{center}', file=write_out)
print(r'{\LARGE ABRicate}\\', file=write_out)
print(r'\vspace{8mm}', file=write_out)
print(r'{\Large Software Version}\\', file=write_out)
print(r'\vspace{3mm}', file=write_out)
print(r'\begin{tabular}{c}', file=write_out)
print(r'\hline', file=write_out)
print(r'%s \\' % ab_version, file=write_out)
print(r'\hline', file=write_out)
print(r'\end{tabular}', file=write_out)
print(r'\vspace{8mm}', file=write_out)
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
print(r'\vspace{8mm}', file=write_out)
print(r'\begin{center}', file=write_out)
print(r'{\Large AMR Gene Parameters}\\', file=write_out)
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

#ABRICATE_DATA_OUTPUT

print(r'\newpage', file=write_out)
print(r'', file=write_out)
print(r'\begin{landscape}', file=write_out)
print(r'\center{\textbf{\Large Abricate - ResFinder}}', file=write_out)
print(r'', file=write_out)
print(r'\vspace{5mm}', file=write_out)
print(r'\resizebox{27cm}{!}{', file=write_out)
print(r'\begin{tabular}{ |l|c|c|c|c|c|c|c|c|c|l| }', file=write_out)
print(r'\hline', file=write_out)
print(r'SEQUENCE&START&END&GENE&COVERAGE&GAPS&\%COVERAGE&\%IDENTITY&ACCESSION&PRODUCT\\', file=write_out)
print(r'\hline', file=write_out)
if os.stat(ab_resfinder_file[0]).st_size == 0:
   print(r'No Antimicrobial Resistance Genes found.\\', file=write_out)
else:
   with open(ab_resfinder_file[0], "r") as f:
      for line in f:
         line = line.strip()
         line = line + "\\\\"
         print(r'%s' % line, file=write_out)
print(r'\hline', file=write_out)
print(r'\end{tabular}}', file=write_out)
print(r'', file=write_out)
print(r'\vspace{10mm}', file=write_out)
print(r'', file=write_out)
print(r'\center{\textbf{\Large Abricate - NCBI}}', file=write_out)
print(r'', file=write_out)
print(r'\vspace{5mm}', file=write_out)
print(r'', file=write_out)
print(r'\resizebox{27cm}{!}{', file=write_out)
print(r'\begin{tabular}{ |l|c|c|c|c|c|c|c|c|c|l| }', file=write_out)
print(r'\hline', file=write_out)
print(r'SEQUENCE&START&END&GENE&COVERAGE&GAPS&\%COVERAGE&\%IDENTITY&ACCESSION&PRODUCT\\', file=write_out)
print(r'\hline', file=write_out)
if os.stat(ab_ncbi_file[0]).st_size ==0:
   print(r'No Antimicrobial Resistance Genes found.\\', file=write_out)
else:
   with open(ab_ncbi_file[0], "r") as f:
      for line in f:
         line = line.strip()
         line = line + "\\\\"
         print(r'%s' % line, file=write_out)
print(r'\hline', file=write_out)
print(r'\end{tabular}}', file=write_out)
print(r'', file=write_out)
print(r'\vspace{5mm}', file=write_out)
print(r'', file=write_out)
print(r'\end{landscape}', file=write_out)


##AMRFinder

print(r'\newpage', file=write_out)
print(r'', file=write_out)
print(r'\vspace{5mm}', file=write_out)
print(r'', file=write_out)
print(r'\begin{center}', file=write_out)
print(r'{\LARGE AMRFinder}\\', file=write_out)
print(r'\vspace{8mm}', file=write_out)
print(r'', file=write_out)
print(r'{\Large Software Versions}\\', file=write_out)
print(r'\vspace{3mm}', file=write_out)
print(r'\begin{tabular}{l|c}', file=write_out)
print(r'\hline', file=write_out)
print(r'%s \\' % cwl_version, file=write_out)
print(r'%s \\' % dock_version, file=write_out)
print(r'\hline', file=write_out)
print(r'\end{tabular}', file=write_out)
print(r'\vspace{8mm}', file=write_out)
print(r'', file=write_out)
print(r'{\Large Database Version}\\', file=write_out)
print(r'\vspace{3mm}', file=write_out)
print(r'', file=write_out)
print(r'\begin{tabular}{l|c}', file=write_out)
print(r'\hline', file=write_out)
print(r'{} \\' .format(amrdb_version), file=write_out)
print(r'\hline', file=write_out)
print(r'\end{tabular}', file=write_out)
print(r'\end{center}', file=write_out)
print(r'', file=write_out)
print(r'\vspace{5mm}', file=write_out)
print(r'\noindent', file=write_out)
print(r'\textbf{\Large Output}\\', file=write_out)
print(r'\noindent\makebox[\linewidth]{\rule{19cm}{0.4pt}}', file=write_out)
print(r'Definitions are from \href{https://github.com/ncbi/pipelines/tree/master/amr\_finder}{https://github.com/ncbi/pipelines/tree/master/amr\_finder}.\\', file=write_out)
print(r'', file=write_out)
print(r'\begin{center}', file=write_out)
print(r'\begin{minipage}{6in}', file=write_out)
print(r'\begin{basedescript}{\desclabelstyle{\nextlinelabel}}', file=write_out)
print(r'\item[Target Identifier] This is from the FASTA defline for the DNA sequence', file=write_out)
print(r'\item[Contig id] Contig name', file=write_out)
print(r'\item[Start] 1-based coordinate of first nucleotide coding from protein in DNA sequence on contig', file=write_out)
print(r'\item[Stop] 1-based coordinate of last nucleotide coding for protein in DNA sequence on contig', file=write_out)
print(r'\item[Gene symbol] Gene or gene-family symbol for protein hit', file=write_out)
print(r'\item[Protein name] Full-text name for the protein', file=write_out)
print(r'\item[Method] Type of hit found by AMRFinder one of five options', file=write_out)
print(r'\indent', file=write_out)
print(r'\item[ALLELE] 100\% sequence match over 100\% of length to a protein named at the allele level in the AMRFinder database', file=write_out)
print(r'\item[EXACT] 100\% sequence match over 100\% of length to a protein in the database that is not a named allele', file=write_out)
print(r'\item[BLAST] BLAST alignment is \textgreater 90\% of length and \textgreater 90\% identity to a protein in a the AMRFinder database', file=write_out)
print(r'\item[PARTIAL] BLAST alignment is \textgreater 50\% of length, but \textless 90\% of length and \textgreater 90\% identity', file=write_out)
print(r'\end{basedescript}', file=write_out)
print(r'\end{minipage}', file=write_out)
print(r'\end{center}', file=write_out)
print(r'', file=write_out)
print(r'\newpage', file=write_out)
print(r'\begin{center}', file=write_out)
print(r'\begin{minipage}{6in}', file=write_out)
print(r'\begin{basedescript}{\desclabelstyle{\nextlinelabel}}', file=write_out)
print(r'\item[HMM] HMM was hit above the cutoff, but there was not a BLAST hit that met standards for BLAST or PARTIAL', file=write_out)
print(r'\noident', file=write_out)
print(r'\item[Target length] The length of the query protein. The length of the BLAST hit for translated-DNA searches', file=write_out)
print(r'\item[Reference protein length] The length of the AMR Protein in the database (NA if HMM-only hit)', file=write_out)
print(r'\item[\% Coverage of reference protein] \% covered by blast hit (NA if HMM-only hit', file=write_out)
print(r'\item[\% Identity to reference protein] \% amino-acid identity to reference protein (NA if HMM-only hit)', file=write_out)
print(r'\item[Alignment length] Length of BLAST alignment in amino acids (NA if HMM-only hit)', file=write_out)
print(r'\item[Accession of closest protein] RefSeq accession for protein hit by BLAST (NA if HMM-only hit)', file=write_out)
print(r'\item[Name of closest protein] Full name assigned to the AMRFinder database protein (NA if HMM-only hit)', file=write_out)
print(r'\item[HMM id] Accession for the HMM', file=write_out)
print(r'\item[HMM description] The family name associated with the HMM', file=write_out)
print(r'\end{basedescript}', file=write_out)
print(r'\end{minipage}', file=write_out)
print(r'\end{center}', file=write_out)
print(r'', file=write_out)

#AMRFinder Results Page
print(r'\newpage', file=write_out)
print(r'', file=write_out)
print(r'\begin{landscape}', file=write_out)
print(r'\center{\textbf{\Large AMRFinder}', file=write_out)
print(r'', file=write_out)
print(r'\vspace{30mm}', file=write_out)
if os.stat(amrfinder_file[0]).st_size == 0:
   print(r'\resizebox{15cm}{!}{', file=write_out)
   print(r'\begin{tabular}{c{15cm}}', file=write_out)
   print(r'\hline', file=write_out)
   print(r'{\large No Antimicrobial Resistance Genes found.\\}', file=write_out)
   print(r'\hline', file=write_out)
   print(r'\end{tabular}}', file=write_out)
else:
   print(r'\resizebox{27cm}{!}{', file=write_out)
   print(r'\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|c|c|}', file=write_out)
   print(r'Contig id&Start&Stop&Gene symbol&Protein name&\begin{rotate}{60}Method\end{rotate}&\begin{rotate}{60}Target length\end{rotate}&\begin{rotate}{60}Reference protein length\end{rotate}& \begin{rotate}{60}\% Coverage of reference protein\end{rotate}&\begin{rotate}{60}\% Identity of reference protein\end{rotate}& \begin{rotate}{60}Alignment length \end{rotate}& \begin{rotate}{60}Accession of closest protein \end{rotate}& \begin{rotate}{60}HMM id\end{rotate}  \\', file=write_out)
   print(r'\hline', file=write_out)
   with open(amrfinder_file[0], "r") as f:
      for line in f:
         line = line.strip()
         line = line + "\\\\"
         print(r'%s' % line, file=write_out)
   print(r'\hline', file=write_out)

   print(r'\end{tabular}}', file=write_out)
print(r'\end{document}', file=write_out)

write_out.close()

os.system("pdflatex -interaction=nonstopmode {}" .format(tex_document))
os.system("pdflatex -interaction=nonstopmode {}" .format(tex_document))

