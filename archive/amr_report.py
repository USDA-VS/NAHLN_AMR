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
import numpy as np
import math
import svgwrite

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
parser.add_argument('-n', '--genome_length', action='store', dest='genlength', type=int, required=False, help='Optional: Genome length if MLST will not identifiy the species.')
args = parser.parse_args()
print ("\nSET ARGUMENTS: ")
print (args)
print("")

assembly_stats = args.assemble_stats
quality_file = args.quality_stats
seqsero_file = args.seqsero
mlst_file = args.mlst
lcontigs_file = args.lcontigs

sample_name = re.sub('_.*', '', quality_file)

ts = time.time()
st = datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H-%M-%S')

####Lookup- MLST Name and Approx Genome Size 

lookup = {"abaumannii": {"Scientific Name": "Acinetobacter baumannii", "Approx Genome Size (Mb)": 4.3}, "abaumannii_2 ": {"Scientific Name": "Acinetobacter baumannii", "Approx Genome Size (Mb)": 4.3}, "achromobacter": {"Scientific Name": "Achromobacter spp.", "Approx Genome Size (Mb)": 6.5}, "aeromonas": {"Scientific Name": "Aeromonas spp.", "Approx Genome Size (Mb)": 5.0}, "afumigatus": {"Scientific Name": "Aspergillus fumigatus", "Approx Genome Size (Mb)": 30.0}, "aphagocytophilum": {"Scientific Name": "Anaplasma phagocytophilum", "Approx Genome Size (Mb)": 1.5}, "arcobacter": {"Scientific Name": "Arcobacter spp.", "Approx Genome Size (Mb)": 2.5}, "bbacilliformis": {"Scientific Name": "Bartonella bacilliformis", "Approx Genome Size (Mb)": 1.5}, "bcc": {"Scientific Name": "Burkholderia cepacia", "Approx Genome Size (Mb)": 8.6}, "bcereus": {"Scientific Name": "Bacillus cereus", "Approx Genome Size (Mb)": 5.4}, "bhampsonii": {"Scientific Name": "Brachyspira hampsonii", "Approx Genome Size (Mb)": 3.0}, "bhenselae": {"Scientific Name": "Bartonella henselae", "Approx Genome Size (Mb)": 2.0}, "bhyodysenteriae": {"Scientific Name": "Brachyspira hyodysenteriae", "Approx Genome Size (Mb)": 3.0}, "bintermedia": {"Scientific Name": "Brachyspira intermedia", "Approx Genome Size (Mb)": 3.0}, "blicheniformis": {"Scientific Name": "Bacillus licheniformis", "Approx Genome Size (Mb)": 4.2}, "bordetella": {"Scientific Name": "Bordetella spp.", "Approx Genome Size (Mb)": 5.0}, "borrelia": {"Scientific Name": "Borrelia spp.", "Approx Genome Size (Mb)": 1.0}, "bpilosicoli": {"Scientific Name": "Brachyspira pilosicoli", "Approx Genome Size (Mb)": 2.5}, "bpseudomallei": {"Scientific Name": "Burkholderia pseudomallei", "Approx Genome Size (Mb)": 7.0}, "brachyspira": {"Scientific Name": "Brachyspira spp.", "Approx Genome Size (Mb)": 3.0}, "brucella": {"Scientific Name": "Brucella spp.", "Approx Genome Size (Mb)": 3.0}, "bsubtilis": {"Scientific Name": "Bacillus subtilis", "Approx Genome Size (Mb)": 4.2}, "calbicans": {"Scientific Name": "Candida albicans", "Approx Genome Size (Mb)": 28.6}, "campylobacter": {"Scientific Name": "Campylobacter spp.", "Approx Genome Size (Mb)": 2.0}, "cbotulinum": {"Scientific Name": "Clostridium botulinum", "Approx Genome Size (Mb)": 4.0}, "cconcisus": {"Scientific Name": "Campylobacter concisus", "Approx Genome Size (Mb)": 2.0}, "cdifficile": {"Scientific Name": "Clostridium difficile", "Approx Genome Size (Mb)": 4.0}, "cdiphtheriae": {"Scientific Name": "Corynebacterium diphtheriae", "Approx Genome Size (Mb)": 2.5}, "cfetus": {"Scientific Name": "Campylobacter fetus", "Approx Genome Size (Mb)": 2.0}, "cfreundii": {"Scientific Name": "Citrobacter freundii", "Approx Genome Size (Mb)": 5.4}, "cglabrata": {"Scientific Name": "Candida glabrata", "Approx Genome Size (Mb)": 12.3}, "chelveticus": {"Scientific Name": "Campylobacter helveticus", "Approx Genome Size (Mb)": 2.0}, "chlamydiales": {"Scientific Name": "Chlamydiales spp.", "Approx Genome Size (Mb)": 1.0}, "chyointestinalis": {"Scientific Name": "Campylobacter hyointestinalis", "Approx Genome Size (Mb)": 2.0}, "cinsulaenigrae": {"Scientific Name": "Campylobacter insulaenigrae", "Approx Genome Size (Mb)": 2.0}, "ckrusei": {"Scientific Name": "Candida krusei", "Approx Genome Size (Mb)": 11.7}, "clanienae": {"Scientific Name": "Campylobacter lanienae", "Approx Genome Size (Mb)": 2.0}, "clari": {"Scientific Name": "Campylobacter lari", "Approx Genome Size (Mb)": 2.0}, "cmaltaromaticum": {"Scientific Name": "Carnobacterium maltaromaticum", "Approx Genome Size (Mb)": 3.7}, "cronobacter": {"Scientific Name": "Cronobacter spp.", "Approx Genome Size (Mb)": 4.5}, "csepticum": {"Scientific Name": "Clostrodium septicum", "Approx Genome Size (Mb)": 4.0}, "csinensis": {"Scientific Name": "Clonorchis sinensis", "Approx Genome Size (Mb)": 550.0}, "csputorum": {"Scientific Name": "Campylobacter sputorum", "Approx Genome Size (Mb)": 2.0}, "ctropicalis": {"Scientific Name": "Candida tropicalis", "Approx Genome Size (Mb)": 15.0}, "cupsaliensis": {"Scientific Name": "Campylobacter upsaliensis", "Approx Genome Size (Mb)": 2.0}, "dnodosus": {"Scientific Name": "Dichelobacter nodosus", "Approx Genome Size (Mb)": 1.4}, "ecloacae": {"Scientific Name": "Enterobacter cloacae", "Approx Genome Size (Mb)": 5.6}, "ecoli": {"Scientific Name": "Escherichia coli", "Approx Genome Size (Mb)": 5.0}, "ecoli_2": {"Scientific Name": "Escherichia coli", "Approx Genome Size (Mb)": 5.0}, "edwardsiella": {"Scientific Name": "Edwardsiella spp.", "Approx Genome Size (Mb)": 4.0}, "efaecalis": {"Scientific Name": "Enterococcus faecalis", "Approx Genome Size (Mb)": 3.5}, "efaecium": {"Scientific Name": "Enterococcus faecium", "Approx Genome Size (Mb)": 3.0}, "fpsychrophilum ": {"Scientific Name": "Flavobacterium psychrophilum", "Approx Genome Size (Mb)": 3.0}, "ganatis": {"Scientific Name": "Gallibacterium anatis", "Approx Genome Size (Mb)": 3.0}, "hcinaedi": {"Scientific Name": "Helicobacter cinaedi", "Approx Genome Size (Mb)": 2.0}, "hinfluenzae": {"Scientific Name": "Haemophilus influenzae", "Approx Genome Size (Mb)": 2.0}, "hparasuis": {"Scientific Name": "Haemophilus parasuis", "Approx Genome Size (Mb)": 2.3}, "hpylori": {"Scientific Name": "Helicobacter pylori", "Approx Genome Size (Mb)": 2.0}, "hsuis": {"Scientific Name": "Helicobacter suis", "Approx Genome Size (Mb)": 1.5}, "kaerogenes": {"Scientific Name": "Klebsiella aerogenes", "Approx Genome Size (Mb)": 5.3}, "kkingae": {"Scientific Name": "Kingella kingae", "Approx Genome Size (Mb)": 2.0}, "koxytoca": {"Scientific Name": "Klebsiella oxytoca", "Approx Genome Size (Mb)": 7.0}, "kpneumoniae": {"Scientific Name": "Klebsiella pneumoniae", "Approx Genome Size (Mb)": 6.0}, "kseptempunctata": {"Scientific Name": "Kudoa septempunctata", "Approx Genome Size (Mb)": 25.0}, "leptospira": {"Scientific Name": "Leptospira spp.", "Approx Genome Size (Mb)": 4.0}, "leptospira_2": {"Scientific Name": "Leptospira spp.", "Approx Genome Size (Mb)": 4.0}, "leptospira_3": {"Scientific Name": "Leptospira spp.", "Approx Genome Size (Mb)": 4.0}, "lmonocytogenes": {"Scientific Name": "Lysteria monocytogenes", "Approx Genome Size (Mb)": 3.0}, "lsalivarius": {"Scientific Name": "Lactobacillus salivarius", "Approx Genome Size (Mb)": 2.0}, "mabscessus": {"Scientific Name": "Mycobacterium abscessus complex", "Approx Genome Size (Mb)": 5.0}, "magalactiae": {"Scientific Name": "Mycoplasma agalactiae", "Approx Genome Size (Mb)": 1.0}, "mbovis": {"Scientific Name": "Mycoplasma bovis", "Approx Genome Size (Mb)": 1.0}, "mcanis": {"Scientific Name": "Macrococcus canis", "Approx Genome Size (Mb)": 2.0}, "mcaseolyticus": {"Scientific Name": "Macrococcus caseolyticus", "Approx Genome Size (Mb)": 2.0}, "mcatarrhalis": {"Scientific Name": "Moraxella catarrhalis", "Approx Genome Size (Mb)": 2.0}, "mhaemolytica": {"Scientific Name": "Mannheimia haemolytica", "Approx Genome Size (Mb)": 3.0}, "mhyopneumoniae ": {"Scientific Name": "Mycoplasma hyopneumoniae", "Approx Genome Size (Mb)": 1.0}, "mhyorhinis ": {"Scientific Name": "Mycoplasma hyorhinis", "Approx Genome Size (Mb)": 1.0}, "miowae": {"Scientific Name": "Mycoplasma iowae", "Approx Genome Size (Mb)": 1.0}, "mmassiliense": {"Scientific Name": "Mycobacterium massiliense", "Approx Genome Size (Mb)": 5.0}, "mplutonius": {"Scientific Name": "Melissococcus plutonius", "Approx Genome Size (Mb)": 2.0}, "mpneumoniae": {"Scientific Name": "Mycoplasma pneumoniae", "Approx Genome Size (Mb)": 1.0}, "msynoviae": {"Scientific Name": "Mycoplasma synoviae", "Approx Genome Size (Mb)": 1.0}, "neisseria": {"Scientific Name": "Neisseria spp.", "Approx Genome Size (Mb)": 2.0}, "orhinotracheale": {"Scientific Name": "Ornithobacterium rhinotracheale", "Approx Genome Size (Mb)": 2.5}, "otsutsugamushi": {"Scientific Name": "Orientia tsutsugamushi", "Approx Genome Size (Mb)": 2.0}, "pacnes": {"Scientific Name": "Propionibacterium acnes", "Approx Genome Size (Mb)": 2.5}, "paeruginosa": {"Scientific Name": "Pseudomonas aeruginosa", "Approx Genome Size (Mb)": 6.0}, "pdamselae": {"Scientific Name": "Photobacterium damselae", "Approx Genome Size (Mb)": 4.5}, "pfluorescens": {"Scientific Name": "Pseudomonas fluorescens", "Approx Genome Size (Mb)": 6.5}, "pgingivalis ": {"Scientific Name": "Porphyromonas gingivalis", "Approx Genome Size (Mb)": 2.5}, "plarvae": {"Scientific Name": "Paenibacillus larvae", "Approx Genome Size (Mb)": 4.0}, "pmultocida_multihost": {"Scientific Name": "Pasteurella multocida", "Approx Genome Size (Mb)": 2.0}, "pmultocida_rirdc": {"Scientific Name": "Pasteurella multocida", "Approx Genome Size (Mb)": 2.0}, "ppentosaceus": {"Scientific Name": "Pediococcus pentosaceus", "Approx Genome Size (Mb)": 2.0}, "ranatipestifer": {"Scientific Name": "Riemerella anatipestifer", "Approx Genome Size (Mb)": 2.0}, "rhodococcus": {"Scientific Name": "Rhodococcus equi", "Approx Genome Size (Mb)": 5.0}, "sagalactiae": {"Scientific Name": "Streptococcus agalactiae", "Approx Genome Size (Mb)": 2.0}, "saureus": {"Scientific Name": "Staphylococcus aureus", "Approx Genome Size (Mb)": 2.5}, "sbsec": {"Scientific Name": "Streptococcus bovis/equinus complex", "Approx Genome Size (Mb)": 2.0}, "scanis": {"Scientific Name": "Streptococcus canis", "Approx Genome Size (Mb)": 2.0}, "sdysgalactiae": {"Scientific Name": "Streptococcus dysgalactiae", "Approx Genome Size (Mb)": 2.0}, "senterica": {"Scientific Name": "Salmonella enterica", "Approx Genome Size (Mb)": 5.0}, "sepidermidis ": {"Scientific Name": "Staphylococcus epidermidis", "Approx Genome Size (Mb)": 2.5}, "sgallolyticus": {"Scientific Name": "Streptococcus gallolyticus", "Approx Genome Size (Mb)": 2.0}, "shaemolyticus ": {"Scientific Name": "Staphylococcus haemolyticus", "Approx Genome Size (Mb)": 2.5}, "shominis": {"Scientific Name": "Staphylococcus hominis", "Approx Genome Size (Mb)": 2.5}, "sinorhizobium": {"Scientific Name": "Sinorhizobium spp.", "Approx Genome Size (Mb)": 7.0}, "slugdunensis": {"Scientific Name": "Staphylococcus lugdunensis", "Approx Genome Size (Mb)": 2.5}, "smaltophilia": {"Scientific Name": "Stenotrophomonas maltophilia", "Approx Genome Size (Mb)": 5.0}, "soralis": {"Scientific Name": "Streptococcus oralis", "Approx Genome Size (Mb)": 2.0}, "spneumoniae": {"Scientific Name": "Streptococcus pneumoniae", "Approx Genome Size (Mb)": 2.0}, "spseudintermedius ": {"Scientific Name": "Staphylococcus pseudintermedius", "Approx Genome Size (Mb)": 2.5}, "spyogenes": {"Scientific Name": "Streptococcus pyogenes", "Approx Genome Size (Mb)": 2.0}, "ssuis": {"Scientific Name": "Streptococcus suis", "Approx Genome Size (Mb)": 2.0}, "sthermophilus": {"Scientific Name": "Streptococcus thermophilus", "Approx Genome Size (Mb)": 2.0}, "sthermophilus_2": {"Scientific Name": "Streptococcus thermophilus", "Approx Genome Size (Mb)": 2.0}, "streptomyces": {"Scientific Name": "Streptomyces spp.", "Approx Genome Size (Mb)": 9.0}, "suberis": {"Scientific Name": "Streptococcus uberis", "Approx Genome Size (Mb)": 2.0}, "szooepidemicus": {"Scientific Name": "Streptococcus zooepidemicus", "Approx Genome Size (Mb)": 2.0}, "taylorella": {"Scientific Name": "Taylorella spp.", "Approx Genome Size (Mb)": 2.0}, "tenacibaculum ": {"Scientific Name": "Tenacibaculum spp.", "Approx Genome Size (Mb)": 3.0}, "tvaginalis": {"Scientific Name": "Trichomonas vaginalis", "Approx Genome Size (Mb)": 170.0}, "vcholerae": {"Scientific Name": "Vibrio cholerae", "Approx Genome Size (Mb)": 4.0}, "vibrio": {"Scientific Name": "Vibrio spp.", "Approx Genome Size (Mb)": 4.0}, "vparahaemolyticus": {"Scientific Name": "Vibrio parahaemolyticus", "Approx Genome Size (Mb)": 4.0}, "vtapetis": {"Scientific Name": "Vibrio tapetis", "Approx Genome Size (Mb)": 4.0}, "vvulnificus": {"Scientific Name": "Vibrio vulnificus", "Approx Genome Size (Mb)": 4.0}, "wolbachia": {"Scientific Name": "Wolbachia spp.", "Approx Genome Size (Mb)": 1.0}, "xfastidiosa": {"Scientific Name": "Xylella fastidiosa", "Approx Genome Size (Mb)": 2.5}, "yersinia": {"Scientific Name": "Yersinia spp.", "Approx Genome Size (Mb)": 4.5}, "ypseudotuberculosis": {"Scientific Name": "Yersinia pseudotuberculosis", "Approx Genome Size (Mb)": 5.0}, "yruckeri": {"Scientific Name": "Yersinia ruckeri", "Approx Genome Size (Mb)": 4.0}}

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
q301 = q30[0:4]
per = " \%"
q30p1 = '%s %s' % (q301, str(per))
q30 = merged_df.get("q30pass2", "n/a")
q302 = q30[0:4]
q30p2 = '%s %s' % (q302, str(per))
numscaff = merged_df.get("Number of scaffolds", "n/a")
sizescaff = "{:,}".format(merged_df.get("Total size of scaffolds", "n/a"))
sscaff = merged_df.get("Total size of scaffolds", "n/a")
longscaff = "{:,}".format(merged_df.get("Longest scaffold", "n/a"))
gt1k = merged_df.get("Number of scaffolds > 1K nt", "n/a")
n50len = merged_df.get("N50 scaffold length", "n/a")
l50cnt = merged_df.get("L50 scaffold count", "n/a")

cov = '%s %s' % (args.coverage, str(per))
ident = '%s %s' % (args.identity, str(per))
glength = "{:}".format(args.genlength)

####SeqSero Data
if seqsero_file is not None:
    seqsero_indx = pd.read_csv(seqsero_file, sep='\t', skiprows=5, index_col=0, header=None)
    predant = seqsero_indx.iloc[0][1]
    predsubsp = seqsero_indx.iloc[1][1]
    predst = seqsero_indx.iloc[2][1]
    seqserocomment = seqsero_indx.tail(1).index.values[0]

####MLST Data
if mlst_file is not None:
    try:
        mlst_data = pd.read_csv(mlst_file, sep='\t', header=None)

        mlst_id = mlst_data.iloc[0,1]
        mlst_id = mlst_id.partition('_')

        mlst_sp = lookup.get(mlst_id[0], {}).get("Scientific Name", "Not Identified")

        if mlst_sp == "Not Identified":
            mlst_st = "Not Identified"
            mlst_sch = "Not Identified"

        mlst_sch = mlst_id[2]
        mlst_st = mlst_data.iloc[0,2]
        mlst_g1 = mlst_data.iloc[0,3]
        mlst_g2 = mlst_data.iloc[0,4]
        mlst_g3 = mlst_data.iloc[0,5]
        mlst_g4 = mlst_data.iloc[0,6]
        mlst_g5 = mlst_data.iloc[0,7]
        mlst_g6 = mlst_data.iloc[0,8]
        mlst_g7 = mlst_data.iloc[0,9]
        #mlst_version = os.popen("mlst --version").readlines()[0]
    except IndexError:
        mlst_sp = "Not Identified"
        mlst_st = "Not Identified"
    mlst_version = os.popen("mlst --version").readlines()[0]

####Percent of Genome on contigs > 1k nt
lcontig_data = pd.read_csv(lcontigs_file, sep='\t', header=None)

gl = (lcontig_data[0][0]/sscaff)*100
rgl = round(gl, 2)
glc = '%s %s' % (str(rgl), str(per))

####Sequencing Scoring

seqcor=[[1.0000000,0.4412253,0.9569653,0.4266978,0.3393063],
[0.4412253,1.0000000,0.4446494,0.9527558,0.2464773],
[0.9569653,0.4446494,1.0000000,0.4330931,0.3315736],
[0.4266978,0.9527558,0.4330931,1.0000000,0.2367824],
[0.3393063,0.2464773,0.3315736,0.2367824,1.0000000]]

if glength != 'None':
   sp_len = glength
   sp_method = 'Genome length was provided for the analysis.'
else:
   sp_len = lookup.get(mlst_id[0], {}).get("Approx Genome Size (Mb)", "5")
   if mlst_id[0] == '-':
      sp_method = 'Default genome length of 5Mb was used.'
   else:
      sp_method = 'Genome length estimated based on MLST identification.'

gen_len = "{:,}".format(int(sp_len) * 1000000)
reads1 = merged_df.get("count1", "n/a")
reads2 = merged_df.get("count2", "n/a")
seqcov = int(reads1)*250/(int(sp_len) * 1000000)
x = [(float(quality1)-30), (float(quality2)-30), (float(q301)-80), (float(q302)-80), seqcov]
mx = np.matrix(x)
mseqcor = np.matrix(seqcor)
mxt = mx.transpose()
seqv = math.sqrt(mx*mseqcor*mxt)
seqval = (seqv*800/90) + 1180

if seqval > 1980 :
   seqval = 1980
elif seqval < 1180:
   seqval = 1180


#Title Header

mysvg = svgwrite.Drawing( size=(2100 ,180)) 


#Create horizontal gradient with red, yellow, and green
vert_grad = svgwrite.gradients.LinearGradient(start=(0,0), end=(0,1), id="vert_lin_grad")
vert_grad.add_stop_color(offset='0%', color='rgb(0,86,67)', opacity=None)
vert_grad.add_stop_color(offset='100%', color='rgb(0,45,114)', opacity=None)
#horz_grad.add_stop_color(offset='100%', color='limegreen', opacity=None)
mysvg.defs.add(vert_grad)


#draw a box and add gradient by about definition #id, adjust 5 so indicator can go outside gradient
#mysvg.add(mysvg.rect( [0,70], [2100,71], stroke='rgb(0,86,67)', fill='rgb(0,45,114)'))

mysvg.add(mysvg.line( start=(80,7), end=(2020,7), stroke_width="4",  stroke='black'))
mysvg.add(mysvg.line( start=(80,170), end=(2020,170), stroke_width="4",  stroke='black'))

#Add title on the x-axis
mysvg.add(mysvg.text('Bacterial Whole Genome Sequencing Report', insert=(80,120), fill='rgb(0,45,114)', font_size='90px', style = "font-family:Arial", font_weight="bold"))


mysvg.saveas('header.svg')
os.popen("inkscape -z --export-png=header.png --file=header.svg")

mysvg = svgwrite.Drawing( size=(2100 ,100))

#Create horizontal gradient with red, yellow, and green
horz_grad = svgwrite.gradients.LinearGradient(start=(0,0), end=(1,0), id="horz_lin_grad")
horz_grad.add_stop_color(offset='0%', color='crimson', opacity=None)
horz_grad.add_stop_color(offset='50%', color='yellow', opacity=None)
horz_grad.add_stop_color(offset='100%', color='limegreen', opacity=None)
mysvg.defs.add(horz_grad)

#draw a box and add gradient by about definition #id, adjust 5 so indicator can go outside gradient
mysvg.add(mysvg.rect( [0,0], [2100,100], stroke='rgb(0,44,118)', fill="rgb(0,44,118)"))

#Add title on the x-axis
mysvg.add(mysvg.text('Quality Scale', insert=(1450,40), fill='white', font_size='45px', style = "font-family:Arial"))

#draw a box and add gradient by about definition #id, adjust 5 so indicator can go outside gradient
mysvg.add(mysvg.rect( [1180,60], [800,10], stroke='black', fill="url(#horz_lin_grad)"))

#Add another indicator box
mysvg.add( mysvg.ellipse( [seqval,65], [10,15], fill='white', stroke_width=1, stroke='midnightblue', opacity=0.85))

#Add low quality label on the x-axis
mysvg.add(mysvg.text('Low', insert=(1100,75), fill='white', font_size='35px', font_weight="bold", style = "font-family:Arial"))

#Add high quality label on the x-axis
mysvg.add(mysvg.text('High', insert=(1995,75), fill='white', font_size='35px', font_weight="bold", style="font-family:Arial"))

#ADD File Stats header
mysvg.add(mysvg.text('Sequence Statistics', insert=(30,70), fill='white', font_size='50px', font_weight="bold"))

mysvg.saveas('seq.svg')

#os.popen("inkscape -D -z --file=seq.svg --export-pdf=seq.pdf --export-latex")
os.popen("inkscape -z --export-png=seq.png --file=seq.svg")


####Assembly Quality Scale

assemblecor=[[1.00000000,0.4640352,0.39099636,-0.703376083,-0.074391092,-0.27033841,-0.26923388,-0.36627871],
[0.46403525,1.0000000,0.75927035,-0.220017666,-0.101805557,-0.56550423,0.27504920,-0.70388508],
[0.39099636,0.7592704,1.00000000,-0.170697061,-0.054252944,-0.72295517,0.29979708,-0.81000716],
[-0.70337608,-0.2200177,-0.17069706,1.000000000,0.003833981,0.09220126,0.46237217,0.16084638],
[-0.07439109,-0.1018056,-0.05425294,0.003833981,1.000000000,0.05500119,-0.04451447,0.01902384],
[-0.27033841,-0.5655042,-0.72295517,0.092201260,0.055001186,1.00000000,-0.28768966,0.55006987],
[-0.26923388,0.2750492,0.29979708,0.462372175,-0.044514471,-0.28768966,1.00000000,-0.27338410],
[-0.36627871,-0.7038851,-0.81000716,0.160846384,0.019023843,0.55006987,-0.27338410,1.00000000]]


massemblecor = np.matrix(assemblecor)

lscaff = merged_df.get("Longest scaffold", "n/a")
perdiffexpect = abs(100-(float(sscaff)/(int(sp_len)*10000)))
perlongexpect = float(lscaff)/(int(sp_len)*1000000)*100
perlongscaff = float(gt1k)/float(numscaff)*100
pern50expect = int(n50len)/(int(sp_len)*1000000)*100

a=[numscaff, gt1k, l50cnt, rgl, perdiffexpect, perlongexpect, perlongscaff, pern50expect] 

ma = np.matrix(a)
mat = ma.transpose()
assemblev = 1000/(math.sqrt(ma*massemblecor*mat))
assembleval = (assemblev*800/7) + 1180

if assembleval > 1980 :
   assembleval = 1980
elif assembleval < 1180:
   assembleval = 1180


mysvg = svgwrite.Drawing( size=(2100 ,100))

#draw a box and add gradient by about definition #id, adjust 5 so indicator can go outside gradient
mysvg.add(mysvg.rect( [0,0], [2100,100], stroke='rgb(0,44,118)', fill="rgb(0,44,118)"))

#Create horizontal gradient with red, yellow, and green
horz_grad = svgwrite.gradients.LinearGradient(start=(0,0), end=(1,0), id="horz_lin_grad")
horz_grad.add_stop_color(offset='0%', color='crimson', opacity=None)
horz_grad.add_stop_color(offset='50%', color='yellow', opacity=None)
horz_grad.add_stop_color(offset='100%', color='limegreen', opacity=None)
mysvg.defs.add(horz_grad)

#draw a box and add gradient by about definition #id, adjust 5 so indicator can go outside gradient
mysvg.add(mysvg.rect( [1180,60], [800,10], stroke='black', fill="url(#horz_lin_grad)"))

#Add another indicator box
mysvg.add( mysvg.ellipse( [assembleval,65], [10,15], fill='white', stroke_width=1, stroke='midnightblue', opacity=0.85))

#ADD Assembly Stats header
mysvg.add(mysvg.text('Assembly Statistics', insert=(30,70), fill='white', font_size='50px', font_weight="bold"))

#Add low quality label on the x-axis
mysvg.add(mysvg.text('Low', insert=(1100,75), fill='white', font_size='35px', font_weight="bold", style = "font-family:Arial"))

#Add high quality label on the x-axis
mysvg.add(mysvg.text('High', insert=(1995,75), fill='white', font_size='35px', font_weight="bold", style="font-family:Arial"))

#Add title on the x-axis
mysvg.add(mysvg.text('Quality Scale', insert=(1450,40), fill='white', font_size='45px', style = "font-family:Arial"))

mysvg.saveas('assemble.svg')

os.popen("inkscape -z --export-png=assemble.png --file=assemble.svg")



####MLST BAR

mysvg = svgwrite.Drawing( size=(2100 ,100))

#Draw a box and add gradient by about definition #id, adjust 5 so indicator can go outside gradient
mysvg.add(mysvg.rect( [0,0], [2100,100], stroke='rgb(0,44,118)', fill="rgb(0,44,118)"))

#ADD File Stats header
mysvg.add(mysvg.text('Multi Locus Sequence Typing (MLST)', insert=(30,70), fill='white', font_size='50px', font_weight="bold"))

mysvg.saveas('mlst.svg')

os.popen("inkscape -z --export-png=mlst.png --file=mlst.svg")


####Salmonella Serotyping Bar

mysvg = svgwrite.Drawing( size=(2100 ,100))

#Draw a box and add gradient by about definition #id, adjust 5 so indicator can go outside gradient
mysvg.add(mysvg.rect( [0,0], [2100,100], stroke='rgb(0,44,118)', fill="rgb(0,44,118)"))

#ADD File Stats header
mysvg.add(mysvg.text('Serotyping for Salmonella Isolates', insert=(30,70), fill='white', font_size='50px', font_weight="bold"))

mysvg.saveas('seqsero.svg')

os.popen("inkscape -z --export-png=seqsero.png --file=seqsero.svg")



#####Abricate Header for ResFinder

mysvg = svgwrite.Drawing( size=(2100 ,140)) 

#draw a box and add gradient by about definition #id, adjust 5 so indicator can go outside gradient
mysvg.add(mysvg.rect( [0,0], [2100,100], stroke='rgb(0,44,118)', fill="rgb(0,44,118)"))

#Add title on the x-axis
mysvg.add(mysvg.text('ABRicate', insert=(30,70), fill='white', font_size='60px', font_weight="bold"))

#Add title on the x-axis
mysvg.add(mysvg.text('on the', insert=(363,70), fill='white', font_size='40px', style = "fill-opacity: .4"))

#Add title on the x-axis
mysvg.add(mysvg.text('ResFinder', insert=(505,70), fill='white', font_size='60px', font_weight="bold"))

#Add title on the x-axis
mysvg.add(mysvg.text('Database', insert=(880,70), fill='white', font_size='40px', style = "fill-opacity: 0.4"))

#draw a box and add gradient by about definition #id, adjust 5 so indicator can go outside gradient
mysvg.add(mysvg.rect( [0,102], [2100,140], fill='rgb(0,86,67)', style="fill-opacity: .2"))

#Add title on the x-axis
#mysvg.add(mysvg.text('SEQUENCE', insert=(50,127), fill='rgb(0,44,118)', font_size='22px', font_weight="bold"))

#Add title on the x-axis
#mysvg.add(mysvg.text('START', insert=(650,127), fill='rgb(0,44,118)', font_size='22px', font_weight="bold"))

#Add title on the x-axis
#mysvg.add(mysvg.text('END', insert=(803,127), fill='rgb(0,44,118)', font_size='22px', font_weight="bold"))

#Add title on the x-axis
#mysvg.add(mysvg.text('GENE', insert=(991,127), fill='rgb(0,44,118)', font_size='22px', font_weight="bold"))

#Add title on the x-axis
#mysvg.add(mysvg.text('COVERAGE', insert=(1165,127), fill='rgb(0,44,118)', font_size='22px', font_weight="bold"))

#Add title on the x-axis
#mysvg.add(mysvg.text('GAPS', insert=(1359,127), fill='rgb(0,44,118)', font_size='22px', font_weight="bold"))

#Add title on the x-axis
#mysvg.add(mysvg.text('%COV', insert=(1452,127), fill='rgb(0,44,118)', font_size='22px', font_weight="bold"))

#Add title on the x-axis
#mysvg.add(mysvg.text('%IDENT', insert=(1578,127), fill='rgb(0,44,118)', font_size='22px', font_weight="bold"))

#Add title on the x-axis
#mysvg.add(mysvg.text('ACCESSION', insert=(1730,127), fill='rgb(0,44,118)', font_size='22px', font_weight="bold"))

#Add title on the x-axis
#mysvg.add(mysvg.text('PRODUCT', insert=(1935,127), fill='rgb(0,44,118)', font_size='22px', font_weight="bold"))

#Save Resfinder Image
mysvg.saveas('resfinder.svg')
os.popen("inkscape -z --export-png=resfinder.png --file=resfinder.svg")


#####Abricate Header for NCBI

mysvg = svgwrite.Drawing( size=(2100 ,140)) 

#draw a box and add gradient by about definition #id, adjust 5 so indicator can go outside gradient
mysvg.add(mysvg.rect( [0,0], [2100,100], stroke='rgb(0,44,118)', fill="rgb(0,44,118)"))

#Add title on the x-axis
mysvg.add(mysvg.text('ABRicate', insert=(30,70), fill='white', font_size='60px', font_weight="bold"))

#Add title on the x-axis
mysvg.add(mysvg.text('on the', insert=(363,70), fill='white', font_size='40px', style = "fill-opacity: .4"))

#Add title on the x-axis
mysvg.add(mysvg.text('NCBI', insert=(505,70), fill='white', font_size='60px', font_weight="bold"))

#Add title on the x-axis
mysvg.add(mysvg.text('Database', insert=(700,70), fill='white', font_size='40px', style = "fill-opacity: 0.4"))

#draw a box and add gradient by about definition #id, adjust 5 so indicator can go outside gradient
mysvg.add(mysvg.rect( [0,102], [2100,140], fill='rgb(0,86,67)', style="fill-opacity: .2"))

#Add title on the x-axis
#mysvg.add(mysvg.text('SEQUENCE', insert=(50,127), fill='rgb(0,44,118)', font_size='22px', font_weight="bold"))

#Add title on the x-axis
#mysvg.add(mysvg.text('START', insert=(650,127), fill='rgb(0,44,118)', font_size='22px', font_weight="bold"))

#Add title on the x-axis
#mysvg.add(mysvg.text('END', insert=(803,127), fill='rgb(0,44,118)', font_size='22px', font_weight="bold"))

#Add title on the x-axis
#mysvg.add(mysvg.text('GENE', insert=(991,127), fill='rgb(0,44,118)', font_size='22px', font_weight="bold"))

#Add title on the x-axis
#mysvg.add(mysvg.text('COVERAGE', insert=(1165,127), fill='rgb(0,44,118)', font_size='22px', font_weight="bold"))

#Add title on the x-axis
#mysvg.add(mysvg.text('GAPS', insert=(1359,127), fill='rgb(0,44,118)', font_size='22px', font_weight="bold"))

#Add title on the x-axis
#mysvg.add(mysvg.text('%COV', insert=(1452,127), fill='rgb(0,44,118)', font_size='22px', font_weight="bold"))

#Add title on the x-axis
#mysvg.add(mysvg.text('%IDENT', insert=(1578,127), fill='rgb(0,44,118)', font_size='22px', font_weight="bold"))

#Add title on the x-axis
#mysvg.add(mysvg.text('ACCESSION', insert=(1730,127), fill='rgb(0,44,118)', font_size='22px', font_weight="bold"))

#Add title on the x-axis
#mysvg.add(mysvg.text('PRODUCT', insert=(1935,127), fill='rgb(0,44,118)', font_size='22px', font_weight="bold"))

#Save Resfinder Image
mysvg.saveas('ncbi.svg')
os.popen("inkscape -z --export-png=ncbi.png --file=ncbi.svg")

###AMRFinder result header

mysvg = svgwrite.Drawing( size=(2100 ,400)) 

#draw a box and add gradient by about definition #id, adjust 5 so indicator can go outside gradient
mysvg.add(mysvg.rect( [0,0], [2100,100], stroke='rgb(0,44,118)', fill='rgb(0,44,118)'))

#Add title on the x-axis
mysvg.add(mysvg.text('AMRFinder', insert=(30,70), fill='white', font_size='60px', style = "font-family:Arial", font_weight="bold"))

#draw a box and add gradient by about definition #id, adjust 5 so indicator can go outside gradient
mysvg.add(mysvg.rect( [0,102], [2100,400], fill='rgb(0,86,67)', style="fill-opacity: .25"))

#Add title on the x-axis
#mysvg.add(mysvg.text('Contig ID', insert=(10,280), fill='rgb(0,44,118)', font_size='22px', font_weight="bold"))

#Add title on the x-axis
#mysvg.add(mysvg.text('Start', insert=(200,280), fill='rgb(0,44,118)', font_size='22px', font_weight="bold"))

#Add title on the x-axis
#mysvg.add(mysvg.text('End', insert=(300,280), fill='rgb(0,44,118)', font_size='22px', font_weight="bold"))

#Add title on the x-axis
#mysvg.add(mysvg.text('Gene Symbol', insert=(400,280), fill='rgb(0,44,118)', font_size='22px', font_weight="bold"))

#Add title on the x-axis
#mysvg.add(mysvg.text('Protein Name', insert=(800,280), fill='rgb(0,44,118)', font_size='22px', font_weight="bold"))

#Add title on the x-axis
#mysvg.add(mysvg.text('Method', insert=(1200,280), fill='rgb(0,44,118)', font_size='22px', font_weight="bold"))

#Add title on the x-axis
#mysvg.add(mysvg.text('Target Length', insert=(1300,280), fill='rgb(0,44,118)', font_size='22px', font_weight="bold", transform = 'rotate(-45, 1350, 240)'))

#Add title on the x-axis
#mysvg.add(mysvg.text('Ref Protein Length', insert=(1400,280), fill='rgb(0,44,118)', font_size='22px', font_weight="bold", transform = 'rotate(-45, 1450, 240)'))

#Add title on the x-axis
#mysvg.add(mysvg.text('%Cov of Ref Protein', insert=(1500,280), fill='rgb(0,44,118)', font_size='22px', font_weight="bold", transform = 'rotate(-45, 1550, 240)'))

#Add title on the x-axis
#mysvg.add(mysvg.text('%Ident of Ref Protein', insert=(1600,280), fill='rgb(0,44,118)', font_size='22px', font_weight="bold", transform = 'rotate(-45, 1650, 240)'))

#Add title on the x-axis
#mysvg.add(mysvg.text('Alignment Length', insert=(1700,280), fill='rgb(0,44,118)', font_size='22px', font_weight="bold", transform = 'rotate(-45, 1750, 240)'))

#Add title on the x-axis
#mysvg.add(mysvg.text('Accession of', insert=(1890,280), fill='rgb(0,44,118)', font_weight="bold", font_size='22px', transform = 'rotate(-45, 1940, 240)'))
#Add title on the x-axis
#mysvg.add(mysvg.text(' Closest Protein', insert=(1910,280), fill='rgb(0,44,118)', font_weight="bold", font_size='22px', transform = 'rotate(-45, 1960, 240)'))

#Add title on the x-axis
#mysvg.add(mysvg.text('HMM ID', insert=(1950,280), fill='rgb(0,44,118)', font_size='22px', font_weight="bold", transform = 'rotate(-45, 2000, 240)'))


mysvg.saveas('amrfinder.svg')
os.popen("inkscape -z --export-png=amrfinder.png --file=amrfinder.svg")


mysvg = svgwrite.Drawing( size=(2100 ,140)) 

#draw a box and add gradient by about definition #id, adjust 5 so indicator can go outside gradient
mysvg.add(mysvg.rect( [0,0], [2100,100], stroke='rgb(0,44,118)', fill='rgb(0,44,118)'))

#Add title on the x-axis
mysvg.add(mysvg.text('AMRFinder', insert=(30,70), fill='white', font_size='60px', style = "font-family:Arial", font_weight="bold"))

#draw a box and add gradient by about definition #id, adjust 5 so indicator can go outside gradient
mysvg.add(mysvg.rect( [0,102], [2100,140], fill='rgb(0,86,67)', style="fill-opacity: .25"))

mysvg.saveas('amrfinder1.svg')
os.popen("inkscape -z --export-png=amrfinder1.png --file=amrfinder1.svg")


####Abricate
ab_version = os.popen("abricate --version | sed 's/abricate//g'").readlines()[0]
ncbi = os.popen("abricate --list | awk 'BEGIN{OFS=\"\t\"} $1 == \"ncbi\" {print $1, $2, $4}'|tr '\t' '&'").readlines()[0]
res = os.popen("abricate --list | awk 'BEGIN{OFS=\"\t\"} $1 == \"resfinder\" {print $1, $2, $4}'|tr '\t' '&'").readlines()[0]

##NCBI
ab_ncbi_file = glob.glob("*-ncbi-report.tab")


##Resfinder
ab_resfinder_file = glob.glob("*-resfinder-report.tab")


####AMRFinder
#cwl_version = os.popen("/usr/local/bin/amrfinder/amrfinder.pl -v| sed '1!d' |sed 's/:/\&/g'").readlines()[0]
#dock_version = os.popen("/usr/local/bin/amrfinder/amrfinder.pl -v| sed '3!d' |sed 's/container:/container\&/g'").readlines()[0]
#amrdb_version = os.popen("/usr/local/bin/amrfinder/amrfinder.pl -v| sed '4!d' | sed 's/: /\&/g'").readlines()[0]
amr_version = os.popen("/usr/local/bin/amrfinder/amrfinder.pl -v| sed '1!d' |sed 's/version $Revision://g'|sed 's/\$//g'").readlines()[0]

###Abricate documentation header
mysvg = svgwrite.Drawing( size=(2100 ,170)) 

#draw a box and add gradient by about definition #id, adjust 5 so indicator can go outside gradient
mysvg.add(mysvg.rect( [0,0], [2100,125], stroke='rgb(0,86,67)', fill='rgb(0,86,67)'))

#Add title on the x-axis
mysvg.add(mysvg.text('ABRicate', insert=(30,90), fill='white', font_size='70px', font_weight="bold"))

#Add title on the x-axis
mysvg.add(mysvg.text('Version %s' % (ab_version), insert=(1600,90), fill='white', font_size='55px', font_style = "italic"))

#draw a box and add gradient by about definition #id, adjust 5 so indicator can go outside gradient
mysvg.add(mysvg.rect( [0,110], [2100,170], stroke='rgb(0,44,118)', fill='rgb(0,44,118)'))

#Add title on the x-axis
mysvg.add(mysvg.text('https://github.com/tseemann/mlst', insert=(70,155), fill='white', font_size='40px'))

mysvg.saveas('abricate_doc.svg')
os.popen("inkscape -z --export-png=abricate_doc.png --file=abricate_doc.svg")


###AMRFinder documentation header
mysvg = svgwrite.Drawing( size=(2100 ,170)) 

#draw a box and add gradient by about definition #id, adjust 5 so indicator can go outside gradient
mysvg.add(mysvg.rect( [0,0], [2100,125], stroke='rgb(0,86,67)', fill='rgb(0,86,67)'))

#Add title on the x-axis
mysvg.add(mysvg.text('AMRFinder', insert=(30,90), fill='white', font_size='70px', font_weight="bold"))

#Add title on the x-axis
mysvg.add(mysvg.text('Revision %s' % (amr_version) , insert=(1600,90), fill='white', font_size='55px', font_style = "italic"))

#draw a box and add gradient by about definition #id, adjust 5 so indicator can go outside gradient
mysvg.add(mysvg.rect( [0,110], [2100,170], stroke='rgb(0,44,118)', fill='rgb(0,44,118)'))

#Add title on the x-axis
mysvg.add(mysvg.text('https://github.com/ncbi/amr/wiki', insert=(70,155), fill='white', font_size='40px'))

mysvg.saveas('amrfind_doc.svg')
os.popen("inkscape -z --export-png=amrfind_doc.png --file=amrfind_doc.svg")



amrfinder_file = glob.glob("*-amrfinder-report.tab")

spadesver = os.popen("spades.py --version").readlines()[0]

tex_document = sample_name + ".tex"

write_out = open(tex_document, 'w')

 #################################
            # LATEX

print(r'\documentclass[a4paper,12pt]{article}', file=write_out)
print(r'\usepackage[margin=0.5in]{geometry}', file=write_out)
print(r'\usepackage{graphicx}', file=write_out)
print(r'\usepackage[table]{xcolor}', file=write_out)
print(r'\usepackage{hyperref}', file=write_out)
print(r'\hypersetup{colorlinks = true, linkcolor = [RGB]{10,10,44}, urlcolor = [RGB]{10,10,44}, citecolor = [RGB]{10,10,44}, anchorcolor = [RGB]{10,10,44}}', file=write_out)
print(r'\usepackage{xcolor}', file=write_out)
print(r'\usepackage{tabularx}', file=write_out)
print(r'\usepackage{float}', file=write_out)
print(r'\usepackage{multirow}', file=write_out)
print(r'\usepackage{charter}', file=write_out)
print(r'\usepackage{mdwlist}', file=write_out)
print(r'\usepackage{fancyhdr}', file=write_out)
print(r'\usepackage{pdflscape}', file=write_out)
print(r'\usepackage{rotating}', file=write_out)
print(r'\usepackage[lastpage,user]{zref}', file=write_out)
print(r'\usepackage{wrapfig}', file=write_out)
print(r'\usepackage{calc}', file=write_out)
print(r'\usepackage{tcolorbox}', file=write_out)
print(r'\usepackage{tikz}', file=write_out)
print(r'\usepackage{multicol}', file=write_out)
print(r'\tcbuselibrary{skins}', file=write_out)
print(r'\setlength{\headheight}{30pt}', file=write_out)
print(r'\setlength{\footskip}{10pt}', file=write_out)
print(r'\fancyhead[L]{\includegraphics[scale=0.15]{/home/shared/usdalogo.png}} \fancyhead[R]{\textbf{Isolate ID:}{%s}}' % (sample_name), file=write_out)
print(r'\lfoot{\today}', file=write_out)
print(r'\cfoot{\thepage\ of \zref[abspage]{LastPage}}', file=write_out) 
print(r'\pagestyle{fancy}', file=write_out)
print(r'\thispagestyle{plain}', file=write_out)
print(r'\renewcommand{\thepage}{Page \arabic{page}}', file=write_out)
print(r'\includegraphics[scale=0.2]{/home/shared/usdalogo.png}', file=write_out)
print(r'\definecolor{midnightblue}{RGB}{0,44,118}', file=write_out)
print(r'\definecolor{usdagreen}{RGB}{0,84,67}', file=write_out)
print(r'\newcommand{\MYhref}[3][\textcolor{midnightblue}]{\href{#2}{\color{#1}{#3}}}%', file=write_out)
print(r'\begin{document}', file=write_out)
print(r'', file=write_out)
print(r'\vspace{7mm}', file=write_out)
#print(r'\centerline{\Huge {\bf{Whole Genome Sequencing Report}}}', file=write_out)
print(r'\includegraphics[scale=0.333]{header.png}', file=write_out)
print(r'', file=write_out)
print(r'\vspace{7mm}', file=write_out)
print(r'', file=write_out)
#print(r'\begin{minipage}[t]{0.5\textwidth}', file=write_out)
#print(r'\begin{multicols}{2}', file=write_out)
print(r'{\large \today}', file=write_out)
print(r'', file=write_out)
print(r'\vspace{4mm}', file=write_out)
print(r'', file=write_out)
print(r'\textbf{Sample ID:} {\large %s}' % (sample_name), file=write_out) 
#print(r'\end{minipage}', file=write_out)
print(r'\vspace{4mm}', file=write_out)
print(r'', file=write_out)
#print(r'\begin{minipage}[t]{0.5\textwidth}', file=write_out)
print(r'\textbf{Sequencing Technology}', file=write_out)
print(r'\vspace{2mm}', file=write_out)
print(r'', file=write_out)
print(r'\begin{tabular}{ p{11cm} }', file=write_out)
print(r'Nextera XT DNA Library Preparation \\', file=write_out)
print(r'MiSeq 2 x 250 Read Generation \\', file=write_out)
print(r'\end{tabular}', file=write_out)
#print(r'\end{minipage}', file=write_out)
#print(r'\end{multicols}', file=write_out)
print(r'', file=write_out)
print(r'\vspace{3mm}', file=write_out)
print(r'', file=write_out)
#print(r'\begin{wrapfigure}{r}{0.5\textwidth}', file=write_out)
#print(r'\begin{center}', file=write_out)
#print(r'\fontsize{10}{12}\selectfont', file=write_out)
#print(r'\fontfamily{cms}\selectfont', file=write_out)
#print(r'\vspace{-15mm}', file=write_out)
#print(r'\def\svgscale{0.5}', file=write_out)
#print(r'\input{seq.pdf_tex}', file=write_out)
#print(r'\end{center}', file=write_out)
#print(r'\end{wrapfigure}', file=write_out)
print(r'', file=write_out)
#print(r'\tcbset{colback=blue!5!white, colframe=blue!5!black}', file=write_out)
#print(r'\begin{tcolorbox}[enhanced, title style image=seq.png]', file=write_out)
#print(r'\tcbox[left=0mm, right=2mm, top=0mm, bottom=0mm, boxsep=0mm, toptitle=1mm, bottomtitle=1mm, title=\textbf{File Stats}]{%', file=write_out)
print(r'\includegraphics[scale=0.333]{seq.png}', file=write_out)
print(r'', file=write_out)
#print(r'\vspace{-1mm}', file=write_out)
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
print(r'\begin{tabular}{ p{3.2cm} | p{3cm} p{11cm} }', file=write_out)
#print(r'\hline', file=write_out)
print(r'Sequence Depth & %sX  & Calculated by Number of reads x 250/Genome Length \\' % (str(seqcov)), file=write_out)
print(r'\hline', file=write_out)
print(r'Genome Length & %sbp & %s \\' % (gen_len, sp_method), file=write_out)
#print(r'\#reads*250/genome length & %s\\ ' % (sp_method), file=write_out)
print(r'\hline', file=write_out)
print(r'\end{tabular}', file=write_out)
print(r'', file=write_out)
print(r'\vspace{5mm}', file=write_out)
print(r'', file=write_out)
#print(r'\textbf{Assembly}', file=write_out)
#print(r'\vspace{2mm}', file=write_out)
print(r'\includegraphics[scale=0.333]{assemble.png}', file=write_out)
print(r'', file=write_out)
print(r'\begin{tabularx}{\columnwidth}{X|X|X|X|X|X|X}', file=write_out)
print(r'\hline', file=write_out)
print(r'Scaffolds & Total length & Longest scaffold & Scaffolds \textgreater 1K nt & Genome \textgreater 1K nt & N50 & L50 \\', file=write_out)
print(r'\hline', file=write_out)
print(r'%s & %s & %s & %s & %s & %s & %s \\' % (numscaff, sizescaff, longscaff, gt1k, glc, n50len, l50cnt), file=write_out)
print(r'\hline', file=write_out)
print(r'\end{tabularx}', file=write_out)
print(r'', file=write_out)
print(r'{\footnotesize\ De novo assembly performed using \href{http://cab.spbu.ru/software/spades/} {\textcolor{midnightblue}{%s}.} }\\,' % spadesver, file=write_out)

if mlst_file is not None:
    print(r'\vspace{3mm}', file=write_out)
    print(r'', file=write_out)
    #print(r'\textbf{Multi Locus Sequence Typing (MLST)}', file=write_out)
    #print(r'\vspace{2mm}', file=write_out)
    #print(r'', file=write_out)
    print(r'\includegraphics[scale=0.333]{mlst.png}', file=write_out)
    print(r'', file=write_out)
    print(r'\vspace{-0.5mm}', file=write_out)
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
        print(r'{\footnotesize\ Data obatined using %s.  Software website: \href{https://github.com/tseemann/mlst} {\textcolor{midnightblue}{https://github.com/tseemann/mlst}}}\\' % mlst_version, file=write_out)
        print(r'{\footnotesize\ Information on MLST schemes and allelic profiles can be found at \href{https://pubmlst.org} {\textcolor{midnightblue}{pubMLST.org}} }\\', file=write_out) 
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
        print(r'{\footnotesize\ Data obtained using %s.  Software website: \href{https://github.com/tseemann/mlst} {\textcolor{midnightblue}{https://github.com/tseemann/mlst}}}\\' % mlst_version, file=write_out)
        print(r'{\footnotesize\ Information on MLST schemes and allelic profiles can be found at \href{https://pubmlst.org} {\textcolor{midnightblue}{pubMLST.org}} }\\', file=write_out) 
        print(r'\end{tabular}', file=write_out)
        print(r'', file=write_out)

    if mlst_id[0] == "senterica":  
        print(r'\vspace{4mm}', file=write_out)
        #print(r'\textbf{Serotyping for Salmonella Isolates}', file=write_out)
        #print(r'\vspace{2mm}', file=write_out)
        print(r'\includegraphics[scale=0.333]{seqsero.png}', file=write_out)
        print(r'', file=write_out)
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
        print(r'{\footnotesize\ Data obtained using SeqSero.  Software website: \href{https://github.com/denglab/SeqSero2} {\textcolor{midnightblue}{https://github.com/denglab/SeqSero2}} }\\', file=write_out)
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
    print(r'{\footnotesize\ Data obtained using SeqSero.  Software website: \href{https://github.com/denglab/SeqSero2} {\textcolor{midnightblue}{https://github.com/denglab/SeqSero2}} }\\', file=write_out)
    print(r'\end{tabular}', file=write_out)
    print(r'', file=write_out)

print(r'', file=write_out)

#AMR START
print(r'\newpage', file=write_out)
print(r'', file=write_out)
print(r'\begin{landscape}', file=write_out)
print(r'\vspace{10mm}', file=write_out)
print(r'', file=write_out)
#print(r'\centerline{\bf{\textcolor{midnightblue}{\fontsize{35}{45}{\selectfont Antimicrobial Resistance Analysis}}}}', file=write_out)
print(r'\centerline{\bf{\textcolor{midnightblue}{\Huge Antimicrobial Resistance Analysis}}}', file=write_out)
print(r'', file=write_out)
print(r'\vspace{5mm}', file=write_out)
print(r'\vspace{5mm}', file=write_out)
print(r'\noindent\makebox[\linewidth]{\rule{16cm}{0.4pt}}', file=write_out)
print(r'\vspace{-3mm}', file=write_out)
print(r'\begin{center}', file=write_out)
print(r'\begin{minipage}{8in}', file=write_out)
print(r'Results were obtained using \href{https://github.com/tseemann/abricate}{\textcolor{midnightblue}{ABRicate}} and \href{https://github.com/ncbi/pipelines/tree/master/amr\_finder}{\textcolor{midnightblue}{AMRFinder}}.  ABRicate locates genes based on nucleotide BLAST searching a database, which is specified in the filename. Searches return results when the minimum percent coverage and percent identity are met.  AMRFinder uses BLASTX to search a hierarchy of gene families with predetermined cutoffs.\\', file=write_out)
print(r'\end{minipage}', file=write_out)
print(r'\end{center}', file=write_out)
print(r'', file=write_out)
print(r'\vspace{-3mm}', file=write_out)
print(r'', file=write_out)
print(r'\noindent\makebox[\linewidth]{\rule{16cm}{0.4pt}}', file=write_out)
print(r'\vspace{15mm}', file=write_out)
print(r'', file=write_out)

#ABRICATE_DATA_OUTPUT
print(r'\includegraphics[scale=0.485]{resfinder.png}', file=write_out)
print(r'', file=write_out)
#print(r'\vspace{-0.5mm}', file=write_out)
print(r'', file=write_out)
if os.stat(ab_resfinder_file[0]).st_size == 0:
   print(r'\begin{tabular}{p{26.5cm}}', file=write_out)
   print(r'\hline', file=write_out)
   print(r'\vspace{1mm}', file=write_out)
   print(r'\centerline{\Large No Antimicrobial Resistance Genes found.}', file=write_out)
   print(r'\vspace{2mm}', file=write_out)
else:
   print(r'\vspace{-5mm}', file=write_out)
   print(r'\resizebox{27cm}{!}{', file=write_out)
   print(r'\begin{tabular}{ l|c|c|c|c|c|c|c|c|l }', file=write_out)
   print(r'\hline', file=write_out)
   print(r'SEQUENCE&START&END&GENE&COVERAGE&GAPS&\%COVERAGE&\%IDENTITY&ACCESSION&PRODUCT\\', file=write_out)
   print(r'\hline', file=write_out)
   with open(ab_resfinder_file[0], "r") as f:
      for line in f:
         line = line.strip()
         line = line + "\\\\"
         print(r'%s' % line, file=write_out)
print(r'\noalign{\global\arrayrulewidth=0.75mm}', file=write_out)
print(r'\arrayrulecolor{midnightblue}\hline', file=write_out)
print(r'\end{tabular}}', file=write_out)
print(r'', file=write_out)
print(r'\noalign{\global\arrayrulewidth=0.1mm}', file=write_out)
print(r'\vspace{20mm}', file=write_out)
print(r'', file=write_out)
print(r'\includegraphics[scale=0.485]{ncbi.png}', file=write_out)
print(r'', file=write_out)
if os.stat(ab_ncbi_file[0]).st_size ==0:
   print(r'\begin{tabular}{p{26.5cm}}', file=write_out)
   print(r'\hline', file=write_out)
   print(r'\vspace{1mm}', file=write_out)
   print(r'\centerline{\Large No Antimicrobial Resistance Genes found.}', file=write_out)
   print(r'\vspace{2mm}', file=write_out)
else:
   print(r'\vspace{-5mm}', file=write_out)
   print(r'\resizebox{27cm}{!}{', file=write_out)
   print(r'\begin{tabular}{ l|c|c|c|c|c|c|c|c|l }', file=write_out)
   print(r'\hline', file=write_out)
   print(r'SEQUENCE&START&END&GENE&COVERAGE&GAPS&\%COVERAGE&\%IDENTITY&ACCESSION&PRODUCT\\', file=write_out)
   print(r'\hline', file=write_out)
   with open(ab_ncbi_file[0], "r") as f:
      for line in f:
         line = line.strip()
         line = line + "\\\\"
         print(r'%s' % line, file=write_out)
print(r'\noalign{\global\arrayrulewidth=0.75mm}', file=write_out)
print(r'\arrayrulecolor{midnightblue}\hline', file=write_out)
print(r'\end{tabular}}', file=write_out)
print(r'', file=write_out)
print(r'\noalign{\global\arrayrulewidth=0.1mm}', file=write_out)
print(r'\vspace{5mm}', file=write_out)
print(r'', file=write_out)
#print(r'\end{landscape}', file=write_out)


#AMRFinder Results Page
print(r'\newpage', file=write_out)
print(r'', file=write_out)
print(r'\vspace{10mm}', file=write_out)
print(r'', file=write_out)#print(r'\begin{landscape}', file=write_out)
#print(r'\includegraphics[scale=0.485]{amrfinder.png}', file=write_out)
#print(r'', file=write_out)
#print(r'\vspace{30mm}', file=write_out)
if os.stat(amrfinder_file[0]).st_size == 0:
   print(r'\includegraphics[scale=0.485]{amrfinder1.png}', file=write_out)
   print(r'', file=write_out)
   print(r'\begin{tabular}{p{26.5cm}}', file=write_out)
   print(r'\hline', file=write_out)
   print(r'\vspace{1mm}', file=write_out)
   print(r'\centerline{\Large No Antimicrobial Resistance Genes found.}', file=write_out)
   print(r'\vspace{2mm}', file=write_out)
   print(r'\noalign{\global\arrayrulewidth=0.75mm}', file=write_out)
   print(r'\arrayrulecolor{midnightblue}\hline', file=write_out)
   print(r'\end{tabular}}', file=write_out)
else:
   print(r'\includegraphics[scale=0.485]{amrfinder.png}', file=write_out)
   print(r'', file=write_out)
   print(r'\vspace{-3mm}', file=write_out)
   print(r'\resizebox{27cm}{!}{', file=write_out)
   print(r'\begin{tabular}{c|c|c|c|c|c|c|c|c|c|c|c|c}', file=write_out)
   print(r'\textcolor{midnightblue}{\bf Contig id}&\textcolor{midnightblue}{\bf Start}&\textcolor{midnightblue}{\bf Stop}&\textcolor{midnightblue}{\bf Gene symbol}&\textcolor{midnightblue}{\bf Protein name}&\begin{rotate}{60}\textcolor{midnightblue}{\bf Method}\end{rotate}&\begin{rotate}{60}\textcolor{midnightblue}{\bf Target length}\end{rotate}&\begin{rotate}{60}\textcolor{midnightblue}{\bf Reference protein length}\end{rotate}& \begin{rotate}{60}\textcolor{midnightblue}{\bf \% Coverage of reference protein}\end{rotate}&\begin{rotate}{60}\textcolor{midnightblue}{\bf \% Identity of reference protein}\end{rotate}& \begin{rotate}{60}\textcolor{midnightblue}{\bf Alignment length}\end{rotate}& \begin{rotate}{60}\textcolor{midnightblue}{\bf Accession of closest protein} \end{rotate}& \begin{rotate}{60}\textcolor{midnightblue}{\bf HMM id}\end{rotate}  \\', file=write_out)
   print(r'\vspace{1mm}', file=write_out)
   print(r'\hline', file=write_out)
   with open(amrfinder_file[0], "r") as f:
      for line in f:
         line = line.strip()
         line = line + "\\\\"
         print(r'%s' % line, file=write_out)
   print(r'\noalign{\global\arrayrulewidth=0.75mm}', file=write_out)
   print(r'\arrayrulecolor{midnightblue}\hline', file=write_out)
   print(r'\end{tabular}}', file=write_out)
print(r'\noalign{\global\arrayrulewidth=0.1mm}', file=write_out)
print(r'\end{landscape}', file=write_out)


##Data Definitions
print(r'\newpage', file=write_out)
print(r'\begin{figure}', file=write_out)
print(r'\centering', file=write_out)
print(r'\vspace{-5mm}', file=write_out)
print(r'', file=write_out)
print(r'\includegraphics[scale=0.333]{abricate_doc.png}', file=write_out)
print(r'', file=write_out)
print(r'\end{figure}', file=write_out)
print(r'', file=write_out)
print(r'\begin{center}', file=write_out)
print(r'', file=write_out)
print(r'', file=write_out)
print(r'{\large Database Versions} \\', file=write_out)
print(r'\vspace{3mm}', file=write_out)
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
#print(r'\vspace{3mm}', file=write_out)
print(r'\begin{center}', file=write_out)
print(r'{\large Parameters}\\', file=write_out)
print(r'\vspace{3mm}', file=write_out)
print(r'\begin{tabular}{l|c}', file=write_out)
print(r'\hline', file=write_out)
print(r'Minimum coverage & %s \\ ' % (cov), file=write_out)
print(r'Minimum identity & %s \\' % (ident), file=write_out)
print(r'\hline', file=write_out)
print(r'\end{tabular}', file=write_out)
#print(r'\vspace{2mm}', file=write_out)
print(r'\end{center}', file=write_out)
print(r'', file=write_out)
#print(r'\vspace{3mm}', file=write_out)
print(r'\vspace{-2mm}', file=write_out)
print(r'', file=write_out)
print(r'\noindent', file=write_out)
#print(r'\begin{figure}', file=write_out)
print(r'\textbf{\large Summary Output}\\', file=write_out)
print(r'\vspace{-3mm}', file=write_out)
print(r'', file=write_out)
#print(r'\noindent\makebox[\linewidth]{\rule{19cm}{0.4pt}}', file=write_out)
print(r'\hrule', file=write_out)
print(r'', file=write_out)
print(r'\vspace{1mm}', file=write_out)
print(r'This abbreviated section appears in the first tab of each ABRicate Excel workbook. \\', file=write_out)
#print(r'\end{figure}', file=write_out)
print(r'', file=write_out)
print(r'\begin{center}', file=write_out)
print(r'\begin{minipage}{7in}', file=write_out)
print(r'\begin{basedescript}{\desclabelstyle{\pushlabel}}', file=write_out)
print(r'\item[NUM\_FOUND-] The number of distinct genes identified in the analysis of the sample.  This is NOT the total number of hits if gene duplicates are identified.', file=write_out)
print(r'\item[Gene List and Percent Coverage-] A gene list with percent coverage is also given by isolate.  If multiple identifications of the same gene were made in a single isolate then the percent coverage information is given in a colon separated list in the order the identifications are listed in the full analysis.', file=write_out)
print(r'\end{basedescript}', file=write_out)
print(r'\end{minipage}', file=write_out)
print(r'\end{center}', file=write_out)
print(r'', file=write_out)
#print(r'\newpage', file=write_out)
#print(r'', file=write_out)
print(r'\vspace{4mm}', file=write_out)
print(r'\noindent', file=write_out)
print(r'\textbf{\large Full Analysis Output}\\', file=write_out)
print(r'\vspace{-3mm}', file=write_out)
print(r'', file=write_out)
#print(r'\noindent\makebox[\linewidth]{\rule{19cm}{0.4pt}}', file=write_out)
print(r'\hrule', file=write_out)
print(r'', file=write_out)
print(r'\vspace{1mm}', file=write_out)
print(r'Details of the full analysis are located on subsequent tabs of the ABRicate workbooks.\\', file=write_out)
print(r'', file=write_out)
print(r'\begin{center}', file=write_out)
print(r'\begin{minipage}{7in}', file=write_out)
print(r'\begin{basedescript}{\desclabelstyle{\pushlabel}}', file=write_out)
print(r'\item[SEQUENCE-] The analysis is ran on scaffolds(also known as nodes or contigs) output by the assembler. Node length and coverage are included with the name.', file=write_out)
print(r'\item[FILE-] This field refers to the original file output by the analysis prior to the workbook creation.', file=write_out)
print(r'\item[START-] This is the position in the scaffold where the alignment with the given reference gene starts.', file=write_out)
print(r'\item[END-] This is the position of the scaffold where the alignment with the given reference gene ends.', file=write_out)
print(r'\item[GENE-] The reference gene to which the alignment with the scaffold is performed.', file=write_out)
print(r'\item[COVERAGE-] The positions (given as a range) of the reference gene that align with the scaffold divided by the length of the gene sequence. (aligned length/gene length)', file=write_out)
print(r'\item[COVERAGE\_MAP-] This gives an overview of the alignment relative to the reference gene. ....-no alignment ====-alignment /-gap in alignment.', file=write_out)
print(r'\item[GAPS-] The number of gaps in the alignment of the reference gene to the scaffold.', file=write_out)
print(r'\item[\%COV-] The percent of the reference gene covered by the alignment with the scaffold. (Percent coverage)', file=write_out)
print(r'\item[\%IDENT-] The percent identity of the reference gene with the scaffold. (Percent Identity)', file=write_out)
print(r'\item[DATABASE-] The gene database that was used for comparison to the scaffolds.', file=write_out)
print(r'\item[ACCESSION-] The accession number of the reference gene in the database used for the comparison.', file=write_out)
print(r'\end{basedescript}', file=write_out)
print(r'\end{minipage}', file=write_out)
print(r'\end{center}', file=write_out)
print(r'', file=write_out)




##AMRFinder

print(r'\newpage', file=write_out)
#print(r'', file=write_out)
print(r'\begin{figure}', file=write_out)
print(r'\centering', file=write_out)
#print(r'\begin{center}', file=write_out)
print(r'\vspace{-5mm}', file=write_out)
print(r'', file=write_out)
print(r'\includegraphics[scale=0.333]{amrfind_doc.png}', file=write_out)
print(r'Definitions were taken from the AMRFinder documentation.\\', file=write_out)
#print(r'\vspace{-3mm}', file=write_out)
#print(r'', file=write_out)
print(r'\noindent\makebox[\linewidth]{\rule{17cm}{0.4pt}}', file=write_out)
print(r'\end{figure}', file=write_out)
print(r'', file=write_out)
print(r'\begin{center}', file=write_out)
print(r'\begin{minipage}{6.5in}', file=write_out)
print(r'\begin{basedescript}{\desclabelstyle{\pushlabel}}', file=write_out)
print(r'\item[Target Identifier-] This is from the FASTA defline for the DNA sequence', file=write_out)
print(r'\item[Contig id-] Contig name', file=write_out)
print(r'\item[Start-] 1-based coordinate of first nucleotide coding from protein in DNA sequence on contig', file=write_out)
print(r'\item[Stop-] 1-based coordinate of last nucleotide coding for protein in DNA sequence on contig', file=write_out)
print(r'\item[Gene symbol-] Gene or gene-family symbol for protein hit', file=write_out)
print(r'\item[Protein name-] Full-text name for the protein', file=write_out)
print(r'\item[Method-] Type of hit found by AMRFinder one of five options', file=write_out)
print(r'\indent', file=write_out)
print(r'\item[ALLELE-] 100\% sequence match over 100\% of length to a protein named at the allele level in the AMRFinder database', file=write_out)
print(r'\item[EXACT-] 100\% sequence match over 100\% of length to a protein in the database that is not a named allele', file=write_out)
print(r'\item[BLAST-] BLAST alignment is \textgreater 90\% of length and \textgreater 90\% identity to a protein in a the AMRFinder database', file=write_out)
print(r'\item[PARTIAL-] BLAST alignment is \textgreater 50\% of length, but \textless 90\% of length and \textgreater 90\% identity', file=write_out)
print(r'\item[HMM-] HMM was hit above the cutoff, but there was not a BLAST hit that met standards for BLAST or PARTIAL', file=write_out)
print(r'\noident', file=write_out)
print(r'\item[Target length-] The length of the query protein. The length of the BLAST hit for translated-DNA searches', file=write_out)
print(r'\item[Reference protein length-] The length of the AMR Protein in the database (NA if HMM-only hit)', file=write_out)
print(r'\item[\% Coverage of reference protein-] \% covered by blast hit (NA if HMM-only hit', file=write_out)
print(r'\item[\% Identity to reference protein-] \% amino-acid identity to reference protein (NA if HMM-only hit)', file=write_out)
print(r'\end{basedescript}', file=write_out)
print(r'\end{minipage}', file=write_out)
print(r'\end{center}', file=write_out)
print(r'', file=write_out)
#print(r'\newpage', file=write_out)
print(r'\begin{center}', file=write_out)
print(r'\begin{minipage}{6.5in}', file=write_out)
print(r'\begin{basedescript}{\desclabelstyle{\pushlabel}}', file=write_out)
print(r'\item[Alignment length-] Length of BLAST alignment in amino acids (NA if HMM-only hit)', file=write_out)
print(r'\item[Accession of closest protein-] RefSeq accession for protein hit by BLAST (NA if HMM-only hit)', file=write_out)
print(r'\item[Name of closest protein-] Full name assigned to the AMRFinder database protein (NA if HMM-only hit)', file=write_out)
print(r'\item[HMM id-] Accession for the HMM', file=write_out)
print(r'\item[HMM description-] The family name associated with the HMM', file=write_out)
print(r'\end{basedescript}', file=write_out)
print(r'\end{minipage}', file=write_out)
print(r'\end{center}', file=write_out)
print(r'', file=write_out)


print(r'\end{document}', file=write_out)

write_out.close()

os.system("pdflatex -interaction=nonstopmode {}" .format(tex_document))
os.system("pdflatex -interaction=nonstopmode {}" .format(tex_document))

print(seqcov)
print(a)
print(assemblev, assembleval)
