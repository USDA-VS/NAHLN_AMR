#!/usr/bin/env python

__version__ = "0.0.1"

import os
import sys
import re
import argparse
import textwrap
import math
import numpy as np
import pandas as pd
import svgwrite

class AMR_Latex_Report:
    ''' 
    '''

    def __init__(self, **kwargs):
        fastq_scaling_variable = kwargs.get('fastq_scaling_variable', None)
        assembly_scaling_variable = kwargs.get('assembly_scaling_variable', None)
        self.genome_size = kwargs.get('genome_size', None)
        self.genome_coverage_depth = kwargs.get('genome_coverage_depth', None) # this is calculated on MLST identification/FASTQ else SPAdes
        self.coverage_method = kwargs.get('coverage_method', None)
        self.size_method = kwargs.get('size_method', None)
        abricate_ab_version = kwargs.get('abricate_ab_version', None)
        amr_version = kwargs.get('amr_version', None)
        read1_fastq = kwargs.get('read1_fastq', None)
        read2_fastq = kwargs.get('read2_fastq', None)

        #Create horizontal gradient with red, yellow, and green
        vert_grad = svgwrite.gradients.LinearGradient(start=(0,0), end=(0,1), id="vert_lin_grad")
        vert_grad.add_stop_color(offset='0%', color='rgb(0,86,67)', opacity=None)
        vert_grad.add_stop_color(offset='100%', color='rgb(0,45,114)', opacity=None)
        #horz_grad.add_stop_color(offset='100%', color='limegreen', opacity=None)
        mysvg = svgwrite.Drawing(size=(2100 ,180)) 
        mysvg.defs.add(vert_grad)

        #draw a box and add gradient by about definition #id, adjust 5 so indicator can go outside gradient
        mysvg.add(mysvg.line( start=(80,7), end=(2020,7), stroke_width="4",  stroke='black'))
        mysvg.add(mysvg.line( start=(80,170), end=(2020,170), stroke_width="4",  stroke='black'))

        #Add title on the x-axis
        mysvg.add(mysvg.text('Bacterial Whole Genome Sequencing Report', insert=(80,120), fill='rgb(0,45,114)', font_size='76px', style = "font-family:Arial", font_weight="bold"))

        mysvg.saveas('header.svg')
        os.popen("inkscape -z --export-png=header.png --file=header.svg")
        mysvg = svgwrite.Drawing( size=(2100 ,100))

        if read1_fastq:
            if fastq_scaling_variable:
                #Create horizontal gradient with red, yellow, and green
                horz_grad = svgwrite.gradients.LinearGradient(start=(0,0), end=(1,0), id="horz_lin_grad")
                horz_grad.add_stop_color(offset='0%', color='crimson', opacity=None)
                horz_grad.add_stop_color(offset='50%', color='yellow', opacity=None)
                horz_grad.add_stop_color(offset='100%', color='limegreen', opacity=None)
                mysvg.defs.add(horz_grad)

            #draw a box and add gradient by about definition #id, adjust 5 so indicator can go outside gradient
            mysvg.add(mysvg.rect( [0,0], [2100,100], stroke='rgb(0,44,118)', fill="rgb(0,44,118)"))
            if fastq_scaling_variable:
                #Add title on the x-axis
                mysvg.add(mysvg.text('Quality Scale', insert=(1450,40), fill='white', font_size='45px', style = "font-family:Arial"))
                #draw a box and add gradient by about definition #id, adjust 5 so indicator can go outside gradient
                mysvg.add(mysvg.rect( [1180,60], [800,10], stroke='black', fill="url(#horz_lin_grad)"))
                #Add another indicator box
                mysvg.add( mysvg.ellipse([fastq_scaling_variable,65], [10,15], fill='white', stroke_width=1, stroke='midnightblue', opacity=0.85))
                #Add low quality label on the x-axis
                mysvg.add(mysvg.text('Low', insert=(1100,75), fill='white', font_size='35px', font_weight="bold", style = "font-family:Arial"))
                #Add high quality label on the x-axis
                mysvg.add(mysvg.text('High', insert=(1995,75), fill='white', font_size='35px', font_weight="bold", style="font-family:Arial"))
            #ADD File Stats header
            mysvg.add(mysvg.text('Sequence Statistics', insert=(30,70), fill='white', font_size='50px', font_weight="bold"))
            mysvg.saveas('seq.svg')
            os.popen("inkscape -z --export-png=seq.png --file=seq.svg")

        ####Assembly Quality Scale
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
        mysvg.add( mysvg.ellipse( [assembly_scaling_variable,65], [10,15], fill='white', stroke_width=1, stroke='midnightblue', opacity=0.85))
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
        ###Abricate documentation header
        mysvg = svgwrite.Drawing( size=(2100 ,170)) 
        #draw a box and add gradient by about definition #id, adjust 5 so indicator can go outside gradient
        mysvg.add(mysvg.rect( [0,0], [2100,125], stroke='rgb(0,86,67)', fill='rgb(0,86,67)'))
        #Add title on the x-axis
        mysvg.add(mysvg.text('ABRicate', insert=(30,90), fill='white', font_size='70px', font_weight="bold"))
        #Add title on the x-axis
        mysvg.add(mysvg.text('Version %s' % (abricate_ab_version), insert=(1600,90), fill='white', font_size='55px', font_style = "italic"))
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

    def latex_document(self, **kwargs):
        LOGO = "/project/bioinformatic_databases/usdalogo.png"
        sample_name = kwargs.get('sample_name', None)
        read1_fastq = kwargs.get('read1_fastq', None)
        read2_fastq = kwargs.get('read2_fastq', None)
        read1_file_size = kwargs.get('read1_file_size', None)
        read2_file_size = kwargs.get('read2_file_size', None)
        read1_read_average = kwargs.get('read1_read_average', None)
        read2_read_average = kwargs.get('read2_read_average', None)
        read1_reads_gt_q30 = kwargs.get('read1_reads_gt_q30', None)
        read2_reads_gt_q30 = kwargs.get('read2_reads_gt_q30', None)
        sampling_size = kwargs.get('sampling_size', None)
        stat_contig_count = kwargs.get('stat_contig_count', None)
        stat_total_contig_lengths = kwargs.get('stat_total_contig_lengths', None)
        stat_longest_contig = kwargs.get('stat_longest_contig', None)
        stat_greater_one_kb_count = kwargs.get('stat_greater_one_kb_count', None)
        spades_version = kwargs.get('spades_version', None)
        stat_n50 = kwargs.get('stat_n50', None)
        stat_l50 = kwargs.get('stat_l50', None)
        mlst_file = kwargs.get('mlst_file', None)
        mlst_scheme = kwargs.get('mlst_scheme', None)
        mlst_st = kwargs.get('mlst_st', None)
        mlst_detail = kwargs.get('mlst_detail', None)
        mlst_species_lookup = kwargs.get('mlst_species_lookup', None)
        mlst_version = kwargs.get('mlst_version', None),
        seqsero2_serotype = kwargs.get('seqsero2_serotype', None)
        seqsero2_antigenic = kwargs.get('seqsero2_antigenic', None)
        seqsero2_subspecies = kwargs.get('seqsero2_subspecies', None)
        seqsero_file = kwargs.get('seqsero_file', None)
        seqserocomment = kwargs.get('seqserocomment', None)
        amrfinder_file = kwargs.get('amrfinder_file', None)
        mininum_report = kwargs.get('mininum_report', None)
        abricate_report = kwargs.get('abricate_report', None)
        abricate_depth = kwargs.get('abricate_depth', None)
        abricate_coverage = kwargs.get('abricate_coverage', None)
        ab_ncbi_file = kwargs.get('ab_ncbi_file', None)
        ab_resfinder_file = kwargs.get('ab_resfinder_file', None)
        rgl = kwargs.get('rgl', None)
        abricate_ncbi_version_date = kwargs.get('abricate_ncbi_version_date', None)
        abricate_res_version_date = kwargs.get('abricate_res_version_date', None)
        abricate_ncbi_seq_number = kwargs.get('abricate_ncbi_seq_number', None)
        abricate_res_seq_number = kwargs.get('abricate_res_seq_number', None)

        tex_document = sample_name + ".tex"
        write_out = open(tex_document, 'w')
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
        print(r'\fancyhead[L]{\includegraphics[scale=0.15]{%s}} \fancyhead[R]{\textbf{Isolate ID:}{%s}}' % (LOGO, sample_name), file=write_out)
        print(r'\lfoot{\today}', file=write_out)
        print(r'\cfoot{\thepage\ of \zref[abspage]{LastPage}}', file=write_out) 
        print(r'\pagestyle{fancy}', file=write_out)
        print(r'\thispagestyle{plain}', file=write_out)
        print(r'\renewcommand{\thepage}{Page \arabic{page}}', file=write_out)
        print(r'\includegraphics[scale=0.2]{%s}' % (LOGO), file=write_out)
        print(r'\definecolor{midnightblue}{RGB}{0,44,118}', file=write_out)
        print(r'\definecolor{usdagreen}{RGB}{0,84,67}', file=write_out)
        print(r'\newcommand{\MYhref}[3][\textcolor{midnightblue}]{\href{#2}{\color{#1}{#3}}}%', file=write_out)
        print(r'\begin{document}', file=write_out)
        print(r'', file=write_out)
        print(r'\vspace{7mm}', file=write_out)
        print(r'\includegraphics[scale=0.333]{header.png}', file=write_out)
        print(r'', file=write_out)
        print(r'\vspace{7mm}', file=write_out)
        print(r'', file=write_out)
        print(r'{\large \today}', file=write_out)
        print(r'', file=write_out)
        print(r'\vspace{4mm}', file=write_out)
        print(r'', file=write_out)
        print(r'\textbf{Sample ID:} {\large %s}' % (sample_name), file=write_out) 
        print(r'\vspace{4mm}', file=write_out)
        print(r'', file=write_out)
        if read1_fastq:
            if read2_fastq:
                q301 = (read1_reads_gt_q30/sampling_size) * 100
                q302 = (read2_reads_gt_q30/sampling_size) * 100
                print(r'\textbf{Sequencing Technology}', file=write_out)
                print(r'\vspace{2mm}', file=write_out)
                print(r'', file=write_out)
                print(r'\begin{tabular}{ p{11cm} }', file=write_out)
                print(r'Nextera XT DNA Library Preparation \\', file=write_out)
                print(r'MiSeq 2 x 250 Read Generation \\', file=write_out)
                print(r'\end{tabular}', file=write_out)
                print(r'', file=write_out)
                print(r'\vspace{3mm}', file=write_out)
                print(r'', file=write_out)
                print(r'\includegraphics[scale=0.333]{seq.png}', file=write_out)
                print(r'', file=write_out)
                print(r'', file=write_out)
                print(r'\begin{tabular}{ l | p{7cm} | p{7cm} }', file=write_out)
                print(r'\hline', file=write_out)
                print(r'Filename &  %s & %s \\' % (read1_fastq.replace('_', '-'), read2_fastq.replace('_', '-')), file=write_out )
                print(r'\hline', file=write_out)
                print(r'File size & %s & %s \\' % (read1_file_size, read2_file_size), file=write_out)
                print(f'Mean Read Score & {read1_read_average:.2f} & {read2_read_average:.2f} \\\\', file=write_out)
                print(f'Q30 Passing & {q301:.1f}\\% & {q302:.1f}\\% \\\\', file=write_out)
                print(r'\hline', file=write_out)
                print(r'\end{tabular}', file=write_out)
                print(r'', file=write_out)
                print(r'\begin{tabular}{ p{3.2cm} | p{3cm} | p{11cm} }', file=write_out)
                print(f'Sequence Depth & {self.genome_coverage_depth:,.1f}X  & {self.coverage_method} \\\\', file=write_out)
                print(r'\hline', file=write_out)
                print(f'Genome Length & {self.genome_size:,}bp & {self.size_method} \\\\', file = write_out)
                print(r'\hline', file=write_out)
                print(r'\end{tabular}', file=write_out)
                print(r'', file=write_out)
            else:
                q301 = (read1_reads_gt_q30/sampling_size) * 100
                print(r'\textbf{Sequencing Technology}', file=write_out)
                print(r'\vspace{2mm}', file=write_out)
                print(r'', file=write_out)
                print(r'\begin{tabular}{ p{11cm} }', file=write_out)
                print(r'Ion Torrent Read Generation \\', file=write_out)
                print(r'\end{tabular}', file=write_out)
                print(r'', file=write_out)
                print(r'\vspace{3mm}', file=write_out)
                print(r'', file=write_out)
                print(r'\includegraphics[scale=0.333]{seq.png}', file=write_out)
                print(r'', file=write_out)
                print(r'', file=write_out)
                print(r'\begin{tabular}{ l | p{14cm} }', file=write_out)
                print(r'\hline', file=write_out)
                print(r'Filename &  %s \\' % (read1_fastq.replace('_', '-')), file=write_out )
                print(r'\hline', file=write_out)
                print(r'File size & %s \\' % (read1_file_size), file=write_out)
                print(f'Mean Read Score & {read1_read_average:.2f} \\\\', file=write_out)
                print(f'Q30 Passing & {q301:.1f}\\% \\\\', file=write_out)
                print(r'\hline', file=write_out)
                print(r'\end{tabular}', file=write_out)
                print(r'', file=write_out)
                print(r'\begin{tabular}{ p{3.2cm} | p{3cm} | p{11cm} }', file=write_out)
                print(f'Sequence Depth & {self.genome_coverage_depth:,.1f}X  & {self.coverage_method} \\\\', file=write_out)
                print(r'\hline', file=write_out)
                print(f'Genome Length & {self.genome_size:,}bp & {self.size_method} \\\\', file = write_out)
                print(r'\hline', file=write_out)
                print(r'\end{tabular}', file=write_out)
                print(r'', file=write_out)
            
        print(r'\vspace{5mm}', file=write_out)
        print(r'', file=write_out)
        print(r'\includegraphics[scale=0.333]{assemble.png}', file=write_out)
        print(r'', file=write_out)
        print(r'\begin{tabularx}{\columnwidth}{X|X|X|X|X|X|X}', file=write_out)
        print(r'\hline', file=write_out)
        print(r'Scaffolds & Total length & Longest scaffold & Scaffolds \textgreater 1K nt & Genome \textgreater 1K nt & N50 & L50 \\', file=write_out)
        print(r'\hline', file=write_out)
        print(f'{stat_contig_count:,} & {stat_total_contig_lengths:,} & {stat_longest_contig:,} & {stat_greater_one_kb_count:,} & {rgl} & {stat_n50:,} & {stat_l50:,} \\\\', file=write_out)
        print(r'\hline', file=write_out)
        print(r'\end{tabularx}', file=write_out)
        if not read1_fastq:
            print(r'', file=write_out)
            print(r'\begin{tabularx}{\columnwidth}{ p{3.2cm} | p{3cm} | p{11cm} }', file=write_out)
            if self.genome_coverage_depth == 1:
                self.coverage_method = "Depth of coverage not provided in FASTA"
            print(f'Sequence Depth & {self.genome_coverage_depth:,.1f}X  & {self.coverage_method} \\\\', file=write_out)
            print(r'\hline', file=write_out)
            print(r'\end{tabularx}', file=write_out)
        print(r'', file=write_out)
        print(r'{\footnotesize\ De novo assembly performed using \href{http://cab.spbu.ru/software/spades/} {\textcolor{midnightblue}{%s}.} }\\' % spades_version, file=write_out)

        if mlst_file is not None:
            print(r'\vspace{3mm}', file=write_out)
            print(r'', file=write_out)
            print(r'\includegraphics[scale=0.333]{mlst.png}', file=write_out)
            print(r'', file=write_out)
            print(r'\vspace{-0.5mm}', file=write_out)
            print(r'', file=write_out)
            print(r'\begin{tabular}{ p{8.2cm} | p{9cm} }', file=write_out)
            print(r'\hline', file=write_out)

            if mlst_scheme == '':
                print(r'Organism ID: \textbf{\large No MLST scheme found} & Sequence type: \textbf{\large n/a} \\', file=write_out)
            else:
                print(r'Organism ID: \textbf{\large %s} & Schema-Sequence type: \textbf{\large %s-%s} \\' % (mlst_species_lookup, mlst_scheme.replace('_', '\_'), mlst_st), file=write_out)

            if mlst_st == 'Not Identified':
                print(r'\hline', file=write_out)
                print(r'\end{tabular}', file=write_out)
                print(r'', file=write_out)
                print(r'\begin{tabular}{ p{18cm} }', file=write_out)
                print(r'{\footnotesize\ Data obtained using %s.  Software website: \href{https://github.com/tseemann/mlst} {\textcolor{midnightblue}{https://github.com/tseemann/mlst}}}\\' % mlst_version, file=write_out)
                print(r'{\footnotesize\ Information on MLST schemes and allelic profiles can be found at \href{https://pubmlst.org} {\textcolor{midnightblue}{pubMLST.org}} }\\', file=write_out) 
                print(r'\end{tabular}', file=write_out)
                print(r'', file=write_out)
            else:
                print(r'\hline', file=write_out)
                print(r'\end{tabular}', file=write_out)
                print(r'', file=write_out)
                print(r'\begin{tabular} { p {17.2cm} }', file=write_out)
                print(r'\ MLST Detail: ', end='', file=write_out)
                for each in mlst_detail:
                    each = each.replace('_', '\_')
                    print(f'{each}', end=' ', file=write_out)
                print(f' \\\\', file=write_out)
                print(r'\hline', file=write_out)
                print(r'\end{tabular}', file=write_out)
                if not read1_fastq:
                    print(r'', file=write_out)
                    print(r'\begin{tabularx}{\columnwidth}{ p{3.2cm} | p{3cm} | p{11cm} }', file=write_out)
                    print(f'Genome Length & {self.genome_size:,}bp & {self.size_method} \\\\', file = write_out)
                    print(r'\hline', file=write_out)
                    print(r'\end{tabularx}', file=write_out)
                print(r'', file=write_out)
                print(r'\begin{tabular}{ p{17.2cm} }', file=write_out)
                print(r'{\footnotesize\ Data obtained using %s.  Software website: \href{https://github.com/tseemann/mlst} {\textcolor{midnightblue}{https://github.com/tseemann/mlst}}}\\' % mlst_version, file=write_out)
                print(r'{\footnotesize\ Information on MLST schemes and allelic profiles can be found at \href{https://pubmlst.org} {\textcolor{midnightblue}{pubMLST.org}} }\\', file=write_out) 
                print(r'\end{tabular}', file=write_out)
                print(r'', file=write_out)

        if seqsero_file != None and seqsero2_serotype != '-:-:-':
            print(r'\vspace{4mm}', file=write_out)
            print(r'\includegraphics[scale=0.333]{seqsero.png}', file=write_out)
            print(r'', file=write_out)
            print(r'', file=write_out)
            print(r'\begin{tabular}{ p{6.2cm} | p{11cm} }', file=write_out)
            print(r'\hline', file=write_out)
            print(r'Predicted serotype(s) & \textbf{\large %s} \\' % seqsero2_serotype, file=write_out)
            print(r'Predicted antigenic profile & %s \\' % seqsero2_antigenic, file=write_out)
            print(r'Predicted subspecies & %s \\' % seqsero2_subspecies, file=write_out)
            print(r'\hline', file=write_out)
            print(r'\end{tabular}', file=write_out)
            print(r'', file=write_out)
            print(r'\begin{tabular}{ p{17.2cm} }', file=write_out)
            if seqserocomment != 'Note:':
                 print(r'%s \\' % seqserocomment, file=write_out)
                 print(r'\hline', file=write_out)
            print(r'{\footnotesize\ Data obtained using SeqSero.  Software website: \href{https://github.com/denglab/SeqSero2} {\textcolor{midnightblue}{https://github.com/denglab/SeqSero2}} }\\', file=write_out)
            print(r'\end{tabular}', file=write_out)
            print(r'', file=write_out)

        print(r'', file=write_out)
        if not mininum_report:
            #AMR START
            print(r'\newpage', file=write_out)
            print(r'', file=write_out)
            print(r'\begin{landscape}', file=write_out)
            print(r'\vspace{10mm}', file=write_out)
            print(r'', file=write_out)
            print(r'\centerline{\bf{\textcolor{midnightblue}{\Huge Antimicrobial Resistance Analysis}}}', file=write_out)
            print(r'', file=write_out)
            print(r'\vspace{5mm}', file=write_out)
            print(r'\vspace{5mm}', file=write_out)
            print(r'\noindent\makebox[\linewidth]{\rule{16cm}{0.4pt}}', file=write_out)
            print(r'\vspace{-3mm}', file=write_out)
            print(r'\begin{center}', file=write_out)
            print(r'\begin{minipage}{8in}', file=write_out)
            if abricate_report:
                print(r'Results were obtained using \href{https://github.com/tseemann/abricate}{\textcolor{midnightblue}{ABRicate}} and \href{https://github.com/ncbi/pipelines/tree/master/amr\_finder}{\textcolor{midnightblue}{AMRFinder}}.  ABRicate locates genes based on nucleotide BLAST searching a database, which is specified in the filename. Searches return results when the minimum percent coverage and percent identity are met.  AMRFinder uses BLASTX to search a hierarchy of gene families with predetermined cutoffs.\\', file=write_out)
            else:
                print(r'Results were obtained using \href{https://github.com/ncbi/pipelines/tree/master/amr\_finder}{\textcolor{midnightblue}{AMRFinder}}.  AMRFinder uses BLASTX to search a hierarchy of gene families with predetermined cutoffs.\\', file=write_out)
            print(r'\end{minipage}', file=write_out)
            print(r'\end{center}', file=write_out)
            print(r'', file=write_out)
            print(r'\vspace{-3mm}', file=write_out)
            print(r'', file=write_out)
            print(r'\noindent\makebox[\linewidth]{\rule{16cm}{0.4pt}}', file=write_out)
            print(r'\vspace{15mm}', file=write_out)
            print(r'', file=write_out)

            if abricate_report:
                #ABRICATE_DATA_OUTPUT
                print(r'\includegraphics[scale=0.485]{resfinder.png}', file=write_out)
                print(r'', file=write_out)
                print(r'', file=write_out)
                ab_resfinder_lines = sum(1 for line in open(ab_resfinder_file))
                if ab_resfinder_lines > 1:
                    print(r'\vspace{-5mm}', file=write_out)
                    print(r'\resizebox{27cm}{!}{', file=write_out)
                    print(r'\begin{tabular}{ l|c|c|c|c|c|c|c|c|l }', file=write_out)
                    print(r'\hline', file=write_out)
                    print(r'SEQUENCE&START&END&GENE&COVERAGE&GAPS&\%COVERAGE&\%IDENTITY&ACCESSION&PRODUCT\\', file=write_out)
                    print(r'\hline', file=write_out)
                    df = pd.read_csv(ab_resfinder_file, sep='\t')
                    df = df[['SEQUENCE', 'START', 'END', 'GENE', 'COVERAGE', 'GAPS', '%COVERAGE', '%IDENTITY', 'ACCESSION', 'PRODUCT']]
                    for index_num in df.index:
                        series = df.T[index_num]
                        string = f'{series[0]}&{series[1]}&{series[2]}&{series[3]}&{series[4]}&{series[5]}&{series[6]}&{series[7]}&{series[8]}&{series[9]}'
                        fix = string.replace('_', '\_')
                        print(f'{fix}', end='\\\\\n', file=write_out)
                else:
                    print(r'\begin{tabular}{p{26.5cm}}', file=write_out)
                    print(r'\hline', file=write_out)
                    print(r'\vspace{1mm}', file=write_out)
                    print(r'\centerline{\Large No Antimicrobial Resistance Genes found.}', file=write_out)
                    print(r'\vspace{2mm}', file=write_out)
                print(r'\noalign{\global\arrayrulewidth=0.75mm}', file=write_out)
                print(r'\arrayrulecolor{midnightblue}\hline', file=write_out)
                print(r'\end{tabular}}', file=write_out)
                print(r'', file=write_out)
                print(r'\noalign{\global\arrayrulewidth=0.1mm}', file=write_out)
                print(r'\vspace{20mm}', file=write_out)
                print(r'', file=write_out)
                print(r'\includegraphics[scale=0.485]{ncbi.png}', file=write_out)
                print(r'', file=write_out)
                ab_ncbi_lines = sum(1 for line in open(ab_ncbi_file))
                if ab_ncbi_lines > 1:
                    print(r'\vspace{-5mm}', file=write_out)
                    print(r'\resizebox{27cm}{!}{', file=write_out)
                    print(r'\begin{tabular}{ l|c|c|c|c|c|c|c|c|l }', file=write_out)
                    print(r'\hline', file=write_out)
                    print(r'SEQUENCE&START&END&GENE&COVERAGE&GAPS&\%COVERAGE&\%IDENTITY&ACCESSION&PRODUCT\\', file=write_out)
                    print(r'\hline', file=write_out)
                    df = pd.read_csv(ab_ncbi_file, sep='\t')
                    df = df[['SEQUENCE', 'START', 'END', 'GENE', 'COVERAGE', 'GAPS', '%COVERAGE', '%IDENTITY', 'ACCESSION', 'PRODUCT']]
                    for index_num in df.index:
                        series = df.T[index_num]
                        string = f'{series[0]}&{series[1]}&{series[2]}&{series[3]}&{series[4]}&{series[5]}&{series[6]}&{series[7]}&{series[8]}&{series[9]}'
                        fix = string.replace('_', '\_')
                        print(f'{fix}', end='\\\\\n', file=write_out)
                else:
                    print(r'\begin{tabular}{p{26.5cm}}', file=write_out)
                    print(r'\hline', file=write_out)
                    print(r'\vspace{1mm}', file=write_out)
                    print(r'\centerline{\Large No Antimicrobial Resistance Genes found.}', file=write_out)
                    print(r'\vspace{2mm}', file=write_out)
                print(r'\noalign{\global\arrayrulewidth=0.75mm}', file=write_out)
                print(r'\arrayrulecolor{midnightblue}\hline', file=write_out)
                print(r'\end{tabular}}', file=write_out)
                print(r'', file=write_out)
                print(r'\noalign{\global\arrayrulewidth=0.1mm}', file=write_out)
                print(r'\vspace{5mm}', file=write_out)
                print(r'', file=write_out)
                print(r'\newpage', file=write_out)
                print(r'', file=write_out)
            #AMRFinder Results Page
            #print(r'\newpage', file=write_out)
            #print(r'', file=write_out)
            print(r'\vspace{10mm}', file=write_out)
            print(r'', file=write_out)
            amrfinder_lines = sum(1 for line in open(amrfinder_file))
            if amrfinder_lines > 1:
                print(r'\includegraphics[scale=0.485]{amrfinder.png}', file=write_out)
                print(r'', file=write_out)
                print(r'\vspace{-3mm}', file=write_out)
                print(r'\resizebox{27cm}{!}{', file=write_out)
                print(r'\begin{tabular}{c|c|c|c|c|c|c|c|c|c|c|c|c}', file=write_out)
                print(r'\textcolor{midnightblue}{\bf Contig id}&\textcolor{midnightblue}{\bf Start}&\textcolor{midnightblue}{\bf Stop}&\textcolor{midnightblue}{\bf Gene symbol}&\textcolor{midnightblue}{\bf Sequence name}&\begin{rotate}{60}\textcolor{midnightblue}{\bf Scope}\end{rotate}&\begin{rotate}{60}\textcolor{midnightblue}{\bf \% Coverage of reference seqeunce}\end{rotate}&\begin{rotate}{60}\textcolor{midnightblue}{\bf \% Identity to reference sequence}\end{rotate}&\textcolor{midnightblue}{\bf Class}&\textcolor{midnightblue}{\bf Subclass} \\', file=write_out)
                print(r'\vspace{1mm}', file=write_out)
                print(r'\hline', file=write_out)
                df = pd.read_csv(amrfinder_file, sep='\t', header=None, skiprows=1)
                for index_num in df.index:
                    series = df.T[index_num]
                    string = f'{series[1]}&{series[2]}&{series[3]}&{series[5]}&{series[6]}&{series[7]}&{series[15]}&{series[16]}&{series[10]}&{series[11]}'
                    fix = string.replace('_', '\_')
                    print(f'{fix}', end='\\\\\n', file=write_out)
            else:
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
            print(r'\noalign{\global\arrayrulewidth=0.1mm}', file=write_out)
            print(r'\end{landscape}', file=write_out)
            ##Data Definitions
            if abricate_report:
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
                 print(f' NCBI & {int(abricate_ncbi_seq_number):,} & {abricate_ncbi_version_date} \\\\', file=write_out)
                 print(f' ResFinder & {int(abricate_res_seq_number):,} & {abricate_res_version_date} \\\\', file=write_out)
                 print(r'\hline', file=write_out)
                 print(r'\end{tabular}', file=write_out)
                 print(r'\end{center}', file=write_out)
                 print(r'', file=write_out)
                 print(r'\begin{center}', file=write_out)
                 print(r'{\large Parameters}\\', file=write_out)
                 print(r'\vspace{3mm}', file=write_out)
                 print(r'\begin{tabular}{l|c}', file=write_out)
                 print(r'\hline', file=write_out)
                 print(f'Minimum coverage & {abricate_depth}\% \\\\', file=write_out)
                 print(f'Minimum identity & {abricate_coverage}\% \\\\', file=write_out)
                 print(r'\hline', file=write_out)
                 print(r'\end{tabular}', file=write_out)
                 print(r'\end{center}', file=write_out)
                 print(r'', file=write_out)
                 print(r'\vspace{-2mm}', file=write_out)
                 print(r'', file=write_out)
                 print(r'\noindent', file=write_out)
                 print(r'\textbf{\large Summary Output}\\', file=write_out)
                 print(r'\vspace{-3mm}', file=write_out)
                 print(r'', file=write_out)
                 print(r'\hrule', file=write_out)
                 print(r'', file=write_out)
                 print(r'\vspace{1mm}', file=write_out)
                 print(r'This abbreviated section appears in the first tab of each ABRicate Excel workbook. \\', file=write_out)
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
                 print(r'\vspace{4mm}', file=write_out)
                 print(r'\noindent', file=write_out)
                 print(r'\textbf{\large Full Analysis Output}\\', file=write_out)
                 print(r'\vspace{-3mm}', file=write_out)
                 print(r'', file=write_out)
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
            print(r'\begin{figure}', file=write_out)
            print(r'\centering', file=write_out)
            print(r'\vspace{-5mm}', file=write_out)
            print(r'', file=write_out)
            print(r'\includegraphics[scale=0.333]{amrfind_doc.png}', file=write_out)
            print(r'Definitions were taken from the AMRFinder documentation.\\', file=write_out)
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
        os.system(f'pdflatex -interaction=nonstopmode {tex_document} > /dev/null 2>&1')
        os.system(f'pdflatex -interaction=nonstopmode {tex_document} > /dev/null 2>&1')    

if __name__ == "__main__": # execute if directly access by the interpreter
    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

    ---------------------------------------------------------
    Place description

    '''), epilog='''---------------------------------------------------------''')
    
    parser.add_argument('-n', '--name', action='store', dest='name', help='description')
    parser.add_argument('-x', '--xname', action='store', dest='name', required=True, default=85, help='Required with default value-- description')
    parser.add_argument('-b', '--myboolean', action='store_true', dest='myboolean', help='description')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.basename(__file__)}: version {__version__}')
    args = parser.parse_args()
    
    print(f'\n{os.path.basename(__file__)} SET ARGUMENTS:')
    print(args)
    print("\n")

# Created March 2020 by Tod Stuber
