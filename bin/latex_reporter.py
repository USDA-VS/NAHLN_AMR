#!/usr/bin/env python

__version__ = "0.0.2"

import os
import sys
import re
import argparse
import textwrap
import math
import numpy as np
import pandas as pd
import svgwrite
from cairosvg import svg2png
from PIL import Image


class AMR_Latex_Report:
    def __init__(self, **kwargs):
        fastq_scaling_variable   = kwargs.get('fastq_scaling_variable', None)
        assembly_scaling_variable = kwargs.get('assembly_scaling_variable', None)
        self.genome_size          = kwargs.get('genome_size', None)
        self.genome_coverage_depth = kwargs.get('genome_coverage_depth', None)
        self.coverage_method      = kwargs.get('coverage_method', None)
        self.size_method          = kwargs.get('size_method', None)
        self.cur_dir              = os.getcwd()
        abricate_ab_version       = kwargs.get('abricate_ab_version', None)
        amr_version               = kwargs.get('amr_version', None)
        read1_fastq               = kwargs.get('read1_fastq', None)
        read2_fastq               = kwargs.get('read2_fastq', None)
        self.read1_read_length    = kwargs.get('read1_read_length', None)
        self.read2_read_length    = kwargs.get('read2_read_length', None)
        self.logo_path            = kwargs.get('logo_path', None)

        # ── Header banner ────────────────────────────────────────────────────
        vert_grad = svgwrite.gradients.LinearGradient(
            start=(0, 0), end=(0, 1), id="vert_lin_grad")
        vert_grad.add_stop_color(offset='0%',   color='rgb(0,86,67)',  opacity=None)
        vert_grad.add_stop_color(offset='100%', color='rgb(0,45,114)', opacity=None)
        mysvg = svgwrite.Drawing(size=(2100, 180))
        mysvg.defs.add(vert_grad)
        mysvg.add(mysvg.line(start=(80, 7),   end=(2020, 7),   stroke_width="4", stroke='black'))
        mysvg.add(mysvg.line(start=(80, 170), end=(2020, 170), stroke_width="4", stroke='black'))
        mysvg.add(mysvg.text('Bacterial Whole Genome Sequencing Report',
                             insert=(80, 120), fill='rgb(0,45,114)',
                             font_size='76px', style="font-family:Arial", font_weight="bold"))
        mysvg.saveas(os.path.join(self.cur_dir, 'header.svg'))
        with open('header.svg', 'r') as f:
            svg2png(bytestring=f.read(), write_to='header.png')

        # ── Sequence Statistics banner (FASTQ only) ───────────────────────
        mysvg = svgwrite.Drawing(size=(2100, 100))
        if read1_fastq:
            if fastq_scaling_variable:
                horz_grad = svgwrite.gradients.LinearGradient(
                    start=(0, 0), end=(1, 0), id="horz_lin_grad")
                horz_grad.add_stop_color(offset='0%',   color='crimson',   opacity=None)
                horz_grad.add_stop_color(offset='50%',  color='yellow',    opacity=None)
                horz_grad.add_stop_color(offset='100%', color='limegreen', opacity=None)
                mysvg.defs.add(horz_grad)
            mysvg.add(mysvg.rect([0, 0], [2100, 100],
                                 stroke='rgb(0,44,118)', fill="rgb(0,44,118)"))
            if fastq_scaling_variable:
                mysvg.add(mysvg.text('Quality Scale', insert=(1450, 40),
                                     fill='white', font_size='45px', style="font-family:Arial"))
                mysvg.add(mysvg.rect([1180, 60], [800, 10],
                                     stroke='black', fill="url(#horz_lin_grad)"))
                mysvg.add(mysvg.ellipse([fastq_scaling_variable, 65], [10, 15],
                                        fill='white', stroke_width=1,
                                        stroke='midnightblue', opacity=0.85))
                mysvg.add(mysvg.text('Low',  insert=(1100, 75), fill='white',
                                     font_size='35px', font_weight="bold", style="font-family:Arial"))
                mysvg.add(mysvg.text('High', insert=(1995, 75), fill='white',
                                     font_size='35px', font_weight="bold", style="font-family:Arial"))
            mysvg.add(mysvg.text('Sequence Statistics', insert=(30, 70),
                                 fill='white', font_size='50px', font_weight="bold"))
            mysvg.saveas('seq.svg')
            with open('seq.svg', 'r') as f:
                svg2png(bytestring=f.read(), write_to='seq.png')

        # ── Assembly Statistics banner ────────────────────────────────────
        mysvg = svgwrite.Drawing(size=(2100, 100))
        mysvg.add(mysvg.rect([0, 0], [2100, 100],
                             stroke='rgb(0,44,118)', fill="rgb(0,44,118)"))
        horz_grad = svgwrite.gradients.LinearGradient(
            start=(0, 0), end=(1, 0), id="horz_lin_grad")
        horz_grad.add_stop_color(offset='0%',   color='crimson',   opacity=None)
        horz_grad.add_stop_color(offset='50%',  color='yellow',    opacity=None)
        horz_grad.add_stop_color(offset='100%', color='limegreen', opacity=None)
        mysvg.defs.add(horz_grad)
        mysvg.add(mysvg.rect([1180, 60], [800, 10],
                             stroke='black', fill="url(#horz_lin_grad)"))
        mysvg.add(mysvg.ellipse([assembly_scaling_variable, 65], [10, 15],
                                fill='white', stroke_width=1,
                                stroke='midnightblue', opacity=0.85))
        mysvg.add(mysvg.text('Assembly Statistics', insert=(30, 70),
                             fill='white', font_size='50px', font_weight="bold"))
        mysvg.add(mysvg.text('Low',  insert=(1100, 75), fill='white',
                             font_size='35px', font_weight="bold", style="font-family:Arial"))
        mysvg.add(mysvg.text('High', insert=(1995, 75), fill='white',
                             font_size='35px', font_weight="bold", style="font-family:Arial"))
        mysvg.add(mysvg.text('Quality Scale', insert=(1450, 40),
                             fill='white', font_size='45px', style="font-family:Arial"))
        mysvg.saveas('assemble.svg')
        with open('assemble.svg', 'r') as f:
            svg2png(bytestring=f.read(), write_to='assemble.png')

        # ── MLST banner ───────────────────────────────────────────────────
        mysvg = svgwrite.Drawing(size=(2100, 100))
        mysvg.add(mysvg.rect([0, 0], [2100, 100],
                             stroke='rgb(0,44,118)', fill="rgb(0,44,118)"))
        mysvg.add(mysvg.text('Multi Locus Sequence Typing (MLST)', insert=(30, 70),
                             fill='white', font_size='50px', font_weight="bold"))
        mysvg.saveas('mlst.svg')
        svg2png(url='mlst.svg', write_to='mlst.png')

        # ── SeqSero banner ────────────────────────────────────────────────
        mysvg = svgwrite.Drawing(size=(2100, 100))
        mysvg.add(mysvg.rect([0, 0], [2100, 100],
                             stroke='rgb(0,44,118)', fill="rgb(0,44,118)"))
        mysvg.add(mysvg.text('Serotyping for Salmonella Isolates', insert=(30, 70),
                             fill='white', font_size='50px', font_weight="bold"))
        mysvg.saveas('seqsero.svg')
        with open('seqsero.svg', 'r') as f:
            svg2png(bytestring=f.read(), write_to='seqsero.png')

        # ── ABRicate / ResFinder banner ───────────────────────────────────
        def _abricate_banner(db_label, db_x, filename_base):
            d = svgwrite.Drawing(size=(2100, 140))
            d.add(d.rect([0, 0],   [2100, 100], stroke='rgb(0,44,118)', fill='rgb(0,44,118)'))
            d.add(d.text('ABRicate', insert=(30, 70),   fill='white', font_size='60px', font_weight="bold"))
            d.add(d.text('on the',   insert=(363, 70),  fill='white', font_size='40px', style="fill-opacity:.4"))
            d.add(d.text(db_label,   insert=(505, 70),  fill='white', font_size='60px', font_weight="bold"))
            d.add(d.text('Database', insert=(db_x, 70), fill='white', font_size='40px', style="fill-opacity:0.4"))
            d.add(d.rect([0, 102], [2100, 140], fill='rgb(0,86,67)', style="fill-opacity:.2"))
            d.saveas(f'{filename_base}.svg')
            with open(f'{filename_base}.svg', 'r') as f:
                svg2png(bytestring=f.read(), write_to=f'{filename_base}.png')

        _abricate_banner('ResFinder', 880, 'resfinder')
        _abricate_banner('NCBI',      700, 'ncbi')

        # ── AMRFinder banners (results + doc-page variant) ────────────────
        for svg_name, save_fn in [('amrfinder.svg',  lambda d: svg2png(bytestring=open('amrfinder.svg').read(),  write_to='amrfinder.png')),
                                   ('amrfinder1.svg', lambda d: svg2png(url='amrfinder1.svg', write_to='amrfinder1.png'))]:
            d = svgwrite.Drawing(size=(2100, 140))
            d.add(d.rect([0, 0],   [2100, 100], stroke='rgb(0,44,118)', fill='rgb(0,44,118)'))
            d.add(d.text('AMRFinder', insert=(30, 70), fill='white',
                         font_size='60px', style="font-family:Arial", font_weight="bold"))
            d.add(d.rect([0, 102], [2100, 140], fill='rgb(0,86,67)', style="fill-opacity:.25"))
            d.saveas(svg_name)
            save_fn(d)

        # ── Documentation banners ─────────────────────────────────────────
        def _doc_banner(title, version_str, url_str, filename_base):
            d = svgwrite.Drawing(size=(2100, 170))
            d.add(d.rect([0, 0],   [2100, 125], stroke='rgb(0,86,67)', fill='rgb(0,86,67)'))
            d.add(d.text(title,       insert=(30, 90),   fill='white', font_size='70px', font_weight="bold"))
            d.add(d.text(version_str, insert=(1600, 90), fill='white', font_size='55px', font_style="italic"))
            d.add(d.rect([0, 110],  [2100, 170], stroke='rgb(0,44,118)', fill='rgb(0,44,118)'))
            d.add(d.text(url_str,     insert=(70, 155),  fill='white', font_size='40px'))
            d.saveas(f'{filename_base}.svg')
            with open(f'{filename_base}.svg', 'r') as f:
                svg2png(bytestring=f.read(), write_to=f'{filename_base}.png')

        _doc_banner('ABRicate',  f'Version {abricate_ab_version}',
                    'https://github.com/tseemann/mlst',      'abricate_doc')
        _doc_banner('AMRFinder', f'Revision {amr_version}',
                    'https://github.com/ncbi/amr/wiki',      'amrfind_doc')

    def escape_latex(self, text):
        """Escape all special LaTeX characters.  \\ must be first."""
        if not isinstance(text, str):
            text = str(text)
        replacements = [
            ('\\', r'\textbackslash{}'),
            ('&',  r'\&'),
            ('%',  r'\%'),
            ('$',  r'\$'),
            ('#',  r'\#'),
            ('_',  r'\_'),
            ('{',  r'\{'),
            ('}',  r'\}'),
            ('~',  r'\textasciitilde{}'),
            ('^',  r'\textasciicircum{}'),
        ]
        for char, replacement in replacements:
            text = text.replace(char, replacement)
        return text

    # helpers
    def _e(self, value):
        """Escape any value for safe LaTeX output."""
        return self.escape_latex(str(value) if value is not None else '')

    def _fmt(self, value, spec=','):
        """Format a number with commas (or other spec) then escape."""
        try:
            return self.escape_latex(format(int(value), spec))
        except (ValueError, TypeError):
            return self._e(value)

    def _p(self, text, f):
        """Single-line helper: print text to file."""
        print(text, file=f)

    # latex_document
    def latex_document(self, logo_path,  sequencing_type=None, **kwargs):
        #LOGO = os.path.join(os.path.dirname(__file__), '../docs/img/usdalogo.png')
        LOGO = logo_path

        # all wargs
        sample_name             = kwargs.get('sample_name', '')
        read1_fastq             = kwargs.get('read1_fastq', None)
        read2_fastq             = kwargs.get('read2_fastq', None)
        read1_file_size         = kwargs.get('read1_file_size', None)
        read2_file_size         = kwargs.get('read2_file_size', None)
        read1_read_average      = kwargs.get('read1_read_average', None)
        read2_read_average      = kwargs.get('read2_read_average', None)
        read1_read_length       = kwargs.get('read1_read_length', None)
        read2_read_length       = kwargs.get('read2_read_length', None)
        read1_reads_gt_q30      = kwargs.get('read1_reads_gt_q30', None)
        read2_reads_gt_q30      = kwargs.get('read2_reads_gt_q30', None)
        sampling_size           = kwargs.get('sampling_size', None)
        stat_contig_count       = kwargs.get('stat_contig_count', None)
        stat_total_contig_lengths = kwargs.get('stat_total_contig_lengths', None)
        stat_longest_contig     = kwargs.get('stat_longest_contig', None)
        stat_greater_one_kb_count = kwargs.get('stat_greater_one_kb_count', None)
        spades_version          = kwargs.get('spades_version', None)
        stat_n50                = kwargs.get('stat_n50', None)
        stat_l50                = kwargs.get('stat_l50', None)
        mlst_file               = kwargs.get('mlst_file', None)
        mlst_scheme             = kwargs.get('mlst_scheme', None)
        mlst_st                 = kwargs.get('mlst_st', None)
        mlst_detail             = kwargs.get('mlst_detail', None)
        mlst_species_lookup     = kwargs.get('mlst_species_lookup', None)
        mlst_version            = kwargs.get('mlst_version', None),
        seqsero2_serotype       = kwargs.get('seqsero2_serotype', None)
        seqsero2_antigenic      = kwargs.get('seqsero2_antigenic', None)
        seqsero2_subspecies     = kwargs.get('seqsero2_subspecies', None)
        seqsero_file            = kwargs.get('seqsero_file', None)
        seqserocomment          = kwargs.get('seqserocomment', None)
        amrfinder_file          = kwargs.get('amrfinder_file', None)
        mininum_report          = kwargs.get('mininum_report', None)
        abricate_report         = kwargs.get('abricate_report', None)
        abricate_mincov         = kwargs.get('abricate_mincov', None)  # Legacy "depth" param
        abricate_minid          = kwargs.get('abricate_minid', None)    # Legacy "coverage" param
        ab_ncbi_file            = kwargs.get('ab_ncbi_file', None)
        ab_resfinder_file       = kwargs.get('ab_resfinder_file', None)
        rgl                     = kwargs.get('rgl', None)
        abricate_ncbi_version_date = kwargs.get('abricate_ncbi_version_date', None)
        abricate_res_version_date  = kwargs.get('abricate_res_version_date', None)
        abricate_ncbi_seq_number   = kwargs.get('abricate_ncbi_seq_number', None)
        abricate_res_seq_number    = kwargs.get('abricate_res_seq_number', None)
        software_versions       = kwargs.get('software_versions', None)
        bracken_pie_file        = kwargs.get('bracken_pie_file', None)

        def clean_filename(path):
            if not path:
                return ''
            fname = os.path.basename(path)
            if fname.startswith('reads_input/'):
                fname = fname.replace('reads_input/', '')
            if fname.startswith('reads input/'):
                fname = fname.replace('reads input/', '')
            return fname

        safe_sample_name    = self._e(sample_name)
        safe_seq_type       = self._e(sequencing_type).replace('\\_', ' ') if sequencing_type else None
        safe_r1_fastq       = self._e(clean_filename(read1_fastq))
        safe_r2_fastq       = self._e(clean_filename(read2_fastq))
        safe_r1_size        = self._e(read1_file_size)  if read1_file_size  else 'N/A'
        safe_r2_size        = self._e(read2_file_size)  if read2_file_size  else 'N/A'
        safe_coverage_method = self._e(self.coverage_method)
        safe_size_method    = self._e(self.size_method)
        safe_spades         = self._e(spades_version or 'N/A')
        safe_mlst_scheme    = self._e(mlst_scheme)    if mlst_scheme    else ''
        safe_mlst_st        = self._e(mlst_st)        if mlst_st        else ''
        safe_mlst_lookup    = self._e(mlst_species_lookup) if mlst_species_lookup else ''
        safe_mlst_version   = self._e(mlst_version)   if mlst_version   else 'N/A'
        safe_serotype       = self._e(seqsero2_serotype)  if seqsero2_serotype  else ''
        safe_antigenic      = self._e(seqsero2_antigenic) if seqsero2_antigenic else ''
        safe_subspecies     = self._e(seqsero2_subspecies) if seqsero2_subspecies else ''
        safe_seqcomment     = self._e(seqserocomment)  if seqserocomment else ''
        safe_ncbi_date      = self._e(abricate_ncbi_version_date) if abricate_ncbi_version_date else 'N/A'
        safe_res_date       = self._e(abricate_res_version_date)  if abricate_res_version_date  else 'N/A'
        safe_rgl            = self._e(rgl)             if rgl is not None else '0'
        safe_depth_str      = f'{self.genome_coverage_depth:,.1f}X'
        safe_genome_bp      = f'{self.genome_size:,}bp'

        tex_document = sample_name + ".tex"
        w = open(tex_document, 'w')

        for line in [
            r'\documentclass[a4paper,12pt]{article}',
            r'\usepackage[top=0.4in,bottom=0.4in,left=0.5in,right=0.5in]{geometry}',
            r'\usepackage{graphicx}',
            r'\usepackage[table]{xcolor}',
            r'\usepackage{hyperref}',
            r'\usepackage{seqsplit}',
            r'\hypersetup{colorlinks=true,linkcolor=[RGB]{10,10,44},urlcolor=[RGB]{10,10,44},'
            r'citecolor=[RGB]{10,10,44},anchorcolor=[RGB]{10,10,44}}',
            r'\usepackage{xcolor}',
            r'\usepackage{tabularx}',
            r'\usepackage{float}',
            r'\usepackage{multirow}',
            r'\usepackage{charter}',
            r'\usepackage{mdwlist}',
            r'\usepackage{fancyhdr}',
            r'\usepackage{pdflscape}',
            r'\usepackage{rotating}',
            r'\usepackage[lastpage,user]{zref}',
            r'\usepackage{wrapfig}',
            r'\usepackage{calc}',
            r'\usepackage{tcolorbox}',
            r'\usepackage{tikz}',
            r'\usepackage{multicol}',
            r'\usepackage{longtable}',
            r'\tcbuselibrary{skins}',
            r'\setlength{\headheight}{40pt}',
            r'\setlength{\footskip}{15pt}',
        ]:
            print(line, file=w)

        print(r'\fancyhead[L]{\includegraphics[scale=0.15]{%s}} \fancyhead[R]{\textbf{Isolate ID:}{%s}}'
              % (LOGO, safe_sample_name), file=w)
        print(r'\lfoot{\today}', file=w)
        print(r'\cfoot{\thepage\ of \zref[abspage]{LastPage}}', file=w)
        for line in [
            r'\pagestyle{fancy}',
            r'\thispagestyle{plain}',
            r'\renewcommand{\thepage}{Page \arabic{page}}',
        ]:
            print(line, file=w)
        print(r'\includegraphics[scale=0.2]{%s}' % LOGO, file=w)
        print(r'\definecolor{midnightblue}{RGB}{0,44,118}', file=w)
        print(r'\definecolor{usdagreen}{RGB}{0,84,67}', file=w)
        print(r'\newcommand{\MYhref}[3][\textcolor{midnightblue}]{\href{#2}{\color{#1}{#3}}}%', file=w)
        print(r'\begin{document}', file=w)
        print(r'\vspace{7mm}', file=w)
        print(r'\includegraphics[scale=0.25]{header.png}', file=w)
        print(r'\vspace{7mm}', file=w)
        print(r'{\large \today}\\', file=w)
        print(r'\vspace{4mm}', file=w)
        print(r'\textbf{Sample ID:} {\large %s}\\' % safe_sample_name, file=w)
        print(r'\vspace{4mm}', file=w)

        if read1_fastq:
            if safe_seq_type:
                print(r'\textbf{Sequencing Technology:}\\', file=w)
                print(r'\vspace{2mm}', file=w)
                print(r'\begin{tabular}{l}', file=w)
                print(f'{safe_seq_type} \\\\', file=w)
                print(r'\end{tabular}', file=w)
                print(r'\\', file=w)
                print(r'\vspace{3mm}', file=w)

            print(r'\includegraphics[scale=0.25]{seq.png}', file=w)

            # Paired-End
            if read2_fastq:
                q301 = (read1_reads_gt_q30 / sampling_size * 100) if sampling_size and sampling_size > 0 else 0
                q302 = (read2_reads_gt_q30 / sampling_size * 100) if sampling_size and sampling_size > 0 else 0
                print(r'\begin{tabular}{ l | p{7cm} | p{7cm} }', file=w)
                print(r'\hline', file=w)
                print(r'\textbf{Filename} & \textbf{\seqsplit{%s}} & \textbf{\seqsplit{%s}} \\[2ex]'
                      % (safe_r1_fastq, safe_r2_fastq), file=w)
                print(r'\hline', file=w)
                print(r'File size & %s & %s \\' % (safe_r1_size, safe_r2_size), file=w)
                print(f'Mean Read Score & {read1_read_average:.2f} & {read2_read_average:.2f} \\\\', file=w)
                print(f'Mean Read Length & {read1_read_length:.1f}bp & {read2_read_length:.1f}bp \\\\', file=w)
                print(f'Q30 Passing & {q301:.1f}\\% & {q302:.1f}\\% \\\\', file=w)
                print(r'\hline', file=w)
                print(r'\end{tabular}', file=w)
            # Single-End
            else:
                q301 = (read1_reads_gt_q30 / sampling_size * 100) if sampling_size and sampling_size > 0 else 0
                print(r'\begin{tabular}{ l | p{7cm} }', file=w)
                print(r'\hline', file=w)
                print(r'\textbf{Filename} & \textbf{\seqsplit{%s}} \\[2ex]' % safe_r1_fastq, file=w)
                print(r'\hline', file=w)
                print(r'File size & %s \\' % safe_r1_size, file=w)
                print(f'Mean Read Score & {read1_read_average:.2f} \\\\', file=w)
                print(f'Mean Read Length & {read1_read_length:.1f}bp \\\\', file=w)
                print(f'Q30 Passing & {q301:.1f}\\% \\\\', file=w)
                print(r'\hline', file=w)
                print(r'\end{tabular}', file=w)

            # Bracken pie chart
            if bracken_pie_file and os.path.exists(bracken_pie_file):
                print(r'\begin{figure}[H]', file=w)
                print(r'\centering', file=w)
                print(r'\includegraphics[width=0.8\textwidth]{%s}' % bracken_pie_file, file=w)
                print(r"\caption{Relative abundance of identified species as estimated by Bracken. "
                      r"Species representing less than 2\% of the total are grouped into the "
                      r"`Other' category.}", file=w)
                print(r'\end{figure}', file=w)

        # Assembly statistics
        print(r'\\', file=w)
        print(r'\vspace{5mm}', file=w)
        print(r'\noindent\begin{minipage}{\columnwidth}', file=w)
        print(r'\includegraphics[scale=0.25]{assemble.png}', file=w)
        print(r'\begin{tabularx}{\columnwidth}{X|X|X|X|X|X|X}', file=w)
        print(r'\hline', file=w)
        print(r'Scaffolds & Total length & Longest scaffold & Scaffolds \textgreater 1K nt '
              r'& Genome \textgreater 1K nt & N50 & L50 \\', file=w)
        print(r'\hline', file=w)
        print('%s & %s & %s & %s & %s & %s & %s \\\\'
              % (self._fmt(stat_contig_count),
                 self._fmt(stat_total_contig_lengths),
                 self._fmt(stat_longest_contig),
                 self._fmt(stat_greater_one_kb_count),
                 safe_rgl,
                 self._fmt(stat_n50),
                 self._fmt(stat_l50)), file=w)
        print(r'\hline', file=w)

        # Sequence Depth + Genome Length rows (merged into assembly table)
        if read1_fastq:
            print(r'\multicolumn{2}{l}{Sequence Depth} & \multicolumn{5}{l}{%s \quad %s} \\'
                  % (safe_depth_str, safe_coverage_method), file=w)
            print(r'\hline', file=w)
            print(r'\multicolumn{2}{l}{Genome Length} & \multicolumn{5}{l}{%s \quad %s} \\'
                  % (safe_genome_bp, safe_size_method), file=w)
            print(r'\hline', file=w)
        else:
            # FASTA-only: depth row only
            depth_method = self._e("Depth of coverage not provided in FASTA") \
                if self.genome_coverage_depth == 1 else safe_coverage_method
            print(r'\multicolumn{2}{l}{Sequence Depth} & \multicolumn{5}{l}{%s \quad %s} \\'
                  % (safe_depth_str, depth_method), file=w)
            print(r'\hline', file=w)

        print(r'\end{tabularx}', file=w)
        print(r'\end{minipage}', file=w)
        print(r'', file=w)
        print(r'{\footnotesize\ De novo assembly performed using '
              r'\href{http://cab.spbu.ru/software/spades/}'
              r'{\textcolor{midnightblue}{%s}.} }\\' % safe_spades, file=w)

        # MLST section
        if mlst_file is not None:
            safe_mlst_org = self._e(mlst_scheme) if mlst_scheme else 'N/A'
            
            print(r'\vspace{4mm}', file=w)
            print(r'\noindent\begin{minipage}{\textwidth}', file=w)
            print(r'\noindent', file=w)
            print(r'\includegraphics[scale=0.25]{mlst.png}', file=w)
            print(r'\vspace{2mm}', file=w)
            print(r'\noindent', file=w)
            print(r'\begin{tabular}{ p{8.5cm} | p{8.5cm} }', file=w)
            print(r'\hline', file=w)
            
            if mlst_scheme == '':
                print(r'Organism ID: \textbf{\large No MLST scheme found} '
                    r'& Sequence type: \textbf{\large n/a} \\', file=w)
            else:
                print(r'Organism ID: \textbf{\large %s} & Schema-Sequence type: \textbf{\large %s-%s} \\'
                    % (safe_mlst_lookup, safe_mlst_scheme, safe_mlst_st), file=w)

            if mlst_st == 'Not Identified':
                print(r'\hline', file=w)
                print(r'\end{tabular}', file=w)
                print(r'\\[2mm]', file=w)
                print(r'\noindent', file=w)
                print(r'\begin{tabular}{ p{17.2cm} }', file=w)
                print(r'{\footnotesize\ Data obtained using %s.  '
                    r'Software website: \href{https://github.com/tseemann/mlst}'
                    r'{\textcolor{midnightblue}{https://github.com/tseemann/mlst}}}\\' % safe_mlst_version, file=w)
                print(r'{\footnotesize\ Information on MLST schemes and allelic profiles can be found at '
                    r'\href{https://pubmlst.org}{\textcolor{midnightblue}{pubMLST.org}} }\\', file=w)
                print(r'\end{tabular}', file=w)
            else:
                print(r'\hline', file=w)
                print(r'\end{tabular}', file=w)
                print(r'\\[2mm]', file=w)
                print(r'\noindent', file=w)
                print(r'\begin{tabular}{ p{17.2cm} }', file=w)
                print(r'\ MLST Detail: ', end='', file=w)
                for token in (mlst_detail or []):
                    print(self._e(token), end=' ', file=w)
                print(r' \\', file=w)
                print(r'\hline', file=w)
                print(r'\end{tabular}', file=w)

                if not read1_fastq:
                    print(r'\\[2mm]', file=w)
                    print(r'\noindent', file=w)
                    print(r'\begin{tabularx}{\columnwidth}{ p{3.2cm} | p{3cm} | p{11cm} }', file=w)
                    print(f'Genome Length & {safe_genome_bp} & {safe_size_method} \\\\', file=w)
                    print(r'\hline', file=w)
                    print(r'\end{tabularx}', file=w)

                print(r'\\[2mm]', file=w)
                print(r'\noindent', file=w)
                print(r'\begin{tabular}{ p{17.2cm} }', file=w)
                print(r'{\footnotesize\ Data obtained using %s.  '
                    r'Software website: \href{https://github.com/tseemann/mlst}'
                    r'{\textcolor{midnightblue}{https://github.com/tseemann/mlst}}}\\' % safe_mlst_version, file=w)
                print(r'{\footnotesize\ Information on MLST schemes and allelic profiles can be found at '
                    r'\href{https://pubmlst.org}{\textcolor{midnightblue}{pubMLST.org}} }\\', file=w)
                print(r'\end{tabular}', file=w)
            
            print(r'\end{minipage}', file=w)

        # SeqSero
        if seqsero_file is not None and seqsero2_serotype != '- -:-:-':
            print(r'\noindent\begin{minipage}{\textwidth}', file=w)
            print(r'\noindent', file=w)
            print(r'\includegraphics[scale=0.25]{seqsero.png}', file=w)
            print(r'\vspace{2mm}', file=w)
            print(r'\noindent', file=w)
            print(r'\begin{tabular}{ p{6.2cm} | p{11cm} }', file=w)
            print(r'\hline', file=w)
            print(r'Predicted serotype(s) & \textbf{\large %s} \\' % safe_serotype, file=w)
            print(r'Predicted antigenic profile & %s \\'              % safe_antigenic, file=w)
            print(r'Predicted subspecies & %s \\'                     % safe_subspecies, file=w)
            print(r'\hline', file=w)
            print(r'\end{tabular}', file=w)
            print(r'\\[2mm]', file=w)
            print(r'\noindent', file=w)
            print(r'\begin{tabular}{ p{17.2cm} }', file=w)
            if seqserocomment != 'Note:':
                print(r'%s \\' % safe_seqcomment, file=w)
                print(r'\hline', file=w)
            print(r'{\footnotesize\ Data obtained using SeqSero.  '
                r'Software website: \href{https://github.com/denglab/SeqSero2}'
                r'{\textcolor{midnightblue}{https://github.com/denglab/SeqSero2}} }\\', file=w)
            print(r'\end{tabular}', file=w)
            print(r'\end{minipage}', file=w)

        # AMR sectio
        if not mininum_report:
            print(r'\newpage', file=w)
            print(r'\fancyhf{}', file=w)
            print(r'\fancyhead[L]{\includegraphics[scale=0.15]{%s}}' % LOGO, file=w)
            print(r'\fancyhead[R]{\textbf{Isolate ID:}{%s}}' % safe_sample_name, file=w)
            print(r'\cfoot{\thepage\ of \zref[abspage]{LastPage}}', file=w)
            print(r'\begin{landscape}', file=w)
            print(r'\setlength{\headheight}{25pt}', file=w)
            print(r'\setlength{\footskip}{10pt}', file=w)
            print(r'\thispagestyle{fancy}', file=w)
            print(r'\vspace{3mm}', file=w)
            print(r'\centerline{\bf{\textcolor{midnightblue}{\Huge Antimicrobial Resistance Analysis}}}', file=w)
            print(r'\vspace{5mm}', file=w)
            print(r'\noindent\makebox[\linewidth]{\rule{16cm}{0.4pt}}', file=w)
            print(r'\vspace{-3mm}', file=w)
            print(r'\begin{center}', file=w)
            print(r'\begin{minipage}{8in}', file=w)
            if abricate_report:
                print(r'Results were obtained using '
                      r'\href{https://github.com/tseemann/abricate}{\textcolor{midnightblue}{ABRicate}} '
                      r'and \href{https://github.com/ncbi/pipelines/tree/master/amr\_finder}'
                      r'{\textcolor{midnightblue}{AMRFinder}}.  ABRicate locates genes based on nucleotide '
                      r'BLAST searching a database, which is specified in the filename. Searches return '
                      r'results when the minimum percent coverage and percent identity are met.  '
                      r'AMRFinder uses BLASTX to search a hierarchy of gene families with predetermined cutoffs.\\', file=w)
            else:
                print(r'Results were obtained using '
                      r'\href{https://github.com/ncbi/pipelines/tree/master/amr\_finder}'
                      r'{\textcolor{midnightblue}{AMRFinder}}.  AMRFinder uses BLASTX to search a '
                      r'hierarchy of gene families with predetermined cutoffs.\\', file=w)
            print(r'\end{minipage}', file=w)
            print(r'\end{center}', file=w)
            print(r'\vspace{-3mm}', file=w)
            print(r'\noindent\makebox[\linewidth]{\rule{16cm}{0.4pt}}', file=w)
            print(r'\vspace{15mm}', file=w)

            # ABRicate tables 
            ABRICATE_COLS = ['SEQUENCE', 'START', 'END', 'GENE', 'COVERAGE',
                             'GAPS', '%COVERAGE', '%IDENTITY', 'ACCESSION', 'PRODUCT']

            def _abricate_table(tab_file, banner_png):
                print(r'\noindent', file=w)
                print(r'\includegraphics[scale=0.25]{%s}' % banner_png, file=w)
                line_count = 0
                if tab_file and os.path.exists(tab_file):
                    line_count = sum(1 for _ in open(tab_file))
                if line_count > 1:
                    print(r'\vspace{2mm}', file=w)
                    print(r'\noindent', file=w)
                    print(r'{\small', file=w)
                    print(r'\renewcommand{\arraystretch}{0.9}', file=w)
                    # Use longtable instead of tabular for automatic page breaks
                    print(r'\begin{longtable}{ p{3.5cm}|c|c|c|c|c|c|c|c|p{5cm} }', file=w)
                    # First page header
                    print(r'\hline', file=w)
                    print(r'SEQUENCE&START&END&GENE&COVERAGE&GAPS&\%COVERAGE&\%IDENTITY&ACCESSION&PRODUCT\\', file=w)
                    print(r'\hline', file=w)
                    print(r'\endfirsthead', file=w)
                    # Continuation header for subsequent pages
                    print(r'\multicolumn{10}{c}{\textit{\small (Continued from previous page)}} \\', file=w)
                    print(r'\hline', file=w)
                    print(r'SEQUENCE&START&END&GENE&COVERAGE&GAPS&\%COVERAGE&\%IDENTITY&ACCESSION&PRODUCT\\', file=w)
                    print(r'\hline', file=w)
                    print(r'\endhead', file=w)
                    # Footer for continuation pages
                    print(r'\hline', file=w)
                    print(r'\multicolumn{10}{r}{\textit{\small (Continued on next page)}} \\', file=w)
                    print(r'\endfoot', file=w)
                    # Footer for last page
                    print(r'\noalign{\global\arrayrulewidth=0.75mm}', file=w)
                    print(r'\arrayrulecolor{midnightblue}\hline', file=w)
                    print(r'\endlastfoot', file=w)
                    # Data rows
                    df = pd.read_csv(tab_file, sep='\t')
                    df = df[ABRICATE_COLS]
                    for idx in df.index:
                        safe_cells = []
                        for i in range(len(df.columns)):
                            cell_value = self.escape_latex(str(df.iloc[idx, i]))
                            if i == 0:
                                cell_value = r'\seqsplit{' + cell_value + r'}'
                            safe_cells.append(cell_value)
                        print(' & '.join(safe_cells) + r' \\', file=w)
                    print(r'\end{longtable}', file=w)
                    print(r'}', file=w)
                else:
                    print(r'\vspace{2mm}', file=w)
                    print(r'\begin{center}', file=w)
                    print(r'\begin{tabular}{|p{24cm}|}', file=w)
                    print(r'\noalign{\global\arrayrulewidth=0.75mm}', file=w)
                    print(r'\arrayrulecolor{midnightblue}', file=w)
                    print(r'\hline', file=w)
                    print(r'\vspace{1mm}', file=w)
                    print(r'{\Large No Antimicrobial Resistance Genes found.}', file=w)
                    print(r'\vspace{2mm}', file=w)
                    print(r'\hline', file=w)
                    print(r'\end{tabular}', file=w)
                    print(r'\end{center}', file=w)

            if abricate_report:
                _abricate_table(ab_resfinder_file, 'resfinder.png')
                print(r'\vspace{8mm}', file=w)
                print(r'\thispagestyle{fancy}', file=w)
                _abricate_table(ab_ncbi_file, 'ncbi.png')
                print(r'\vspace{2mm}', file=w)
                print(r'\end{landscape}', file=w)
                print(r'\fancyhf{}', file=w)
                print(r'\fancyhead[L]{\includegraphics[scale=0.15]{%s}}' % LOGO, file=w)
                print(r'\fancyhead[R]{\textbf{Isolate ID:}{%s}}' % safe_sample_name, file=w)
                print(r'\lfoot{\today}', file=w)
                print(r'\cfoot{\thepage\ of \zref[abspage]{LastPage}}', file=w)
                print(r'\rfoot{}', file=w)
                print(r'\newpage', file=w)

            # ── AMRFinder table ────────────────────────────────────────────
            print(r'\newpage', file=w)
            print(r'\fancyhf{}', file=w)
            print(r'\fancyhead[L]{\includegraphics[scale=0.15]{%s}}' % LOGO, file=w)
            print(r'\fancyhead[R]{\textbf{Isolate ID:}{%s}}' % safe_sample_name, file=w)
            print(r'\cfoot{\thepage\ of \zref[abspage]{LastPage}}', file=w)
            print(r'\begin{landscape}', file=w)
            print(r'\setlength{\headheight}{25pt}', file=w)
            print(r'\setlength{\footskip}{10pt}', file=w)
            print(r'\thispagestyle{fancy}', file=w)
            print(r'\vspace{3mm}', file=w)

            amrfinder_lines = 0
            if amrfinder_file and os.path.exists(amrfinder_file):
                amrfinder_lines = sum(1 for _ in open(amrfinder_file))

            if amrfinder_lines > 1:
                print(r'\noindent', file=w)
                print(r'\includegraphics[scale=0.25]{amrfinder.png}', file=w)
                print(r'\vspace{2mm}', file=w)
                print(r'\noindent', file=w)
                print(r'{\footnotesize', file=w)
                print(r'\renewcommand{\arraystretch}{0.9}', file=w)
                # Use longtable instead of tabular for automatic page breaks
                print(r'\begin{longtable}{p{2.5cm}|c|c|c|p{3cm}|c|c|c|c|c|c}', file=w)
                # First page header
                print(r'\textcolor{midnightblue}{\bf Contig id}&'
                      r'\textcolor{midnightblue}{\bf Start}&'
                      r'\textcolor{midnightblue}{\bf Stop}&'
                      r'\textcolor{midnightblue}{\bf Gene symbol}&'
                      r'\textcolor{midnightblue}{\bf Sequence name}&'
                      r'\textcolor{midnightblue}{\bf Scope}&'
                      r'\textcolor{midnightblue}{\bf Element Subtype}&'
                      r'\textcolor{midnightblue}{\bf Class}&'
                      r'\textcolor{midnightblue}{\bf Subclass}&'
                      r'\textcolor{midnightblue}{\bf \% Coverage}&'
                      r'\textcolor{midnightblue}{\bf \% Identity} \\', file=w)
                print(r'\hline', file=w)
                print(r'\endfirsthead', file=w)
                # Continuation header for subsequent pages
                print(r'\multicolumn{11}{c}{\textit{\small (Continued from previous page)}} \\', file=w)
                print(r'\hline', file=w)
                print(r'\textcolor{midnightblue}{\bf Contig id}&'
                      r'\textcolor{midnightblue}{\bf Start}&'
                      r'\textcolor{midnightblue}{\bf Stop}&'
                      r'\textcolor{midnightblue}{\bf Gene symbol}&'
                      r'\textcolor{midnightblue}{\bf Sequence name}&'
                      r'\textcolor{midnightblue}{\bf Scope}&'
                      r'\textcolor{midnightblue}{\bf Element Subtype}&'
                      r'\textcolor{midnightblue}{\bf Class}&'
                      r'\textcolor{midnightblue}{\bf Subclass}&'
                      r'\textcolor{midnightblue}{\bf \% Coverage}&'
                      r'\textcolor{midnightblue}{\bf \% Identity} \\', file=w)
                print(r'\hline', file=w)
                print(r'\endhead', file=w)
                # Footer for continuation pages
                print(r'\hline', file=w)
                print(r'\multicolumn{11}{r}{\textit{\small (Continued on next page)}} \\', file=w)
                print(r'\endfoot', file=w)
                # Footer for last page
                print(r'\noalign{\global\arrayrulewidth=0.75mm}', file=w)
                print(r'\arrayrulecolor{midnightblue}\hline', file=w)
                print(r'\endlastfoot', file=w)
                # Data rows
                df = pd.read_csv(amrfinder_file, sep='\t', header=None, skiprows=1)
                for idx in df.index:
                    s = df.iloc[idx]
                    if (s.iloc[9] in ("AMR", "POINT")) and s.iloc[10] != "EFFLUX":
                        cols = [s.iloc[i] for i in [1, 2, 3, 5, 6, 7, 9, 10, 11, 15, 16]]
                        safe_cells = []
                        for i, c in enumerate(cols):
                            cell_value = self.escape_latex(str(c))
                            # Wrap Contig id (index 0) and Sequence name (index 4) with \seqsplit
                            if i in [0, 4]:
                                cell_value = r'\seqsplit{' + cell_value + r'}'
                            safe_cells.append(cell_value)
                        print(' & '.join(safe_cells) + r' \\', file=w)
                print(r'\end{longtable}', file=w)
                print(r'}', file=w)  # ← Close footnotesize
            else:
                print(r'\noindent', file=w)
                print(r'\includegraphics[scale=0.35]{amrfinder1.png}', file=w)
                print(r'\vspace{2mm}', file=w)
                print(r'\begin{center}', file=w)  # ← ADD: Explicit centering environment
                print(r'\begin{tabular}{|p{24cm}|}', file=w)  # ← CHANGE: Add borders, reduce width
                print(r'\noalign{\global\arrayrulewidth=0.75mm}', file=w)  # ← MOVE: Set width before first hline
                print(r'\arrayrulecolor{midnightblue}', file=w)  # ← MOVE: Set color before first hline
                print(r'\hline', file=w)
                print(r'\vspace{1mm}', file=w)
                print(r'{\Large No Antimicrobial Resistance Genes found.}', file=w)  # ← REMOVE \centerline
                print(r'\vspace{2mm}', file=w)
                print(r'\hline', file=w)
                print(r'\end{tabular}', file=w)
                print(r'\end{center}', file=w)  # ← ADD: Close centering environment

            print(r'\end{landscape}', file=w)
            print(r'\setlength{\headheight}{40pt}', file=w)
            print(r'\setlength{\footskip}{15pt}', file=w)
            print(r'\fancyhf{}', file=w)
            print(r'\fancyhead[L]{\includegraphics[scale=0.15]{%s}}' % LOGO, file=w)
            print(r'\fancyhead[R]{\textbf{Isolate ID:}{%s}}' % safe_sample_name, file=w)
            print(r'\lfoot{\today}', file=w)
            print(r'\cfoot{\thepage\ of \zref[abspage]{LastPage}}', file=w)
            print(r'\rfoot{}', file=w)

            # ── ABRicate definitions page ──────────────────────────────────
            if abricate_report:
                print(r'\newpage', file=w)
                print(r'\begin{figure}', file=w)
                print(r'\centering', file=w)
                print(r'\vspace{2mm}', file=w)
                print(r'\includegraphics[scale=0.25]{abricate_doc.png}', file=w)
                print(r'\end{figure}', file=w)
                print(r'\begin{center}', file=w)
                print(r'{\large Database Versions} \\', file=w)
                print(r'\vspace{3mm}', file=w)
                print(r'\begin{tabular}{l|c|c}', file=w)
                print(r'\hline', file=w)
                print(r'\textbf{Database} & \textbf{Sequences} & \textbf{Date Updated} \\', file=w)
                print(r'\hline', file=w)
                print(f' NCBI & {int(abricate_ncbi_seq_number):,} & {safe_ncbi_date} \\\\', file=w)
                print(f' ResFinder & {int(abricate_res_seq_number):,} & {safe_res_date} \\\\', file=w)
                print(r'\hline', file=w)
                print(r'\end{tabular}', file=w)
                print(r'\end{center}', file=w)
                print(r'\begin{center}', file=w)
                print(r'{\large Parameters}\\', file=w)
                print(r'\vspace{3mm}', file=w)
                print(r'\begin{tabular}{l|c}', file=w)
                print(r'\hline', file=w)
                print(f'Minimum Coverage & {abricate_mincov}\\% \\\\', file=w)
                print(f'Minimum Identity & {abricate_minid}\\% \\\\', file=w)
                print(r'\hline', file=w)
                print(r'\end{tabular}', file=w)
                print(r'\end{center}', file=w)
                print(r'\vspace{-2mm}', file=w)
                print(r'\noindent', file=w)
                print(r'\textbf{\large Summary Output}\\', file=w)
                print(r'\vspace{-3mm}', file=w)
                print(r'\hrule', file=w)
                print(r'\vspace{1mm}', file=w)
                print(r'This abbreviated section appears in the first tab of each ABRicate Excel workbook. \\', file=w)
                print(r'\begin{center}', file=w)
                print(r'\begin{minipage}{7in}', file=w)
                print(r'\begin{basedescript}{\desclabelstyle{\pushlabel}}', file=w)
                print(r'\item[NUM\_FOUND-] The number of distinct genes identified in the analysis of the sample.  '
                      r'This is NOT the total number of hits if gene duplicates are identified.', file=w)
                print(r'\item[Gene List and Percent Coverage-] A gene list with percent coverage is also given by '
                      r'isolate.  If multiple identifications of the same gene were made in a single isolate then '
                      r'the percent coverage information is given in a colon separated list in the order the '
                      r'identifications are listed in the full analysis.', file=w)
                print(r'\end{basedescript}', file=w)
                print(r'\end{minipage}', file=w)
                print(r'\end{center}', file=w)
                print(r'\vspace{4mm}', file=w)
                print(r'\noindent', file=w)
                print(r'\textbf{\large Full Analysis Output}\\', file=w)
                print(r'\vspace{-3mm}', file=w)
                print(r'\hrule', file=w)
                print(r'\vspace{1mm}', file=w)
                print(r'Details of the full analysis are located on subsequent tabs of the ABRicate workbooks.\\', file=w)
                for item, desc in [
                    (r'SEQUENCE-',     r'The analysis is ran on scaffolds (also known as nodes or contigs) output by the assembler. Node length and coverage are included with the name.'),
                    (r'FILE-',         r'This field refers to the original file output by the analysis prior to the workbook creation.'),
                    (r'START-',        r'This is the position in the scaffold where the alignment with the given reference gene starts.'),
                    (r'END-',          r'This is the position of the scaffold where the alignment with the given reference gene ends.'),
                    (r'GENE-',         r'The reference gene to which the alignment with the scaffold is performed.'),
                    (r'COVERAGE-',     r'The positions (given as a range) of the reference gene that align with the scaffold divided by the length of the gene sequence. (aligned length/gene length)'),
                    (r'COVERAGE\_MAP-',r'This gives an overview of the alignment relative to the reference gene. ....-no alignment ====-alignment /-gap in alignment.'),
                    (r'GAPS-',         r'The number of gaps in the alignment of the reference gene to the scaffold.'),
                    (r'\%COV-',        r'The percent of the reference gene covered by the alignment with the scaffold. (Percent coverage)'),
                    (r'\%IDENT-',      r'The percent identity of the reference gene with the scaffold. (Percent Identity)'),
                    (r'DATABASE-',     r'The gene database that was used for comparison to the scaffolds.'),
                    (r'ACCESSION-',    r'The accession number of the reference gene in the database used for the comparison.'),
                ]:
                    print(r'\begin{center}\begin{minipage}{7in}\begin{basedescript}{\desclabelstyle{\pushlabel}}', file=w)
                    print(r'\item[%s] %s' % (item, desc), file=w)
                    print(r'\end{basedescript}\end{minipage}\end{center}', file=w)

            # AMRFinder definitions page
            print(r'\newpage', file=w)
            print(r'\begin{figure}', file=w)
            print(r'\centering', file=w)
            print(r'\vspace{2mm}', file=w)
            print(r'\includegraphics[scale=0.25]{amrfind_doc.png}', file=w)
            print(r'Definitions were taken from the AMRFinder documentation.\\', file=w)
            print(r'\noindent\makebox[\linewidth]{\rule{17cm}{0.4pt}}', file=w)
            print(r'\end{figure}', file=w)
            print(r'\begin{center}', file=w)
            print(r'\begin{minipage}{6.5in}', file=w)
            print(r'\begin{basedescript}{\desclabelstyle{\pushlabel}}', file=w)
            for item, desc in [
                (r'Target Identifier-', r'This is from the FASTA defline for the DNA sequence'),
                (r'Contig id-',         r'Contig name'),
                (r'Start-',             r'1-based coordinate of first nucleotide coding from protein in DNA sequence on contig'),
                (r'Stop-',              r'1-based coordinate of last nucleotide coding for protein in DNA sequence on contig'),
                (r'Gene symbol-',       r'Gene or gene-family symbol for protein hit'),
                (r'Protein name-',      r'Full-text name for the protein'),
                (r'Method-',            r'Type of hit found by AMRFinder one of five options'),
                (r'ALLELE-',            r'100\% sequence match over 100\% of length to a protein named at the allele level in the AMRFinder database'),
                (r'EXACT-',             r'100\% sequence match over 100\% of length to a protein in the database that is not a named allele'),
                (r'BLAST-',             r'BLAST alignment is \textgreater 90\% of length and \textgreater 90\% identity to a protein in the AMRFinder database'),
                (r'PARTIAL-',           r'BLAST alignment is \textgreater 50\% of length, but \textless 90\% of length and \textgreater 90\% identity'),
                (r'HMM-',              r'HMM was hit above the cutoff, but there was not a BLAST hit that met standards for BLAST or PARTIAL'),
                (r'Target length-',     r'The length of the query protein. The length of the BLAST hit for translated-DNA searches'),
                (r'Reference protein length-', r'The length of the AMR Protein in the database (NA if HMM-only hit)'),
                (r'Scope-',             r'The AMRFinderPlus database is split into core AMR proteins that are expected to have an effect on resistance and plus proteins of interest added with less stringent inclusion criteria.'),
                (r'Element subtype-',   r'Further elaboration of functional category into (ANTIGEN, BIOCIDE, HEAT, METAL, PORIN) if more specific category is available, otherwise the element is repeated'),
                (r'\% Coverage of reference protein-', r'\% covered by blast hit (NA if HMM-only hit)'),
                (r'\% Identity to reference protein-', r'\% amino-acid identity to reference protein (NA if HMM-only hit)'),
            ]:
                print(r'\item[%s] %s' % (item, desc), file=w)
            print(r'\end{basedescript}', file=w)
            print(r'\end{minipage}', file=w)
            print(r'\end{center}', file=w)
            print(r'\begin{center}', file=w)
            print(r'\begin{minipage}{6.5in}', file=w)
            print(r'\begin{basedescript}{\desclabelstyle{\pushlabel}}', file=w)
            for item, desc in [
                (r'Alignment length-',           r'Length of BLAST alignment in amino acids (NA if HMM-only hit)'),
                (r'Accession of closest protein-',r'RefSeq accession for protein hit by BLAST (NA if HMM-only hit)'),
                (r'Name of closest protein-',    r'Full name assigned to the AMRFinder database protein (NA if HMM-only hit)'),
                (r'HMM id-',                     r'Accession for the HMM'),
                (r'HMM description-',            r'The family name associated with the HMM'),
            ]:
                print(r'\item[%s] %s' % (item, desc), file=w)
            print(r'\end{basedescript}', file=w)
            print(r'\end{minipage}', file=w)
            print(r'\end{center}', file=w)

            # Software Versions page
            if software_versions:
                print(r'\newpage', file=w)
                print(r'\section*{Software Versions}', file=w)
                print(r'\begin{small}', file=w)

                current_env = None
                env_data = {}
                for line in software_versions.strip().split('\n'):
                    line = line.strip()
                    if line.startswith('---') and line.endswith('---'):
                        current_env = line.strip('-').strip()
                        env_data.setdefault(current_env, [])
                        continue
                    if (not line or line.startswith('#') or line.startswith('==')
                            or line.startswith('Name')
                            or line == 'AMR Pipeline - Software Versions Report'
                            or current_env is None):
                        continue
                    if ':' in line:
                        sw, ver = line.split(':', 1)
                        env_data[current_env].append((sw.strip().strip('"'), ver.strip().strip('"')))
                    elif ' ' in line:
                        parts = line.split(None, 1)
                        if len(parts) == 2:
                            env_data[current_env].append((parts[0].strip(), parts[1].strip()))

                for env_name, packages in env_data.items():
                    if not packages:
                        continue
                    print(r'\subsection*{' + self._e(env_name) + r'}', file=w)
                    print(r'\begin{longtable}{|p{4cm}|p{10cm}|}', file=w)
                    print(r'\hline', file=w)
                    print(r'\textbf{Software} & \textbf{Version} \\', file=w)
                    print(r'\hline', file=w)
                    print(r'\endfirsthead', file=w)
                    print(r'\hline', file=w)
                    print(r'\textbf{Software} & \textbf{Version} \\', file=w)
                    print(r'\hline', file=w)
                    print(r'\endhead', file=w)
                    print(r'\hline', file=w)
                    print(r'\endfoot', file=w)
                    for sw, ver in packages:
                        print(f'{self._e(sw)} & {self._e(ver)} \\\\', file=w)
                        print(r'\hline', file=w)
                    print(r'\end{longtable}', file=w)
                    print(r'\vspace{0.3cm}', file=w)

                print(r'\end{small}', file=w)

            # Close document 
            print(r'\end{document}', file=w)
            w.close()

            # Run pdflatex
            print("--- Running pdflatex pass 1 ---", flush=True)
            ret1 = os.system(f'pdflatex -interaction=nonstopmode {tex_document}')
            exit_code1 = os.WEXITSTATUS(ret1) if os.WIFEXITED(ret1) else ret1
            print(f"--- pdflatex pass 1 exit code: {exit_code1} ---", flush=True)

            print("--- Running pdflatex pass 2 ---", flush=True)
            ret2 = os.system(f'pdflatex -interaction=nonstopmode {tex_document}')
            exit_code2 = os.WEXITSTATUS(ret2) if os.WIFEXITED(ret2) else ret2
            print(f"--- pdflatex pass 2 exit code: {exit_code2} ---", flush=True)

            pdf_file = tex_document.replace('.tex', '.pdf')
            if not os.path.exists(pdf_file):
                log_file = tex_document.replace('.tex', '.log')
                if os.path.exists(log_file):
                    print("=== PDFLATEX LOG (last 100 lines) ===", flush=True)
                    with open(log_file, 'r', errors='replace') as lf:
                        lines = lf.readlines()
                    print(''.join(lines[-100:]), flush=True)
                raise RuntimeError(
                    f"pdflatex failed to create PDF (pass1={exit_code1}, pass2={exit_code2}). "
                    f"Check log above for details."
                )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='PROG',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''
    ---------------------------------------------------------
    AMR LaTeX Report Generator
    ---------------------------------------------------------
    '''),
        epilog='---------------------------------------------------------'
    )
    parser.add_argument('-v', '--version', action='version',
                        version=f'{os.path.basename(__file__)}: version {__version__}')
    args = parser.parse_args()

# Created March 2020 by Tod Stuber, modified by RW in 2026