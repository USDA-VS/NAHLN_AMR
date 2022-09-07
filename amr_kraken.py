#!/usr/bin/env python

__version__ = "0.0.1"

import os
import sys
import shutil
import glob
import argparse
import textwrap
import pandas as pd
import multiprocessing
multiprocessing.set_start_method('spawn', True)
import matplotlib.pyplot as plt
plt.style.use('seaborn-colorblind')
from matplotlib import cm
cmap = cm.get_cmap('coolwarm')
import traceback
import urllib3
import xmltodict
from Bio import Entrez
import matplotlib.pyplot as plt
from numpy.random import rand
from numpy import arange

from amr_file_setup import Setup
from amr_file_setup import bcolors
from amr_file_setup import Banner
from amr_file_setup import Latex_Report
from amr_file_setup import Excel_Stats

class Kraken_Identification(Setup):
    ''' 
    Assemble reads using Spades assembler.
    Paired or single reads
    '''

    def __init__(self, FASTA=None, FASTQ_R1=None, FASTQ_R2=None, directory='kraken', influenza=None, debug=None):

        Setup.__init__(self, FASTA=FASTA, FASTQ_R1=FASTQ_R1, FASTQ_R2=FASTQ_R2, debug=debug)
        self.directory = directory
        self.influenza = influenza
        if influenza:
            self.kraken_db = "/project/bioinformatic_databases/databases/kraken/flu_jhu"
        
    def run(self,):
        self.print_run_time('Kraken')
        kraken_db = self.kraken_db
        cpus = self.cpus
        sample_name = self.sample_name
        FASTQ_list = self.FASTQ_list
        FASTA =  self.FASTA
        cwd = self.cwd
        if self.influenza: #Need to run JHU database using Kraken1
            if len(FASTQ_list) == 2:
                os.system(f'kraken --db {self.kraken_db} --paired {FASTQ_list[0]} {FASTQ_list[1]} > {sample_name}_outputkraken.txt')
            elif len(FASTQ_list) == 1:
                os.system(f'kraken --db {self.kraken_db} {FASTQ_list[0]} > {sample_name}_outputkraken.txt')
            else:
                os.system(f'kraken --db {self.kraken_db} {FASTA} > {sample_name}_outputkraken.txt')
            os.system(f'kraken-report --db {self.kraken_db} {sample_name}_outputkraken.txt > {sample_name}_reportkraken.txt')
            os.system(f'dvl_krakenreport2krona.sh -i {sample_name}_reportkraken.txt -k {self.kraken_db} -t {sample_name}-jhu-output.txt -o {sample_name}-jhu-Krona_id_graphic.html')
            if not os.path.exists(self.directory):
                os.mkdir(self.directory)
            shutil.move(f'{sample_name}-jhu-Krona_id_graphic.html', self.directory)
            os.remove(f'{sample_name}-jhu-output.txt')
        else:
            if len(FASTQ_list) == 2:
                os.system(f'kraken2 --db {kraken_db} --threads {cpus} --paired {FASTQ_list[0]} {FASTQ_list[1]} --output {sample_name}_outputkraken.txt --report {sample_name}_reportkraken.txt')
            elif len(FASTQ_list) == 1:
                os.system(f'kraken2 --db {kraken_db} --threads {cpus} {FASTQ_list[0]} --output {sample_name}_outputkraken.txt --report {sample_name}_reportkraken.txt')
            else:
                os.system(f'kraken2 --db {kraken_db} --threads {cpus} {FASTA} --output {sample_name}_outputkraken.txt --report {sample_name}_reportkraken.txt')

        if os.path.exists(f'{cwd}/{sample_name}_outputkraken.txt'):
            output = f'{cwd}/{sample_name}_outputkraken.txt'
        else:
            print(f'\n### Error: Kraken report did not complete')
            sys.exit(0)
        if os.path.exists(f'{cwd}/{sample_name}_reportkraken.txt'):
            report = f'{cwd}/{sample_name}_reportkraken.txt'
        else:
            print(f'\n### Error: Kraken report did not complete')
            sys.exit(0)

        if self.directory:
            if not os.path.exists(self.directory):
                os.mkdir(self.directory)
            shutil.move(report, self.directory)
            shutil.move(output, self.directory)
            self.report = f'{cwd}/{self.directory}/{sample_name}_reportkraken.txt'
            self.output = f'{cwd}/{self.directory}/{sample_name}_outputkraken.txt'
            log_file = open("kraken_log.txt", "a")
            try:
                log_file.write(f'DB used: {os.readlink(self.kraken_db)}')
            except OSError:
                log_file.write(f'DB used: {self.kraken_db}')
            log_file.close()
            shutil.move("kraken_log.txt", self.directory)

    def krona_make_graph(self, report):
        '''
        kreport2krona.py -r 20-037580-001s_reportkraken.txt -o sample.krona 
        kreport2krona.py --intermediate-ranks -r 20-037580-001s_reportkraken.txt -o sample.krona
        '''
        os.system(f'kreport2krona.py --intermediate-ranks -r {report} -o {self.sample_name}.krona')
        os.system(f'ktImportText {self.sample_name}.krona -o {self.sample_name}_{self.date_stamp}_krona.html')
        os.remove(f'{self.sample_name}.krona')

        if os.path.exists(f'{self.cwd}/{self.sample_name}_{self.date_stamp}_krona.html'):
            self.krona_html = f'{self.cwd}/{self.sample_name}_{self.date_stamp}_krona.html'
        else:
            print(f'\n### Error: Krona HTML did not complete')
            sys.exit(0)
        if self.directory:
            shutil.move(f'{self.sample_name}_{self.date_stamp}_krona.html', self.directory)
            self.krona_html = f'{self.cwd}/{self.directory}/{self.sample_name}_{self.date_stamp}_krona.html'
            
    def bracken(self, report, output):
        os.system(f'bracken -d {self.kraken_db} -i {report} -o {self.sample_name}-bracken.txt -r 250')
        df = pd.read_csv(f'{self.sample_name}-bracken.txt', sep='\t')
        df.to_excel(f'{self.sample_name}-bracken.xlsx', index=False)
        os.remove(f'{self.sample_name}-bracken.txt')
        self.bracken_excel = f'{os.getcwd()}/{self.sample_name}-bracken.xlsx'
        if self.directory:
            shutil.move(f'{self.sample_name}-bracken.xlsx', self.directory)
            self.bracken_excel = f'{os.getcwd()}/{self.directory}/{self.sample_name}-bracken.xlsx'

class Bracken_Pie_Charts:

    def __init__(self, FASTA=False):
        self.FASTA = FASTA

    def run(self, bracken_excel,):
        df = pd.read_excel(bracken_excel)
        df = df[df['fraction_total_reads'] > 0.01 ]
        df2 = pd.DataFrame([['unclassified', 0, 'na', 0, 0, 0, 1 - df['fraction_total_reads'].sum()]], columns=['name','taxonomy_id','taxonomy_lvl', 'kraken_assigned_reads', 'added_reads', 'new_est_reads', 'fraction_total_reads'])
        df3 = pd.concat([df, df2])
        df3 = df3.set_index('name')
        if self.FASTA:
            plot = df3.plot.pie(y='fraction_total_reads', title='Identification of Assembled Scaffolds', figsize=(9, 5), cmap=cmap, labeldistance=None, legend=True, autopct='%1.1f%%')
        else:  #default FASTQ
            plot = df3.plot.pie(y='fraction_total_reads', title='FASTQ Read Identification', figsize=(9, 5), cmap=cmap, labeldistance=None, legend=True, autopct='%1.1f%%')
        plot.axis('off')
        plot.legend(bbox_to_anchor=(0.9, 0.9))
        plot.yaxis.label.set_visible(False)
        plot.get_figure().savefig(f'{os.getcwd()}/bracken_pie.png', format='png', bbox_inches='tight')
        self.pie_chart = f'{os.getcwd()}/bracken_pie.png'

    def latex(self, build_latex=False):
            tex = build_latex
            print(r'\begin{table}[ht!]', file=tex)
            print(r'\begin{adjustbox}{width=1\textwidth}', file=tex)
            if self.FASTA:
                bracken_pie_banner = Banner("FASTA Identifications")
            else:
                bracken_pie_banner = Banner("FASTQ Identifications")
            print(r'\begin{center}', file=tex)
            print(r'\includegraphics[scale=1]{' + bracken_pie_banner.banner + r'}', file=tex)
            print(r'\end{center}', file=tex)
            print(r'\end{adjustbox}', file=tex)
            print(r'\begin{center}', file=tex)
            print(r'\includegraphics[scale=0.8]{' + self.pie_chart + r'}', file=tex)
            print(r'\end{center}', file=tex)
            print(r'Identified using: \href{https://ccb.jhu.edu/software/kraken2/}{Kraken} and \href{https://ccb.jhu.edu/software/bracken/}{Bracken}', file=tex)
            print(r'\end{table}', file=tex)

class Bracken_Size_Graph:

    def __init__(self, sample_name, bracken_excel, total_contig_lengths, predicted_size=None, target_name=None):
        self.sample_name = sample_name
        self.bracken_excel = bracken_excel
        self.total_contig_lengths = total_contig_lengths
        self.predicted_size = predicted_size
        self.target_name = target_name

    def run(self, ):
        predicted_size = self.predicted_size
        def getxml(taxid):
            url=f'https://api.ncbi.nlm.nih.gov/genome/v0/expected_genome_size?species_taxid={taxid}'
            http = urllib3.PoolManager()
            response = http.request('GET', url)
            try:
                data = xmltodict.parse(response.data)
            except:
                print("Failed to parse xml from response (%s)" % traceback.format_exc())
            return data

        df = pd.read_excel(self.bracken_excel)
        df = df[df['fraction_total_reads'] > 0.01 ]
        top_bracken_id = df['name'].to_list()[0]
        Entrez.email = "tod.p.stuber@usda.gov"     # Always tell NCBI who you are
        handle = Entrez.esearch(db="taxonomy", term=top_bracken_id)
        record = Entrez.read(handle)
        taxid = record['IdList']

        data = getxml(taxid[0])
        error = True
        count = 0

        if predicted_size is None:
            #Iterate Bracken IDs until expected_ungapped_length found
            self.true_top_bracken_id = df['name'].to_list()[0]
            while error == True:
                try:
                    bracken_predicted_size = data['genome_size_response']['expected_ungapped_length']
                    error = False
                except KeyError:
                    count += 1
                    try:
                        top_bracken_id = df['name'].to_list()[count]
                    except IndexError:
                        #if after iterating through Bracken ID list no suitable taxid provides genome length
                        bracken_predicted_size = 5000000
                        error = False
                    handle = Entrez.esearch(db="taxonomy", term=top_bracken_id)
                    record = Entrez.read(handle)
                    try:
                        taxid = record['IdList']
                        data = getxml(taxid[0])
                    except IndexError: #pass if esearch comes up empty, then try the next best Bracken ID
                        pass
            unknown_sample = self.sample_name
            unknown_assemble_size = int(str(self.total_contig_lengths).replace(',', ''))
            self.predicted_sample = top_bracken_id
            predicted_size = int(bracken_predicted_size)
        elif self.target_name:
            self.predicted_sample = self.target_name
            unknown_sample = self.sample_name
            unknown_assemble_size = int(str(self.total_contig_lengths).replace(',', ''))
        else:
            self.predicted_sample = 'predicted target'
            unknown_sample = self.sample_name
            unknown_assemble_size = int(str(self.total_contig_lengths).replace(',', ''))
        ######
        if predicted_size == 5000000:
            count = [f'assembled: {unknown_assemble_size:,} nt', f'expected: appx guess ~{predicted_size:,} nt']
            self.predicted_sample = 'no lookup found'
        else:
            count = [f'assembled: {unknown_assemble_size:,} nt', f'expected: {predicted_size:,} nt']
        platform = [unknown_sample, self.predicted_sample]

        def autolabel(rects):
            biggest_xloc1 = [] # the biggest of the 2 bars needs to be used to align the text as needed if predicted and sample are very different in size
            for each, rect in enumerate(rects):
                ## iterate each index (each)
                width = int(rect.get_width())
                height = rect.get_height()
                yloc1=rect.get_y() + height /2.0
                yloc2=rect.get_y() + height /2.0
                if (width <= 5):
                    # Shift the text to the right side of the right edge
                    xloc1 = width + 1
                    yloc2=yloc2+0.3
                    # Black against white background
                    clr = 'black'
                    align = 'left'
                else:
                    # Shift the text to the left side of the right edge
                    xloc1 = 0.98*width
                    # White on blue
                    clr = 'black'
                    align = 'right'
                yloc1=rect.get_y() + height /2.0
                biggest_xloc1.append(xloc1)
            max_xloc1 = max(biggest_xloc1) #will be just to numbers in list
            min_xloc1 = min(biggest_xloc1)
            ratio = min_xloc1/max_xloc1
            if ratio > 0.75: # if ratio > then put length value in bar line, do this so sample name and length don't squish together
                xloc1 = min_xloc1
            else: # print lenght outside bar line
                xloc1 = max_xloc1

            # Make horizontal bar graphs with largest xloc1
            for each, rect in enumerate(rects):
                ## iterate each index (each)
                width = int(rect.get_width())
                height = rect.get_height()
                yloc1=rect.get_y() + height /2.0
                yloc2=rect.get_y() + height /2.0
                if (width <= 5):
                    # Shift the text to the right side of the right edge
                    yloc2=yloc2+0.3
                    # Black against white background
                    clr = 'black'
                    align = 'left'
                else:
                    # Shift the text to the left side of the right edge
                    # White on blue
                    clr = 'black'
                    align = 'right'
                yloc1=rect.get_y() + height /2.0
                ax.text(xloc1, yloc1, count[each] ,horizontalalignment=align,
                                verticalalignment='center',color=clr,
                                clip_on=True)
                ax.text(1, yloc2, platform[each],horizontalalignment='left',
                                verticalalignment='center',color=clr,
                                clip_on=True)
        # Make horizontal bars with labels
        fig = plt.figure(figsize=(10,.6))
        ax = fig.add_subplot(111)
        rects = ax.barh([1, 1.1], [unknown_assemble_size, predicted_size], align='center', color='#5abfa1', edgecolor='black', height=0.1)
        autolabel(rects)
        ax.grid(False)
        plt.axis('off')
        plt.savefig(f'{os.getcwd()}/size_comparison.png')
        self.size_comparison = f'{os.getcwd()}/size_comparison.png'

    def latex(self, build_latex=False):
        tex = build_latex
        print(r'\begin{table}[ht!]', file=tex)
        print(r'\begin{adjustbox}{width=1\textwidth}', file=tex)
        size_comparison_banner = Banner("Genome Size Comparisons")
        print(r'\begin{center}', file=tex)
        print('\includegraphics[scale=1]{' + size_comparison_banner.banner + '}', file=tex)
        print(r'\end{center}', file=tex)
        print(r'\end{adjustbox}', file=tex)
        print(r'\begin{adjustbox}{width=1\textwidth}', file=tex)
        print('\includegraphics[scale=1.2]{size_comparison.png}', file=tex)
        print(r'\end{adjustbox}', file=tex)
        print(r'\end{table}', file=tex)


if __name__ == "__main__": # execute if directly access by the interpreter

    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

        ---------------------------------------------------------
        Provide either a single FASTA file, single FASTQ or Paired files.
        Usage:
            kraken_run.py -r1 *_R1*fastg.gz
            kraken_run.py -r1 *_R1*fastg.gz -r2 *_R2*fastq.gz -d
            kraken_run.py -f *fasta

        '''), epilog='''---------------------------------------------------------''')

    parser.add_argument('-f', '--FASTA', action='store', dest='FASTA', required=False, help='Provide FASTA file')
    parser.add_argument('-r1', '--FASTQ_R1', action='store', dest='FASTQ_R1', required=False, help='Provide R1 FASTQ gz file, or single read')
    parser.add_argument('-r2', '--FASTQ_R2', action='store', dest='FASTQ_R2', required=False, default=None, help='Provide R2 FASTQ gz file')
    parser.add_argument('-y', '--directory', action='store', dest='directory', required=False, default="kraken", help='Put output to directory')
    parser.add_argument('-i', '--influenza', action='store_true', dest='influenza', default=False, help='Use JHU influenza specific database')
    parser.add_argument('-d', '--debug', action='store_true', dest='debug', default=False, help='keep temp file')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.basename(__file__)}: version {__version__}')
    args = parser.parse_args()

    print(f'\n{os.path.basename(__file__)} SET ARGUMENTS:')
    print(args)

    kraken = Kraken_Identification(FASTA=args.FASTA, FASTQ_R1=args.FASTQ_R1, FASTQ_R2=args.FASTQ_R2, directory=args.directory, influenza=args.influenza, debug=args.debug)
    kraken.run()
    if args.influenza is False:
        krona_html = kraken.krona_make_graph(kraken.report)
        kraken.bracken(kraken.report, kraken.output)

    print('done')
# Created February 2021 by Tod Stuber