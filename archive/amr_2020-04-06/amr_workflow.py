#!/usr/bin/env python

import os
import sys
import shutil
import glob
import json
import time
import traceback
import re
import math
import numpy as np
import pandas as pd
import xlsxwriter
import argparse
import textwrap
import humanize
import smtplib
import multiprocessing
from multiprocessing import Pool
from prettytable import PrettyTable
from dask import delayed
from itertools import repeat as itertools_repeat
from collections import Iterable
from numpy import mean
from email.utils import formatdate
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.mime.base import MIMEBase
from email import encoders
from distutils.dir_util import copy_tree
from datetime import datetime
from concurrent import futures
from collections import OrderedDict
from collections import Counter
from collections import defaultdict
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio import SeqIO
import logging
import inspect
import svgwrite

from simple_spades import trim_reads
from simple_spades import simple_spades
from simple_kraken1 import simple_kraken1
from simple_kraken1 import make_krona

LOGO = os.environ['LOGO'] #grab environment variable

#logging.basicConfig(format='%(levelname)s: %(message)s', filemode='w', filename='debug.log', level=logging.DEBUG)
#logging.getLogger().addHandler(logging.StreamHandler()) #print to console

root_dir = str(os.getcwd())

def warning_log(ex, inspect_getframeinfo, *args):
    logging.warning(f'\nException occured, file: {inspect_getframeinfo.filename}\nfuction: {inspect.stack()[0][3]}, near line in script {inspect_getframeinfo.lineno} --> {type(ex).__name__, ex.args}\nAdditional args: {args}\n\n')
    np.save('parameter_dict.npy', parameter_dict)

def debug_log(ex, inspect_getframeinfo, *args):
    with open ('traceback.txt', 'w') as exception_trace:
        traceback.print_exc(file=exception_trace)
    traceback.print_exc()
    logging.debug(f'\nException occured, file: {inspect_getframeinfo.filename}\nfuction: {inspect.stack()[0][3]}, near line in script {inspect_getframeinfo.lineno} --> {type(ex).__name__, ex.args}\nAdditional args: {args}\n\n')
    np.save('parameter_dict.npy', parameter_dict)

def kickoff(parameter_dict):

    #check available
    try:
        parameter_dict['mlst_version'] = os.popen("mlst --version").readlines()[0]
        os.popen("abricate --version | sed 's/abricate//g'").readlines()[0]
        os.popen("amrfinder.pl -v| sed '1!d' |sed 's/version $Revision://g'|sed 's/\$//g'").readlines()[0]
    except IndexError:
        print("### check program dependencies")
        sys.exit(0)

    email_dict = {}
    email_dict["tod"] =  "tod.p.stuber@aphis.usda.gov"
    email_dict["jess"] =  "Jessica.A.Hicks@aphis.usda.gov"
    email_dict["all"] =  "kristina.lantz@aphis.usda.gov, tod.p.stuber@aphis.usda.gov, jessica.a.hicks@aphis.usda.gov, suelee.robbe-austerman@aphis.usda.gov, patrick.m.camp@aphis.usda.gov"
    parameter_dict['email_list'] = email_dict.get(args.email, None)

    cpu_count = int(multiprocessing.cpu_count())
    limited_cpu_count = int(cpu_count / 8)
    if limited_cpu_count == 0:
        limited_cpu_count = 1
    parameter_dict['cpu_count'] = cpu_count
    parameter_dict['limited_cpu_count'] = limited_cpu_count

    #check that there are only FASTQ in working directory
    if parameter_dict["test_fasta"]:
        print("### FASTA file for testing was provided.  The SPAdes assembly will not be done")
    else:
        all_file_types = glob.glob('*.*')
        all_file_types_count = len([x for x in all_file_types if not re.match(r'.*log', x) and not re.match(r'.*out', x)]) #don't include immediately made .log and slurm .out files
        fastq_check = len(glob.glob('*fastq.gz'))
        if all_file_types_count != fastq_check:
            print("\n### Something more than just FASTQs are in working directory\n\n")
            sys.exit(0)

    #Ensure files are paired
    R1 = glob.glob('*_R1*fastq.gz')
    R2 = glob.glob('*_R2*fastq.gz')
    if len(R1) == len(R2):
        run_loop(parameter_dict)
    else:
        print("\n### Check paired files.  Unpaired files seen by odd number of counted FASTQs\n\n")
        sys.exit(0)

def send_email(email_list, runtime, summary_file, st):
    text = "See attached:  "
    send_from = "tod.p.stuber@aphis.usda.gov"
    send_to = email_list
    msg = MIMEMultipart()
    msg['From'] = send_from
    msg['To'] = send_to
    msg['Date'] = formatdate(localtime=True)
    msg['Subject'] = f'AMR script complete, runtime: {runtime}'
    msg.attach(MIMEText(text))
    part = MIMEBase('application', "octet-stream")
    part.set_payload(open(summary_file, "rb").read())
    encoders.encode_base64(part)
    part.add_header('Content-Disposition', 'attachment; filename="amr_stat_summary_{}.xlsx"' .format(st))
    msg.attach(part)
    smtp = smtplib.SMTP('10.10.8.12')
    smtp.send_message(msg)
    smtp.quit()

def tab_to_excel(tab_file):
    excel_file = tab_file.replace('tab', 'xlsx')
    df = pd.read_csv(tab_file, sep='\t')
    df.to_excel(excel_file)

def run_loop(parameter_dict):

    root_dir = parameter_dict['root_dir']
    limited_cpu_count = parameter_dict['limited_cpu_count']

    startTime = datetime.now()
    ts = time.time()
    st = datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H-%M-%S')
    print("*** Start time: %s" % st)

    list_of_files = glob.glob('*gz')
    for afile in list_of_files:
        prefix_name = re.sub('_.*', '', afile)
        prefix_name = re.sub('\..*', '', prefix_name)
        print(prefix_name)
        if not os.path.exists(prefix_name):
            os.makedirs(prefix_name)
        shutil.move(afile, prefix_name)

    directory_list = [ddd for ddd in os.listdir(os.getcwd()) if os.path.isdir(ddd)]
    directory_size = {}
    for folder in directory_list: #run files by size, smallest to largest
        size = sum(os.path.getsize(os.path.join(dirpath,filename)) for dirpath, dirnames, filenames in os.walk(folder) for filename in filenames)
        directory_size[folder] = size
    directory_size = {k: v for k, v in sorted(directory_size.items(), key=lambda x: x[1])}
    directory_list = [*directory_size] #ordered list
    total_samples = len(directory_list)
    lower_count = 0
    upper_count = 1
    while lower_count < total_samples:
        upper_count = lower_count + limited_cpu_count
        run_list = directory_list[lower_count:upper_count] #create a run list
        for i in run_list:
            directory_list.remove(i)
        total_samples = len(directory_list)
        print(run_list)
        print("*** Iterating directories")
        frames = []
        if parameter_dict['debug_call']: #run just one sample at a time to debug
            for sample_name in run_list:
                print("*** DEBUGGING, SAMPLES RAN INDIVIDUALLY")
                parameter_dict = assembly_pipe(sample_name, parameter_dict)
                parameter_df = pd.DataFrame.from_dict(parameter_dict, orient='index').T #convert parameter_dict to df
                parameter_df_cut = parameter_df[['sample_name', 'File_R1', 'Size_R1', 'mean_quality1', 'q30pass1', 'count1', 'File_R2', 'Size_R2', 'mean_quality2', 'q30pass2', 'count2', 'Number of scaffolds', 'Total size of scaffolds', 'Longest scaffold', 'Number of scaffolds > 1K nt', 'N50 scaffold length', 'L50 scaffold count']]
                frames.append(parameter_df_cut) #frames to concatenate
                parameter_df_cut.to_excel(parameter_dict["sample_directory"] + "/" + parameter_dict['sample_name'] + "_stats.xlsx")
        else:  # run all in run_list in parallel
            print("*** SAMPLES RAN IN PARALLEL")
            # itertools allows additional arguments to pass
            # Need to speed test which why is faster
            with futures.ProcessPoolExecutor(max_workers=limited_cpu_count) as pool:
                for parameter_dict in pool.map(assembly_pipe, run_list, itertools_repeat(parameter_dict)):
                    parameter_df = pd.DataFrame.from_dict(parameter_dict, orient='index').T #convert parameter_dict to df
                    parameter_df_cut = parameter_df[['sample_name', 'File_R1', 'Size_R1', 'mean_quality1', 'q30pass1', 'count1', 'File_R2', 'Size_R2', 'mean_quality2', 'q30pass2', 'count2', 'Number of scaffolds', 'Total size of scaffolds', 'Longest scaffold', 'Number of scaffolds > 1K nt', 'N50 scaffold length', 'L50 scaffold count']]
                    frames.append(parameter_df_cut) #frames to concatenate
                    parameter_df_cut.to_excel(parameter_dict["sample_directory"] + "/" + parameter_dict['sample_name'] + "_stats.xlsx")
    # iterate each directory, run kraken and create krona graph
    os.chdir(root_dir)
    directory_list = [ddd for ddd in os.listdir(os.getcwd()) if os.path.isdir(ddd)]
    if parameter_dict['no_kraken']:
        print("\n*** Kraken not ran")
    else:
        for directory in directory_list:
            os.chdir(directory)
            read1 = glob.glob('*_R1*fastq.gz')[0]
            read2 = glob.glob('*_R1*fastq.gz')[0]
            print(f'\n\n{directory} Kraken running...')
            try:
                krona_path = simple_kraken1(directory, read1, read2, 'host')
            except:
                print(f'\n\n###### {sample_name} assembly not found.  Assume Kraken failed, check FASTQ quality\n\n')
                with open('kraken_failed', 'w') as message:
                    print(f'{sample_name} assembly not found.  Assume Kraken failed, check FASTQ quality\n', file=message)
                parameter_dict['Total size of scaffolds'] = "Kraken failed"
                return parameter_dict
            os.chdir(root_dir)

    #clean up directories
    if parameter_dict['debug_call']:
        print("\n*** Debugging -- all files kept")

    else:
        os.chdir(root_dir)
        directory_list = [ddd for ddd in os.listdir(os.getcwd()) if os.path.isdir(ddd)]
        for directory in directory_list:
            os.chdir(directory)
            os.mkdir("original_reads")
            os.mkdir("kraken")
            os.mkdir("latex")
            try:
                shutil.rmtree('./assembly/spades_output')
            except FileNotFoundError:
                print(f'No ./assembly/spades_output... likely testing with -t FASTA')
                pass
            fastqs = glob.glob('*.gz')
            for fastq in fastqs:
                shutil.move(fastq, "./original_reads")
            kraken_files = glob.glob('*kraken.txt')
            for kraken_file in kraken_files:
                shutil.move(kraken_file, "./kraken")
            try:
                shutil.move(glob.glob('*.tex')[0], "latex")
            except IndexError:
                pass
            files_grabbed = [glob.glob(e) for e in ['*.aux', '*.log', '*.out', '*.svg', '*.png', '*.txt']]
            for file_list in files_grabbed:
                for afile in file_list:
                    os.remove(afile)
            try:
                ss = glob.glob('./seqsero/seqsero_files/*')
                for afile in ss:
                    shutil.move(afile, './seqsero')
                shutil.rmtree('./seqsero/seqsero_files')
                shutil.rmtree('./taxonomy.krona.html.files')
            except FileNotFoundError:
                pass
            os.chdir(root_dir)

    runtime = (datetime.now() - startTime)
    os.chdir(parameter_dict['root_dir'])
    summary_file = root_dir + '/cumulative_stats_' + st + '.xlsx'
    df_concat = pd.concat(frames)
    df_concat.to_excel(summary_file, index=False)

    if parameter_dict['email_list']:
        try:
            send_email(parameter_dict['email_list'], runtime, summary_file, st)
            print(f'*** Email sent to: {parameter_dict["email_list"]}')
        except TimeoutError:
            #debug_log(ex, inspect.getframeinfo(inspect.currentframe()), "Unable to send email with current smtp setting")
            pass

    print("\n\t*** arm_workflow.py has completed ***")
    print("\t   *** runtime: %s:  \n" % runtime)

def assembly_pipe(sample_name, parameter_dict):
    print(f'###### --> {sample_name}')
    sample_directory = parameter_dict['root_dir'] + "/" + sample_name
    parameter_dict['sample_name'] = sample_name
    print(f'###### --> {sample_directory}')
    parameter_dict["sample_directory"] = sample_directory
    os.chdir(sample_directory)
    read1 = glob.glob('*_R1*fastq.gz')
    read2 = glob.glob('*_R2*fastq.gz')
    if len(read1) > 1:
        print("*** More than 1 R1 read")
        parameter_dict["File_R1"] = "More than 1 R1 read"
        return parameter_dict
    if len(read1) != len(read2):
        print("*** Unmatched reads")
        parameter_dict["File_R1"] = "Unmatched reads"
        return parameter_dict
    read1 = f'{sample_directory}/{read1[0]}'
    read2 = f'{sample_directory}/{read2[0]}'
    parameter_dict["read1"] = read1
    parameter_dict["read2"] = read2

    os.mkdir("assembly")
    os.chdir("assembly")
    if parameter_dict["test_fasta"]:
        assembled_fasta = parameter_dict["test_fasta"]
        shutil.copy(f'{root_dir}/{assembled_fasta}', f'{sample_directory}')
        os.chdir(sample_directory)
    else:
        print("\n*** Running trim reads...")
        read1_trimmed, read2_trimmed = trim_reads(sample_name, read1, read2)
        assembled_fasta = simple_spades(sample_name, read1_trimmed, read2_trimmed)
        trimmed_reads = glob.glob('*gz')
        for read in trimmed_reads:
            os.remove(read)
        try:
            SeqIO.write(assembled_fasta, f'{sample_name}_scaffolds.fasta', "fasta")
        except FileNotFoundError:
            print(f'\n\n###### {sample_name} assembly not found.  Assume SPAdes failed, check FASTQ quality\n\n')
            with open('assembly_failed', 'w') as message:
                print(f'{sample_name} assembly not found.  Assume SPAdes failed, check FASTQ quality\n', file=message)
            parameter_dict['Number of scaffolds'] = "SPAdes failed"
            return parameter_dict
        os.chdir(sample_directory)
        shutil.copy(f'./assembly/{sample_name}_scaffolds.fasta', 'spades_scaffolds.fasta')

    spades_scaffolds = f'{sample_directory}/spades_scaffolds.fasta'
    
    os.mkdir("mlst")
    os.chdir("mlst")
    print("\n*** Running MLST...")
    os.system(f'mlst {spades_scaffolds} > {sample_name}_mlst.txt')
    mlst_file = f'{sample_directory}/mlst/{sample_name}_mlst.txt'
    os.chdir(sample_directory)

    os.mkdir("seqsero")
    os.chdir("seqsero")
    print("\n*** Running SeqSero2...")
    os.system(f'SeqSero2_package.py -d seqsero_files -t 2 -i {read1} {read2}')
    if os.path.isfile(f'{sample_directory}/seqsero/seqsero_files/SeqSero_result.txt'):
        seqsero_file = f'{sample_directory}/seqsero/seqsero_files/SeqSero_result.txt'
    else:
        print(f'\t####### {sample_directory} #########\n\tSEQSERO2 FAILED\n\t####### {sample_directory} #########')
        seqsero_file = None
    os.chdir(sample_directory)

    os.mkdir("quality_stats")
    os.chdir("quality_stats")
    print("\n*** Getting FASTQ stats...")
    os.system(f'fastq_quality_amr.py -r1 {read1} -r2 {read2}')
    quality_stats = f'{sample_directory}/quality_stats/{sample_name}_stat_quality30.txt'
    os.system(f'assemblathon_reformat_stats.sh {spades_scaffolds}')
    assemblathon_stats = f'{sample_directory}/spades_scaffolds-assemblathon_reformat_stats.txt'
    os.chdir(sample_directory)
    stat_list = glob.glob('*assemblathon*')
    for assemblathon_stat in stat_list:
        shutil.move(assemblathon_stat, "quality_stats")
    assemblathon_stats = f'{sample_directory}/quality_stats/spades_scaffolds-assemblathon_reformat_stats.txt'
    
    os.mkdir("abricate")
    os.chdir("abricate")
    depth = parameter_dict["abricate_depth"]
    coverage = parameter_dict["abricate_coverage"]
    print("\n*** Running abricate resfinder (default)...")
    os.system(f'abricate --db resfinder --mincov {depth} --minid {coverage} {spades_scaffolds} > {sample_name}-resfinder.tab')
    tab_to_excel(f'{sample_name}-resfinder.tab')
    print("\n*** Running abricate NCBI...")
    os.system(f'abricate --db ncbi --mincov {depth} --minid {coverage} {spades_scaffolds} > {sample_name}-ncbi.tab')
    tab_to_excel(f'{sample_name}-ncbi.tab')
    print("\n*** Running abricate plasmidfinder...")
    os.system(f'abricate --db plasmidfinder --mincov {depth} --minid {coverage} {spades_scaffolds} > {sample_name}-plasmidfinder.tab')
    tab_to_excel(f'{sample_name}-plasmidfinder.tab')
    ab_resfinder_file = f'{sample_directory}/abricate/{sample_name}-resfinder.tab'
    ab_ncbi_file = f'{sample_directory}/abricate/{sample_name}-ncbi.tab'
    ab_plasmid_file = f'{sample_directory}/abricate/{sample_name}-plasmidfinder.tab'
    
    os.chdir(sample_directory)

    os.mkdir("amrfinder")
    os.chdir("amrfinder")
    print("\n*** Running amrfinder...")
    os.system(f'amrfinder.pl -n {spades_scaffolds} >> {sample_name}-amrfinder.tab')
    tab_to_excel(f'{sample_name}-amrfinder.tab')
    os.chdir(sample_directory)
    amrfinder_file = f'{sample_directory}/amrfinder/{sample_name}-amrfinder.tab' 

    parameter_dict["File_R1"] = os.path.basename(read1)
    parameter_dict["File_R2"] = os.path.basename(read2)
    
    parameter_dict = amr_report(parameter_dict, spades_scaffolds, assemblathon_stats, quality_stats, amrfinder_file, seqsero=seqsero_file, mlst=mlst_file, ab_resfinder_file=ab_resfinder_file, ab_ncbi_file=ab_ncbi_file, ab_plasmid_file=ab_plasmid_file)

    if parameter_dict['debug_call']:
        np.save('parameter_dict.npy', parameter_dict)
        #import with: parameter_dict = np.load('parameter_dict.npy').item()
    return parameter_dict

def merge_stats(parameter_dict):
    qf = pd.read_csv(parameter_dict["quality_file"], sep='\t')
    af = pd.read_csv(parameter_dict["assembly_stats"], sep='\t')
    af = af.drop(['File_R1', 'File_R2'], axis=1) # these columns conflict with keys in parameter_dict
    df = pd.concat([af, qf], axis=1, sort=False)
    merged_dict = df.to_dict(orient='index')[0]
    #assemblathon script doesn't pull read size when FASTQ not in working directory, needs updating
    merged_dict['Size_R1'] = humanize.naturalsize(os.path.getsize(parameter_dict["read1"]))
    merged_dict['Size_R2'] = humanize.naturalsize(os.path.getsize(parameter_dict["read2"]))
    return merged_dict

def seqsero_call(parameter_dict):
    seqsero_file = parameter_dict["seqsero_file"]
    if seqsero_file is not None:
        try:
            seqsero_indx = pd.read_csv(seqsero_file, sep='\t', skiprows=5, index_col=0, header=None)
            parameter_dict['predant'] = seqsero_indx.loc['Predicted antigenic profile:', :][1]
            parameter_dict['predsubsp'] = seqsero_indx.loc['Predicted subspecies:', :][1]
            parameter_dict['predst'] = seqsero_indx.iloc[2][1]
            parameter_dict['seqserocomment'] = seqsero_indx.tail(1).index.values[0]
        except (AttributeError, IndexError) as e:
            parameter_dict['predant'] = 'no data'
            parameter_dict['predsubsp'] = 'no data'
            parameter_dict['predst'] = 'no data'
            parameter_dict['seqserocomment'] = 'no data'
        return(parameter_dict)

def mlst_call(parameter_dict):
    parameter_dict['mlst_sp'] = "Not Identified"
    mlst_file = parameter_dict['mlst_file']
    look_up = parameter_dict['look_up']
    
    if mlst_file is not None:
        try:
            mlst_data = pd.read_csv(mlst_file, sep='\t', header=None)
            mlst_id = mlst_data.iloc[0,1]
            parameter_dict['mlst_id'] = mlst_id
            mlst_id = mlst_id.partition('_')
            mlst_sp = look_up.get(mlst_id[0], {}).get("Scientific Name", "Not Identified")
            parameter_dict['mlst_sp'] = mlst_sp
            parameter_dict['mlst_sch'] = mlst_id[2]
            parameter_dict['mlst_st'] = mlst_data.iloc[0,2]
            parameter_dict['mlst_g1'] = mlst_data.iloc[0,3]
            parameter_dict['mlst_g2'] = mlst_data.iloc[0,4]
            parameter_dict['mlst_g3'] = mlst_data.iloc[0,5]
            parameter_dict['mlst_g4'] = mlst_data.iloc[0,6]
            parameter_dict['mlst_g5'] = mlst_data.iloc[0,7]
            parameter_dict['mlst_g6'] = mlst_data.iloc[0,8]
            parameter_dict['mlst_g7'] = mlst_data.iloc[0,9]
        except IndexError:
            #debug_log(ex, inspect.getframeinfo(inspect.currentframe()), "MLST error")
            parameter_dict['mlst_st'] = "Not Identified"
            parameter_dict['mlst_id'] = None
            parameter_dict['mlst_sch'] = None
            parameter_dict['mlst_st'] = None
            parameter_dict['mlst_g1'] = None
            parameter_dict['mlst_g2'] = None
            parameter_dict['mlst_g3'] = None
            parameter_dict['mlst_g4'] = None
            parameter_dict['mlst_g5'] = None
            parameter_dict['mlst_g6'] = None
            parameter_dict['mlst_g7'] = None

    return parameter_dict

def genome_gt_1kb(parameter_dict):
    assembled_fasta = SeqIO.parse(parameter_dict["assemblied_fasta_path"], "fasta")
    measure_length = 0
    for contig in assembled_fasta:
        if len(contig) > 1000:
            measure_length = measure_length + len(contig)
    gl = (measure_length/parameter_dict['Total size of scaffolds'])*100
    rgl = round(gl, 2)
    glc = '%s %s' % (str(rgl), str(" \%"))
    parameter_dict['glc'] = glc
    parameter_dict['rgl'] = rgl
    return parameter_dict

def sequence_score(parameter_dict):
    look_up = parameter_dict['look_up']

    seqcor=[[1.0000000,0.4412253,0.9569653,0.4266978,0.3393063],
    [0.4412253,1.0000000,0.4446494,0.9527558,0.2464773],
    [0.9569653,0.4446494,1.0000000,0.4330931,0.3315736],
    [0.4266978,0.9527558,0.4330931,1.0000000,0.2367824],
    [0.3393063,0.2464773,0.3315736,0.2367824,1.0000000]]

    mlst_id = parameter_dict['mlst_id']
    if parameter_dict['genlength']:
        sp_len = parameter_dict['genlength']
        sp_method = 'Genome length was provided for the analysis.'
    else:
        pre_sp_len = look_up.get(mlst_id, {})
        sp_len = pre_sp_len.get("Approx Genome Size (Mb)", "5") #default to size of 5
    try:
        if mlst_id[0] == '-':
            sp_method = 'Default genome length of 5Mb was used.'
        else:
            sp_method = 'Genome length estimated based on MLST identification.'
    except TypeError:
        sp_method = 'Default genome length of 5Mb was used.'

    reads1 = parameter_dict["count1"]
    reads2 = parameter_dict["count2"]
    seqcov = int(reads1+reads2)*240/(int(sp_len) * 1000000) 
    #seqcov = int(reads1)*250/(int(sp_len) * 1000000)
    q30 = parameter_dict["q30pass1"]
    q301 = q30[0:4]
    q30 = parameter_dict["q30pass2"]
    q302 = q30[0:4]
    x = [(float(parameter_dict["mean_quality1"])-30), (float(parameter_dict["mean_quality2"])-30), (float(q301)-80), (float(q302)-80), seqcov]
    mx = np.matrix(x)
    mseqcor = np.matrix(seqcor)
    mxt = mx.transpose()
    seqv = math.sqrt(mx * mseqcor * mxt)
    seqval = (seqv*800/90) + 1180

    if seqval > 1980 :
        seqval = 1980
    elif seqval < 1180:
        seqval = 1180

    parameter_dict['sp_method'] = sp_method
    parameter_dict['seqcov'] = seqcov
    parameter_dict['seqval'] = seqval
    parameter_dict['sp_len'] = int(sp_len)
    return(parameter_dict)

def report_setup(parameter_dict):
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
    seqval = parameter_dict['seqval']
    mysvg.add( mysvg.ellipse( [seqval,65], [10,15], fill='white', stroke_width=1, stroke='midnightblue', opacity=0.85))
    #Add low quality label on the x-axis
    mysvg.add(mysvg.text('Low', insert=(1100,75), fill='white', font_size='35px', font_weight="bold", style = "font-family:Arial"))
    #Add high quality label on the x-axis
    mysvg.add(mysvg.text('High', insert=(1995,75), fill='white', font_size='35px', font_weight="bold", style="font-family:Arial"))
    #ADD File Stats header
    mysvg.add(mysvg.text('Sequence Statistics', insert=(30,70), fill='white', font_size='50px', font_weight="bold"))
    mysvg.saveas('seq.svg')
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
    sp_len = parameter_dict["sp_len"]
    gt1k = parameter_dict["Number of scaffolds > 1K nt"]
    numscaff = parameter_dict["Number of scaffolds"]
    n50len = parameter_dict["N50 scaffold length"]
    l50cnt = parameter_dict["L50 scaffold count"]
    rgl = parameter_dict["rgl"]

    lscaff = parameter_dict["Longest scaffold"]
    perdiffexpect = abs(100-(float(parameter_dict['Total size of scaffolds'])/(int(sp_len)*10000)))
    perlongexpect = float(lscaff)/(int(sp_len)*1000000)*100
    perlongscaff = float(gt1k)/float(numscaff)*100
    pern50expect = int(n50len)/(int(sp_len)*1000000)*100

    quality_scale_list = [numscaff, gt1k, l50cnt, rgl, perdiffexpect, perlongexpect, perlongscaff, pern50expect] 

    ma = np.matrix(quality_scale_list)
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
    ab_version = os.popen("abricate --version | sed 's/abricate//g'").readlines()[0]
    parameter_dict["ab_version"] = ab_version
    ncbi = os.popen("abricate --list | awk 'BEGIN{OFS=\"\t\"} $1 == \"ncbi\" {print $1, $2, $4}'|tr '\t' '&'").readlines()[0]
    parameter_dict["ncbi"] = ncbi
    res = os.popen("abricate --list | awk 'BEGIN{OFS=\"\t\"} $1 == \"resfinder\" {print $1, $2, $4}'|tr '\t' '&'").readlines()[0]
    parameter_dict["res"] = res
    ####AMRFinder
    amr_version = os.popen("amrfinder.pl -v| sed '1!d' |sed 's/version $Revision://g'|sed 's/\$//g'").readlines()[0]

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
    spadesver = os.popen("spades.py --version").readlines()[0]
    parameter_dict['spadesver'] = spadesver
    parameter_dict['quality_scale_list'] = quality_scale_list
    parameter_dict['assemblev'] = assemblev
    parameter_dict['assembleval'] = assembleval
    return(parameter_dict)

def latex_document(parameter_dict):
    try:
        #dump keys as variables (does not work)
        # for key, item in parameter_dict.items():
        #     key = key.replace(' ', '_')
        #     key = key.replace('-', '_')
        #     key = key.replace('>', 'gt')
        #     print(key)
        #     exec(key + '=item')
        assembly_stats = parameter_dict.get("assembly_stats", None)
        quality_file = parameter_dict.get("quality_file", None)
        seqsero_file = parameter_dict.get("seqsero_file", None)
        mlst_file = parameter_dict.get("mlst_file", None)
        ab_resfinder_file = parameter_dict.get("ab_resfinder_file", None)
        ab_ncbi_file = parameter_dict.get("ab_ncbi_file", None)
        ab_plasmid_file = parameter_dict.get("ab_plasmid_file", None)
        assemblied_fasta_path = parameter_dict.get("assemblied_fasta_path", None)
        amrfinder_file = parameter_dict.get("amrfinder_file", None)
        abricate_depth = parameter_dict.get("abricate_depth", None)
        abricate_coverage = parameter_dict.get("abricate_coverage", None)
        genlength = parameter_dict.get("genlength", None)
        sample_name = parameter_dict.get("sample_name", None)
        Sample = parameter_dict.get("Sample", None)
        File_R1 = parameter_dict.get("File_R1", None)
        Size_R1 = parameter_dict.get("Size_R1", None)
        File_R2 = parameter_dict.get("File_R2", None)
        Size_R2 = parameter_dict.get("Size_R2", None)
        Number_of_scaffolds = parameter_dict.get("Number of scaffolds", None)
        Total_size_of_scaffolds = parameter_dict.get("Total size of scaffolds", None)
        Longest_scaffold = parameter_dict.get("Longest scaffold", None)
        Number_of_scaffolds_gt_1K_nt = parameter_dict.get("Number of scaffolds > 1K nt", None)
        N50_scaffold_length = parameter_dict.get("N50 scaffold length", None)
        L50_scaffold_count = parameter_dict.get("L50 scaffold count", None)
        mean_quality1 = parameter_dict.get("mean_quality1", None)
        mean_quality2 = parameter_dict.get("mean_quality2", None)
        q30pass1 = parameter_dict.get("q30pass1", None)
        q30pass2 = parameter_dict.get("q30pass2", None)
        count1 = parameter_dict.get("count1", None)
        count2 = parameter_dict.get("count2", None)
        mlst_sp = parameter_dict.get('mlst_sp', None)
        mlst_id = parameter_dict.get("mlst_id", None)
        mlst_sch = parameter_dict.get("mlst_sch", None)
        mlst_st = parameter_dict.get("mlst_st", None)
        mlst_g1 = parameter_dict.get("mlst_g1", None)
        mlst_g2 = parameter_dict.get("mlst_g2", None)
        mlst_g3 = parameter_dict.get("mlst_g3", None)
        mlst_g4 = parameter_dict.get("mlst_g4", None)
        mlst_g5 = parameter_dict.get("mlst_g5", None)
        mlst_g6 = parameter_dict.get("mlst_g6", None)
        mlst_g7 = parameter_dict.get("mlst_g7", None)
        mlst_version = parameter_dict.get("mlst_version", None)
        sp_method = parameter_dict.get("sp_method", None)
        seqcov = parameter_dict.get("seqcov", None)
        seqval = parameter_dict.get("seqval", None)
        sp_len = parameter_dict.get("sp_len", None)
        glc = parameter_dict.get("glc", None)
        rgl = parameter_dict.get("rgl", None)
        ab_version = parameter_dict.get("ab_version", None)
        ncbi = parameter_dict.get("ncbi", None)
        res = parameter_dict.get("res", None)
        spadesver = parameter_dict.get("spadesver", None)
        quality_scale_list = parameter_dict.get("quality_scale_list", None)
        assemblev = parameter_dict.get("assemblev", None)
        assembleval = parameter_dict.get("assembleval", None)
        predant = parameter_dict.get("predant", None)
        predsubsp = parameter_dict.get("predsubsp", None)
        predst = parameter_dict.get("predst", None)
        seqserocomment = parameter_dict.get("seqserocomment", None)

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
        print(r'Filename &  %s & %s \\' % (File_R1.replace('_', '\_'), File_R2.replace('_', '\_')), file=write_out )
        print(r'\hline', file=write_out)
        print(r'File size & %s & %s \\' % (Size_R1, Size_R2), file=write_out)
        print(r'Mean Read Score & %s & %s \\' % (mean_quality1, mean_quality2), file=write_out)
        print(r'Q30 Passing & %s & %s \\' % (q30pass1.replace('%', ' \%'), q30pass2.replace('%', ' \%')), file=write_out)
        print(r'\hline', file=write_out)
        print(r'\end{tabular}', file=write_out)
        print(r'', file=write_out)
        print(r'\begin{tabular}{ p{3.2cm} | p{3cm} p{11cm} }', file=write_out)
        print(r'Sequence Depth & %sX  & Calculated by Number of reads x 250/Genome Length \\' % (str(seqcov)), file=write_out)
        print(r'\hline', file=write_out)
        print(r'Genome Length & %sbp & %s \\' % ("{:,}".format(int(sp_len) * 1000000), sp_method), file=write_out)
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
        print(r'%s & %s & %s & %s & %s & %s & %s \\' % (Number_of_scaffolds, "{:,}".format(Total_size_of_scaffolds), "{:,}".format(Longest_scaffold), Number_of_scaffolds_gt_1K_nt, glc, N50_scaffold_length, L50_scaffold_count), file=write_out)
        print(r'\hline', file=write_out)
        print(r'\end{tabularx}', file=write_out)
        print(r'', file=write_out)
        print(r'{\footnotesize\ De novo assembly performed using \href{http://cab.spbu.ru/software/spades/} {\textcolor{midnightblue}{%s}.} }\\,' % spadesver, file=write_out)

        if mlst_file is not None:
            print(r'\vspace{3mm}', file=write_out)
            print(r'', file=write_out)
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

            if mlst_id == "senterica":  
                print(r'\vspace{4mm}', file=write_out)
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

        #AMRFinder Results Page
        print(r'\newpage', file=write_out)
        print(r'', file=write_out)
        print(r'\vspace{10mm}', file=write_out)
        print(r'', file=write_out)
        amrfinder_lines = sum(1 for line in open(amrfinder_file))
        if amrfinder_lines > 1:
            print(r'\includegraphics[scale=0.485]{amrfinder.png}', file=write_out)
            print(r'', file=write_out)
            print(r'\vspace{-3mm}', file=write_out)
            print(r'\resizebox{27cm}{!}{', file=write_out)
            print(r'\begin{tabular}{c|c|c|c|c|c|c|c|c|c|c|c|c}', file=write_out)
            print(r'\textcolor{midnightblue}{\bf Contig id}&\textcolor{midnightblue}{\bf Start}&\textcolor{midnightblue}{\bf Stop}&\textcolor{midnightblue}{\bf Gene symbol}&\textcolor{midnightblue}{\bf Protein name}&\begin{rotate}{60}\textcolor{midnightblue}{\bf Method}\end{rotate}&\begin{rotate}{60}\textcolor{midnightblue}{\bf Target length}\end{rotate}&\begin{rotate}{60}\textcolor{midnightblue}{\bf Reference protein length}\end{rotate}& \begin{rotate}{60}\textcolor{midnightblue}{\bf \% Coverage of reference protein}\end{rotate}&\begin{rotate}{60}\textcolor{midnightblue}{\bf \% Identity of reference protein}\end{rotate}& \begin{rotate}{60}\textcolor{midnightblue}{\bf Alignment length}\end{rotate}& \begin{rotate}{60}\textcolor{midnightblue}{\bf Accession of closest protein} \end{rotate}& \begin{rotate}{60}\textcolor{midnightblue}{\bf HMM id}\end{rotate}  \\', file=write_out)
            print(r'\vspace{1mm}', file=write_out)
            print(r'\hline', file=write_out)
            df = pd.read_csv(amrfinder_file, sep='\t', header=None, skiprows=1)
            for index_num in df.index:
                series = df.T[index_num]
                string = f'{series[1]}&{series[2]}&{series[3]}&{series[5]}&{series[6]}&{series[7]}&{series[8]}&{series[9]}&{series[10]}&{series[11]}&{series[12]}&{series[13]}&{series[15]}'
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
        os.system("pdflatex -interaction=nonstopmode {}".format(tex_document))
        os.system("pdflatex -interaction=nonstopmode {}".format(tex_document))

    except KeyError:
        #debug_log(ex, inspect.getframeinfo(inspect.currentframe()), "See Example 2 from help option --> parameter_dict = np.load('parameter_dict.npy').item()")
        pass
    
def amr_report(parameter_dict, contigs, assemble_stats, quality_file, amrfinder_file, seqsero=None, mlst=None, ab_resfinder_file=None, ab_ncbi_file=None, ab_plasmid_file=None, genlength=None, **kwargs):

    parameter_dict["assemblied_fasta_path"] = contigs
    parameter_dict["assembly_stats"] = assemble_stats
    parameter_dict["quality_file"] = quality_file
    parameter_dict["amrfinder_file"] = amrfinder_file
    parameter_dict["seqsero_file"] = seqsero
    parameter_dict["mlst_file"] = mlst
    parameter_dict["ab_resfinder_file"] = ab_resfinder_file
    parameter_dict["ab_ncbi_file"] = ab_ncbi_file
    parameter_dict["ab_plasmid_file"] = ab_plasmid_file
    parameter_dict["genlength"] = genlength

    ####Lookup- MLST Name and Approx Genome Size 
    real_path = os.path.dirname(os.path.realpath(__file__))
    lookup_genome_size = real_path + "/lookup_genome_size.json"
    with open(lookup_genome_size) as infile:
        look_up = json.load(infile)
    parameter_dict["look_up"] = look_up

    print("\n*** Merging stats...")
    merged_dict = merge_stats(parameter_dict)
    parameter_dict = {**parameter_dict, **merged_dict} #merge dictionaries
    print("\n*** Gathering MLST stats...")
    parameter_dict = mlst_call(parameter_dict)
    print("\n*** Getting sequence score...")
    parameter_dict = sequence_score(parameter_dict)
    print("\n*** Getting contigs larger than 1Kb...")
    parameter_dict = genome_gt_1kb(parameter_dict)
    print("\n*** Setting up the report...")
    parameter_dict = report_setup(parameter_dict)
    print("\n*** Preparing SeqSero...")
    if os.path.isfile(f'{parameter_dict["sample_directory"]}/seqsero/seqsero_files/SeqSero_result.txt'):
        parameter_dict = seqsero_call(parameter_dict)
    latex_document(parameter_dict)

    return parameter_dict


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

    ---------------------------------------------------------
    Example1:  ## MUST BE "spades_scaffolds.fasta" named ##
    ipython: %run -d -b 1150 /home/tstuber/git/gitlab/stuber/amr_python/amr_workflow.py -d -t spades_scaffolds.fasta
    
    Example2:
    ipython
    import numpy as np
    parameter_dict = np.load('parameter_dict.npy').item()
    from amr_workflow import latex_document #test report generation, must have export PYTHONPATH="${PYTHONPATH}:${HOME}/git/gitlab/stuber/amr_python" in shell PATH
    latex_document(parameter_dict) #put into jupyter notebook for helpful debugging

    Usage:
    arm_workflow.py is called on a working directory containing paired FASTQs

    Additional:
    If NAHLM samples... separate after running amr_workflow.sh on all samples.  Separate as EC (ecoli), MH, SIG (staph intermedius group) then run nahln_amr_updates.sh on each group.  Before running make sure the assemblies (.fasta) are in the sample folders.  May need to copy them from the assembly folder to the root of the sample folder.  Same with amrfind.tab.
    for i in *; do cp -v $i/assembly/*scaffolds.fasta $i; done
    for i in *; do cp -v $i/amrfinder/*amrfinder.tab $i; done

    '''), epilog='''---------------------------------------------------------''')

    parser.add_argument('-i', '--isescan', action='store_true', dest='isescan', help='OPTIONAL: run ISEScan')
    parser.add_argument('-a', '--abricate_depth', action='store', dest='abricate_depth', default=0, help='OPTIONAL: percent average depth cutoff for abricate, aka: --mincov, cvb use -a 50')
    parser.add_argument('-b', '--abricate_coverage', action='store', dest='abricate_coverage', default=75, help='OPTIONAL: percent genome coverage cutoff for abricate, aka: --minid, cvb use -b 90')
    parser.add_argument('-d', '--debug', action='store_true', dest='debug_call', help='debug, run without pool.map')
    parser.add_argument('-n', '--no_kraken', action='store_true', dest='no_kraken', help='skip kraken/krona')
    parser.add_argument('-t', '--test_fasta', action='store', dest='test_fasta', help='OPTIONAL: when debugging provide FASTA skip the long assembly process.  When using FASTA working directory must also have FASTQ files. Recommend using with -n option.')
    parser.add_argument('-m', '--email', action='store', dest='email', help='Send email at completion')
    parser.add_argument('-c', '--cluster', action='store_true', dest='cluster', help='Run on HPC, cluster.scale(5)')
    args = parser.parse_args()

    parameter_dict = {
        "isescan": args.isescan,
        "abricate_depth": args.abricate_depth,
        "abricate_coverage": args.abricate_coverage,
        "debug_call": args.debug_call,
        "no_kraken": args.no_kraken,
        "test_fasta": args.test_fasta,
        "cluster": args.cluster,
        # default values
        'sample_name': 'na', 
        'File_R1': 'na', 
        'Size_R1': 'na', 
        'mean_quality1': 'na', 
        'q30pass1': 'na', 
        'count1': 'na', 
        'File_R2': 'na', 
        'Size_R2': 'na', 
        'mean_quality2': 'na', 
        'q30pass2': 'na', 
        'count2': 'na', 
        'Number of scaffolds': 'na', 
        'Total size of scaffolds': 'na', 
        'Longest scaffold': 'na', 
        'Number of scaffolds > 1K nt': 'na', 
        'N50 scaffold length': 'na', 
        'L50 scaffold count': 'na',
    }

    parameter_dict['root_dir'] = os.getcwd()
    kickoff(parameter_dict)
