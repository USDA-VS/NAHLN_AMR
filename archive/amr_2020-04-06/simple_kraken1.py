#!/usr/bin/env python

import os
import sys
import argparse
import textwrap
import multiprocessing
import pandas as pd
import smtplib
from email.utils import formatdate
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.mime.base import MIMEBase
from email import encoders

KDB = os.environ['KDB'] #grab environment variable

def send_email(email_list, name, database, krona_path):
    text = "See attached Krona graph."
    msg = MIMEMultipart()
    msg['From'] = "tod.p.stuber@usda.gov"
    msg['To'] = email_list
    msg['Date'] = formatdate(localtime=True)
    msg['Subject'] = f'Kraken1 {database} {name}'
    msg.attach(MIMEText(text))
    part = MIMEBase('application', "octet-stream")
    part.set_payload(open(krona_path, "rb").read())
    encoders.encode_base64(part)
    part.add_header('Content-Disposition', f'attachment; filename={name}-krona_graph.html')
    msg.attach(part)
    smtp = smtplib.SMTP('10.10.8.12')
    smtp.send_message(msg)
    smtp.quit()

def select_database(arg_in):
    if arg_in == "flu":
        database = KDB + "/kraken/flu_jhu/fludb_20150820_with_hosts"
    elif arg_in == "host":
        database = KDB + "/kraken/kraken1"
    else:
        print("###Exiting --> \n provide database as flu or host, if calling from cmdline use -f or -h")
    return database

def make_krona(name, db_name, report, output, database):
    kraken_path = os.getcwd()
    if db_name == "flu":
        os.system(f'krakenreport2krona.sh -i {report} -k {database} -t {output} -o {name}-jhu-Krona_id_graphic.html')
    elif db_name == "host":
        output = pd.read_csv(output, sep='\t')
        output.drop(output.columns[[0, 3, 4]], axis=1, inplace=True)
        output.to_csv("cut_output.txt", sep='\t', index=False)
        os.system(f'ktImportTaxonomy cut_output.txt')
        os.rename('taxonomy.krona.html', f'{name}-taxonomy.krona.html')
        absolute_file_path = kraken_path + f'/{name}-taxonomy.krona.html'
        return(absolute_file_path)

def simple_kraken1(name, read1, read2, *args):
    db_name = args[0]
    database = select_database(db_name)
    cpus = int(multiprocessing.cpu_count()/2)
    os.system(f'kraken --db {database} --threads {cpus} --paired {read1} {read2} > {name}-outputkraken.txt')
    os.system(f'kraken-report --db {database} {name}-outputkraken.txt > {name}-reportkraken.txt')
    krona_path = make_krona(name, db_name, f'{name}-reportkraken.txt', f'{name}-outputkraken.txt', database)
    return krona_path

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

    ---------------------------------------------------------

    '''), epilog='''---------------------------------------------------------''')

    parser.add_argument('-n', '--name', action='store', dest='name', required=True, help='sample name')
    parser.add_argument('-r1', '--read1', action='store', dest='read1', required=True, help='R1 read')
    parser.add_argument('-r2', '--read2', action='store', dest='read2', required=True, help='R2 read')
    parser.add_argument('-f', '--flu', action='store_true', dest='flu', help='Use flu database')
    parser.add_argument('-o', '--host', action='store_true', dest='host', help='Use host database')
    parser.add_argument('-m', '--email', action='store', dest='email', help='Send email at completion')
    args = parser.parse_args()

    email_dict = {}
    email_dict["tod"] =  "tod.p.stuber@aphis.usda.gov"
    email_dict["jess"] =  "Jessica.A.Hicks@aphis.usda.gov"
    email={}
    email['email_list'] = email_dict.get(args.email, None)

    if args.flu:
        database = "flu"
    if args.host:
        database = "host"

    krona_path = simple_kraken1(args.name, args.read1, args.read2, database)
    if email['email_list']:
        send_email(email['email_list'], args.name, database, krona_path)
