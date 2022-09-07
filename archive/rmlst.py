#!/usr/bin/env python3
# Upload contigs file to PubMLST rMLST species identifier via RESTful API
# Written by Keith Jolley
# Copyright (c) 2018, University of Oxford
# Licence: GPL3
# updated by Tod Stuber 2020-03-25

import os
import sys
import requests
import argparse
import textwrap
import base64

class rMLST:

    def __init__(self, FASTA):
        self.rank = None
        self.taxon = None
        self.support = None
        self.taxonomy = None
        uri = 'http://rest.pubmlst.org/db/pubmlst_rmlst_seqdef_kiosk/schemes/1/sequence'
        with open(FASTA, 'r') as open_fasta: 
            fasta = open_fasta.read()
        payload = '{"base64":true,"details":true,"sequence":"' + base64.b64encode(fasta.encode()).decode() + '"}'
        response = requests.post(uri, data=payload)
        if response.status_code == requests.codes.ok:
            data = response.json()
            try: 
                data['taxon_prediction']
            except KeyError:
                print("No match provided by rMLST")
                sys.exit(0)
            for match in data['taxon_prediction']:
                self.rank = match['rank']
                self.taxon = match['taxon']
                self.support = match['support']
                self.taxonomy = match['taxonomy']
        else:
            print(response.text)

if __name__ == "__main__": # execute if directly access by the interpreter
    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

    ---------------------------------------------------------
    Place description

    ''' ), epilog='''---------------------------------------------------------''')

    parser.add_argument('-f', '--FASTA', action='store', dest='FASTA', required=True, default=None, help='Assembly FASTA file')
    args = parser.parse_args()
    print(f'\n{os.path.basename(__file__)} SET ARGUMENTS:')
    print(args)

    rmlst = rMLST(args.FASTA)
    print(f'Rank: {rmlst.rank}')
    print(f'Taxon: {rmlst.taxon}')
    print(f'Suppart: {rmlst.support}')
    print(f'Taxonomy: {rmlst.taxonomy}')
