#!/usr/bin/env python

import os
import sys
import logging
from argparse import ArgumentParser
from Bio import SeqIO

logger = logging.getLogger('CDSXl8')
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s %(levelname)s: %(funcName)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)

def main():
    parser = ArgumentParser(description="Translate CDS sequence into protein sequence")
    parser.add_argument('--codon-table', type=int, default=1,
                        help='Codon table to use for translation')
    parser.add_argument('--sequence-id', type=str, required=True,
                        help='Sequence identifier to translate')
    parser.add_argument('--cds-fasta-file', type=str, required=True,
                        help='FASTA formatted file containing the sequence')
    parser.add_argument('--description', type=str,
                        help='Description to give the translated sequence')
    args = parser.parse_args()
    seq = None
    for rec in SeqIO.parse(args.cds_fasta_file, 'fasta'):
       if rec.id == args.sequence_id:
           seq = rec
           break
    description = args.description
    if not description:
        description = seq.description
    logger.info('Found sequence:\n%s', seq.translate(table=1, id=seq.id, description=description).format('fasta'))

if __name__ == '__main__':
    sys.exit(main())
