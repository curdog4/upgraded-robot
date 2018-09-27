#!/usr/bin/env python

import os
import sys
import logging
from argparse import ArgumentParser
from Bio import SeqIO, BiopythonWarning

logger = logging.getLogger('CDS2Prot')
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s %(levelname)s: %(funcName)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)

def main():
    parser = ArgumentParser(description='Translate nucleotide to protein sequences')
    parser.add_argument('--informat', type=str, default='fasta',
                        help='Format of input')
    parser.add_argument('--outformat', type=str, default='fasta',
                        help='Format of output')
    parser.add_argument('--xtable', type=int, default=1,
                        help='Codon translation table')
    args = parser.parse_args()
    logger.info('Piping sys.stdin through SeqIO to output CDS -> protein sequences...')
    logger.info('Input format is %s', args.informat)
    logger.info('Output format is %s', args.outformat)
    logger.info('Codon table is %d', args.xtable)
    for record in SeqIO.parse(sys.stdin, args.informat):
        try:
            sys.stdout.write(record.translate(table=args.xtable, id=record.id, description='description removed').format('fasta') + '\n')
        except BiopythonWarning as err:
            logger.warn('Caught BiopythonWarning about sequence:\n%s\n\n%s', record, err)

if __name__ == '__main__':
    sys.exit(main())
