#!/usr/bin/env python

import os
import sys
import logging
from argparse import ArgumentParser
from Bio import SeqIO, BiopythonWarning
from Bio.Alphabet import generic_dna
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger('CDS2Prot')
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s %(levelname)s: %(funcName)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)


class ORFFinder:
    """Find the longest ORF in a given sequence
    "seq" is a string, if "start" is not provided any codon can be the start of
    and ORF. If muliple ORFs have the longest length the first one encountered
    is printed

    Code writted to answer this challenge at Biostar:
    http://biostar.stackexchange.com/questions/5902/

    (This code includes improvements from Brad Chapman)

    Source: https://gist.github.com/dwinter/933737
    """
    def __init__(self, seq, start=[], stop=["TAG", "TAA", "TGA"]):
        self.seq = seq
        self.start = start
        self.stop = stop
        self.result = ("+",0,0,0,0)
        self.winner = 0

    def _reverse_comp(self):
        swap = {"A":"T", "T":"A", "C":"G", "G":"C", "N":"N"}
        return "".join(swap[b] for b in self.seq)

    def _print_current(self):
        print "frame %s%s position %s:%s (%s nucleotides)" % self.result

    def codons(self, frame):
        """ A generator that yields DNA in one codon blocks
        "frame" counts for 0. This function yelids a tuple (triplet, index) with
        index relative to the original DNA sequence
        """
        start = frame
        while start + 3 <= len(self.seq):
            yield (self.seq[start:start+3], start)
            start += 3

    def run_one(self, frame_number, direction):
        """ Search in one reading frame """
        orf_start = None
        for c, index in self.codons(frame_number):
            if (c not in self.stop and (c in self.start or not self.start)
                and orf_start is None):
                orf_start = index + 1 # we'll return the result as 1-indexed
            elif c in self.stop and orf_start:
                self._update_winner(orf_start, index, direction, frame_number)
                orf_start = None
        if orf_start:
            self._update_winner(orf_start, index, direction, frame_number)

    def _update_winner(self, orf_start, index, direction, frame_number):
        orf_end = index + 3 # because index is realitve to start of codons
        L = (orf_end - orf_start) + 1
        if L > self.winner:
            self.winner = L
            self.result = (direction, frame_number+1, orf_start, orf_end, L)

    def run(self):
        direction = "+"
        for frame in range(3):
            self.run_one(frame, direction)
        direction = "-"
        self.seq = self._reverse_comp()
        for frame in range(3):
            self.run_one(frame, direction)
        #self._print_current()

def find_ORF(seq_string):
    finder = ORFFinder(seq_string)
    finder.run()
    return finder

def seq_is_mod_three(sequence):
    if type(sequence) is str:
        return len(sequence) % 3
    elif type(sequence) is Seq:
        return len(str(sequence)) % 3
    elif type(sequence) is SeqRecord:
        return len(str(sequence.seq)) % 3
    else:
        raise TypeError('Unrecognized type for sequence: %s', type(sequence))

def begin_start_codon(sequence, ctable):
    #logger.info('Start codons: %s', ctable.start_codons)
    #logger.info('Stop codons: %s', ctable.stop_codons)
    first_codon = ''
    if type(sequence) is str:
        first_codon = sequence[0:3]
    elif type(sequence) is Seq:
        first_codon = str(sequence)[0:3]
    elif type(sequence) is SeqRecord:
        first_codon = str(sequence.seq)[0:3]
    else:
        raise TypeError('Unrecognized type for sequence: %s', type(sequence))
    #logger.info('First codon is %s. Is it a start codon? %s', first_codon, bool(first_codon in ctable.start_codons))
    return bool(first_codon in ctable.start_codons)

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
    if type(args.xtable) is int:
        ctable = CodonTable.unambiguous_dna_by_id[args.xtable]
    else:
        ctable = CodonTable.unambiguous_dna_by_name[args.xtable]
    ridx = 0
    for record in SeqIO.parse(sys.stdin, args.informat):
        ridx += 1
        try:
            orf = ORFFinder(str(record.seq), ctable.start_codons, ctable.stop_codons)
            orf.run()
            #logger.info('Find ORF returned %s', orf.result)
            (d, f, s, e, l) = orf.result
            rec = SeqRecord(Seq(str(record.seq[s-1:e])), id=record.id, description='description removed')
            #logger.info('ORF Record:\n%s', rec.translate(args.xtable, id=record.id, description='description removed').format('fasta'))
            if seq_is_mod_three(record):
                logger.warning('Sequence %d length is not divisible by three. Using sequence from ORF search.', ridx)
                sys.stdout.write(rec.translate(args.xtable, id=record.id, description='description removed').format('fasta'))
            elif not begin_start_codon(record, ctable):
                logger.warning('Sequence %d does not begin with a start codon. Using sequence from ORF search.', ridx)
                sys.stdout.write(rec.translate(args.xtable, id=record.id, description='description removed').format('fasta'))
            else:
                sys.stdout.write(record.translate(table=args.xtable, id=record.id, description='description removed').format('fasta') + '\n')
        except BiopythonWarning as err:
            logger.warn('Caught BiopythonWarning about sequence %d:\n%s\n\n%s', ridx, record, err)

if __name__ == '__main__':
    sys.exit(main())
