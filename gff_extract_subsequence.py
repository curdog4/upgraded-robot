'''
Extract genomic sequence range from provided genome FASTA and GFF files for the provided sequence ID and range
'''

import os
import sys
import logging
import logging.config
import json
import math
import argparse

import gffutils
from Bio import SeqIO

PADDING_FACTOR=0.25

logConfFile = 'logging.json'
with open(logConfFile, 'r') as fd:
    logConfContent = fd.read()
logConfData = json.loads(logConfContent)
logging.config.dictConfig(logConfData)
logger = logging.getLogger('SplitChimeras')


class Error(Exception):
    pass


class GffExtractionError(Error):
    pass


class SequenceMappingError(Error):
    pass


def adjust_coords(seq_range, mrna, gene, exons, pad_factor):
    start_pos = stop_pos = None
    offset = seq_range[0]
    remaining = feature_len = abs(seq_range[1] - seq_range[0])
    padding = int(math.ceil(pad_factor * float(feature_len)))
    logger.debug('Feature length is %d, padding is %d', feature_len, padding)
    exons = sorted(exons)
    if gene.strand == '-':
        exons.reverse()
    logger.debug('Exons: %s', exons)
    logger.debug('Exonic length: %d', sum([x[1] - x[0] for x in exons]))
    for exon in exons:
        logger.debug('Exon %s, length %d', exon, exon[1] - exon[0])
        exon_range = range(exon[0], exon[1]+1)
        if gene.strand == '-':
            if not stop_pos:
                if exon[1] - offset < exon[0]:
                    logger.debug('Feature %s, strand %s, offset %d, not within this exon...',
                                 mrna.id, gene.strand, offset)
                    offset -= exon[1] - exon[0] + 1
                    continue
                stop_pos = exon[1] - offset
                logger.debug('Feature %s, strand %s, stop position %d', mrna.id,
                             gene.strand, stop_pos)
                if stop_pos - exon[0] > feature_len:
                    remaining = 0
                    start_pos = stop_pos - feature_len
                    logger.debug('Feature %s, strand %s, start position %d',
                                 mrna.id, gene.strand, start_pos)
                    break
                else:
                    remaining -= stop_pos - exon[0] + 1
                    logger.debug('Feature %s, strand %s, remaining %d', mrna.id,
                                 gene.strand, remaining)
                    continue
            else:
                if exon[1] - exon[0] > remaining:
                    start_pos = exon[1] - remaining
                    logger.debug('Feature %s, strand %s, start position %d',
                                 mrna.id, gene.strand, start_pos)
                    break
                else:
                    remaining -= exon[1] - exon[0] + 1
                    logger.debug('Feature %s, strand %s, remaining %d', mrna.id,
                                 gene.strand, remaining)
                    continue
        else:
            if not start_pos:
                if exon[0] + offset > exon[1]:
                    logger.debug('Feature %s, strand %s, offset %d, not within this exon...',
                                 mrna.id, gene.strand, offset)
                    offset -= exon[1] -exon[0] + 1
                    continue
                start_pos = exon[0] + offset
                logger.debug('Feature %s, strand %s, start position %d', mrna.id,
                             gene.strand, start_pos)
                if exon[1] - start_pos > feature_len:
                    remaining = 0
                    stop_pos = start_pos + feature_len
                    logger.debug('Feature %s, strand %s, stop position %d',
                                 mrna.id, gene.strand, stop_pos)
                    break
                else:
                    remaining -= exon[1] - start_pos + 1
                    logger.debug('Feature %s, strand %s, remaining %d', mrna.id,
                                 gene.strand, remaining)
                    continue
            else:
                if exon[1] - exon[0] > remaining:
                    stop_pos = exon[0] + remaining
                    logger.debug('Feature %s, strand %s, stop position %d',
                                 mrna.id, gene.strand, stop_pos)
                    break
                else:
                    remaining -= exon[1] - exon[0] + 1
                    logger.debug('Feature %s, strand %s, remaining %d', mrna.id,
                                 gene.strand, remaining)
                    continue
    if not (start_pos and stop_pos):
        raise SequenceMappingError('Unable to determine start/end coordinates: start:%s end:%s' % (start_pos, stop_pos))
    start_pos -= padding
    stop_pos += padding
    return start_pos, stop_pos

def get_gff_info(db, seq_id, seq_range, pad_factor):
    gffinfo = {}
    mrna = db[seq_id]
    logger.debug('Feature mRNA %s: start=%s, stop=%s', seq_id, mrna.start,
                 mrna.stop)
    feature_id = seq_id
    logger.debug('Finding gene feature %s in the GFF DB', feature_id)
    errors = False
    while db[feature_id].featuretype != 'gene':
        logger.debug('recurse to find parent of %s', feature_id)
        if db[feature_id].attributes.get('Parent'):
            feature_id = db[feature_id].attributes.get('Parent')[0]
        else:
            logger.error('no parent for %s', feature_id)
            errors = True
            break
    if errors:
        logger.error('Problems processing feature %s', feature_id)
        return None
    gene = db[feature_id]
    gffinfo.update({'chr': gene.chrom})
    logger.info('Feature gene %s, chr:%s strand:%s start:%d stop:%d', feature_id, gene.chrom,
                gene.strand, gene.start, gene.stop)
    exons = []
    for f in db.children(feature_id):
        if f.featuretype == 'exon':
            exons.append((f.start, f.stop))
    if len(exons) == 0:
        logger.error('No exons found for feature %s', feature_id)
        return None
    exons = sorted(exons)
    logger.debug('Exon coordinates for feature %s: %s', feature_id, exons)
    try:
        start, stop = adjust_coords(seq_range, mrna, gene, exons, pad_factor)
    except SequenceMappingError as err:
        logger.error('Failed to adjust coordinates to genome: %s', err)
        return None
    logger.debug('Genome adusted coordinates for %s: %s:%s => %s:%s', seq_id, seq_range[0], seq_range[1], start, stop)
    gffinfo.update({'coords': (start, stop)})
    return gffinfo


def main():
    parser = argparse.ArgumentParser(description='GFF Extract')
    loglevel_group = parser.add_mutually_exclusive_group()
    loglevel_group.add_argument('--debug', action='store_true',
                                help='Produce debugging output')
    loglevel_group.add_argument('--verbose', action='store_true',
                                help='Produce more verbose output')
    input_group = parser.add_argument_group(title='Input Files',
                                            description='Required input files')
    input_group.add_argument('--gff', required=True,
                             help='Input GFF file')
    input_group.add_argument('--genome', required=True,
                             help='Input genome FASTA file')
    parser.add_argument('--seqid', required=True, type=str,
                        help='Sequence ID to map back into genome for extraction')
    parser.add_argument('--range', required=True,
                        help='start:end range from the sequence to map back into the genome for extraction')
    parser.add_argument('--outdir', default=os.getcwd(),
                        help='Output directory')
    parser.add_argument('--padding', default=PADDING_FACTOR, type=float,
                        help='Fraction of the alignment length to use as padding')
    args, extra = parser.parse_known_args()
    seq_range = list(set(args.range.split(':', 1)))
    if len(seq_range) != 2:
        logger.error('Invalid sequence range parameter value %s', args.range)
        sys.exit(1)
    for i in range(len(seq_range)):
        try:
            seq_range[i] = int(seq_range[i])
        except TypeError:
            logger.error('Invalid sequence range component value %s. Must be an integer.', e)
            sys.exit(1)
    seq_range = tuple(set(seq_range))
    if args.debug:
        logger.setLevel(logging.DEBUG)
    elif args.verbose:
        logger.setLevel(logging.INFO)
    logger.info('Starting')
    logger.debug('Args=%s extra=%s', args, extra)
    dbfile = os.path.splitext(os.path.basename(args.gff))[0] + '.db'
    dbfile = os.path.join(args.outdir, dbfile)
    if not os.path.exists(dbfile):
        logger.info('Initialize the database...')
        db = gffutils.create_db(args.gff, dbfile, keep_order=False, merge_strategy='create_unique', id_spec=['ID'], sort_attribute_values=False)
        db.execute('CREATE INDEX coordinate_index ON features(seqid, start, end)')
    else:
        logger.info('Loading genome features from %s', dbfile)
        db = gffutils.FeatureDB(dbfile)
    if not os.path.exists(args.genome):
        logger.error('Genome FASTA %s not found', args.genome)
        return 1
    if not db[args.seqid]:
        logger.error('Sequence ID %s not found in the GFF data', args.seqid)
    seq_gff_info = get_gff_info(db, args.seqid, seq_range, args.padding)
    if not seq_gff_info:
        logger.error('Failed to get GFF info for sequence ID %s', args.seqid)
        return None
    logger.info('Loading genome FASTA data')
    genomeSeqs = SeqIO.index(args.genome, 'fasta')
    chromRec = genomeSeqs[seq_gff_info['chr']]
    subseq = chromRec[seq_gff_info['coords'][0]:seq_gff_info['coords'][1]]
    subseq.id = '%s.A' % args.seqid
    logger.debug('%s subsequence [%s]:\n%s', args.seqid, seq_gff_info['coords'], subseq.format('fasta'))


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        logger.warning('Caught keyboard interrupt. Cleaning up and exiting early.')
        sys.exit(1)
