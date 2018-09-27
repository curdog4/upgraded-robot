#!/usr/bin/env python3

import os
import sys
import json
import argparse
import tempfile
import textwrap
import subprocess

from Bio import SeqIO

import logging 
logger = logging.getLogger('BlastFilter')
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s %(levelname)s: %(funcName)s: %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)


def check_args(args):
    errors = 0
    if args.bitscore:
        if args.bitscore < 0.0 or args.bitscore > 100.0:
            logger.error('Bit score value must be a valid percent between 0.0 and 100.0')
            errors += 1
    if args.coverage:
        if args.coverage < 0.0 or args.coverage > 100.0:
            logger.error('Coverage must be a valid percent between 0.0 and 100.0')
            errors += 1
    if args.evalue:
        if args.evalue < 0.0 or args.evalue > 100.0:
            logger.error('E-value mest be a valid percent between 0.0 and 100.0')
            errors += 1
    if args.identity:
        if args.identity < 0.0 or args.identity > 100.0:
            logger.error('Identity must be a valid percent between 0.0 and 100.0')
            errors += 1
    if errors:
        return True
    return False

def check_bitscore(alignment, max_bitscore, max_variance):
    min_bitscore = max_bitscore * (max_variance / 100.0)
    if alignment.get('bit_score') < min_bitscore:
        return True
    return False
 
def check_coverage(alignment, query, min_coverage):
    # multiply alignment length, which is protein sequence length, by 3 to get codon sequence length
    alen = float(alignment.get('alignment_length') * 3)
    qlen = float(len(query.get('sequence')))
    logger.info('Alignment length:%d, query length:%d', alen, qlen)
    coverage = alen / qlen * 100.0
    coverage2 = alen / float(alignment.get('query_end')) * 100.0
    logger.info('Coverage comparison: %.1f vs. %.1f', coverage, coverage2)
    alignment.update({'coverage': coverage, 'coverage2': coverage2})
    if coverage < min_coverage:
        return coverage
    return False

def check_count(count, max_count):
    if count > max_count:
        return True
    return False

def check_evalue(alignment, min_evalue, max_variance):
    max_evalue = min_evalue * (1.0 + max_variance / 100.0)
    if alignment.get('evalue') > max_evalue:
        return True
    return False

def check_identity(alignment, pct_identity):
    if alignment.get('pct_identity') < pct_identity:
        return True
    return False

def get_sequence_data(dbh, seqid=None):
    data = []
    dbh.seek(0)
    for record in SeqIO.parse(dbh, 'fasta'):
        if seqid:
            if record.id == seqid:
                data.append({
                    'identity': record.id,
                    'sequence': '%s' % record.seq,
                    'text': record.format('fasta')
                })
                break
        else:
            data.append({
                'identity': record.id,
                'sequence': record.seq,
                'text': record.format('fasta')
            })
    return data

def make_query(seqid, qdb):
    query = get_sequence_data(qdb, seqid)[0]
    return query

def parse_result(rtext):
    data = []
    for line in rtext.splitlines():
        if line.startswith('#'):
            continue
        fields = [x.strip() for x in line.split('\t')]
        if len(fields) != 12:
            logger.error('Wrong number of fields (%d != 12) for line:\n%s', len(fields), line)
            continue
        data.append(
          {
            'query_id': fields[0],
            'subject_id': fields[1],
            'pct_identity': float(fields[2]),
            'alignment_length': int(fields[3]),
            'mismatches': int(fields[4]),
            'gaps': int(fields[5]),
            'query_start': int(fields[6]),
            'query_end': int(fields[7]),
            'subject_start': int(fields[8]),
            'subject_end': int(fields[9]),
            'evalue': float(fields[10]),
            'bit_score': float(fields[11])
          })

    return data

def run_blast():
    '''Now, just reading from STDIN
    '''
    rdata = sys.stdin.read()
    logger.info('BLAST results data:\n%s', rdata)
    qresult = parse_result(rdata)
    return qresult

def write_query_file(query):
    fd, qtemp = tempfile.mkstemp()
    with open(qtemp, 'w') as fh:
        fh.write(query.get('text'))
    return qtemp

def main():
    """Open the file argv[0], presume to be blast output format 7:
    Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
    """
    logger.info('Starting')
    data = {}
    parser = argparse.ArgumentParser(description='Filter BLAST+ results')
    parser.add_argument('-b', '--blast', default='blastx', type=str,
                        help='BLAST+ binary to use, can be absolute or relative path, or resolved via PATH')
    parser.add_argument('-i', '--id', required=False, type=str,
                        help='Query sequence ID')
    parser.add_argument('-f', '--query-file', required=False, type=argparse.FileType('r'),
                        help='FASTA formatted query sequence, in lieu of the query ID')
    parser.add_argument('-q', '--query-cds', required=True, type=argparse.FileType('r'),
                        help='FASTA formatted CDS DB containing the query sequence')
    parser.add_argument('-Q', '--query-prot', required=False, type=argparse.FileType('r'),
                        help='FASTA formatted PROT DB containing the query sequence')
    parser.add_argument('-s', '--subject-cds', required=True, type=argparse.FileType('r'),
                        help='FASTA formatted CDS DB containing subject sequences.')
    parser.add_argument('-S', '--subject-prot', required=True, type=argparse.FileType('r'),
                        help='FASTA formatted protein DB to compare the query against.')
    parser.add_argument('--bitscore', required=False, type=float,
                        help='Maximum variance of a bit score value as a percent of the best score')
    parser.add_argument('--count', required=False, type=int,
                        help='Maximum count a sequence ID can match')
    parser.add_argument('--coverage', required=False, type=float,
                        help='Minimum percent of query length an alignment must match')
    parser.add_argument('--evalue', required=False, type=float,
                        help='Maximum variance of an E-value score as a percent of the best score')
    parser.add_argument('--identity', required=False, type=float,
                        help='Minimum percent identity an alignment must match')

    args = parser.parse_args()
    logger.debug('Args: %s', args)
    if check_args(args):
        logger.error('There were errors processing command-line arguments. Aborting.')
        return 1
    clean_qfile = False
    if args.id:
        qfile = open(write_query_file(make_query(args.id, args.query_cds)), 'r')
        clean_qfile = True
        logger.debug('Query: %s', query)
    elif args.query_file:
        qfile = args.query_file
    else:
        logger.info('Running against entire query CDS database.')
        qfile = args.query_cds

    qdata = get_sequence_data(qfile)

    result = run_blast()
    for qry in qdata:
        for aln in result:
            logger.debug('Alignment data: %s', aln)
            sid = aln.get('subject_id')
            qid = aln.get('query_id')
            if qry['identity'] != qid:
                 logger.debug('Skip result not related to this query sequence ID: %s != %s', qry['identity'], qid)
                 continue
            result_subset = [x for x in result if x['query_id'] == qid]
            subject = get_sequence_data(args.subject_cds, sid)
            #if not 'query' in globals():
            #    query = get_sequence_data(args.query_cds, qid)
            logger.info('Query length is %d', len(qry.get('sequence')))
            if args.coverage:
                coverage_pct = check_coverage(aln, qry, args.coverage)
                logger.info('Coverage for %s is %3f', sid, aln.get('coverage'))
                if coverage_pct:
                    logger.warning('Dropping %s due to poor coverage (<%.1f%%): %.1f', sid, args.coverage, coverage_pct)
                    continue
            if args.count:
                subject_list = [x.get('subject_id') for x in result_subset]
                count = subject_list.count(sid)
                if check_count(count, args.count):
                    logger.warning('Dropping %s due to high incident count (>2): %d', sid, count)
                    continue
            if args.identity:
                if check_identity(aln, args.identity):
                    logger.warning('Dropping %s due to low identity (<60.0%%): %3f', sid, aln.get('pct_identity'))
                    continue
            if args.evalue:
                evalue_scores = list(set([x.get('evalue') for x in result_subset]))
                evalue_scores.sort()
                if evalue_scores[0] == 0:
                    evalue_scores.pop(0)
                logger.debug('E-value scores: %s', evalue_scores)
                min_evalue = evalue_scores[0]
                if check_evalue(aln, min_evalue, args.evalue):
                    logger.warning('Dropping %s due to poor e-value (>5%% of highest e-value score %.1e): %.1e', sid, min_evalue, aln.get('evalue'))
                    continue
            if args.bitscore:
                bit_scores = list(set([x.get('bit_score') for x in result_subset]))
                bit_scores.sort()
                max_bitscore = bit_scores[-1]
                if check_bitscore(aln, max_bitscore, args.bitscore):
                    logger.warning('Dropping %s due to low bit score (<5%% of highest bit score %d): %d', sid, max_bitscore, aln.get('bit_score'))
                    continue

            if data.get(sid):
                data[sid].append(aln)
            else:
                data.update({sid: [aln]})

    if clean_qfile:
        os.unlink(qfile)
    logger.info('Filtered data, %d records', len(data))

    # need to output as blast outfmt 6
    #AL8G30230.t1	LOC_Os08g08840.1	83.59	323	49	1	187	1155	68	386	9e-160	459
    for seq in data.values():
        for rec in seq:
            logger.info('Record: %s', json.dumps(rec))
            try:
                outline = '{query_id:s}\t{subject_id:s}\t{pct_identity:.1f}\t{alignment_length:d}\t{mismatches:d}\t{gaps:d}\t{query_start:d}\t{query_end:d}\t{subject_start:d}\t{subject_end:d}\t{evalue:.0e}\t{bit_score:g}\n'.format(**rec)
            except ValueError:
                outline = '{query_id:s}\t{subject_id:d}\t{pct_identity:.1f}\t{alignment_length:d}\t{mismatches:d}\t{gaps:d}\t{query_start:d}\t{query_end:d}\t{subject_start:d}\t{subject_end:d}\t{evalue:.0e}\t{bit_score:g}\n'.format(**rec)
            sys.stdout.write(outline)
    logger.info('Complete')
    return 0


if __name__ == '__main__':
    sys.exit(main())
