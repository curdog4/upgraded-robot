#!/usr/bin/env python3

import os
import sys
import json
import argparse
import tempfile
import textwrap
import subprocess

import logging 
logger = logging.getLogger('BlastFilter')
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s %(levelname)s: %(funcName)s: %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)


def get_sequence_data(seqid, dbh):
    data = {}
    dbh.seek(0)
    fdata = dbh.read()
    records = fdata.split('>')
    for record in records:
        if record.startswith(seqid):
            lines = record.splitlines()
            data = {
                'identity': lines[0],
                'sequence': ''.join(lines[1:]),
                'text': '>' + record,
            }
            break
    return data

def make_query(seqid, qdb):
    query = get_sequence_data(seqid, qdb)
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
            'bit_score': int(fields[11])
          })

    return data

def run_blast(blast, qfile, sdb):
    opts = [blast, '-evalue', '1e-50', '-num_threads', '15', '-outfmt', '7',
            '-query', qfile, '-db', sdb.name]
    '''
    blast -query results/PPC3.cds -db data/Alyrata/Alyrata_384_v2.1.protein.fa -out results/PPC-Alyrata_family.align -evalue 1e-50 -num_threads 4 -outfmt 7
    '''
    logger.debug('Running: %s', ' '.join(opts))
    result = subprocess.run(opts, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0:
        logger.error('BLAST %s returned non-zero %d: %s', blast, result.returncode, result.stderr.decode())
        sys.exit(1)
    rdata = result.stdout.decode()
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
    args = parser.parse_args()
    logger.debug('Args: %s', args)
    clean_qfile = False
    query = None
    if args.id:
        query = make_query(args.id, args.query_cds)
        qfile = write_query_file(query)
        clean_qfile = True
        logger.debug('Query: %s', query)
    elif args.query_file:
        qfile = args.query_file.name
        logger.warning('Run using a query file not yet implemented...')
        #return 1
    else:
        logger.info('Running against entire query CDS database.')
        qfile = args.query_cds.name
        logger.warning('Run against entire query CDS database not yet implemented...')
        #return 1

    result = run_blast(args.blast, qfile, args.subject_prot)
    if clean_qfile:
        os.unlink(qfile)
    for aln in result:
        logger.debug('Alignment data: %s', aln)
        sid = aln.get('subject_id')
        if not sid.endswith('1'):
            logger.warning('Dropping %s due to non-primary subject sequence', sid)
            continue
        qid = aln.get('query_id')
        if not qid.endswith('1'):
            logger.warning('Dropping %s due to non-primary query sequence', qid)
            continue
        result_subset = [x for x in result if x['query_id'] == qid]
        subject_list = [x.get('subject_id') for x in result_subset]
        evalue_scores = list(set([x.get('evalue') for x in result_subset]))
        evalue_scores.sort()
        bit_scores = list(set([x.get('bit_score') for x in result_subset]))
        bit_scores.sort()
        if evalue_scores[0] == 0:
            evalue_scores.pop(0)
        logger.debug('E-value scores: %s', evalue_scores)
        subject = get_sequence_data(sid, args.subject_cds)
        if not 'query' in globals():
            query = get_sequence_data(qid, args.query_cds)
        qlen = len(query.get('sequence'))
        slen = len(subject.get('sequence'))
        coverage = float(slen) / float(qlen) * 100.0
        aln.update({'coverage': coverage})
        count = subject_list.count(sid)
        logger.info('Coverage for %s is %3f', sid, coverage)
        if coverage < 80.0:
            logger.warning('Dropping %s due to low coverage (<80.0%%): %3f', sid, coverage)
            continue
        if count > 2:
            logger.warning('Dropping %s due to high incident count (>2): %d', sid, count)
            continue
        if aln.get('pct_identity') < 60.0:
            logger.warning('Dropping %s due to low identity (<60.0%%): %3f', sid, aln.get('pct_identity'))
            continue
        max_evalue = evalue_scores[0]
        if aln.get('evalue') > 0.95 * max_evalue:
            logger.warning('Dropping %s due to low e-value (>5%% of highest e-value score %3f): %3f', sid, max_evalue, aln.get('evalue'))
            continue
        max_bitscore = bit_scores[-1]
        if aln.get('bit_score') < 0.95 * max_bitscore:
            logger.warning('Dropping %s due to low bit score (<5%% of highest bit score %d): %d', sid, max_bitscore, aln.get('bit_score'))
            continue

        if data.get(sid):
            data[sid].append(aln)
        else:
            data.update({sid: [aln]})

    logger.info('Filtered data:\n%s', json.dumps(data, indent=4))

    logger.info('Complete')
    return 0




if __name__ == '__main__':
    sys.exit(main())
