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
handler.setLevel(logging.DEBUG)
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
    logger.debug('Sequence length: %d', len(data.get('sequence')))
    return data

def make_query(seqid, qdb):
    query = get_sequence_data(seqid, qdb)
    return query

def parse_result(rtext):
    logger.debug('Got %d characters of data', len(rtext))
    data = []
    for line in rtext.splitlines():
        if line.startswith('#'):
            continue
        fields = [x.strip() for x in line.split('\t')]
        logger.debug('Got %d fields from line: %s', len(fields), line)
        if len(fields) != 12:
            logger.error('Wrong number of fields (%d != 12) for line:\n%s', len(fields), line)
            continue
        data.append(
          {
            'query_id': fields[0],
            'subject_id': fields[1],
            'pct_identity': fields[2],
            'alignment_length': fields[3],
            'mismatches': fields[4],
            'gaps': fields[5],
            'query_start': fields[6],
            'query_end': fields[7],
            'subject_start': fields[8],
            'subject_end': fields[9],
            'evalue': fields[10],
            'bit_score': fields[11]
          })

    logger.info('Final data:\n%s', json.dumps(data, indent=4))
    return data

def run_blast(blast, query, sdb):
    qfile = write_query_file(query)
    opts = [blast, '-evalue', '1e-50', '-num_threads', '4', '-outfmt', '7',
            '-query', qfile, '-db', sdb.name]
    '''
    blast -query results/PPC3.cds -db data/Alyrata/Alyrata_384_v2.1.protein.fa -out results/PPC-Alyrata_family.align -evalue 1e-50 -num_threads 4 -outfmt 7
    '''
    logger.debug('Running: %s', ' '.join(opts))
    result = subprocess.run(opts, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    os.unlink(qfile)
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
    data = {}
    parser = argparse.ArgumentParser(description='Filter BLAST+ results')
    parser.add_argument('-b', '--blast', default='blastx', type=str,
                        help='BLAST+ binary to use')
    parser.add_argument('-i', '--id', required=True, type=str,
                        help='Query sequence ID')
    parser.add_argument('-q', '--query-cds', required=True, type=argparse.FileType('r'),
                        help='FASTA formatted CDS DB containing the query sequence')
    parser.add_argument('-s', '--subject-cds', required=True, type=argparse.FileType('r'),
                        help='FASTA formatted CDS DB containing subject sequences.')
    parser.add_argument('-S', '--subject-prot', required=True, type=argparse.FileType('r'),
                        help='FASTA formatted protein DB to compare the query against.')
    args = parser.parse_args()
    query = make_query(args.id, args.query_cds)
    logger.debug('Query: %s', query)

    result = run_blast(args.blast, query, args.subject_prot)
    for aln in result:
        sid = aln.get('subject_id')
        subject = get_sequence_data(sid, args.subject_cds)
        qlen = len(query.get('sequence'))
        slen = len(subject.get('sequence'))
        coverage = float(slen) / float(qlen) * 100.0
        logger.info('Coverage for %s is %3f', subject.get('subject_id'), coverage)

    return 0




if __name__ == '__main__':
    sys.exit(main())
