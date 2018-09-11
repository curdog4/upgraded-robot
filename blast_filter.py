#!/usr/bin/env python3

import os
import sys
import json

import textwrap

import logging 
logger = logging.getLogger('BlastFilter')
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s %(levelname)s: %(funcName)s: %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)

def usage():
    print(textwrap.dedent("""
    Usage: %s <filename>

    Filename is the output from blast, using outformat 6 or 7
    """ % sys.argv[0]))
    return 0

def main():
    """Open the file argv[0], presume to be blast output format 7:
    Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
    """
    data = {}
    if len(sys.argv) != 2:
        usage()
        return 1
    if not os.path.exists(sys.argv[1]):
        logger.error('File %s does not exist', sys.argv[1])
        return 1
    fdata = None
    with open(sys.argv[1], 'r') as fd:
        fdata = fd.read()
    if not fdata:
        logger.error('No data read from file %s', sys.argv[1])
        return 1
    logger.debug('Got %d characters of data', len(fdata))
    for line in fdata.splitlines():
        if line.startswith('#'):
            continue
        fields = [x.strip() for x in line.split('\t')]
        logger.debug('Got %d fields from line: %s', len(fields), line)
        if not data.get(fields[1]):
            data.update({
              fields[1]: {
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
              }
            })
        else:
            logger.warning('Skipping entry for %s, already known: %s <-> %s', fields[1], fields, data.get(fields[1]))

    logger.info('Final data:\n%s', json.dumps(data, indent=4))


if __name__ == '__main__':
    sys.exit(main())
