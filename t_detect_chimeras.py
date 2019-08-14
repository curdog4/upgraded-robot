'''Try and detect self and multigene chimeric sequences
'''

import os
import sys
import logging

logConfFile = 'logging.json'
with open(logConfFile, 'r') as fd:
    logConfContent = fd.read()
logConfData = json.loads(logConfContent)
logging.config.dictConfig(logConfData)
logger = logging.getLogger('DetectChimeric')

DATA_SRC = 'seqs/Mecry_04G114660.1_blastx.tbl'
COVERAGE_THRESHOLD = 0.60
IDENTITY_THRESHOLD = 0.40
MINIMUM_LENGTH = 100

def qcov(hsp):
    return abs(int(hsp[11]) - int(hsp[10])) + 1

def main():
    with open(DATA_SRC, 'r') as fd:
        data = fd.read()
    blast_map = {}
    for line in data.splitlines():
        line = line.strip()
        hsp = line.split('\t')
        for i in (5, 10, 11):
            hsp[i] = float(hsp[i])
        query, hit = hsp[0], hsp[2]

        blast_map.set_defaults(query, {})
        blast_map[query].set_defaults(hit, {})
        if blast_map[query][hit]:

    return 0


if __name__ == '__main__':
    sys.exit(main())
