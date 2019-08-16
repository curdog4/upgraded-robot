'''Try and detect self and multigene chimeric sequences
'''

import os
import sys
import glob
import logging
import logging.config
import json

logConfFile = 'logging.json'
with open(logConfFile, 'r') as fd:
    logConfContent = fd.read()
logConfData = json.loads(logConfContent)
logging.config.dictConfig(logConfData)
logger = logging.getLogger('DetectChimeras')

DATA_SRC = glob.glob('seqs/*_blastx.tbl')
COVERAGE_THRESHOLD = 0.30
IDENTITY_THRESHOLD = 0.20
MINIMUM_LENGTH = 100

def qcov(hsp):
    return abs(hsp[1] - hsp[0]) + 1

def separated(s1, s2):
    ''' given two hsps, return True if
    overlap les than 20% of the shorter and overlap less than 60 bp
    '''
    length1 = qcov(s1)
    length2 = qcov(s2)
    start = min(s1[0], s1[1], s2[0], s2[1])
    end = max(s1[0], s1[1], s2[0], s2[1])
    overlap = length1 + length2 - (end - start) + 1
    # value of overlap can < 0 but only the upper limit maters
    if overlap < min(60, 0.2 * min(length1, length2)):
        return True
    return False

def main():
    '''Fields in each row:
    0-qseqid 1-qlen 2-sseqid 3-slen 4-qframe 5-pident 6-nident 7-length 8-mismatch
    9-gapopen 10-qstart 11-qend 12-sstart 13-send 14-evalue 16-bitscore
    '''
    for fname in DATA_SRC:
        with open(fname, 'r') as fd:
            data = fd.read()
        blast_map = {}
        for line in data.splitlines():
            line = line.strip()
            hsp = line.split('\t')
            for i in (5, 10, 11):
                hsp[i] = float(hsp[i])

            start_end = (hsp[10], hsp[11])
            if hsp[10] > hsp[11]:
                start_end = (hsp[11], hsp[10])
            if hsp[5] < (IDENTITY_THRESHOLD * 100.0):
                logger.warning('Dropping HSP for %s due to low identity: %f < %f',
                               subject, hsp[5], IDENTITY_THRESHOLD * 100.0)
                continue
            elif qcov(start_end) < MINIMUM_LENGTH:
                logger.warning('Dropping HSP for %s due to insufficient length: %d < %d',
                               subject, qcov(start_end), MINIMUM_LENGTH)
                continue  # ignore low similarity or short HSPs
            query, subject = hsp[0], hsp[2]

            blast_map.setdefault(query, {'multigene': []})
            blast_map[query].setdefault(subject, [])
            if blast_map[query][subject]:
                blast_map[query][subject].append(start_end)
            else:
                blast_map[query][subject] = [start_end]
            blast_map[query]['multigene'].append(start_end)

        for qid in list(blast_map.keys())[:]:
            for sid in list(blast_map[qid].keys())[:]:
                if qid == 'multigene':
                    continue
                if len(blast_map[qid][sid]) < 2:
                    blast_map[qid].pop(sid)
            if len(blast_map[qid]['multigene']) < 2:
                if len(blast_map[qid].keys()) < 3:
                    blast_map.pop(qid)
                else:
                    blast_map[qid].pop('multigene')

        #logger.debug('BLAST Map:\n%s', json.dumps(blast_map, indent=4))
        for query in list(blast_map.keys())[:]:
            subjects = list(set(list(blast_map[query].keys())[:]) - set(['multigene']))
            for subject in subjects:
                newlist = []
                logger.info('Looking at %s <-> %s HSPs for possible self-chimeric sequences', query, subject)
                rangesets = blast_map[query][subject][:]
                for i in range(len(rangesets)):
                    distinct = True
                    for j in range(i+1, len(rangesets)):
                        if not separated(rangesets[i], rangesets[j]):
                            start = min(rangesets[i][0], rangesets[j][0])
                            end = max(rangesets[i][1], rangesets[j][1])
                            distinct = False
                            rangesets[i] = (start, end)
                    if distinct:
                        newlist.append(rangesets[i])
                    else:
                        for j in range(len(newlist)):
                            if not separated(rangesets[i], newlist[j]):
                                start = min(rangesets[i][0], newlist[j][0])
                                end = max(rangesets[i][1], newlist[j][1])
                                newlist[j] = (start, end)
                blast_map[query][subject] = newlist
            logger.info('Looking at possible %s multigene chimeric sequences', query)
            rangesets = blast_map[query]['multigene'][:]
            newlist = []
            for i in range(len(rangesets)):
                distinct = True
                for j in range(i+1, len(rangesets)):
                    if not separated(rangesets[i], rangesets[j]):
                        start = min(rangesets[i][0], rangesets[j][0])
                        end = max(rangesets[i][1], rangesets[j][1])
                        distinct = False
                        rangesets[i] = (start, end)
                if not distinct:
                    new = True
                    for j in range(len(newlist)):
                        if not separated(rangesets[i], newlist[j]):
                            new = False
                            start = min(rangesets[i][0], newlist[j][0])
                            end = max(rangesets[i][1], newlist[j][1])
                            newlist[j] = (start, end)
                    if new:
                        newlist.append(rangesets[i])
            blast_map[query]['multigene'] = newlist
            if len(blast_map[query]['multigene']) < 2:
                blast_map.pop(query)
        logger.debug('Reduced BLAST Map:\n%s', json.dumps(blast_map, indent=4))
    return 0


if __name__ == '__main__':
    sys.exit(main())
