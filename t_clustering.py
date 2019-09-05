'''Test use of kmeans clustering to determine chimeras in blastx table results
Non-standard BLAST+ Table Field Values:
0-qseqid 1-qlen 2-sseqid 3-slen 4-qframe 5-pident 6-nident 7-length
8-mismatch 9-gapopen 10-qstart 11-qend 12-sstart 13-send 14-evalue
15-bitscore

LAST Tab Table Field Values:
0-score 1-name1 2-start1 3-alnSize1 4-strand1 5-seqSize1 6-name2
7-start2 8-alnSize2 9-strand2 10-seqSize2 11-blocks

LAST BlastTab Table Field Values:
0-query id 1-subject id 2-% identity 3-alignment length 4-mismatches
5- gap opens 6-q. start 7-q. end 8-s. start 9-s. end 10-evalue 11-bit score

LAST BlastTab+ Table Field Values:
0-query id 1-subject id 2-% identity 3-alignment length 4-mismatches
5-gap opens 6-q. start 7-q. end 8-s. start 9-s. end 10-evalue 11-bit score
12-query length 13-subject length

Logic borrowed from:
https://blog.cambridgespark.com/how-to-determine-the-optimal-number-of-clusters-for-k-means-clustering-14f27070048f
https://jtemporal.com/kmeans-and-elbow-method/

'''

import concurrent.futures
import os
import sys
import argparse
import glob
import io
import logging
import logging.config
import json
import math
import time

import pandas as pd
from sklearn.preprocessing import MinMaxScaler
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt

logConfFile = 'logging.json'
with open(logConfFile, 'r') as fd: 
    logConfContent = fd.read()
logConfData = json.loads(logConfContent)
logging.config.dictConfig(logConfData)
logger = logging.getLogger('DetectChimeras')

BITSCORE_RATIO = 2.0
CCOV_THRESHOLD = 0.60
LENGTH_THRESHOLD = 100
PIDENT_THRESHOLD = 20.0
QCOV_THRESHOLD = 0.30
SCORE_RATIO = 1.0
SCOV_THRESHOLD = 0.70

def calculate_data(data):
    wcss = []
    for n in range(1, min(len(data), 15)):
        kmeans = KMeans(n_clusters=n)
        kmeans.fit(X=data)
        wcss.append(kmeans.inertia_)
    return wcss

def get_blastx_table_data(fname, ftype, args):
    with open(fname, 'r') as fdesc:
        fdata = fdesc.read()
    coordinates = []
    lines = fdata.splitlines()
    hsp_cnt = len(lines)
    for line in lines:
        if line.startswith('#'):
            logger.debug('Skipping non-HSP line [%d/%d]: %s', lines.index(line) + 1, hsp_cnt, line)
            continue
        fields = line.split('\t')
        # Field value conversions to int or float
        int_fields = [1, 3, 7, 10, 11, 12, 13]
        float_fields = [5, 15]
        if ftype == 1:  # LAST TAB
            int_fields = [0, 2, 3, 5, 7, 8, 10]
            float_fields = []
        elif ftype == 2:  # LAST BlastTab
            int_fields = [3, 4, 5, 6, 7, 8, 9]
            float_fields = [2, 11]
        elif ftype == 3:  # LAST BlastTab+
            int_fields = [3, 4, 5, 6, 7, 8, 9, 12, 13]
            float_fields = [2, 11]
        for idx in int_fields:
            fields[idx] = int(fields[idx])
        for idx in float_fields:
            fields[idx] = float(fields[idx])

        if is_fit(fields, ftype, args):
            coordinates.append(fields)
        else:
            logger.debug('Discarding unfit HSP [%d/%d]: %s', lines.index(line) + 1, hsp_cnt, line)
    return coordinates


def is_fit(fields, ftype, args):
    '''Checks for 'fitness' of the result
    '''
    fit = True
    _scov = scov(fields, ftype)
    _slen = slen(fields, ftype)
    _qcov = qcov(fields, ftype)
    _qlen = qlen(fields, ftype)
    if ftype == 0:
        # This is if the file format is BLAST+ Custom
        if fields[1] < 3 * fields[3] :
            # subject longer than query
            logger.debug('HSP subject longer than query: %s: %d < 3 * %d', fields[2], fields[1], fields[3])
            fit = False
        if fields[5] < args.pident:
            # insufficient identity
            logger.debug('HSP has insufficient percent identity: %.3f < %.3f', fields[5], args.pident)
            fit = False
        if float(fields[3]) / fields[15] < args.bitscore_ratio:
            # insufficient bitscore value
            logger.debug('HSP has insufficient bitscore ratio: %.3f / %.3f < %.3f', float(fields[3]), fields[15], args.bitscore_ratio)
            fit = False
    elif ftype == 1:
        # This is if the file format is 1 (LAST TAB)
        if fields[10] < 3 * fields[5]:
            # subject longer than query
            logger.debug('HSP subject longer than query: %s: %d < 3 * %d', fields[1], fields[10], fields[5])
            fit = False
        # Cannot perform this check
        #if fields[] < args.pident:
        #    # insufficent identity
        #    fit = False
        if float(fields[0]) / float(fields[5]) < args.score_ratio:
            # insufficient score value
            logger.debug('HSP has insufficient score ratio: %.3f / %.3f < %.3f', float(fields[0]), fields[5], args.score_ratio)
            fit = False
    elif ftype == 2:
        # This is if the file format is 2 (LAST BlastTab)
        # Cannot perform this check
        #if fields[] < 3 * fields[]:
        #    # subject longer than query
        #    fit = False
        if fields[2] < args.pident:
            # insufficient identity
            logger.debug('HSP has insufficient percent identity: %.3f < %.3f', fields[2], args.pident)
            fit = False
        # Cannot perform this check
        #if float(fields[]) / fields[11] < args.bitscore_ratio:
        #    # insufficient bitscore value
        #    logger.debug('HSP has insufficient bitscore ratio: %.3f / %.3f < %.3f', float(fields[]), fields[11], args.bitscore_ratio)
        #    fit = False
    elif ftype == 3:
        # This is if the file format is 3 (LAST BlastTab+)
        if fields[12] < 3 * fields[13]:
            # subject longer than query
            logger.debug('HSP subject longer than query: %s: %d < 3 * %d', fields[1], fields[12], fields[13])
            fit = False
        if fields[2] < args.pident:
            # insufficient identity
            logger.debug('HSP has insufficient percent identity: %.3f < %.3f', fields[2], args.pident)
            fit = False
        #if float(fields[13]) / fields[11] < args.bitscore_ratio:
        #    # insufficient bitscore value
        #    logger.debug('HSP has insufficient bitscore ratio: %.3f / %.3f < %.3f', float(fields[13]), fields[11], args.bitscore_ratio)
        #    fit = False
        if float(fields[14]) / float(_slen) < args.score_ratio:
            # insufficient score value
            logger.debug('HSP has insufficient score ratio: %.3f / %.3f < %.3f', float(fields[14]), _slen, args.score_ratio)
            fit = False
    else:
        logger.error('Cannot determine fitness of HSP from table of type %s', ftype)
        fit = False

    if _qlen < args.qlen:
        # insufficient alignment length
        logger.debug('HSP has insufficient query alignment length: %d < %d', _qlen, args.qlen)
        fit = False
    if _qcov < args.qcov:
        # insufficient query coverage
        logger.debug('HSP has insufficient query coverage: %.3f < %.3f', _qcov, args.qcov)
        fit = False
    if _scov < args.scov:
        # insufficient subject coverage
        logger.debug('HSP has insufficient subject coverage: %.3f < %.3f', _scov, args.scov)
        fit = False
    return fit

def optimal_number_of_clusters(wcss):
    x1, y1 = 2, wcss[0]
    x2, y2 = 20, wcss[len(wcss)-1]

    distances = []
    for i in range(1, len(wcss) - 1):
        x0 = i + 1
        y0 = wcss[i]
        numerator = abs((y2-y1)*x0 - (x2-x1)*y0 + x2*y1 - y2*x1)
        denominator = math.sqrt((y2 - y1)**2 + (x2 - x1)**2)
        d = numerator / denominator
        #logger.info('Distance at (%d, %d): %f', x0, y0, d)
        distances.append(d)
    if distances:
        return distances.index(max(distances)) + 2
    return None

def overlap(start1, end1, start2, end2):
    o_start = max(start1, start2)
    o_end = min(end1, end2)
    o_len = abs(o_end - o_start)
    t_len = max(start1, start2, end1, end2) - min(start1, start2, end1, end2)
    o_fraction = o_len / t_len
    return o_fraction

def qcov(hsp, ftype):
    if ftype == 0:  # Blast+ Custom
        return float(qlen(hsp, ftype)) / float(hsp[1])
    elif ftype == 1:  # LAST TAB
        return float(qlen(hsp, ftype)) / float(hsp[10])
    elif ftype == 2:  # LAST BlastTab
        # Cannot perform this check
        return 0.0
    elif ftype == 3:  # LAST BlastTab+
        return float(qlen(hsp, ftype)) / float(hsp[12])
    return None

def qlen(hsp, ftype):
    if ftype == 0:  # Blast+ Custom
        return abs(hsp[11] - hsp[10]) + 1
    elif ftype == 1:  # LAST TAB
        return hsp[8]
    elif ftype in (2, 3):  # LAST BlastTab and BlastTab+
        return abs(hsp[7] - hsp[6]) + 1
    return None

def scov(hsp, ftype):
    if ftype == 0:  # BLAST+ Custom
        return float(slen(hsp, ftype)) / float(hsp[3])
    elif ftype == 1:  # Last TAB
        return float(slen(hsp, ftype)) / float(hsp[5])
    elif ftype == 2:  # LAST BlastTab
        # Cannot perform this check
        return 0.0
    elif ftype == 3:  # LAST BlastTab+
        return float(slen(hsp, ftype)) / float(hsp[13])
    return None

def slen(hsp, ftype):
    if ftype == 0:  # BLAST+ Custom
        return abs(hsp[13] - hsp[12]) + 1
    elif ftype == 1:  # LAST TAB
        return hsp[3]
    elif ftype in (2, 3):  # LAST BlastTab and BlastTab+
        return abs(hsp[9] - hsp[8]) + 1
    return None

def separated(hsp1, hsp2):
    ''' given two hsps, return True if
    overlap less than the smaller of:
    * 20% of the shorter length
    * 60 bp
    '''
    length1 = abs(hsp1[1] - hsp1[0])
    length2 = abs(hsp2[1] - hsp2[0])
    start = min(int(hsp1[0]), int(hsp1[1]), int(hsp2[0]), int(hsp2[1]))
    end = max(int(hsp1[0]), int(hsp1[1]), int(hsp2[0]), int(hsp2[1]))
    overlap = length1 + length2 - (end - start) + 1
    # value of overlap can < 0 but only the upper limit maters
    if overlap < min(60, 0.2 * min(length1, length2)):
        return True
    return False

def main():
    ftype_list = ['custom', 'lasttab', 'blasttab', 'blasttab+']
    parser = argparse.ArgumentParser(description='Test Clustering')
    loglevel_group = parser.add_mutually_exclusive_group()
    loglevel_group.add_argument('--debug', action='store_true',
                                help='Produce debugging output')
    loglevel_group.add_argument('--verbose', action='store_true',
                                help='Produce more verbose output')
    input_group = parser.add_mutually_exclusive_group()
    input_group.add_argument('--file', type=str,
                             help='Process single input file')
    input_group.add_argument('--directory', type=str,
                             help='Process *.tbl in this directory')
    filter_group = parser.add_argument_group(title='Filtering options',
                                             description='Options controlling HSP result filtering')
    filter_group.add_argument('--bitscore-ratio', default=BITSCORE_RATIO, type=float,
                              help='Minimum ratio of subject sequence length to bitscore (for BLAST* format results)')
    filter_group.add_argument('--ccov', default=CCOV_THRESHOLD, type=float,
                              help='Maximum ratio of chimeric range length to query sequence length')
    filter_group.add_argument('--ftype', choices=ftype_list, required=True,
                         help='Input file/s table format')
    filter_group.add_argument('--pident', default=PIDENT_THRESHOLD, type=float,
                              help='Minimum identity threshold for HSP')
    filter_group.add_argument('--qcov', default=QCOV_THRESHOLD, type=float,
                              help='Minimum query coverage threshold for HSP')
    filter_group.add_argument('--qlen', default=LENGTH_THRESHOLD, type=int,
                              help='Minimum query alignment sequence length for HSP')
    filter_group.add_argument('--score-ratio', default=SCORE_RATIO, type=float,
                              help='Minimum ratio of score to subject sequence length (for LAST TAB format results)')
    filter_group.add_argument('--scov', default=SCOV_THRESHOLD, type=float,
                              help='Minimum subject coverage threshold for HSP')
    args, extra = parser.parse_known_args()
    if args.debug:
        logger.setLevel(logging.DEBUG)
    elif args.verbose:
        logger.setLevel(logging.INFO)
    file_list = []
    if args.file:
        if not os.path.exists(args.file):
            raise OSError('File %s not found' % args.file)
        file_list = [args.file]
    elif args.directory:
        if not os.path.exists(args.directory):
            raise OSError('Directory %s not found' % args.directory)
        file_list = glob.glob(os.path.join(args.directory, '*.tbl'))
    if not file_list:
        raise OSError('No files found to process')
    ftype = ftype_list.index(args.ftype)
    t_start = time.time()
    logger.info('Begin')
    coordinate_map = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=12) as executor:
        future_map = {executor.submit(get_blastx_table_data, fname, ftype, args): fname for fname in file_list}
        for future in concurrent.futures.as_completed(future_map):
            fname = future_map[future]
            try:
                data = future.result()
            except Exception as exc:
                logger.error('%r generated an exception: %s', fname, exc)
            else:
                if data:
                    logger.info('%r table returned %d coordinates', fname, len(data))
                    coordinate_map[fname] = data
    logger.info('Processing %d results', len(coordinate_map))

    headers = 'start\tend'
    cntr = 0
    for label, coordinates in coordinate_map.items():
        cntr += 1
        logger.info('Percent completion: %.3f', float(cntr) / float(len(coordinate_map)) * 100.0)
        if len(coordinates) < 3:
            logger.error('Insufficient records to continue for %s', label)
            continue

        sio = io.StringIO()
        sio.write(headers + '\n')
        for coords in coordinates:
            if ftype == 0:
                sio.write('%s\t%s\n' % (coords[10], coords[11]))
            elif ftype == 1:
                sio.write('%s\t%s\n' % (coords[7], coords[7] + coords[8] - 1))
            elif ftype in (2, 3):
                sio.write('%s\t%s\n' % (coords[6], coords[7]))
        sio.seek(0)
        #logger.debug('Raw table data:\n%s', sio.read())
        sio.seek(0)
        data = pd.read_csv(sio, sep='\t')
        #logger.debug('Data:\n%s', data)

        mms = MinMaxScaler()
        mms.fit(data)
        data_transformed = mms.transform(data)
        #logger.debug('Data Xform:\n%s', data_transformed)

        wcss = []
        K = range(1, min(len(coordinates), 15))
        for k in K:
            km = KMeans(n_clusters=k)
            km = km.fit(data_transformed)
            wcss.append(km.inertia_)

        n = optimal_number_of_clusters(wcss)
        logger.info('Predicted optimal number of clusters for %s: %s', label, n)
        if not n:
            logger.error('Number of clusters could not be calculated for %s', label)
            continue

        '''
        plt.plot(K, wcss, 'bx-')
        plt.xlabel('k')
        plt.ylabel('wcss')
        plt.title('Elbow Method For Optimal k (%s)' % label)
        plt.show()
        '''
        km = KMeans(n_clusters=n)
        clusters = km.fit_predict(data_transformed)
        #logger.debug('Clusters:\n%s', clusters)
        centroids = []
        for i in range(n):
            centroids.append([])
        for i in range(len(data.index)):
            row = data.iloc[i]
            cn =  clusters[i]
            #logger.debug('For row %d (%s, %s) the cluster number is %d',
            #             i, row[0], row[1], cn)
            centroids[cn].append((row[0], row[1]))
        #logger.debug('Cluster sorted rows:\n%s', cm)
        for i in range(len(centroids)):
            centroid = centroids[i]
            if not centroid:
                logger.warning('Empty cluster for %s: %s', label, centroid)
                continue
            start = min([x[0] for x in centroid])
            end = max([x[1] for x in centroid])
            logger.info('Possible chimeric gene in region from %s to %s',
                        start, end)
            clen = abs(end - start) + 1
            if ftype == 0:
                qlen = coordinates[0][1]
            elif ftype == 1:
                qlen = coordinates[0][10]
            elif ftype == 3:
                qlen = coordinates[0][12]

            if qlen:
                ccov = float(clen) / float(qlen)
                logger.info('Chimeric sequence covers %.3f percent of the query sequence %s', ccov * 100.0, label)
                if ccov > args.ccov:
                    logger.warning('Dropping range (%d, %d), too long to be likely chimeric', start, end)
                    continue
            for j in range(i+1, len(centroids)):
                _centroid = centroids[j]
                if not _centroid:
                    logger.warning('Empty cluster for %s: %s. Cannot determine overlap.',
                                   label, _centroid)
                    continue
                _start = min([x[0] for x in _centroid])
                _end = max([x[1] for x in _centroid])
                logger.info('Overlap for (%d, %d) and (%d, %d): %.3f', start, end, _start, _end,
                            overlap(start, end, _start, _end))
                if not separated((start, end), (_start, _end)):
                    logger.info('Ranges (%s, %s) and (%s, %s) not separated, should merge',
                                start, end, _start, _end)

    t_end = time.time()
    logger.info('Complete. Elapsed %.3f seconds.', t_end - t_start)
    return 0

if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        logger.info('Received keyboard interrupt. Aborting.')
        sys.exit(1)
