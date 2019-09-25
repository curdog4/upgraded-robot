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

Exonerate Custom formatting:
%ti\t%qi\t%pi\t%em\t%tab\t%tae\t%tal\t%qab\t%qae\t%qal\t%tl\t%ql\t%s
0-target id 1-query id 2-% identity 3-mismatches 4-t.start 5-t.end 6-t.length
7-q.start 8-q.end 9-q.length 10-target length 11-query length 12-raw score

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
import pprint
import time

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastxCommandline
from ckmeans import ckmeans
import gffutils
import numpy as np
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
CCOV_THRESHOLD = 0.80
REFDB_FILE = '/data/gpfs/home/dcurdie/scratch/UniProtKB/uniprot_sprot_plants.fa'
LENGTH_THRESHOLD = 100
NUM_THREADS = 12
PIDENT_THRESHOLD = 20.0
QCOV_THRESHOLD = 0.30
SCORE_RATIO = 1.0
SCOV_THRESHOLD = 0.70


def calculate_clusters_ckmeans(label, ftype, aln_hsps):
    if len(coordinates) < 3:
        logger.error('Insufficient records to continue for %s', label)
        return None
    _start_idx = 10
    _end_idx = 11 
    if ftype == 3:
        _start_idx = 6
        _end_idx = 7
    elif ftype == 4:
        _start_idx = 4
        _end_idx = 5
    aln_ends = []
    aln_weights = []
    logger.info('Parsing %s data into lists', label)
    for hsp in aln_hsps:
        _end = hsp[_end_idx]
        _weight = float(qlen(hsp, ftype)) / float(hsp[1])
        if ftype == 3:  # LAST BlastTab+
            _weight = float(hsp[3]) / float(hsp[12])
        elif ftype == 4:  # Exonerate custom
            _weight = hsp[6] / hsp[10]
        aln_ends.append(_end)
        aln_weights.append(_weight)
    logger.info('Instantiating numpy arrays for %s', label)
    array_ends = np.array(aln_ends)
    array_weights = np.array(aln_weights)
    logger.info('Calculating ckmeans 1d clustering for %S', label)
    #clusters = ckmeans(array_ends, k=(2, 7), weights=array_weights)
    _min_k = len(list(set(aln_ends)))
    if _min_k < 2:
        logger.error('Unable to cluster %s with %d unique values', label, _min_k)
        return None
    elif _min_k > 3:
        k = (2, 3)
    else:
        k = 2
    clusters = ckmeans(array_ends, k=k)
    logger.info('Ckmeans result for %s: number of clusters:%d, centroids:%s', label, len(clusters.centers), clusters.centers)
    logger.info('Processing results for %s', label)
    centroids = []
    for c_idx in range(len(clusters.centers)):
        centroids.append([])
        _members = []
        for m_idx in range(len(clusters.clustering)):
            if clusters.clustering[m_idx] == c_idx:
                _members.append(aln_hsps[m_idx])
        _start = min([int(x[_start_idx]) for x in _members])
        _end = max([int(x[_end_idx]) for x in _members])
    return centroids


def calculate_clusters_kmeans(label, ftype, coordinates):
    if len(coordinates) < 3:
        logger.error('Insufficient records to continue for %s', label)
        return None

    headers = 'start\tend'
    sio = io.StringIO()
    sio.write(headers + '\n')
    for coords in coordinates:
        if ftype == 0:
            sio.write('%s\t%s\n' % (coords[10], coords[11]))
        elif ftype == 1:
            sio.write('%s\t%s\n' % (coords[7], coords[7] + coords[8] - 1))
        elif ftype in (2, 3):
            sio.write('%s\t%s\n' % (coords[6], coords[7]))
        elif ftype == 4:
            sio.write('%s\t%s\n' % (coords[4], coords[5]))
    sio.seek(0)
    #logger.debug('Raw table data for %s:\n%s', label, sio.read())
    sio.seek(0)
    data = pd.read_csv(sio, sep='\t')
    #logger.debug('Data for %s:\n%s', label, data)

    mms = MinMaxScaler()
    mms.fit(data)
    data_transformed = mms.transform(data)
    #logger.debug('Data Xform for %s:\n%s', label, data_transformed)

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
        return None

    km = KMeans(n_clusters=n)
    clusters = km.fit_predict(data_transformed)
    #logger.debug('Clusters for %s:\n%s', label, clusters)
    centroids = []
    for i in range(n):
        centroids.append([])
    for i in range(len(data.index)):
        row = data.iloc[i]
        cn =  clusters[i]
        #logger.debug('For %s row %d (%s, %s) the cluster number is %d',
        #             label, i, row[0], row[1], cn)
        centroids[cn].append((row[0], row[1]))
    #logger.debug('Cluster sorted rows for %s:\n%s', label, cm)
    for i in range(len(centroids)):
        centroids[i] = (min([x[0] for x in centroids[i]]),
                        max([x[1] for x in centroids[i]]))
    return centroids

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
        elif ftype == 4: # Exonerate ryo output
            int_fields = [x for x in range(3,13)]
            float_fields = [2]
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
    elif ftype == 4:
        if fields[10] < 3 * fields[11]:
            # subject (query) longer than query (target)
            logger.debug('HSP subject longer than query: %s: %d < 3 * %d', fields[1], fields[10], fields[11])
            fit = False
        if fields[2] < args.pident:
            logger.debug('HSP has insufficient percent identity: %.3f < %.3f', fields[2], args.pident)
            fit = False
        # Cannot perform this check
        #if float(fields[]) / fields[11] < args.bitscore_ratio:
        #    # insufficient bitscore value
        #    logger.debug('HSP has insufficient bitscore ratio: %.3f / %.3f < %.3f', float(fields[]), fields[11], args.bitscore_ratio)
        #    fit = False
        # Not sure of values for this, yet...
        #if float(fields[12]) / float(_slen) < args.score_ratio:
        #    # insufficient score value
        #    logger.debug('HSP has insufficient score ratio: %.3f / %.3f < %.3f', float(fields[14]), _slen, args.score_ratio)
        #    fit = False
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
    elif ftype == 4:  # Exonerate ryo output
        return float(qlen(hsp, ftype)) / float(hsp[10])
    return None

def qlen(hsp, ftype):
    if ftype == 0:  # Blast+ Custom
        return abs(hsp[11] - hsp[10]) + 1
    elif ftype == 1:  # LAST TAB
        return hsp[8]
    elif ftype in (2, 3):  # LAST BlastTab and BlastTab+
        return abs(hsp[7] - hsp[6]) + 1
    elif ftype == 4:  # Exonerate ryo output
        return hsp[6]
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
    elif ftype == 4:  # Exonerate ryo output
        return float(slen(hsp, ftype)) / float(hsp[11])
    return None

def slen(hsp, ftype):
    if ftype == 0:  # BLAST+ Custom
        return abs(hsp[13] - hsp[12]) + 1
    elif ftype == 1:  # LAST TAB
        return hsp[3]
    elif ftype in (2, 3):  # LAST BlastTab and BlastTab+
        return abs(hsp[9] - hsp[8]) + 1
    elif ftype == 4:  # Exonerate ryo output
        return hsp[9]
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
    ftype_list = ['custom', 'lasttab', 'blasttab', 'blasttab+', 'exonerate']
    cluster_algorithms = ['kmeans', 'ckmeans']
    parser = argparse.ArgumentParser(description='Test Clustering')
    parser.add_argument('--num-threads', type=int, default=NUM_THREADS,
                        help='Number of CPU threads to use')
    clustering_group = parser.add_argument_group(title='Clustering options',
                                                 description='Options for clustering/partitioning the data')
    clustering_group.add_argument('--algorithm', choices=cluster_algorithms, default='kmeans',
                                  help='Clustering/partitioning algorithm to use')
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


    data_group = parser.add_argument_group(title='Data Files',
                                           description='Required data files')
    data_group.add_argument('--gff', required=True,
                             help='Input GFF file')
    data_group.add_argument('--genome', required=True,
                             help='Input genome FASTA file')
    data_group.add_argument('--refdb', default=REFDB_FILE,
                             help='Reference protein FASTA DB file for BLAST')
    data_group.add_argument('--transcriptome', required=True,
                             help='Input transcriptome FASTA file')


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
    if args.algorithm == 'ckmeans' and ftype in [1, 2]:
        logger.error('Algorithm ckmeans not supported for type %s', args.ftype)
        return 1
    t_start = time.time()
    logger.info('Begin')
    ##
    # Ingest input files
    coordinate_map = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.num_threads) as executor:
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

    ##
    # perform k-means clustering to find optimized centroids
    logger.info('Processing %d extracted results', len(coordinate_map))
    cntr = 0
    cluster_map = {}
    cluster_func = eval('calculate_clusters_%s' % args.algorithm)
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.num_threads) as executor:
        future_map = {executor.submit(cluster_func, label, ftype, coordinates): label for label, coordinates in coordinate_map.items()}
        for future in concurrent.futures.as_completed(future_map):
            cntr += 1
            logger.info('Percent completion: %.3f', float(cntr) / float(len(coordinate_map)) * 100.0)
            label = future_map[future]
            try:
                data = future.result()
            except Exception as exc:
                logger.error('%r generated on exception: %s', label, exc)
            else:
                if data:
                    logger.info('%r return %d clusters', label, len(data))
                    cluster_map[label] = data

    t_end = time.time()
    logger.info('Complete. Elapsed %.3f seconds.', t_end - t_start)
    #logger.info('Final chimeric map:\n%s', pprint.pformat(merge_map))
    logger.info('Final chimeric map:\n%s', pprint.pformat(cluster_map))
    '''
    dbfile = os.path.splitext(os.path.basename(args.gff))[0] + '.db'
    dbfile = os.path.join(args.outdir, dbfile)
    if not os.path.exists(dbfile):
        logger.info('Initialize the database...')
        db = gffutils.create_db(args.gff, dbfile, keep_order=False, merge_strategy='create_unique', id_spec=['ID'], sort_attribute_values=False)
        db.execute('CREATE INDEX coordinate_index ON features(seqid, start, end)')
    else:
        logger.info('Loading genome features from %s', dbfile)
        db = gffutils.FeatureDB(dbfile)
    '''

    return 0

if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        logger.info('Received keyboard interrupt. Aborting.')
        sys.exit(1)
