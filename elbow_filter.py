#!/usr/bin/env python

import os
import sys
import logging
import argparse
import numpy as np
import pandas
from sklearn.cluster import KMeans
#from scipy.spacial.distance import cdist

logger = logging.getLogger('ElbowFilter')
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s %(levelname)s: %(funcName)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)

def main():
    parser = argparse.ArgumentParser(description='KMeans Elbow filtering of blast results')
    names = ['query_id', 'subject_id', 'pct_identity', 'alignment_length',
             'mismatches', 'gap_opens', 'query_start', 'query_end',
             'subject_start', 'subject_end', 'evalue', 'bitscore']
    parser.add_argument('--clusters', type=int, default=3,
                        help='Number of clusters to use in KMeans')
    parser.add_argument('--column', type=str, default='bitscore', choices=names,
                        help='Column to calculate KMeans from')
    args = parser.parse_args()
    indata = pandas.read_csv(sys.stdin, sep='\t', header=None, names=names)
    loci = indata[['subject_id']]
    scores = indata[[args.column]]

    kmeans = KMeans(n_clusters=args.clusters).fit(scores)
    logger.info('Centers: %s', kmeans.cluster_centers_.tolist())
    logger.info('Labels: %s', kmeans.labels_.tolist())
    for i, l in enumerate(kmeans.labels_):
        if l != 0:
            # only interested in first cluster
            continue
        #logger.info('Loci: %s, Score: %s', loci.values[i].item(), scores.values[i].item())
        outline = '{0:s}\t{1:s}\t{2:.1f}\t{3:d}\t{4:d}\t{5:d}\t{6:d}\t{7:d}\t{8:d}\t{9:d}\t{10:.0e}\t{11:g}\n'.format(*indata.values[i].tolist())
        sys.stdout.write(outline)


if __name__ == '__main__':
    sys.exit(main())
