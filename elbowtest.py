#!/usr/bin/env python

import os
import sys
import logging
import math
from sklearn.cluster import KMeans
from scipy.spatial.distance import cdist
import numpy as np

logger = logging.getLogger('ElbowTest')
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s %(levelname)s: %(funcName)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)

loci = [["AT5G17630.1" "AT5G17630.1" "AT5G46110.4" "AT5G46110.4" "AT5G54800.1" "AT5G54800.1" "AT1G61800.1" "AT1G61800.1" "AT3G01550.1" "AT3G01550.1" "AT5G33320.1" "AT5G33320.1" "AT4G03950.1" "AT4G03950.1"]]
bitscores = [[337],[322],[256],[242],[787],[523],[555],[534],[218],[206],[221],[216],[291],[278]]


def main():
    min_k = 0
    max_k = len(bitscores)
    narry = np.array(bitscores)
        
    distance_variances = []
    for i in range(max_k):
        k = i + 1
        logger.info('Number of clusters: %d', k)
        kmeans = KMeans(n_clusters=k).fit(narry)
        kmeans.fit(narry)
        logger.info('Centroids: %s', kmeans.cluster_centers_.tolist())
        logger.info('Inertia: %s', kmeans.inertia_)
        variance = kmeans.inertia_ / k
        stdev = math.sqrt(variance)
        logger.info('Variance: %.1f, Stdev: %.1f', variance, stdev)
        distance_variances.append(kmeans.inertia_)
        #logger.info('Shape: %s', narry.shape)
        #cd = cdist(narry, kmeans.cluster_centers_, 'euclidean')
        #logger.info('CDist: %s', cd)
        #logger.info('Min: %s', np.min(cd, axis=1))
        #logger.info('Distortion: %s', sum(np.min(cd, axis=1) / narry.shape[0]))
        logger.info('%s', '='*80)

    firsts = []
    for i in range(len(distance_variances[1:])):
        change = (distance_variances[i+1] - distance_variances[i])
        firsts.append(abs(change))
        chgpct = change * 100.0 / distance_variances[i]
        logger.info('Clusters:%d delta:%.1f pct:%.1f', i+2, change, chgpct)

    logger.info('%s', '='*80)
    seconds = []
    for i in range(len(firsts[1:])):
        change = firsts[i+1] - firsts[i]
        seconds.append(abs(change))
        chgpct = change * 100.0 / firsts[i]
        logger.info('Clusters:%d delta:%.1f pct:%.1f', i+3, change, chgpct)

    return 0

if __name__ == '__main__':
    sys.exit(main())
