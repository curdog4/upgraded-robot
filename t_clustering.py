'''Test use of kmeans clustering to determine chimeras in blastx table results
0-qseqid 1-qlen 2-sseqid 3-slen 4-qframe 5-pident 6-nident 7-length
8-mismatch 9-gapopen 10-qstart 11-qend 12-sstart 13-send 14-evalue
15-bitscore

Logic borrowed from:
https://blog.cambridgespark.com/how-to-determine-the-optimal-number-of-clusters-for-k-means-clustering-14f27070048f
https://jtemporal.com/kmeans-and-elbow-method/

'''

import concurrent.futures
import os
import sys
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

DATA_DIR = 'seqs'
DATA_FILES = glob.glob(os.path.join(DATA_DIR, '*_blastx.tbl'))
LENGTH_THRESHOLD = 100
QCOV_THRESHOLD = 0.30
SCOV_THRESHOLD = 0.70
PIDENT_THRESHOLD = 20.0

SUBLIST = ['seqs/Mecry_04G102350.1_blastx.tbl', 'seqs/Mecry_05G142250.1_blastx.tbl', 'seqs/Mecry_09G233490.1_blastx.tbl',
           'seqs/Mecry_01G024890.1_blastx.tbl', 'seqs/Mecry_04G120760.1_blastx.tbl', 'seqs/Mecry_09G230850.1_blastx.tbl',
           'seqs/Mecry_04G119610.2_blastx.tbl', 'seqs/Mecry_07G182900.1_blastx.tbl', 'seqs/Mecry_05G144460.1_blastx.tbl',
           'seqs/Mecry_01G001370.1_blastx.tbl', 'seqs/Mecry_09G236430.1_blastx.tbl', 'seqs/Mecry_04G119610.1_blastx.tbl',
           'seqs/Mecry_08G212340.1_blastx.tbl', 'seqs/Mecry_02G064460.1_blastx.tbl', 'seqs/Mecry_05G147260.1_blastx.tbl',
           'seqs/Mecry_09G247560.1_blastx.tbl', 'seqs/Mecry_04G120760.2_blastx.tbl', 'seqs/Mecry_04G119610.3_blastx.tbl',
           'seqs/Mecry_05G127730.1_blastx.tbl']

def calculate_data(data):
    wcss = []
    for n in range(1, min(len(data), 15)):
        kmeans = KMeans(n_clusters=n)
        kmeans.fit(X=data)
        wcss.append(kmeans.inertia_)
    return wcss

def get_blastx_table_data(fname):
    with open(fname, 'r') as fdesc:
        fdata = fdesc.read()
    coordinates = []
    for line in fdata.splitlines():
        fields = line.split('\t')
        for idx in [1, 3, 7, 10, 11, 12, 13]:
            fields[idx] = int(fields[idx])
        for idx in [5, 15]:
            fields[idx] = float(fields[idx])
        '''Checks for 'fitness' of the result
        '''
        if fields[1] < 3 * fields[3] :
            # subject longer than query
            continue
        if qlen(fields) < LENGTH_THRESHOLD:
            # insufficient alignment length
            continue
        if fields[5] < PIDENT_THRESHOLD:
            # insufficient identity
            continue
        if qcov(fields) < QCOV_THRESHOLD:
            # insufficient query coverage
            continue
        if scov(fields) < SCOV_THRESHOLD:
            # insufficient subject coverage
            continue
        if float(fields[3]) / 2 < fields[15]:
            # insufficient bitscore value
            continue

        coordinates.append(fields)
    return coordinates

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

def qcov(hsp):
    return float(qlen(hsp)) / float(hsp[1])

def qlen(hsp):
    return abs(hsp[11] - hsp[10]) + 1

def scov(hsp):
    return  float(slen(hsp)) / float(hsp[3]) 

def slen(hsp):
    return abs(hsp[13] - hsp[12]) + 1

def main():
    start = time.time()
    logger.info('Begin')
    coordinate_map = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=12) as executor:
        #future_map = {executor.submit(get_blastx_table_data, fname): fname for fname in ['seqs/Mecry_07G198110.1_blastx.tbl']}
        future_map = {executor.submit(get_blastx_table_data, fname): fname for fname in SUBLIST}
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

    #headers = 'qseqid,qlen,sseqid,slen,qframe,pident,nident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore'
    headers = 'start\tend'
    for label, coordinates in coordinate_map.items():
        if len(coordinates) < 3:
            logger.error('Insufficient records to continue for %s', label)
            continue

        sio = io.StringIO()
        sio.write(headers + '\n')
        for coords in coordinates:
            sio.write('%s\t%s\n' % (coords[10], coords[11]))
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
        cm = []
        for i in range(n):
            cm.append([])
        for i in range(len(data.index)):
            row = data.iloc[i]
            cn =  clusters[i]
            #logger.debug('For row %d (%s, %s) the cluster number is %d',
            #             i, row[0], row[1], cn)
            cm[cn].append((row[0], row[1]))
        #logger.debug('Cluster sorted rows:\n%s', cm)
        for cntr in cm:
            start = min([x[0] for x in cntr])
            end = max([x[1] for x in cntr])
            logger.info('Possible chimeric gene in region from %s to %s',
                        start, end)

    end = time.time()
    logger.info('Complete. Elapsed %.3f seconds.', end - start)
    return 0

if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        logger.info('Received keyboard interrupt. Aborting.')
        sys.exit(1)
