'''
Given input GFF, genome FASTA, and transcriptome FASTA files
* use blastx of each gene feature in the transcriptome against UniProt/SwissProt
  DB to attempt to detect chimeric genes
* use the GFF file and genome FASTA to generate updated GFF, genome FASTA, and
  transcriptome FASTA files which split the chimeric genes into separate
  features

Substantially influenced by / borrow from:
* Defusion: https://github.com/wjidea/defusion
* Optimize Assembler: https://bitbucket.org/yangya/optimize_assembler/
'''

import os
import sys
import logging
import logging.config
import argparse
import json
import glob
import math

'''
Top-level third-party / external dependencies: gffutils, biopython, python2.7
Package dependency tree:
- gffutils=0.9=py_1 (channel:bioconda)
  - argcomplete=1.9.5=py27_0
    - python >=2.7,<2.8.0a0
  - argh=0.26.1=py27_0
    - python 2.7*
  - pyfaidx=0.5.5.2=py_1 (channel:bioconda)
    - python
    - six=1.12.0=py27_0
      - python >=2.7,<2.8.0a0
  - python
  - simplejson=3.8.1=py27_0 (channel:bioconda)
    - python 2.7*
  - six=1.12.0=py27_0
    - python >=2.7,<2.8.0a0
- biopython=1.70=np112py27_1 (Bio) (channel:bioconda)
  - mmtf-python=1.0.2=py27_0 (channel:bioconda)
    - msgpack-python=0.6.1=py27hfd86e86_1
      - libgcc-ng=8.2.0=hdf63c60_1g
        - _libgcc_mutex * main
      - libstdcxx-ng=8.2.0=hdf63c60_1
      - python >=2.7,<2.8.0a0
    - python 2.7*
  - numpy=1.12.1=py27h9378851_1
    - libgcc-ng=8.2.0=hdf63c60_1
      - _libgcc_mutex * main
    - libgfortran-ng=7.3.0=hdf63c60_0
    - python >=2.7,<2.8.0a0
    - mkl=2018.0.3=1
      - intel-openmp=2019.1=144
    - blas=1.0=mkl
  - python=2.7.15=h9bab390_6
    - libffi >=3.2.1,<4.0a0
    - libgcc-ng=8.2.0=hdf63c60_1
      - _libgcc_mutex * main
    - libstdcxx-ng=8.2.0=hdf63c60_1
    - ncurses=6.1=he6710b0_1
      - libgcc-ng=8.2.0=hdf63c60_1
        - _libgcc_mutex * main
      - libstdcxx-ng=8.2.0=hdf63c60_1
    - openssl=1.1.1c=h7b6447c_1
      - ca-certificates=2019.5.15=0
      - libgcc-ng=8.2.0=hdf63c60_1
        - _libgcc_mutex * main
    - readline=7.0=h7b6447c_5
      - libgcc-ng=8.2.0=hdf63c60_1
        - _libgcc_mutex * main
      - ncurses=6.1=he6710b0_1
        - libgcc-ng=8.2.0=hdf63c60_1
          - _libgcc_mutex * main
        - libstdcxx-ng=8.2.0=hdf63c60_1
    - sqlite=3.26.0=h7b6447c_0
      - libedit=3.1.20181209=hc058e9b_0
        - libgcc-ng=8.2.0=hdf63c60_1
          - _libgcc_mutex * main
        - ncurses=6.1=he6710b0_1
          - libgcc-ng=8.2.0=hdf63c60_1
            - _libgcc_mutex * main
          - libstdcxx-ng=8.2.0=hdf63c60_1
      - libgcc-ng=8.2.0=hdf63c60_1
        - _libgcc_mutex * main
    - tk=8.6.8=hbc83047_0
      - libgcc-ng=8.2.0=hdf63c60_1
        - _libgcc_mutex * main
      - zlib=1.2.11=h7b6447c_3
        - libgcc-ng=8.2.0=hdf63c60_1
          - _libgcc_mutex * main
    - zlib=1.2.11=h7b6447c_3
      - libgcc-ng=8.2.0=hdf63c60_1
        - _libgcc_mutex * main
    - pip=19.0.3=py27_0
      - python >=2.7,<2.8.0a0
      - setuptools=40.8.0=py27_0
        - certifi=2019.6.16=py27_0
          - python >=2.7,<2.8.0a0
        - python >=2.7,<2.8.0a0
      - wheel=0.33.1=py27_0
        - python >=2.7,<2.8.0a0
        - setuptools=40.8.0=py27_0
          - certifi=2019.6.16=py27_0
            - python >=2.7,<2.8.0a0
          - python >=2.7,<2.8.0a0
  - reportlab=3.5.19=py27he686d34_0
    - freetype=2.9.1=h8a8886c_1
      - libgcc-ng=8.2.0=hdf63c60_1
        - _libgcc_mutex * main
      - libpng=1.6.37=hbc83047_0
        - libgcc-ng=8.2.0=hdf63c60_1
          - _libgcc_mutex * main
        - zlib=1.2.11=h7b6447c_3
          - libgcc-ng=8.2.0=hdf63c60_1
            - _libgcc_mutex * main
      - zlib=1.2.11=h7b6447c_3
        - libgcc-ng=8.2.0=hdf63c60_1
          - _libgcc_mutex * main
    - libgcc-ng=8.2.0=hdf63c60_1
      - _libgcc_mutex * main
    - pillow=6.0.0=py27h34e0f95_0
      - freetype=2.9.1=h8a8886c_1
        - libgcc-ng=8.2.0=hdf63c60_1
          - _libgcc_mutex * main
        - libpng=1.6.37=hbc83047_0
          - libgcc-ng=8.2.0=hdf63c60_1
            - _libgcc_mutex * main
          - zlib=1.2.11=h7b6447c_3
            - libgcc-ng=8.2.0=hdf63c60_1
              - _libgcc_mutex * main
      - jpeg=9b=h024ee3a_2
        - libgcc-ng=8.2.0=hdf63c60_1
          - _libgcc_mutex * main
      - libgcc-ng=8.2.0=hdf63c60_1
        - _libgcc_mutex * main
      - libtiff=4.0.10=h2733197_2
        - jpeg=9b=h024ee3a_2
          - libgcc-ng=8.2.0=hdf63c60_1
            - _libgcc_mutex * main
        - libgcc-ng=8.2.0=hdf63c60_1
          - _libgcc_mutex * main
        - libstdcxx-ng=8.2.0=hdf63c60_1
        - xz=5.2.4=h14c3975_4
          - libgcc-ng=8.2.0=hdf63c60_1
            - _libgcc_mutex * main
        - zlib=1.2.11=h7b6447c_3
          - libgcc-ng=8.2.0=hdf63c60_1
            - _libgcc_mutex * main
        - zstd=1.3.7=h0b5b093_0
          - libgcc-ng=8.2.0=hdf63c60_1
            - _libgcc_mutex * main
          - libstdcxx-ng=8.2.0=hdf63c60_1
          - xz=5.2.4=h14c3975_4
            - libgcc-ng=8.2.0=hdf63c60_1
              - _libgcc_mutex * main
          - zlib=1.2.11=h7b6447c_3
            - libgcc-ng=8.2.0=hdf63c60_1
              - _libgcc_mutex * main
      - olefile=0.46=py27_0
        - python >=2.7,<2.8.0a0
      - python >=2.7,<2.8.0a0
      - tk=8.6.8=hbc83047_0
        - libgcc-ng=8.2.0=hdf63c60_1
          - _libgcc_mutex * main
        - zlib=1.2.11=h7b6447c_3
          - libgcc-ng=8.2.0=hdf63c60_1
            - _libgcc_mutex * main
      - zlib=1.2.11=h7b6447c_3
        - libgcc-ng=8.2.0=hdf63c60_1
          - _libgcc_mutex * main
    - python >=2.7,<2.8.0a0

'''
import gffutils
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastxCommandline


COVERAGE_THRESHOLD = 0.60
IDENTITY_THRESHOLD = 0.40
LENGTH_CUTOFF = 100
PADDING_FACTOR=0.25
#UNIPROT_DBFILE = '/data/gpfs/home/dcurdie/scratch/UniProtKB/uniprot_sprot.20190627-151946.fasta'
UNIPROT_DBFILE = '/data/gpfs/home/dcurdie/scratch/VectorBase/mosquito.fa'

logConfFile = 'logging.json'
with open(logConfFile, 'r') as fd:
    logConfContent = fd.read()
logConfData = json.loads(logConfContent)
logging.config.dictConfig(logConfData)
logger = logging.getLogger('SplitChimeras')


def adjust_coords(coords, mrna, gene, exons, pad_factor):
    pos_ptr = mrna.start
    start_pos = stop_pos = None
    feature_len = coords[1] - coords[0]
    padding = int(math.ceil(pad_factor * float(feature_len)))
    logger.debug('Feature length is %d, padding is %d', feature_len, padding)
    offset = 0
    exons = sorted(exons)
    logger.debug('Proceeding through %s mRNA starting from %d', mrna.id, pos_ptr)
    logger.debug('Exons: %s', exons)
    while pos_ptr <= mrna.stop:
        #logger.debug('Position pointer: %d', pos_ptr)
        for exon in exons:
            if pos_ptr >= exon[0] and pos_ptr <= exon[1]:
                #logger.debug('Within exon %s at position %d', exon, pos_ptr)
                if start_pos:
                    feature_len -= 1
                elif offset == coords[0]:
                    start_pos = pos_ptr
                    feature_len -= 1
                offset += 1
                if not feature_len:
                    stop_pos = pos_ptr
                    break
            #logger.debug('Not within exon %s at position %d', exon, pos_ptr)
        if start_pos and stop_pos:
            start_pos -= padding
            stop_pos += padding
            break
        pos_ptr += 1
    return start_pos, stop_pos


def qcov(hsp):
    return abs(int(hsp[11]) - int(hsp[10])) + 1


def separated(hsp1, hsp2):
    ''' given two hsps, return True if
    overlap les than 20% of the shorter and overlap less than 60 bp
    '''
    length1 = qcov(hsp1)
    length2 = qcov(hsp2)
    start = min(int(hsp1[10]), int(hsp1[11]), int(hsp2[10]), int(hsp2[11]))
    end = max(int(hsp1[10]), int(hsp1[11]), int(hsp2[10]), int(hsp2[11]))
    overlap = length1 + length2 - (end - start) + 1
    # value of overlap can < 0 but only the upper limit maters
    if overlap < min(60, 0.2 * min(length1, length2)):
        return True
    return False


def blast_sequence_files(args, suffix='fa'):
    outdir = os.path.join(args.outdir, 'seqs')
    logger.info('Blast files in %s against reference protein DB', outdir)
    logger.debug('Reference protein db=%s', args.refdb)
    inputFastaFiles = glob.glob(os.path.join(outdir, '*.%s' % suffix))
    logger.debug('Processing %d sequence files', len(inputFastaFiles))
    blastfmt = '"6 qseqid qlen sseqid slen qframe pident nident length mismatch gapopen qstart qend sstart send evalue bitscore"'
    blastmap = {}
    for inputFastaFile in [x for x in inputFastaFiles]: # if x.endswith('seqs/Mecry_04G114660.1.fasta')]:
        label = os.path.splitext(os.path.basename(inputFastaFile))[0]
        outfname = os.path.join(os.path.dirname(inputFastaFile), '%s_blastx.tbl' % label)
        blastCmd = NcbiblastxCommandline(cmd='blastx', query=inputFastaFile,
                                         db=args.refdb, task='blastx',
                                         evalue='1e-30', outfmt=blastfmt,
                                         max_hsps='10', num_threads='4')  # maybe also/instead max_target_seqs='100'

        logger.debug('BLASTing %s', os.path.basename(inputFastaFile))
        blastResults = blastCmd()
        if len(blastResults[0]) == 0:
            logger.warn('Zero-length results BLASTing %s', inputFastaFile)
            continue
        if len(blastResults[1]) > 0:
            logger.warn('Error output from BLASTing %s:\n%s', inputFastaFile,
                        blastResults[1])
        logger.debug('Writing blast output of %s to %s', inputFastaFile, outfname)
        with open(outfname, 'w') as fd:
            fd.write(blastResults[0])
        records = blastResults[0].splitlines()
        '''
        if not os.path.exists(outfname):
            continue
        with open(outfname, 'r') as fd:
            records = fd.readlines()
        '''
        logger.info('Got %d hits for sequence %s', len(records), label)
        for record in records:
            fields = record.split('\t')
            featurename = fields[0]
            subject_id = fields[2]
            matchkey = '%s\0%s' % (featurename, subject_id)
            ##
            # filtering criteria: query len, percent identity, coverage, bitscore
            if qcov(fields) < LENGTH_CUTOFF:
                logger.debug('skip record for %s with insufficent query length: %s < %s',
                             fields[2], qcov(fields), LENGTH_CUTOFF)
            if float(fields[5]) / 100.0 < args.identity:
                logger.debug('skip record for %s with low identity: %s < %s',
                             fields[2], float(fields[5]) / 100.0, args.identity)
                continue
            cov = float(fields[7]) / float(fields[3])
            if cov < args.coverage:
                logger.debug('skip record for %s with low coverage: %s < %s',
                             fields[2], cov, args.coverage)
                continue
            if float(fields[3]) / 2.0 < float(fields[15]):
                logger.debug('skip record for %s with low bitscore value: %s / 2 < %s',
                             fields[2], fields[3], fields[15])
            blastmap.setdefault(featurename, {})
            blastmap[featurename].setdefault(matchkey, [])
            blastmap[featurename][matchkey].append(fields)
    return blastmap
        

def separate_transcripts_for_blast(args, suffix='fa'):
    outdir = os.path.join(args.outdir, 'seqs')
    logger.debug('Output individual sequence record FASTA files from %s in %s',
                 args.transcriptome, outdir)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    seqData = SeqIO.parse(args.transcriptome, 'fasta')
    for seqRecord in seqData:
        fname = '%s.%s' % (seqRecord.id, suffix)
        logger.debug('Output sequence record for %s to %s in %s', seqRecord.id,
                     fname, outdir) 
        if not os.path.exists(os.path.join(outdir, fname)):
            SeqIO.write(seqRecord, os.path.join(outdir, fname), 'fasta')

def list_chimeric_features(blastmap):
    chimeric_features = []
    for f in list(blastmap.keys())[:]:
        feature_pairs = []
        for k, v in list(blastmap[f].items())[:]:
            if len(v[:]) < 2:
                # not likely chimeric
                blastmap[f].pop(k)
                continue
            for i in range(len(v[:])):
                separate = True
                for j in range(i+1, len(v[:])):
                    separate = min(separate, separated(v[i], v[j]))
                if separate:
                    logger.debug('Likely chimeric feature %s: key:%s, hits:%d', f, k.replace(f, ''), len(v))
                    logger.debug('HSP: %s', v[i])
                    feature_pairs.append(v[i])
        for i in range(len(feature_pairs)):
            separate = True
            for j in range(i+1, len(feature_pairs)):
                separate = min(separate, separated(feature_pairs[i], feature_pairs[j]))
            if separate:
                chimeric_features.append(f)
    chimeric_features = list(set(chimeric_features))
    return chimeric_features

def map_feature_coords(features, blastmap, db, pad_factor):
    coords = {}
    for label in features:
        mrna = db[label]
        logger.debug('Feature mRNA %s: start=%s, stop=%s', label, mrna.start,
                     mrna.stop)
        feature_id = label
        logger.debug('Finding gene feature %s in the GFF DB', feature_id)
        errors = False
        while db[feature_id].featuretype != 'gene':
            logger.debug('recurse to find parent of %s', feature_id)
            if db[feature_id].attributes.get('Parent'):
                feature_id = db[feature_id].attributes.get('Parent')[0]
            else:
                logger.error('no parent for %s', feature_id)
                errors = True
                break
        if errors:
            logger.error('Problems processing feature %s', feature_id)
            continue
        gene = db[feature_id]
        logger.info('Feature gene %s, chr:%s start:%d stop:%d', feature_id, gene.chrom,
                    gene.start, gene.stop)
        exons = []
        for f in db.children(feature_id):
            if f.featuretype == 'exon':
                exons.append((f.start, f.stop))
        if len(exons) == 0:
            logger.error('No exons found for feature %s', feature_id)
            continue
        exons = sorted(exons)
        logger.debug('Exon coordinates for feature %s: %s', feature_id, exons)
        coords.setdefault(feature_id, {'chr': gene.chrom, 'coords': []})
        for k, v in blastmap[label].items():
            if len(v) >= 2:
                for hit in v:
                    start, stop = int(hit[10]), int(hit[11])
                    if start > stop:
                        start, stop = stop, start
                    if len(coords[feature_id]['coords']) == 0:
                        coords[feature_id]['coords'].append((start, stop))
                    else:
                        found = False
                        for coord in coords[feature_id]['coords']:
                            if start < coord[1] and stop > coord[0]:
                                if start > 0.8 * coord[0] and stop < 1.2 * coord[1]:
                                    found = True
                                    start, stop = min(start, coord[0]), max(stop, coord[1])
                                    coords[feature_id]['coords'][coords[feature_id]['coords'].index(coord)] = (start, stop)
                                    break
                        if not found:
                            coords[feature_id]['coords'].append((start, stop))
        for idx in range(len(coords[feature_id]['coords'])):
            start, stop = coords[feature_id]['coords'][idx]
            logger.debug('Adjust coordinates %s to genome', (start, stop))
            start, stop = adjust_coords((start, stop), mrna, gene, exons, pad_factor)
            logger.debug('Genome adjusted coordinates are %s', (start, stop))
            coords[feature_id]['coords'][idx] = (start, stop)
            
    return coords

def main():
    parser = argparse.ArgumentParser(description='Split chimeric genes')
    loglevel_group = parser.add_mutually_exclusive_group()
    loglevel_group.add_argument('--debug', action='store_true',
                                help='Produce debugging output')
    loglevel_group.add_argument('--verbose', action='store_true',
                                help='Produce more verbose output')
    input_group = parser.add_argument_group(title='Input Files',
                                            description='Required input files')
    input_group.add_argument('--gff', required=True,
                             help='Input GFF file')
    input_group.add_argument('--genome', required=True,
                             help='Input genome FASTA file')
    input_group.add_argument('--refdb', default=UNIPROT_DBFILE,
                             help='Reference protein FASTA DB file for BLAST')
    input_group.add_argument('--transcriptome', required=True,
                             help='Input transcriptome FASTA file')
    blast_group = parser.add_argument_group(title='BLAST filtering',
                                            description='Options controlling BLAST result filtering')
    blast_group.add_argument('--identity', default=IDENTITY_THRESHOLD, type=float,
                             help='Minimum identity threshold for BLAST HSP')
    blast_group.add_argument('--coverage', default=COVERAGE_THRESHOLD, type=float,
                             help='Minimum coverage threshold for BLAST HSP')
    parser.add_argument('--outdir', default=os.getcwd(),
                        help='Output directory')
    parser.add_argument('--padding', default=PADDING_FACTOR, type=float,
                        help='Fraction of the alignment length to use as padding')
    args, extra = parser.parse_known_args()
    if args.debug:
        logger.setLevel(logging.DEBUG)
    elif args.verbose:
        logger.setLevel(logging.INFO)
    logger.info('Starting')
    logger.debug('Args=%s extra=%s', args, extra)
    dbfile = os.path.splitext(os.path.basename(args.gff))[0] + '.db'
    dbfile = os.path.join(args.outdir, dbfile)
    if not os.path.exists(dbfile):
        logger.info('Initialize the database...')
        db = gffutils.create_db(args.gff, dbfile, keep_order=False, merge_strategy='create_unique', id_spec=['ID'], sort_attribute_values=False)
        db.execute('CREATE INDEX coordinate_index ON features(seqid, start, end)')
    else:
        logger.info('Loading genome features from %s', dbfile)
        db = gffutils.FeatureDB(dbfile)
    if not os.path.exists(args.genome):
        logger.error('Genome FASTA %s not found', args.genome)
        return 1
    logger.info('Splitting out individual sequence FASTA files')
    separate_transcripts_for_blast(args, suffix='fasta')
    logger.info('BLASTing each sequence file against UniProt/SwissProt DB')
    blastmap = blast_sequence_files(args, suffix='fasta')
    logger.debug('Return BLAST map:\n%s', json.dumps(blastmap))
    chimeric_features = list_chimeric_features(blastmap)
    logger.info('Got %d likely chimeric features', len(chimeric_features))
    logger.debug('Likely chimeric features: %s', chimeric_features)
    logger.info('Mapping likely chimeric features to genome')
    coords = map_feature_coords(chimeric_features, blastmap, db, args.padding)

    logger.info('Loading genome FASTA data')
    genomeSeqs = SeqIO.index(args.genome, 'fasta')
    logger.info('Write out likely chimeric sequences')
    for l, m in coords.items():
        chromRec = genomeSeqs[m['chr']]
        idx = 96
        for coord in m['coords']:
            idx += 1
            if coord[0] < coord[1]:
                subseq = chromRec[coord[0]:coord[1]]
            else:
                subseq = chromRec[coord[1]:coord[0]]
            subseq.id = l + '.%s' % chr(idx)
            logger.debug('%s subsequence [%s]', l, coord)
            outf = os.path.join(args.outdir, '%s.fasta' % subseq.id)
            with open(outf, 'w') as fd:
                fd.write(subseq.format('fasta'))
            logger.info('Wrote subsequence %s to %s', subseq.id, outf)
    logger.info('Complete')
    return 0

if __name__ == '__main__':
    sys.exit(main())
