'''
Given input GFF, genome FASTA, and transcriptome FASTA files
* use blastx of each gene feature in the transcriptome against UniProt/SwissProt
  DB to attempt to detect chimeric genes
* use the GFF file and genome FASTA to generate updated GFF, genome FASTA, and
  transcriptome FASTA files which split the chimeric genes into separate
  features

TBD: Tell us how!!
'''

import os
import sys
import logging
import logging.config
import argparse
import json
import glob
from pprint import pformat

import gffutils
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastxCommandline


#UNIPROT_DBFILE = '/data/gpfs/home/dcurdie/scratch/UniProtKB/uniprot_sprot.20190627-151946.fasta'
UNIPROT_DBFILE = '/data/gpfs/home/dcurdie/scratch/VectorBase/mosquito.fa'
IDENTITY_THRESHOLD = 0.40
COVERAGE_THRESHOLD = 0.60

logConfFile = 'logging.json'
with open(logConfFile, 'r') as fd:
    logConfContent = fd.read()
logConfData = json.loads(logConfContent)
logging.config.dictConfig(logConfData)
logger = logging.getLogger('SplitChimeras')


def adjust_coords(coords, mrna, gene, exons):
    pos_ptr = mrna.start
    start_pos = stop_pos = None
    feature_len = coords[1] - coords[0]
    offset = 0
    exons = sorted(exons)
    logger.debug('Proceeding through %s mRNA starting from %d', mrna.id, pos_ptr)
    logger.debug('Exons: %s', exons)
    while pos_ptr <= mrna.stop:
        for exon in exons:
            if pos_ptr >= exon[0] and pos_ptr <= exon[1]:
                #logger.debug('Within exon %s at position %d', exon, pos_ptr)
                if start_pos:
                    feature_len -= 1
                if offset == coords[0]:
                    logger.debug('Found start position at %d', pos_ptr)
                    start_pos = pos_ptr
                    offset += 1
                    break
                if not feature_len:
                    logger.debug('Found stop position at %d', pos_ptr)
                    stop_pos = pos_ptr
                    offset += 1
                    break
                offset += 1
            #logger.debug('Not within exon %s at position %d', exon, pos_ptr)
        if start_pos and stop_pos:
            break
        pos_ptr += 1
    return start_pos, stop_pos

def blast_sequence_files(args, suffix='fa'):
    outdir = os.path.join(args.outdir, 'seqs')
    logger.info('Blast files in %s against UniProt/SwissProt DB', outdir)
    inputFastaFiles = glob.glob(os.path.join(outdir, '*.%s' % suffix))
    logger.debug('Processing %d sequence files', len(inputFastaFiles))
    blastfmt = '"6 qseqid qlen sseqid slen qframe pident nident length mismatch gapopen qstart qend sstart send evalue bitscore"'
    blastmap = {}
    for inputFastaFile in inputFastaFiles[0:10]:
        label = os.path.splitext(os.path.basename(inputFastaFile))[0]
        outfname = os.path.join(os.path.dirname(inputFastaFile), '%s_blastx.tbl' % label)
        blastCmd = NcbiblastxCommandline(cmd='blastx', query=inputFastaFile,
                                         db=UNIPROT_DBFILE, task='blastx',
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
            # filtering criteria: percent identity, coverage
            if float(fields[5]) / 100.0 < IDENTITY_THRESHOLD:
                logger.debug('skip record with low identity: %s < %s',
                             float(fields[5]) / 100.0, IDENTITY_THRESHOLD)
                continue
            cov = float(fields[7]) / float(fields[1])
            if cov < COVERAGE_THRESHOLD:
                logger.debug('skip record with low coverage: %s < %s',
                             cov, COVERAGE_THRESHOLD)
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
    for f in blastmap.keys():
        for k, v in blastmap[f].items():
            if len(v) < 2:
                # not likely chimeric
                blastmap[f].pop(k)
                continue
            logger.debug('Likely chimeric feature %s: key:%s, hits:%d', f, k.replace(f, ''), len(v))
            chimeric_features.append(f)
    chimeric_features = list(set(chimeric_features))
    return chimeric_features

def map_feature_coords(features, blastmap, db):
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
            start, stop = adjust_coords((start, stop), mrna, gene, exons)
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
    input_group.add_argument('--transcriptome', required=True,
                             help='Input transcriptome FASTA file')
    parser.add_argument('--outdir', default=os.getcwd(),
                        help='Output directory')
    args, extra = parser.parse_known_args()
    if args.debug:
        logger.setLevel(logging.DEBUG)
    elif args.verbose:
        logger.setLevel(logging.INFO)
    logger.info('Starting')
    logger.debug('Args=%s extra=%s', args, extra)
    dbfile = os.path.splitext(args.gff)[0] + '.db'
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
    logger.debug('Return BLAST map:\n%s', pformat(blastmap))
    chimeric_features = list_chimeric_features(blastmap)
    logger.info('Got %d likely chimeric features', len(chimeric_features))
    logger.debug('Likely chimeric features: %s', chimeric_features)
    logger.info('Mapping likely chimeric features to genome')
    coords = map_feature_coords(chimeric_features, blastmap, db)

    logger.info('Loading genome FASTA data')
    genomeSeqs = SeqIO.index(args.genome, 'fasta')
    logger.info('Write out likely chimeric sequences')
    for l, m in coords.items():
        chromRec = genomeSeqs[m['chr']]
        idx = 0
        for coord in m['coords']:
            idx += 1
            if coord[0] < coord[1]:
                subseq = chromRec[coord[0]:coord[1]]
            else:
                subseq = chromRec[coord[1]:coord[0]]
            subseq.id = l + '.%d' % idx
            logger.debug('%s subsequence [%s]', l, coord)
            outf = os.path.join(args.outdir, '%s.fasta' % subseq.id)
            with open(outf, 'w') as fd:
                fd.write(subseq.format('fasta'))
            logger.info('Wrote subsequence %s to %s', subseq.id, outf)
    logger.info('Complete')
    return 0

if __name__ == '__main__':
    sys.exit(main())
