'''
TRIm, MAp, and MErge sequence reads
Supported trimming tools:
- trimgalore
- trimmomatic
- cutadapt

Mapping is via BWA. Supported algorhythms:
- mem (Maximal Expected Matches), for reads >= 100bp
- backtrack/short, for reads < 100bp

Merging is through BWA. Final output is SAM file.

Input is expect to be in '_R1' and '_R2' named paired-end FASTQ sequence files in the specified top-level directory.

Each identified pair will be processed from the specified input directory. Reference sequence file is also expected to
be within the specified input directory.

Merged output files will be written to the specified output directory.
'''

import os
import sys
import glob
import logging
import argparse
import subprocess

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s %(levelname)s: %(filename)s[%(lineno)d](%(funcName)s) - %(message)s',
                              datefmt='%F %T')
handler.setFormatter(formatter)
logger.addHandler(handler)


class Trimmer():
    def __init__(self, filepairs, adapters=None):
        for f in list(sum(filepairs, ())):
            if not os.path.exists(f):
                raise IOError('File %s does not exist' % r1)
        self.filepairs = filepairs
        if adapters:
            if not os.path.exists(adapters):
                raise IOError('File %s does not exist' % adapters)
        self.adapters = adapters
        self.command = None
        self.options = []
        self.sbatch_opts = []

    def trim(self, options=[]):
        logger.debug('Trimming files with options: %s', options)
        if not options:
            options = self.options
        for r1, r2 in self.filepairs:
            command_options = [self.command]
            command_options.extend(options)
            command_options.extend([r1, r2])
            wrapcmd = ' '.join(command_options)
            sbatch = Sbatch(self.sbatch_opts, wrapcmd)
            result = sbatch.run()
        return None
                                    

class TrimGalore(Trimmer):
    def __init__(self, argc, **kwargs):
        ''' '''
        super(TrimGalore, self).__init__(**kwargs)
        if command in kwargs:
            if os.path.exists(command):
                self.command = command
            else:
                raise IOError('Command %s not found' % command)
        else:
            self.command = 'trim_galore'
        self.options = ['--paired', '--retain_unpaired', '--cores', '4', '--max_n', '40', '--gzip']
        self.sbatch_opts = ['-N', '1', '-c', '16', '--mem=64g']


class Trimomatic(Trimmer):
    def __init__(self, args, **kwargs):
        ''' '''

    def trim(self, opts=[]):
        ''' '''
        return None


class Cutadapt(Trimmer):
    def __init__(self, args, **kwargs):
        ''' '''

    def trim(self, opts=[]):
        ''' '''
        return None


class Aligner():
    def __init__(self, inseqs, refseq, method, options=[]):
        flist = [refseq]
        flist.extend(inseqs)
        for f in flist:
            if not os.path.exists(f):
                raise IOError('File %s does not exist' % f)
        self.inseqs = inseqs
        self.refseq = refseq
        if method not in ['aln', 'mem']:
            raise OSError('Method %s not allowed for alignment' % method)
        self.method = method
        if not options:
            if method == 'aln':
                options = ['-q', 15, '-t', NTHREADS, '-f', outfile]
            elif method == 'mem':
                options = ['@RG\tID:{input}\tPL:ILLUMINA\tLB:{input}\tSM:{input}',
                           '-t', NTHREADS, '-o', outfile]
        self.options = options

    def align(self, options=[]):
        if not options:
            options = self.options
        command_options = [self.command, self.method]
        command_options.extend(options)
        command_options.append(self.refseq)
        command_options.extend(self.inseqs)
        wrapcmd = ' '.join(command_options)
        sbatch_opts = ['-N', 1, '-c', 16, '--mem=64g']
        batch = Sbatch(sbatch_opts, wrapcmd)
        result = batch.run()
        return result


class Merger():
    def __init__(self):
        ''' '''

    def merge(self, opts={}):
        ''' '''
        return None


class JavaJar():
    def __init__(self, jarfile, javaopts=[], jaropts=[]):
        if not os.path.exists(jarfile):
            raise IOError('JAR file %s does not exist' % jarfile)
        self.jarfile = jarfile
        self.javaopts = javaopts
        self.jaropts = jaropts

    def execute(self, javaopts=[], jaropts=[]):
        if not javaopts:
            javaopts = self.javaopts
        if not jaropts:
            jaropts = self.jaropts


class Sbatch():
    def __init__(self, options=[], wrapcmd=None):
        self.command = 'sbatch'
        if not options:
            options = ['-N', 1, '-c', 16, '--mem=64g']
        self.options = options
        self.wrapcmd = wrapcmd

    def run(options=[], wrapcmd=None):
        command_options = [self.command]
        if not options:
            options = self.options
        command_options.extend(options)
        if not wrapcmd:
            wrapcmd = self.wrapcmd
        command_options.extend(['--wrap="%s"' % wrapcmd])

        logger.info('Launching sbatch with: %s', command_options)
        return True
        try:
            result = subprocess.run(command_options,
                                    stderr=subprocess.PIPE,
                                    stdout=subprocess.PIPE,
                                    check=True)
        except subprocess.CalledProcessError as err:
            logger.error('An error occurred queuing sbatch: %s', err)
            return None
        return result


class Samtools():
    def __init__(self):
        ''' '''

    def flagstat(self, opts=[]):
        ''' '''
        return None

    def index(self, opts=[]):
        ''' '''
        return None

    def sort(self, opts=[]):
        ''' '''
        return None

    def view(self, opts=[]):
        ''' '''
        return None

    def run(self, command, opts=[]):
        opts.insert(0, command)

        return None

def findFiles(dname):
    files = glob.glob(os.path.join(dname, '*_R*.fastq'))
    pairdata = {}
    for f in files:
        base, _ = os.path.basename(f).split('_R')
        if pairdata.get(base):
            pairdata[base]['files'].append(f)
            pairdata[base]['files'].sort()
        else:
            pairdata.update({base: {'files': [f]}})
    return pairdata


def main():
    logger.info('Starting')
    trimopts = {
        'trimgalore': TrimGalore,
        'trimmomatic': Trimomatic,
        'cutadapt': Cutadapt
    }
    parser = argparse.ArgumentParser()
    parser.add_argument('--topdir', type=str, required=True,
                        help='Top-level directory containing paired-end sequence data input files')
    parser.add_argument('--refseq', type=str, required=True,
                        help='Reference sequence data file')
    parser.add_argument('--outdir', type=str, required=True,
                        help='Directory to write output files into')
    parser.add_argument('--trimmer', type=str, choices=trimopts.keys(), default='trimmomatic',
                        help='Tool to use for sequence trimming')

    logger.info('Proccessing arguments')
    args = parser.parse_args()

    logger.info('Verifying arguments')
    if not os.path.exists(args.topdir):
        logger.error('Top-level directory %s does not exist', args.topdir)
        return 1
    if not os.path.exists(args.refseq):
        logger.error('Reference sequence file %s does not exist', args.refseq)
        return 1
    if os.path.exists(args.outdir):
        logger.warning('Output directory %s exists, files may get overwritten', args.outdir)
    else:
        os.mkdir(args.outdir)

    logger.info('Finding input data')
    pairdata = findFiles(args.topdir)
    if not pairdata:
        logger.error('No paired-end sequence data files found in %s', args.topdir)
        return 1
    for s in pairdata:
        if not pairdata[s].get('files'):
            pairdata.pop(s)
    if not pairdata:
        logger.error('No paired-end sequence data files found in %s', args.topdir)
        return 1
    logger.info('Found %d sets of files', len(pairdata))

    trimmer = trimopts[args.trimmer]([pairdata[x].get('files') for x in pairdata])

    trimmer.trim()

    logger.info('Complete')
    return 0


if __name__ == "__main__":
    sys.exit(main())
