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

logger = logging.getLogger('trimame')
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s %(levelname)s: %(filename)s[%(lineno)d](%funcName)s) - %(message)s',
                              datefmt='%F %T')
handler.setFormatter(formatter)
logger.addHandler(handler)


class Trimmer():
    def __init__(self, r1, r2, adapters=None):
        for f in (r1, r2):
            if not os.path.exists(f):
                raise IOError('File %s does not exist' % f)
        self.r1 = r1
        self.r2 = r2
        if adapters:
            if not os.path.exists(adapters):
                raise IOError('File %s does not exist' % adapters)
        self.adapters = adapters
        self.command = None
        self.options = []

    def trim(self, options=[]):
        if not options:
            options = self.options
        command_options = [self.command]
        command_options.extend(options)
        try:
            result = subprocess.run(command_options,
                                    stderr=subprocess.PIPE,
                                    stdout=subprocess.PIPE,
                                    check=True)
        except subprocess.CalledProcessError as err:
            logger.error('An error occurred running %s: %s', self.command, err)
            return None
                                    


class TrimGalore(Trimmer):
    def __init__(self, argc, **kwargs):


    def trim(self, opts=[]):


class Trimomatic(Trimmer):
    def __init__(self, args, **kwargs):


    def trim(self, opts=[]):


class Cutadapt(Trimmer):
    def __init__(self, args, **kwargs):


    def trim(self, opts=[]):



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
                options = ['-q', 15, '-t', THREADS, '-f', outfile]
            elif method == 'mem':
                options = ['@RG\tID:{input}\tPL:ILLUMINA\tLB:{input}\tSM:{input}',
                           '-t', THREADS, '-o', outfile]
        self.options = options

    def align(self), options=[]):
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

    def merge(self, opts={}):


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

    def batch(options=[], wrapcmd=None):
        command_options = [self.command]
        if not options:
            options = self.options
        command_options.extend(options)
        if not wrapcmd:
            wrapcmd = self.wrapcmd
        command_options.extend(['--wrap="%s"' % wrapcmd])

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


    def flagstat(self, opts=[]):


    def index(self, opts=[]):


    def sort(self, opts=[]):


    def view(self, opts=[]):


    def run(self, command, opts=[]):
        opts.insert(0, command)


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
    parser = argparse.ArgumentParser()
    parser.add_argument('--topdir', type=str,
                        help='Top-level directory containing paired-end sequence data input files')
    parser.add_argument('--refseq', type=str,
                        help='Reference sequence data file')
    parser.add_argument('--outdir', type=str,
                        help='Directory to write output files into')

    args = parser.parse()

    if not os.path.exists(args.topdir):
        logger.error('Top-level directory %s does not exist', args.topdir)
        return 1
    if not os.path.exists(args.refseq):
        logger.error('Reference sequence file %s does not exist', args.refseq)
        return 1
    if os.path.exists(args.outdir):
        logger.warning('Output directory %s exists, files may get overwritten', args.outdir)
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

    return 0


if __name__ == "__main__":
    sys.exit(main())
