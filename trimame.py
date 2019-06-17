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
import copy

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s %(levelname)s: %(filename)s[%(lineno)d](%(funcName)s) - %(message)s',
                              datefmt='%F %T')
handler.setFormatter(formatter)
logger.addHandler(handler)


class Trimmer(object):
    def __init__(self, filemeta, args, extra):
        for f in list(sum([filemeta[x].get('files') for x in filemeta], [])):
            if not os.path.exists(f):
                raise IOError('File %s does not exist' % r1)
        self.filemeta = filemeta
        self.adapters = None
        if '--adapters' in extra:
            adapters = extra[extra.index('--adapters')+1]
            if not os.path.exists(adapters):
                raise IOError('File %s does not exist' % adapters)
            self.adapters = adapters
        self.command = None
        self.options = []
        self.sbatch_opts = []

    def trim(self, opts=[]):
        ''' '''
        return None
                                    

class TrimGalore(Trimmer):
    def __init__(self, filemeta, args, extra):
        ''' '''
        super(TrimGalore, self).__init__(filemeta, args, extra)
        if '--trimmer_path' in extra:
            command = extra[extra.index('--trimmer_path')+1]
            command = os.path.join(command, 'trim_galore')
            if os.path.exists(command):
                self.command = command
            else:
                raise IOError('Command %s not found' % command)
        else:
            self.command = 'trim_galore'
        self.options = ['--paired', '--retain_unpaired', '--cores', '4', '--max_n', '40', '--gzip']
        if args.outdir:
            self.options.extend(['--outdir', args.outdir])
        self.sbatch_opts = ['-N', '1', '-c', '16', '--mem=64g']

    def trim(self, options=[]):
        if not options:
            options = self.options
        logger.debug('Trimming files with options: %s', options)
        for label in self.filemeta:
            r1, r2 = self.filemeta[label].get('files')
            command_options = [self.command]
            command_options.extend(options)
            command_options.extend([r1, r2])
            wrapcmd = ' '.join(command_options)
            sbatch = Sbatch(self.sbatch_opts, wrapcmd)
            result = sbatch.run()
            if result:
                outf1, outf2 = r1, r2
                if '--basename' in self.options:
                    basename = self.options[self.options.index('--basename')+1]
                    basename = os.path.basename(basename)
                    outf1 = basename + '_val_1.fq'
                    outf2 = basename + '_val_2.fq'
                if '--outdir' in self.options:
                    outdir = self.options[self.options.index('--outdir')+1]
                    outf1 = os.path.join(outdir, os.path.basename(outf1))
                    outf2 = os.path.join(outdir, os.path.basename(outf2))
                self.filemeta[label].update({'trimfiles': [outf1, outf2]})
        return None


class Trimomatic(Trimmer):
    def __init__(self, filemeta, args, extra):
        ''' '''
        super(Trimomatic, self).__init__(filemeta, args, extra)
        command = 'java'
        if '--java_path' in extra:
            javapath = extra[extra.index('--java_path')+1]
            command = os.path.join(javapath, command)
            if os.path.exists(command):
                self.command = command
            else:
                raise IOError('Command %s not found' % command)
        else:
            self.command = command
        jarfile = 'trimmomatic.jar'
        if '--trimmer_path' in extra:
            jarpath = extra[extra.index('--trimmer_path')+1]
            jarfile = os.path.join(jarpath, 'trimmomatic*.jar')
            jarfile = glob.glob(jarfile)[0]  # find first match
            if not jarfile:
                raise IOError('No trimmomatic*.jar found in %s' % jarpath)
        self.jarfile = jarfile
        java_opts = ['-XX:ParallelGCThreads=8', '-Xms24576M', '-Xmx24576M', '-XX:NewSize=6144M',
                     ' -XX:MaxNewSize=6144M']
        if '--java_opts' in extra:
            java_opts = extra[extra.index('--java_opts')+1].split()
        self.java_opts = java_opts
        self.options = ['PE', '-threads', '8', '-trimlog', '{input}_trimm.log', '{r1}', '{r2}',
                        '{outdir}/{input}_trim_R1.fastq.gz', '{outdir}/{input}_unpaired_R1.fastq.gz',
                        '{outdir}/{input}_trim_R2.fastq.gz', '{outdir}/{input}_unpaired_R2.fastq.gz',
                        '#ILLUMINACLIP:{adapters}:2:30:10', 'LEADING:3', 'TRAILING:3', 'SLIDINGWINDOW:4:15',
                        'MINLEN:36']
        if args.outdir:
            self.outdir = args.outdir

    def trim(self, options=[]):
        ''' '''
        if not options:
            options = self.options
        logger.debug('Trimming files with options: %s', options)
        for label in self.filemeta:
            r1, r2 = self.filemeta[label].get('files')
            argmap = {'adapters': self.adapters, 'input': label, 'r1': r1, 'r2': r2}
            outdir = os.path.dirname(r1)
            if self.outdir:
                outdir = self.outdir
            argmap.update({'outdir': outdir})
            _opts = copy.deepcopy(options)
            for opt in _opts[:]:
                idx = options.index(opt)
                opt = opt.format(**argmap)
                _opts[idx] = opt
            command_options = [self.command, self.jarfile]
            command_options.extend(self.java_opts)
            command_options.extend(_opts)
            wrapcmd = ' '.join(command_options)
            sbatch = Sbatch(self.sbatch_opts, wrapcmd)
            result = sbatch.run()
            if result:
                outf1 = _opts[7]
                outf2 = _opts[9]
                self.filemeta[label].update({'trimfiles': [outf1, outf2]})
        return None


class Cutadapt(Trimmer):
    def __init__(self, filemeta, args, extra):
        ''' '''
        super(Cutadapt, self).__init__(filemeta, args, extra)

    def trim(self, opts=[]):
        ''' '''
        return None


class Aligner():
    def __init__(self, filemeta, args, extra):
        if '--bwa_path' in extra:
            command = extra[extra.index('--bwa_path')+1]
            command = os.path.join(command, 'bwa')
            if os.path.exists(command):
                self.command = command
            else:
                raise IOError('Command %s not found' % command)
        else:
            self.command = 'bwa'
        flist = [args.refseq]
        flist.extend(list(sum([filemeta[x].get('trimfiles') for x in filemeta], [])))
        for f in flist:
            if not os.path.exists(f):
                raise IOError('File %s does not exist' % f)
        self.filemeta = filemeta
        self.refseq = args.refseq
        method = 'aln'
        if '--align_method' in extra:
            method = extra[extra.index('--align_method')+1]
            if method not in ['aln', 'mem']:
                raise OSError('Method %s not allowed for alignment' % method)
        self.method = method
        self.options = []
        if method == 'aln':
            self.options = ['-q', '15', '-t', '16']
            self.outflag = '-f'
        elif method == 'mem':
            self.options = ['@RG\tID:{input}\tPL:ILLUMINA\tLB:{input}\tSM:{input}',
                            '-t', '16']
            self.outflag = '-o'
        if args.outdir:
            self.outdir = args.outdir

    def align(self, options=[]):
        if not options:
            options = self.options[:]
        for label in self.filemeta:
            r1, r2 = self.filemeta[label].get('trimfiles')
            outfile = '%s.bam' % label
            #outf2 = '%s_trim_R2.fastq.gz' % label
            if self.outdir:
                outfile = os.path.join(self.outdir, outfile)
                #outf2 = os.path.join(self.outdir, outf2)
            command_options = [self.command, self.method]
            command_options.extend(options)
            if self.method == 'mem':
                command_options[2] = command_options[2].format(input=label)
            command_options.append(self.refseq)
            command_options.extend([r1, r2])
            command_options.extend([self.outflag, outfile])
            wrapcmd = ' '.join(command_options)
            sbatch_opts = ['-N', '1', '-c', '16', '--mem=64g']
            batch = Sbatch(sbatch_opts, wrapcmd)
            result = batch.run()
        return None


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

    def run(self, options=[], wrapcmd=None):
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
    logger.info('Finding paired-end FASTQ sequence files in %s', dname)
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
    schedulers = {
        'local': 'Local execution',
        'sge': 'Sun Grid Engine',
        'slurm': 'SLURM'
    }
    trimopts = {
        'trimgalore': TrimGalore,
        'trimmomatic': Trimomatic,
        'cutadapt': Cutadapt
    }
    parser = argparse.ArgumentParser()
    parser.add_argument('--scheduler', type=str, choices=schedulers.keys(), default='local',
                        help='Scheduler to use')
    parser.add_argument('--topdir', type=str, required=True,
                        help='Top-level directory containing paired-end sequence data input files')
    parser.add_argument('--refseq', type=str, required=True,
                        help='Reference sequence data file')
    parser.add_argument('--outdir', type=str, required=True,
                        help='Directory to write output files into')
    parser.add_argument('--trimmer', type=str, choices=trimopts.keys(), default='trimmomatic',
                        help='Tool to use for sequence trimming')

    logger.info('Processing arguments')
    args, extra = parser.parse_known_args()
    logger.info('Extra arguments: %s', extra)

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

    ##
    # Find file pairs
    logger.info('Finding input data')
    filemeta = findFiles(args.topdir)
    if not filemeta:
        logger.error('No paired-end sequence data files found in %s', args.topdir)
        return 1
    for s in filemeta:
        if not filemeta[s].get('files'):
            filemeta.pop(s)
    if not filemeta:
        logger.error('No paired-end sequence data files found in %s', args.topdir)
        return 1
    logger.info('Found %d sets of files', len(filemeta))

    ##
    # Trim each file in each pair
    logger.info('Instantiating trimmer...')
    trimmer = trimopts[args.trimmer](filemeta, args, extra)
    logger.info('Instantiated %s trimmer object', trimmer)

    trimmer.trim()

    logger.info('Updated file meta:\n%s', filemeta)

    ##
    # Map/align the files to the reference sequence
    logger.info('Instantiating aligner...')
    aligner = Aligner(filemeta, args, extra)
    aligner.align()

    ##
    # Merge the file pairs    
    #logger.info('Instantiating merger...')
    #merger = Merger()

    logger.info('Complete')
    return 0


if __name__ == "__main__":
    sys.exit(main())
