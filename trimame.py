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
import shlex
import shutil
import tempfile

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s %(levelname)s: %(filename)s[%(lineno)d](%(funcName)s) - %(message)s',
                              datefmt='%F %T')
handler.setFormatter(formatter)
logger.addHandler(handler)

NUM_THREADS = 1


class Trimmer(object):
    def __init__(self, filemeta, args, extra):
        self.args = args
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
        '''NOTE: TrimGalore depends on cutadpt (python package)
        '''
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
        self.options = ['--paired', '--retain_unpaired', '--cores', '%s' % NUM_THREADS, '--max_n', '40', '--gzip']
        if args.outdir:
            self.options.extend(['--output_dir', args.outdir])
        self.sbatch_opts = ['-N', '1', '-c', '%s' % NUM_THREADS, '--mem=64g']

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
            job = Scheduler(self.args, self.sbatch_opts, wrapcmd)
            result = job.run()
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
        java_opts = ['-XX:ParallelGCThreads=%s' % NUM_THREADS, '-Xms24576M', '-Xmx24576M', '-XX:NewSize=6144M',
                     ' -XX:MaxNewSize=6144M']
        if '--java_opts' in extra:
            java_opts = extra[extra.index('--java_opts')+1].split()
        self.java_opts = java_opts
        self.options = ['PE', '-threads', '%s' % NUM_THREADS, '-trimlog', '{outdir}/{input}_trimm.log', '{r1}', '{r2}',
                        '{outdir}/{input}_trim_R1.fastq.gz', '{outdir}/{input}_unpaired_R1.fastq.gz',
                        '{outdir}/{input}_trim_R2.fastq.gz', '{outdir}/{input}_unpaired_R2.fastq.gz',
                        'ILLUMINACLIP:{adapters}:2:30:10', 'LEADING:3', 'TRAILING:3', 'SLIDINGWINDOW:4:15',
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
            command_options = [self.command]
            command_options.extend(self.java_opts)
            command_options.extend(['-jar', self.jarfile])
            command_options.extend(_opts)
            wrapcmd = ' '.join(command_options)
            job = Scheduler(self.args, self.sbatch_opts, wrapcmd)
            result = job.run()
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
        self.args = args
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
            self.options = ['-q', '15', '-t', '%s' % NUM_THREADS]
            self.outflag = '-f'
        elif method == 'mem':
            self.options = ['-t', '%s' % NUM_THREADS, '-R',
                            '@RG\\tID:{input}\\tPL:ILLUMINA\\tLB:{input}\\tSM:{input}']
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
                command_options[-1] = """'%s'""" % command_options[-1].format(input=label)
            command_options.append(self.refseq)
            command_options.extend([r1, r2])
            command_options.extend([self.outflag, outfile])
            command_options = ' '.join(command_options)
            sbatch_opts = ['-N', '1', '-c', '16', '--mem=64g']
            job = Scheduler(self.args, sbatch_opts, command_options)
            result = job.run()
            if result:
                self.filemeta[label].update({'alignfiles': [outfile]})
        return None


class JavaJar():
    def __init__(self, jarfile, args, extra, java_opts=[], jar_opts=[]):
        self.args = args
        self.java_opts = java_opts
        self.jar_opts = jar_opts
        self.jarfile = jarfile
        self.command = 'java'
        if '--java_path' in extra:
            java_path = extra[extra.index('--java_path')+1]
            if not os.path.isdir(java_path):
                java_path = os.path.dirname(java_path)
            self.command = os.path.join(java_path, self.command)
            if not os.path.exists(self.command):
                raise OSError('Java not found at %s', self.command)

    def run(self, java_opts=[], jar_opts=[]):
        if not java_opts:
            java_opts = self.java_opts
        if not jar_opts:
            jar_opts = self.jar_opts
        cmd_list = [self.command]
        cmd_list.extend(java_opts)
        cmd_list.extend(['-jar', self.jarfile])
        cmd_list.extend(jar_opts)
        job = Scheduler(self.args, command=cmd_list)
        result = job.run()
        return result

class Local():
    def __init__(self, options=[], command=None):
        if options:
            self.options = options
        if command:
            self.command = command

    def run(self, options=None, command=None):
        if not command:
            if not self.command:
                raise OSError('No command provided')
            command = self.command
        logger.info('Launching local command: %s', command)
        try:
            if type(command) is not list:
                #command = shlex.split(command)
                result = subprocess.run(command,
                                        shell=True,
                                        stderr=subprocess.PIPE,
                                        stdout=subprocess.PIPE,
                                        #capture_output=True,
                                        check=True)
            else:
                result = subprocess.run(command,
                                        stderr=subprocess.PIPE,
                                        stdout=subprocess.PIPE,
                                        #capture_output=True,
                                        check=True)
        except subprocess.CalledProcessError as err:
            logger.error('An error occurred running the command: %s', err)
            return None
        return result

class Sbatch():
    def __init__(self, options=[], wrapcmd=None):
        self.command = 'sbatch'
        if not options:
            options = ['-N', 1, '-c', '%s' % NUM_THREADS, '--mem=64g']
        self.options = options
        self.wrapcmd = wrapcmd

    def run(self, options=[], wrapcmd=None):
        command_options = [self.command]
        if not options:
            options = self.options
        command_options.extend(options)
        if not wrapcmd:
            wrapcmd = self.wrapcmd
        if type(wrapcmd) is list:
            wrapcmd = ' '.join(wrapcmd)
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


class Scheduler():
    def __init__(self, args, options=[], command=None):
        if args.scheduler:
            if args.scheduler == 'slurm':
                self._scheduler = Sbatch(options, command)
            elif args.scheduler == 'sge':
                raise OSError('Scheduler type grid not yet implemented')
            elif args.scheduler == 'local':
                self._scheduler = Local(options, command)
            else:
                raise OSError('Scheduler type %s not understood' % args.scheduler)
        else:
            self._scheduler = Local

    def run(self, options=None, command=None):
        return self._scheduler.run(options, command)


class Samtools():
    def __init__(self, args, extra):
        command = 'samtools'
        if '--samtools_path' in extra:
            samtools_path = extra[extra.index('--samtools_path')+1]
            if os.path.is_file(samtools_path):
                samtools_path = os.path.dirname(samtools_path)
            if not os.path.exists(os.path.join(samtools_path, command)):
                raise OSError('Samtools not found in provided path')
            command = os.path.join(samtools_path, command)
        self.command = command
        self.args = args
        ''' '''

    def faidx(self, opts=[]):
        '''Run samtools faidx with provided options
        '''
        command = 'faidx'
        return self.run(command, opts)

    def flagstat(self, opts=[]):
        '''Run samtools flagstat with provided options
        '''
        command = 'flagstat'
        return self.run(command, opts)

    def index(self, opts=[]):
        '''Run samtools index with provided options
        '''
        command = 'index'
        return self.run(command, opts)

    def sort(self, opts=[]):
        '''Run samtools sort with provided options
        '''
        command = 'sort'
        return self.run(command, opts)

    def view(self, opts=[]):
        '''Run samtools view with provided options
        '''
        command = 'view'
        return self.run(command, opts)

    def run(self, command, opts=[]):
        cmd_list = [self.command, command]
        cmd_list.extend(opts)
        job = Scheduler(self.args, command=cmd_list)
        result = job.run()
        return result


def findFiles(dname):
    logger.info('Finding paired-end FASTQ sequence files in %s', dname)
    files = glob.glob(os.path.join(dname, '*_R[12].fastq'))
    files.extend(glob.glob(os.path.join(dname, '*_R[12].fastq.gz')))
    pairdata = {}
    for f in files:
        base, _ = os.path.basename(f).split('_R')
        pairdata.setdefault(base, {})
        pairdata[base].setdefault('files', [])
        pairdata[base]['files'].append(f)
        pairdata[base]['files'].sort()
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
    parser.add_argument('--threads', type=int, default=4,
                        help='Number of CPU threads/cores to use')

    logger.info('Processing arguments')
    args, extra = parser.parse_known_args()
    logger.info('Args: %s', args)
    logger.info('Extra arguments: %s', extra)

    global NUM_THREADS
    if args.threads:
        NUM_THREADS = args.threads

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

    ##
    # Map/align the files to the reference sequence
    logger.info('Instantiating aligner...')
    aligner = Aligner(filemeta, args, extra)
    aligner.align()

    samtools = Samtools(args, extra)
    for label in filemeta:
        ##
        # Convert to BAM
        tmpfd, tmpfile = tempfile.mkstemp()
        command_opts = ['-b', '-S', '-o', tmpfile, filemeta[label]['alignfiles'][0]]
        logger.info('Converting file %s to BAM using temp file %s', filemeta[label]['alignfiles'][0], tmpfile)
        try:
            ret = samtools.view(command_opts)
        except OSError as err:
            logger.error('Error running samtools view: %s', err)
            sys.exit(1)
        if not os.path.getsize(tmpfile) > 0:
            logger.error('Error running samtools view: zero-size file generated')
            sys.exit(1)
        ##
        # Sort the BAM file
        sortfd, sortfile = tempfile.mkstemp()
        command_opts = ['-@', '%s' % NUM_THREADS, '-o', sortfile, tmpfile]
        logger.info('Sorting the BAM file %s into file %s', tmpfile, sortfile)
        try:
            ret = samtools.sort(command_opts)
        except OSError as err:
            logger.error('Error running samtools sort: %s', err)
            sys.exit(1)
        if not os.path.getsize(sortfile) > 0:
            logger.error('Error running samtools sort: zero-size file generated')
            sys.exit(1)
        sortfd = tmpfd = None
        os.unlink(tmpfile)
        outfile = label + '.bam'
        if args.outdir:
            outfile = os.path.join(args.outdir, outfile)
        shutil.move(sortfile, outfile)
        sortfile = outfile
        ##
        # Run samtools index on sorted BAM file
        command_opts = [sortfile]
        logger.info('Indexing sorted BAM file %s', sortfile)
        try:
            ret = samtools.index(command_opts)
        except OSError as err:
            logger.error('Error running samtools index: %s', err)
            sys.exit(1)
        ##
        # Run samtools flagstat on the sorted BAM file
        outfile = label + '_mapping.stat'
        if args.outdir:
            outfile = os.path.join(args.outdir, outfile)
        command_opts = ['-@', '%s' % NUM_THREADS, sortfile]
        logger.info('Generating mapping stats for sorted, indexed BAM file %s', sortfile)
        try:
            ret = samtools.flagstat(command_opts)
        except OSError as err:
            logger.info('Error running samtools flagstat: %s', err)
            sys.exit(1)
        if ret.stdout:
            with open(outfile, 'wb') as fd:
                fd.write(ret.stdout)
        filemeta[label].update({'sortindex': [sortfile]})

    logger.info('Updated file meta:\n%s', filemeta)

    ##
    # Use picard tools to MarkDuplicates
    java_opts = ['-XX:ParallelGCThreads=%s' % NUM_THREADS, '-XX:-UseGCOverheadLimit', '-Xms24576M',
                 '-Xmx24576M', '-XX:NewSize=6144M', '-XX:MaxNewSize=6144M']
    jarfile = 'picard.jar'
    if '--picard_path' in extra:
        picard_path = extra[extra.index('--picard_path')+1]
        if not os.path.isdir(picard_path):
            picard_path = os.path.dirname(picard_path)
        jarfile = os.path.join(picard_path, jarfile)
    if not os.path.exists(jarfile):
        raise OSError('Picard jar file not found: %s', jarfile)
    for label in filemeta:
        input_filename = filemeta[label]['sortindex'][0]
        output_filename = label + '_dedup.bam'
        metrics_filename = label + '_dedup.metrics'
        if args.outdir:
            output_filename = os.path.join(args.outdir, output_filename)
            metrics_filename = os.path.join(args.outdir, metrics_filename)
        picard_opts = ['MarkDuplicates', 'I=%s' % input_filename, 'O=%s' % output_filename, 'M=%s' % metrics_filename, 'AS=true']
        logger.debug('Mark duplicates with: java %s -jar %s %s', ' '.join(java_opts), jarfile, ' '.join(picard_opts))
        jvm = JavaJar(jarfile, args, extra, java_opts=java_opts, jar_opts=picard_opts)
        try:
            ret = jvm.run()
        except OSError as err:
            logger.error('Error running java jar: %s', err)
            sys.exit(1)

        ##
        # Use sametools index on resultant BAM
        command_opts = [output_filename]
        logger.info('Indexing de-duped BAM file %s', output_filename)
        try:
            ret = samtools.index(command_opts)
        except OSError as err:
            logger.error('Error running samtools index: %s', err)
            sys.exit(1)
        filemeta[label].setdefault('dedup', [output_filename])

    ##
    # Use GenomeToolkit HaplotypeCaller with genome reference
    if not os.path.exists(args.refseq + 'fai'):
        logger.info('Indexing refenence FASTA...')
        command_opts = [args.refseq]
        try:
            ret = samtools.faidx(command_opts)
        except OSError as err:
            logger.error('Error running samtolos faidx: %s', err)
            sys.exit(1)
    if not os.path.exists(os.path.splitext(args.refseq)[0] + '.dict'):
        logger.info('Generating sequence dictionary for reference FASTA...')
        output_filename = os.path.splitext(args.refseq)[0] + '.dict'
        picard_opts = ['CreateSequenceDictionary', 'R=%s' % args.refseq, 'O=%s' % output_filename]
        logger.debug('Create sequence dictionary with: java %s -jar %s %s', ' '.join(java_opts), jarfile, ' '.join(picard_opts))
        jvm = JavaJar(jarfile, args, extra, java_opts=java_opts, jar_opts=picard_opts)
        try:
            ret = jvm.run()
        except OSError as err:
            logger.error('Error running java jar: %s', err)
            sys.exit(1)
        
    jarfile = 'GenomeAnalysisTK.jar'
    if '--gatk_path' in extra:
        gatk_path = extra[extra.index('--gatk_path')+1]
        if not os.path.isdir(gatk_path):
            gatk_path = os.path.dirname(gatk_path)
        jarfile = os.path.join(gatk_path, jarfile)
    if not os.path.exists(jarfile):
        raise OSError('GenomeAnalysisToolkit jar not found at %s', jarfile)
    for label in filemeta:
        input_filename = filemeta[label]['dedup'][0]
        output_filename = label + '_snps.indels.vcf'
        if args.outdir:
            output_filename = os.path.join(args.outdir, output_filename)
        gatk_opts = ['-T', 'HaplotypeCaller', '-R', args.refseq, '-stand_call_conf', '30',
                     '-I', input_filename, '-o', output_filename]
        logger.debug('HaplotypeCaller: java %s -jar %s %s', ' '.join(java_opts), jarfile, ' '.join(gatk_opts))
        jvm = JavaJar(jarfile, args, extra, java_opts=java_opts, jar_opts=gatk_opts)
        try:
            ret = jvm.run()
        except OSError as err:
            logger.error('Error running java jar: %s', err)
            sys.exit(1)

    logger.info('Complete')
    return 0


if __name__ == "__main__":
    sys.exit(main())
