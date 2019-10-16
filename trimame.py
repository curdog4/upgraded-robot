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
import datetime
import re
import shlex
import shutil
import stat
import tempfile
import textwrap
import time

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s %(levelname)s: %(filename)s[%(lineno)d](%(funcName)s) - %(message)s',
                              datefmt='%F %T')
handler.setFormatter(formatter)
logger.addHandler(handler)

NUM_THREADS = 1
JOBID_RE = re.compile(r'^\d+(?:[_.]\d+)?$')
DEFAULT_SBATCH_OPTS = ['--parsable', '-N', '1', '-c', '%s' % NUM_THREADS, '--mem=64g']

class Trimmer(object):
    def __init__(self, filemeta, args, extra):
        self.args = args
        for f in list(sum([filemeta[x].get('files') for x in filemeta], [])):
            if not os.path.exists(f):
                raise IOError('File %s does not exist' % f)
        self.filemeta = filemeta
        self.adapters = None
        if '--adapters' in extra:
            adapters = extra[extra.index('--adapters')+1]
            if not os.path.exists(adapters):
                raise IOError('File %s does not exist' % adapters)
            self.adapters = adapters
        self.command = None
        self.options = []
        self.sbatch_opts = DEFAULT_SBATCH_OPTS

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
        if '--cutadapt_path' in extra:
            cutadapt_path = extra[extra.index('--cutadapt_path')+1]
            if os.path.isdir(cutadapt_path):
                cutadapt_path = os.path.join(cutadapt_path, 'cutadapt')
            if not os.path.exists(cutadapt_path):
                raise IOError('Cutadapt not found at %s' % cutadapt_path)
            self.options.extend(['--path_to_cutadapt', cutadapt_path])
        if args.outdir:
            self.options.extend(['--output_dir', args.outdir])

    def trim(self, options=[]):
        if not options:
            options = self.options
        logger.debug('Trimming files with options: %s', options)
        joblist = []
        scheduler = Scheduler(self.args, self.sbatch_opts, '')
        for label in self.filemeta:
            r1, r2 = self.filemeta[label].get('files')
            command_options = [self.command]
            command_options.extend(options)
            command_options.extend([r1, r2])
            wrapcmd = ' '.join(command_options)
            joblist.append({'label': label, 'command': wrapcmd})
            #job = Scheduler(self.args, self.sbatch_opts, wrapcmd)
            #result = job.run()
        results = scheduler.run(joblist=joblist)
        if not results:
            raise OSError('Failed to trim files')
        rlabels = list(set([x.get('label') for x in results]))
        for label in self.filemeta:
            if not label in rlabels:
                raise OSError('Failed to trim files for %s' % label)
        for result in results:
            label = result.get('label')
            _result = result.get('result')
            if _result.returncode:
                raise OSError('Trim of %s failed with %s: %s' % (label, _result.returncode, _result.stderr or _result.stdout))
            if _result.stderr:
                logger.warning('Trim of %s returned with errors: %s', label, _result.stderr)
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
        joblist = []
        scheduler = Scheduler(self.args, self.sbatch_opts, '')
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
            joblist.append({'label': label, 'command': wrapcmd})
            #job = Scheduler(self.args, self.sbatch_opts, wrapcmd)
            #result = job.run()
        results = scheduler.run(joblist=joblist)
        if not results:
            raise OSError('Failed to trim files')
        rlabels = list(set([x.get('label') for x in results]))
        for label in self.filemeta:
            if not label in rlabels:
                raise OSError('Failed to trim files for %s' % label)
        for result in results:
            label = result.get('label')
            _result = result.get('result')
            if _result.returncode:
                raise OSError('Trim of %s failed with %s: %s' % (label, _result.returncode, _result.stderr or _result.stdout))
            if _result.stderr:
                raise OSError('Trim of %s returned with errors: %s' % (label, _result.stderr))
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
        self.sbatch_opts = DEFAULT_SBATCH_OPTS

    def align(self, options=[]):
        if not options:
            options = self.options[:]
        joblist = []
        scheduler = Scheduler(self.args, self.sbatch_opts, '')
        command_options = [self.command, self.method]
        command_options.extend(options)
        if self.method == 'mem':
            for label in self.filemeta:
                r1, r2 = self.filemeta[label].get('trimfiles')
                outfile = '%s.bam' % label
                if self.outdir:
                    outfile = os.path.join(self.outdir, outfile)
                cmd_opts = command_options[:]
                cmd_opts[-1] = """'%s'""" % cmd_opts[-1].format(input=label)
                cmd_opts.append(self.refseq)
                cmd_opts.extend([r1, r2])
                cmd_opts.extend([self.outflag, outfile])
                cmd_opts = ' '.join(cmd_opts)
                joblist.append({'label': label, 'command': cmd_opts})
                #job = Scheduler(self.args, sbatch_opts, command_options)
                #result = job.run()
                #if result:
                #    self.filemeta[label].update({'alignfiles': [outfile]})
            results = scheduler.run(joblist=joblist)
            if not results:
                raise OSError('Failed to align files')
            rlabels = list(set([x.get('label') for x in results]))
            for label in self.filemeta:
                if not label in rlabels:
                    raise OSError('Failed to align files for %s' % label)
            for result in results:
                label = result.get('label')
                _result = result.get('result')
                if _result.returncode:
                    raise OSError('Align of %s failed with %s: %s' % (label, _result.returncode, _result.stderr or _result.stdout))
                if _result.stderr:
                    raise OSError('Align of %s returned with errors: %s' % (label, _result.stderr))
                outfile = shlex.split(result.get('command'))[-1]
                self.filemeta[label].update({'alignfiles': [outfile]})
        else:
            for label in self.filemeta:
                for trimfile in self.filemeta[label].get('trimfiles'):
                    outfile = '%s_R%d.sai' % (label, self.filemeta[label]['trimfiles'].index(trimfile) + 1)
                    if self.outdir:
                        outfile = os.path.join(self.outdir, outfile)
                    cmdopts = command_options[:]
                    cmdopts.extend([self.outflag, outfile, self.refseq, trimfile])
                    joblist.append({'label': label, 'command': cmdopts})
                    #job = Scheduler(self.args, sbatch_opts, cmdopts)
                    #result = job.run()
                    #if not result:
                    #    raise OSError('Unable to run bwa aln to align the sequence file %s' % trimfile)
                    #self.filemeta[label].setdefault('alignfiles', [])
                    #self.filemeta[label]['alignfiles'].append(outfile)
            results = scheduler.run(joblist=joblist)
            if not results:
                raise OSError('Failed to align files')
            rlabels = list(set([x.get('label') for x in results]))
            for label in self.filemeta:
                if not label in rlabels:
                    raise OSError('Failed to align files for %s' % label)
            for result in results:
                label = result.get('label')
                _result = result.get('result')
                if _result.returncode:
                    raise OSError('Align of %s failed with %s: %s' % (label, _result.returncode, _result.stderr or _result.stdout))
                if _result.stderr:
                    logger.warning('Align of %s returned with errors: %s', label, _result.stderr)
                outfile = shlex.split(result.get('command'))[-3]
                self.filemeta[label].setdefault('alignfiles', [])
                self.filemeta[label]['alignfiles'].append(outfile)

            joblist = []
            for label in self.filemeta:
                outfile = '%s.bam' % label
                if self.outdir:
                    outfile = os.path.join(self.outdir, outfile)
                command_options = [self.command, 'sampe', '-f', outfile, '-r',
                                   '@RG\\tID:{input}\\tPL:ILLUMINA\\tLB:{input}\\tSM:{input}'.format(input=label),
                                   self.refseq]
                for sai in self.filemeta[label].get('alignfiles'):
                    command_options.append(sai)
                for fq in self.filemeta[label].get('trimfiles'):
                    command_options.append(fq)
                joblist.append({'label': label, 'command': command_options})
                #job = Scheduler(self.args, sbatch_opts, command_options)
                #result = job.run()
                #if not result:
                #    raise OSError('Unable to run bwa sampe to merge the alignments')
                #self.filemeta[label]['alignfiles'] = [outfile]
            results = scheduler.run(joblist=joblist)
            if not results:
                raise OSError('Failed to merge align files')
            rlabels = list(set([x.get('label') for x in results]))
            for label in self.filemeta:
                if not label in rlabels:
                    raise OSError('Failed to merge align files for %s' % label)
            for result in results:
                label = result.get('label')
                outfile = shlex.split(result.get('command'))[3]
                self.filemeta[label]['alignfiles'] = [outfile]
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
                raise OSError('Java not found at %s' % self.command)
        self.sbatch_opts = DEFAULT_SBATCH_OPTS

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

    def run(self, joblist=[]):
        if not joblist:
                raise OSError('No jobs provided')
        resultlist = []
        for jobcommand in joblist:
            label, command = jobcommand.values()
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
                result = None
            resultlist.append({'label': label, 'result': result})
        return resultlist

class Sbatch():
    def __init__(self, options=[], wrapcmd=None):
        self.command = 'sbatch'
        if not options:
            #options = ['--parsable', '-N', 1, '-c', '%s' % NUM_THREADS, '--mem=64g']
            options = DEFAULT_SBATCH_OPTS
        else:
            if '--parsable' not in options:
                if type(options) is list:
                    options.append('--parsable')
                else:
                    options += ' --parsable'
        self.options = options
        self.wrapcmd = wrapcmd

    def run(self, joblist=[]):
        if not joblist:
            raise OSError('No jobs provided')
        resultslist = []
        jobids = []
        for jobcommand in joblist:
            label, wrapcmd = jobcommand.values()    
            command_options = [self.command]
            sbatch_opts = self.options[:]
            sbatch_opts.extend(['-o', label + '_out.log', '-e', label + '_err.log'])
            command_options.extend(sbatch_opts)
            command_options.extend(['-J', label])
            if type(wrapcmd) is list:
                wrapcmd = ' '.join(wrapcmd)
            command_options.extend(['--wrap=%s' % wrapcmd])

            logger.info('Launching sbatch with: %s', command_options)
            result = self.wrapsubprocess(command_options)
            if not result:
                raise OSError('No result from sbatch')
            #logger.info('Sbatch process returned %d', result.returncode)
            if result.stdout:
                #logger.info('Stdout from sbatch: %s', result.stdout.decode('utf-8'))
                jobid = result.stdout.decode('utf-8').strip().split(':')[0]
            if result.stderr:
                logger.info('Stderr from sbatch: %s', result.stderr.decode('utf-8'))

            if not jobid:
                raise OSError('No job ID found')
            if not JOBID_RE.match(jobid):
                raise OSError('Unable to match job ID from sbatch output')
            logger.info('Started job %s', jobid)
            jobids.append({'jobid': jobid, 'label': label, 'retries': 5,
                           't_start': datetime.datetime.now(), 'wrapcmd': wrapcmd})

        logger.info('Checking job status, waiting for completion...')
        max_retries = 5
        failed = []
        complete = []
        while jobids:
            time.sleep(1)
            for job in jobids[:]:
                command = ['sacct', '-l', '-n', '-P', '-j', job['jobid']]
                result = self.wrapsubprocess(command)
                status = None
                if not result:
                    logger.error('No result from sbatch for job ID %s', job['jobid'])
                    failed.append(job)
                    jobids.remove(job)
                    continue
                #logger.info('Sacct process returned %d for job ID %s', result.returncode,
                #            job['jobid'])
                if result.stdout:
                    #logger.info('Stdout from sbatch for job ID %s: %s', job['jobid'],
                    #            result.stdout.decode('utf-8'))
                    status = result.stdout.decode('utf-8').strip().splitlines()
                if result.stderr:
                    logger.info('Stderr from sbatch for job ID %s: %s', job['jobid'], 
                                result.stderr.decode('utf-8'))
                if not status:
                    logger.error('Status not available for job %s', job['jobid'])
                    if job['retries']:
                        jobids[jobids.index(job)]['retries'] -= 1
                        continue
                    logger.error('Exceeded retries for job ID %s', job['jobid'])
                    failed.append(job)
                    jobids.remove(job)
                    continue
                if len(status) > 1:
                    logger.warning('Found %d status lines for %s. Using the first...',
                                   len(status), jobid)
                jstatus = status[0].split('|')[23]
                if jstatus in ('COMPLETED', 'FAILED'):
                    t_end = datetime.datetime.now()
                    logger.info('Job %s finished after %s', job['jobid'], t_end - job['t_start'])
                    returncode = int(status[0].split('|')[24].split(':')[0])
                    cmdopts = shlex.split(job['wrapcmd'])[1:]
                    stdout = b''
                    stderr = b''
                    for opt in self.options:
                        stderr = open(job['label'] + '_err.log', 'rb').read()
                        stdout = open(job['label'] + '_out.log', 'rb').read()
                    result = subprocess.CompletedProcess(args=cmdopts, returncode=returncode)
                    result.stderr = stderr
                    result.stdout = stdout
                    resultslist.append({'label': job['label'], 'result': result, 'command': job['wrapcmd']})
                    complete.append(job)
                    jobids.remove(job)
        
        return resultslist

    @classmethod
    def wrapsubprocess(cls, command):
        #logger.info('Attempt to launch command: %s', ' '.join(command))
        try:
            result = subprocess.run(command,
                                    stderr=subprocess.PIPE,
                                    stdout=subprocess.PIPE,
                                    #capture_output=True,
                                    check=True)
        except subprocess.CalledProcessError as err:
            logger.error('Subprocess returned a CalledProcessError: %s', err)
            return None
        if not result:
            logger.error('No result returned')
            return None
        return result


class Scheduler():
    def __init__(self, args, options=[], command=None):
        self.command_queue = []
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

    def add(self, command):
        self.command_queue.append(command)

    def clear(self):
        self.command_queue = []

    def remove(self, command):
        try:
            self.command_queue.remove(command)
        except ValueError:
            pass

    def run(self, joblist=[]):
        if not joblist and self.command_queue:
            return self._scheduler.run(self.command_queue)
        return self._scheduler.run(joblist)


class Samtools():
    def __init__(self, args, extra):
        command = 'samtools'
        if '--samtools_path' in extra:
            samtools_path = extra[extra.index('--samtools_path')+1]
            if os.path.isfile(samtools_path):
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
    global DEFAULT_SBATCH_OPTS
    if args.threads:
        NUM_THREADS = args.threads
        DEFAULT_SBATCH_OPTS[4] = str(args.threads)

    logger.info('Verifying arguments')
    if not os.path.exists(args.topdir):
        logger.error('Top-level directory %s does not exist', args.topdir)
        return 1
    if not os.path.exists(args.refseq):
        logger.error('Reference sequence file %s does not exist', args.refseq)
        return 1
    if args.outdir:
        args.outdir = os.path.abspath(args.outdir)
    if os.path.exists(args.outdir):
        logger.warning('Output directory %s exists, files may get overwritten', args.outdir)
    else:
        os.mkdir(args.outdir)

    ##
    # Find file pairs
    logger.info('Finding input data')
    filemeta = findFiles(os.path.abspath(args.topdir))
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

    logger.debug('Updated filemeta:\n%s', filemeta)

    scheduler = Scheduler(args, DEFAULT_SBATCH_OPTS, '')

    samtools = 'samtools'
    if '--samtools_path' in extra:
        samtools_path = extra[extra.index('--samtools_path')+1]
        if os.path.isfile(samtools_path):
            samtools_path = os.path.dirname(samtools_path)
        if not os.path.exists(os.path.join(samtools_path, samtools)):
            raise OSError('Samtools not found in provided path')
        samtools = os.path.join(samtools_path, samtools)
    joblist = []
    for label in filemeta:
        command_opts = [samtools]
        alignfile = filemeta[label]['alignfiles'][0]
        submap = {'alignfile': alignfile, 'bamfile': alignfile, 'samtools': samtools}
        script_data = '''#!/usr/bin/bash'''
        if not alignfile.endswith('.bam'):
            ##
            # Convert to BAM
            bamfile = os.path.splitext(alignfile)[0] + '.bam'
            submap.update({'bamfile': bamfile})
            script_data += '''
            {samtools} view -b -S -o {bamfile} {alignfile}
            RC=$?
            if [ "$RC" != 0 -o ! -f {bamfile} ]; then
                echo "ERROR: failed to convert alignment file to BAM" >&2
                exit $RC
            elif [ ! -s {bamfile} ]; then
                echo "ERROR: zero-sized BAM file generated" >&2
                exit 1
            fi
            '''.format(**submap)

        ##
        # Sort the BAM file
        sortfile = label + '.bam'
        if args.outdir:
            sortfile = os.path.join(args.outdir, sortfile)
        submap.update({'sortfile': sortfile, 'num_threads': '%s' % NUM_THREADS})
        script_data += '''
        {samtools} sort -@ {num_threads} -o {sortfile} {bamfile}
        RC=$?
        if [ "$RC" != 0 -o ! -f {sortfile} ]; then
            echo "ERROR: failed to sort the BAM file" >&2
            exit $RC
        elif [ ! -s {sortfile} ]; then
            echo "ERROR: zero-sized sorted BAM file generated" >&2
            exit 1
        fi
        {samtools} index {sortfile}
        RC=$?
        if [ "$RC" != 0 ]; then
            echo "ERROR: failed to index the sorted BAM file" >&2
            exit $RC
        fi
        '''.format(**submap)
        ##
        # Run samtools flagstat on the sorted BAM file
        outfile = label + '_mapping.stat'
        if args.outdir:
            outfile = os.path.join(args.outdir, outfile)
        submap.update({'outfile': outfile})
        script_data += '''
        {samtools} flagstat -@ {num_threads} {sortfile} > {outfile}
        RC=$?
        if [ "$RC" != 0 -o ! -s {outfile} ]; then
            echo "ERROR: failed running samtools flagstat" >&2
            exit $RC
        fi
        '''.format(**submap)
        script_file = label + '_samtools_convert.sh'
        if args.outdir:
            script_file = os.path.join(args.outdir, script_file)
        with open(script_file, 'w') as fd:
            fd.write(textwrap.dedent(script_data))
        os.chmod(script_file, stat.S_IREAD | stat.S_IWRITE | stat.S_IEXEC)
        joblist.append({'label': label, 'command': script_file})

    results = scheduler.run(joblist=joblist)
    if not results:
        raise OSError('Failed to merge align files')
    rlabels = list(set([x.get('label') for x in results]))
    for label in filemeta:
        if not label in rlabels:
           raise OSError('Failed to merge align files for %s' % label)
    for result in results:
        label = result.get('label')
        _result = result.get('result')
        if _result.returncode:
            raise OSError('Merge align of %s failed with %s: %s' % (label, _result.returncode, _result.stderr or _result.stdout))
        if _result.stderr:
            raise OSError('Merge align of %s returned with errors: %s' % (label, _result.stderr))
        sortfile = label + '.bam'
        if args.outdir:
            sortfile = os.path.join(args.outdir, sortfile)
        if not os.path.exists(sortfile):
            raise OSError('Failed to find sorted BAM file %s for %s' % (sortfile, label))
        filemeta[label].update({'sortindex': [sortfile]})

    ##
    # Use picard tools to MarkDuplicates
    java = 'java'
    if '--java_path' in extra:
        java_path = extra[extra.index('--java_path')+1]
        if not os.path.isdir(java_path):
            java_path = os.path.dirname(java_path)
        java = os.path.join(java_path, java)
        if not os.path.exists(self.command):
            raise OSError('Java not found at %s' % java)
    java_opts = ['-XX:ParallelGCThreads=%s' % NUM_THREADS, '-XX:-UseGCOverheadLimit', '-Xms24576M',
                 '-Xmx24576M', '-XX:NewSize=6144M', '-XX:MaxNewSize=6144M']
    jarfile = 'picard.jar'
    if '--picard_path' in extra:
        picard_path = extra[extra.index('--picard_path')+1]
        if not os.path.isdir(picard_path):
            picard_path = os.path.dirname(picard_path)
        jarfile = os.path.join(picard_path, jarfile)
    if not os.path.exists(jarfile):
        raise OSError('Picard jar file not found: %s' % jarfile)
    joblist = []
    command = [java]
    command.extend(java_opts)
    command.extend(['-jar', jarfile])
    for label in filemeta:
        input_filename = filemeta[label]['sortindex'][0]
        ##
        # first, need to clean potentially unmapped reads
        cleaned_filename = os.path.splitext(input_filename)[0] + '_cleaned.bam'
        picard_opts = ['CleanSam', 'I=%s' % input_filename, 'O=%s' % cleaned_filename]
        cmd_opts = command[:]
        cmd_opts.extend(picard_opts)
        joblist.append({'label': label, 'command': cmd_opts})
    results = scheduler.run(joblist=joblist)
    if not results:
        raise OSError('Failed to clean align files')
    rlabels = list(set([x.get('label') for x in results]))
    for label in filemeta:
        if not label in rlabels:
           raise OSError('Failed to clean align files for %s' % label)
    for result in results:
        label = result.get('label')
        _result = result.get('result')
        if _result.returncode:
            raise OSError('Clean align of %s failed with %s: %s' % (label, _result.returncode, _result.stderr or _result.stdout))
        if _result.stderr:
            logger.warning('Clean align of %s returned with errors: %s', label, _result.stderr)
        input_filename = filemeta[label]['sortindex'][0]
        cleaned_filename = os.path.splitext(input_filename)[0] + '_cleaned.bam'
        if not os.path.exists(cleaned_filename):
            raise OSError('Failed to find cleaned align file %s for %s' % (cleaned_filename, label))

    joblist = []
    for label in filemeta:
        input_filename = filemeta[label]['sortindex'][0]
        cleaned_filename = os.path.splitext(input_filename)[0] + '_cleaned.bam'
        output_filename = label + '_dedup.bam'
        metrics_filename = label + '_dedup.metrics'
        if args.outdir:
            output_filename = os.path.join(args.outdir, output_filename)
            metrics_filename = os.path.join(args.outdir, metrics_filename)
        picard_opts = ['MarkDuplicates', 'I=%s' % cleaned_filename, 'O=%s' % output_filename, 'M=%s' % metrics_filename, 'AS=true']
        cmd_opts = command[:]
        cmd_opts.extend(picard_opts)
        joblist.append({'label': label, 'command': cmd_opts})
    results = scheduler.run(joblist=joblist)
    if not results:
        raise OSError('Failed to mark duplicates')
    rlabels = list(set([x.get('label') for x in results]))
    for label in filemeta:
        if not label in rlabels:
           raise OSError('Failed to mark duplicates for %s' % label)
    for result in results:
        label = result.get('label')
        _result = result.get('result')
        if _result.returncode:
            raise OSError('Mark duplicates of %s failed with %s: %s' % (label, _result.returncode, _result.stderr or _result.stdout))
        if _result.stderr:
            logger.warning('Mark duplicates of %s returned with errors: %s', label, _result.stderr)
        output_filename = label + '_dedup.bam'
        if args.outdir:
            output_filename = os.path.join(args.outdir, output_filename)
        if not os.path.exists(output_filename):
            raise OSError('Failed to find cleaned align file %s for %s' % (output_filename, label))

    joblist = []
    for label in filemeta:
        output_filename = label + '_dedup.bam'
        if args.outdir:
            output_filename = os.path.join(args.outdir, output_filename)
        ##
        # Use sametools index on resultant BAM
        command_opts = [samtools, 'index', output_filename]
        joblist.append({'label': label, 'command': command_opts})
    results = scheduler.run(joblist=joblist)
    if not results:
        raise OSError('Failed to mark duplicates')
    rlabels = list(set([x.get('label') for x in results]))
    for label in filemeta:
        if not label in rlabels:
           raise OSError('Failed to mark duplicates for %s' % label)
    for result in results:
        label = result.get('label')
        _result = result.get('result')
        if _result.returncode:
            raise OSError('Samtools index of %s failed with %s: %s' % (label, _result.returncode, _result.stderr or _result.stdout))
        if _result.stderr:
            raise OSError('Samtools index of %s returned with errors: %s' % (label, _result.stderr))
        output_filename = label + '_dedup.bam'
        if args.outdir:
            output_filename = os.path.join(args.outdir, output_filename)
        if not os.path.exists(output_filename):
            raise OSError('Failed to find Samtools indexed file %s for %s' % (output_filename, label))
        filemeta[label].setdefault('dedup', [output_filename])

    ##
    # Use GenomeToolkit HaplotypeCaller with genome reference
    if not os.path.exists(args.refseq + '.fai'):
        logger.info('Indexing reference FASTA...')
        results = scheduler.run(joblist=[{'label': 'refindex', 'command': [samtools, 'faidx', args.refseq]}])
        if not results:
            raise OSError('Failed to index reference sequences')
        rlabels = list(set([x.get('label') for x in results]))
        if not 'refindex' in rlabels:
            raise OSError('Failed to index reference sequences')
        if not os.path.exists(args.refseq + '.fai'):
            raise OSError('Failed to index refenence sequences')

    if not os.path.exists(os.path.splitext(args.refseq)[0] + '.dict'):
        logger.info('Generating sequence dictionary for reference FASTA...')
        output_filename = os.path.splitext(args.refseq)[0] + '.dict'
        picard_opts = ['CreateSequenceDictionary', 'R=%s' % args.refseq, 'O=%s' % output_filename]
        logger.debug('Create sequence dictionary with: java %s -jar %s %s', ' '.join(java_opts), jarfile, ' '.join(picard_opts))
        cmd_opts = command[:]
        cmd_opts.extend(picard_opts)
        results = scheduler.run(joblist=[{'label': 'refdict', 'command': cmd_opts}])

        if not results:
            raise OSError('Failed to generate reference sequence dictionary')
        rlabels = list(set([x.get('label') for x in results]))
        if not 'refdict' in rlabels:
            raise OSError('Failed to generate reference sequence dictionary')
        if not os.path.exists(os.path.splitext(args.refseq)[0] + '.dict'):
            raise OSError('Failed to generate reference sequence dictionary')

    jarfile = 'GenomeAnalysisTK.jar'
    if '--gatk_path' in extra:
        gatk_path = extra[extra.index('--gatk_path')+1]
        if not os.path.isdir(gatk_path):
            gatk_path = os.path.dirname(gatk_path)
        jarfile = os.path.join(gatk_path, jarfile)
    if not os.path.exists(jarfile):
        raise OSError('GenomeAnalysisToolkit jar not found at %s' % jarfile)
    command = [java]
    command.extend(java_opts)
    command.extend(['-jar', jarfile])
    joblist = []
    for label in filemeta:
        input_filename = filemeta[label]['dedup'][0]
        output_filename = label + '_snps.indels.vcf'
        if args.outdir:
            output_filename = os.path.join(args.outdir, output_filename)
        gatk_opts = ['-T', 'HaplotypeCaller', '-R', args.refseq, '-stand_call_conf', '30',
                     '-I', input_filename, '-o', output_filename]
        logger.debug('HaplotypeCaller: java %s -jar %s %s', ' '.join(java_opts), jarfile, ' '.join(gatk_opts))
        cmd_opts = command[:]
        cmd_opts.extend(gatk_opts)
        joblist.append({'label': label, 'command': cmd_opts})
    results = scheduler.run(joblist=joblist)
    if not results:
        raise OSError('Failed haplotype caller')
    rlabels = list(set([x.get('label') for x in results]))
    for label in filemeta:
        if not label in rlabels:
           raise OSError('Failed haplotype caller for %s' % label)
    for result in results:
        label = result.get('label')
        _result = result.get('result')
        if _result.returncode:
            raise OSError('Haplotype caller on %s failed with %s: %s' % (label, _result.returncode, _result.stderr or _result.stdout))
        if _result.stderr:
            raise OSError('Haplotype caller on %s returned with errors: %s' % (label, _result.stderr))
        output_filename = label + '_snps.indels.vcf'
        if args.outdir:
            output_filename = os.path.join(args.outdir, output_filename)
        if not os.path.exists(output_filename):
            raise OSError('Failed to find snps/indels file %s for %s' % (output_filename, label))

    logger.info('Complete')
    return 0


if __name__ == "__main__":
    sys.exit(main())
