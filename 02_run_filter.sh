#!/bin/bash
#SBATCH --job-name=_your_job_name_here_
#SBATCH --account=_your_account_name_here_
#SBATCH --partition=_your_partition_name_here_
#SBATCH --chdir=_your_working_directory_here_
#SBATCH --ntasks=1
#SBATCH --mem=65536M 
#SBATCH --cpus-per-task=16
#SBATCH --output=sbatch.out
#SBATCH --error=sbatch.err

OUTFILE="run_filter.log"
REFSEQ="canu.contigs.fasta"

SCRATCH=/data/gpfs/home/dcurdie/scratch
PERL=$SCRATCH/miniconda-3.7/bin/perl
ALLHIC=$SCRATCH/ALLHiC

export PATH=$PATH:$ALLHIC/bin:$ALLHIC/scripts

echo > $OUTFILE < /dev/null

ALN=Ki13_S13-bwa_aln.sam
BASE=$(basename -s.sam $ALN)
echo "Preprocess merged alignments (PreprocessSAMs.pl)" >>$OUTFILE
$PERL $ALLHIC/scripts/PreprocessSAMs.pl $ALN $REFSEQ MBOI >>$OUTFILE 2>&1
echo "Filter merged alignments (filterBAM_forHiC.pl)" >>$OUTFILE
$PERL $ALLHIC/scripts/filterBAM_forHiC.pl ${BASE}.REduced.paired_only.bam ${BASE}.clean.sam >>$OUTFILE 2>&1
echo "Convert to BAM (samtools view)" >>$OUTFILE
samtools view -o ${BASE}.clean.bam -bt ${REFSEQ}.fai ${BASE}.clean.sam >>$OUTFILE 2>&1

touch run.done
exit 0

