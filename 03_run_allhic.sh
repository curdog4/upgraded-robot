#!/bin/bash
#SBATCH --job-name=_your_job_name_here_
#SBATCH --account=_your_account_name_here_
#SBATCH --partition=_your_partition_name_here_
#SBATCH --chdir=_your_working_directory_name_here_
#SBATCH --ntasks=1
#SBATCH --mem=65536M 
#SBATCH --cpus-per-task=16
#SBATCH --output=sbatch.out
#SBATCH --error=sbatch.err

OUTFILE="run_allhic.log"
REFSEQ="canu.contigs.fasta"

SCRATCH=/data/gpfs/home/dcurdie/scratch
PERL=$SCRATCH/miniconda-3.7/bin/perl
ALLHIC=$SCRATCH/ALLHiC

export PATH=$PATH:$ALLHIC/bin:$ALLHIC/scripts

echo > $OUTFILE < /dev/null

##
# Prune step
BAM=Ki13_S13-bwa_aln.clean.bam
ALLTBL=Allele.ctg.table
echo "Pruning $BAM..." >&2 |& tee -a $OUTFILE
$ALLHIC/bin/ALLHiC_prune -i $ALLTBL -b $BAM -r $REFSEQ 2>&1 | tee -a $OUTFILE
if [ $? != 0 ]; then
    echo "ALLHiC_prune exited non-zero[$?]" >&2 |& tee -a $OUTFILE
    exit 1
fi
if [ ! -s prunning.bam ]; then
    echo "ALLHiC_prune output file prunning.bam does not exist or is zero-length" >&2 |& tee -a $OUTFILE
    exit 1
fi

##
# Partition step
RE_SITES="GATCGATC,GAATGATC,GATTGATC,GACTGATC,GAGTGATC,GAATAATC,GAATATTC,GAATACTC,GAATAGTC,GTATAATC,GTATATTC,GTATACTC,GTATAGTC,GCATAATC,GCATATTC,GCATACTC,GCATAGTC,GGATAATC,GGATATTC,GGATACTC,GGATAGTC,GATCAATC,GATCATTC,GATCACTC,GATCAGTC"
NOCTGS=10
echo "Partitioning pruned BAM " 2>&1 |& tee -a $OUTFILE
$ALLHIC/bin/ALLHiC_partition -b prunning.bam -r $REFSEQ -e $RE_SITES -k $NOCTGS 2>&1 |& tee -a $OUTFILE
if [ $? != 0 ]; then
    echo "ALLHiC_partition exited non-zero[$?]" >&2 | tee -a $OUTFILE
    exit 1
fi
for FILE_NAME in prunning.{clm,distribution.txt,pairs.txt,counts_${RE_SITES//,/_}.txt} ; do
    if [ ! -s $FILE_NAME ]; then
        echo "ALLHiC_partition output file $FILE_NAME does not exist or is zero-length" >&2 |& tee -a $OUTFILE
        exit 1
    fi
done

##
# Rescue step (Skipping)
#echo "Rescuing $BAM" 2>&1 |& tee -a $OUTFILE
#CLUSTERS=
#COUNTS=
#$ALLHIC/bin/ALLHiC_rescue -i $BAM -r $REFSEQ -c $CLUSTERS -i $COUNTS 2>&1 |& tee -a $OUTFILE

##
# Optimize step
echo "Extracting BAM using restriction enzyme sites $RE_SITES" 2>&1 |& tee -a $OUTFILE
$ALLHIC/bin/allhic extract $BAM $REFSEQ --RE $RE_SITES
if [ $? != 0 ]; then
    echo "allhic extract exited non-zero[$?]" >&2 |& tee -a $OUTFILE
    exit 1
fi

echo "End of known outputs. Terminating script." >&2 |& tee -a $OUTFILE
exit 0

if [ ! -s SOME_FILE_NAME_HERE ]; then
    echo "allhic extract ouput file SOME_FILE_NAME_HERE does not exist or is zero-length" >&2 |& tee -a $OUTFILE
    exit 1
fi
IDX=0
CLM=${BAM/bam/clm}
while [ $IDX -le $NOCTGS ]; do
    IDX=$(($IDX+1))
    echo "Running optimize pass $IDX..." 2>&1 |& tee -a $OUTFILE
    $ALLHIC/bin/allhic optimize group${IDX}.txt $CLM 2>&1 |& tee -a $OUTFILE
    if [ $? != 0 ]; then
        echo "allhic optimize exited non-zero[$?]" >&2 |& tee -a $OUTFILE
        exit 1
    fi
    if [ ! -s SOME_FILE_NAME_HERE ]; then
        echo "allhic optimize output file SOME_FILE_NAME_HERE does not exist or is zero-length" >&2 |& tee -a $OUTFILE
        exit 1
    fi
done

##
# Build step
echo "Building group locations" >&2 |& tee -a $OUTFILE
$ALLHIC/bin/ALLHiC_build $REFSEQ 2>&1 |& tee -a $OUTFILE
if [ $? != 0 ]; then
    echo "ALLHiC_build exited non-zero[$?]" >&2 |& tee -a $OUTFILE
    exit 1
fi
if [ ! -s SOME_FILE_NAME_HERE ]; then
    echo "ALLHiC_build output file SOME_FILE_NAME_HERE does not exist or is zero-length" >&2 |& tee -a $OUTFILE
    exit 1
fi

##
# Plot step
GROUPS=
CHRN=
SIZE=500k
FORMAT=pdf
echo "Generating the chromatin contact matrix to evaluate genome scaffolding" >&2 |& tee -a $OUTFILE
$ALLHIC/bin/ALLHiC_plot $BAM $GROUPS $CHRN $SIZE $FORMAT
if [ $? != 0 ]; then
    echo "ALLHiC_plot exited non-zero[$?]" >&2 |& tee -a $OUTFILE
    exit 1
fi
if [ ! -s SOME_FILE_NAME_HERE ]; then
    echo "ALLHiC_plot output file SOME_FILE_NAME_HERE does not exist or is zero-length" >&2 |& tee -a $OUTFILE
    exit 1
fi

echo "Run complete" >&2 |& tee -a $OUTFILE

touch run.done
exit 0

