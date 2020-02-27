#!/bin/bash
#SBATCH --job-name=_your_job_name_here
#SBATCH --account=_your_account_name_here_
#SBATCH --partition=_your_partition_name_here
#SBATCH --chdir=_your_job_working_directory_here_
#SBATCH --ntasks=1
#SBATCH --mem=65536M 
#SBATCH --cpus-per-task=16
#SBATCH --output=sbatch.out
#SBATCH --error=sbatch.err

OUTFILE="run_map.log"
REFSEQ="canu.contigs.fasta"
PGLBIN=/data/gpfs/assoc/pgl/bin
TRIMGALORE=$PGLBIN/TrimGalore/trim_galore

echo > $OUTFILE < /dev/null

RFILE=Ki11_S13_L002_R1_001.fastq.gz
echo "Trimming reads file $RFILE with TrimGalore"
$TRIMGALORE --paired --retain_unpaired --cores 16 --max_n 40 --gzip $RFILE ${RFILE/R1/R2}

echo "Index reference (bwa index)" >>$OUTFILE
bwa index -a bwtsw $REFSEQ >>$OUTFILE 2>&1
echo "Generate reference FAIDX (samtools faidx)" >>$OUTFILE
samtools faidx $REFSEQ >>$OUTFILE 2>&1

for RFILE in Ki11_S13_L002_R[12]_001_val_[12].fq.gz; do
    echo "Aligning reads from $RFILE to reference $REFSEQ" >>$OUTFILE
    ALIGN_OUT="$(basename -s.fq.gz $RFILE).sai"
    bwa aln -t 24 $REFSEQ $RFILE >$ALIGN_OUT 2>>$OUTFILE
done

echo "Merging results (bwa sampe)" >>$OUTFILE
R1=Ki11_S13_L002_R1_001_val_1.fq.gz
R2=${R1/_R1_001_val_1/_R2_001_val_2}
A1="$(basename -s.fq.gz $R1).sai"
A2="$(basename -s.fq.gz $R2).sai"
bwa sampe -f Ki13_S13-bwa_aln.sam $REFSEQ $A1 $A2 $R1 $R2 >>$OUTFILE 2>&1

touch run.done
exit 0

