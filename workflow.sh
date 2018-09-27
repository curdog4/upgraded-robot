#!/bin/bash

check_file_exists()
{
    CHECK_FILE=$1
    if [ ! -f $CHECK_FILE ]; then
        echo "ERROR: File $CHECK_FILE does not exist" >&2
        exit 1
    fi
    return 0
}

if [ -z $VIRTUAL_ENV ]; then
    . repo/upgraded-robot/.env/bin/activate
fi

QUERY_FILE=query/GPT1/Athaliana_Osativa_GPT1_protein.fasta
check_file_exists $QUERY_FILE

QUERY_PROT_DB=data/Athaliana/Athaliana_167_TAIR10.protein_primaryTranscriptOnly.fa
check_file_exists $QUERY_PROT_DB

PREFIX="Athaliana-Osativa"

#for DPATH in data/[A-Z]*; do
SAMPLE_SPECIES="Athaliana Atrichopoda Acomosus Osativa Stuberosum"
for DPATH in $SAMPLE_SPECIES; do
    DPATH=data/$DPATH
    SUBJECT=$(basename $DPATH)
    echo "BLAST'ing ${SUBJECT} for GPT1 match..."
    SUBJECT_CDS=( $DPATH/${SUBJECT}*cds_primaryTranscriptOnly.fa )
    if [ ! -e $SUBJECT_CDS ]; then
        echo "ERROR: No subject CDS found for ${SUBJECT}" >&2
        exit 1
    fi
    echo "Subject CDS DB: ${SUBJECT_CDS}"
    SUBJECT_PROT_DB=( $DPATH/${SUBJECT}*protein_primaryTranscriptOnly.fa )
    if [ ! -e $SUBJECT_PROT_DB ]; then
        echo "ERROR: No subject protein DB found for ${SUBJECT}" >&2
        exit 1
    fi
    echo "Subject protein DB: ${SUBJECT_PROT_DB}"
    echo "Running blastx..."
    ##
    # first, unfiltered results
    blastp -query $QUERY_FILE -db $SUBJECT_PROT_DB -evalue 1e-50 -num_threads 16 -outfmt 6 -out results/${PREFIX}_${SUBJECT}_alignment.txt 2>tmp/errors_blast_${SUBJECT}.txt

    ##
    # second, filter results to remove low coverage (<60%) and low identity (<30%)
    cat results/${PREFIX}_${SUBJECT}_alignment.txt | python repo/upgraded-robot/blast_filter.py --query-file $QUERY_FILE --query-cds $QUERY_FILE --subject-cds $SUBJECT_CDS --subject-prot $SUBJECT_PROT_DB --coverage 60 --identity 30 2>tmp/errors_blast_filter_${SUBJECT}.txt | tee results/${PREFIX}_${SUBJECT}_blastfiltered_alignment.txt

    ##
    # finally, filter results through elbow filter, 3 clusters using bitscore value, if we have >= 3 alignments
    A_COUNT=$(wc -l results/${PREFIX}_${SUBJECT}_blastfiltered_alignment.txt | awk '{print $1}')
    if [[ $A_CNT -ge 3 ]]; then
        cat results/${PREFIX}_${SUBJECT}_blastfiltered_alignment.txt | python repo/upgraded-robot/elbow_filter.py 2>tmp/errors_elbow_filter_${SUBJECT}.txt | tee results/${PREFIX}_${SUBJECT}_elbowfiltered_alignment.txt
    else
        cp results/${PREFIX}_${SUBJECT}_blastfiltered_alignment.txt results/${PREFIX}_${SUBJECT}_elbowfiltered_alignment.txt
    fi

    [ -f results/${SUBJECT}_GPT1_candidates.fasta ] && rm -f results/${SUBJECT}_GPT1_candidates_cds.fasta
    [ -f results/${SUBJECT}_GPT1_candidates.fasta ] && rm -f results/${SUBJECT}_GPT1_candidates_prot.fasta
    for SEQID in $(awk -F '\t' '{print $2}' results/${PREFIX}_${SUBJECT}_elbowfiltered_alignment.txt); do
        echo "Getting ${SUBJECT} CDS sequence ID ${SEQID}..."
        ./find_seqid.sh ${SEQID} ${SUBJECT_CDS} | tee -a results/${SUBJECT}_GPT1_candidates_cds.fasta
        echo "Getting ${SUBJECT} protein sequence ID ${SEQID}..."
        #python repo/upgraded-robot/cds_xlate.py --sequence-id $SEQID --description "${SUBJECT} GPT1" --cds-fasta-file ${SUBJECT_CDS} | tee -a results/${SUBJECT}_GPT1_candidates_prot.fasta
        ./find_seqid.sh ${SEQID} ${SUBJECT_PROT_DB} | tee -a results/${SUBJECT}_GPT1_candidates_prot.fasta
    done
    # some cleanup
    for F in tmp/errors_*_${SUBJECT}.txt; do
        if [ ! -s $F ]; then
            rm -f $F
        fi
    done
done

deactivate
