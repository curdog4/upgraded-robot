#!/bin/bash

usage()
{
    echo <<EOS
Usage: $0 <sequence ID> <FASTA DB file>
EOS
    return 0
}

if [ $# != 2 ]; then
    usage
    exit 1
fi

SEQID=$1
DBFILE=$2

if [ ! -e $DBFILE ]; then
    echo "Database $DBFILE does not appear to exist" >&2
    exit 1
fi

cat $DBFILE | tr '\n' '|' | tr '>' '\n' | grep $SEQID | tr '|' '\n' | sed -e "s/^\(.*$SEQID.*\)/>\1/"

exit 0
