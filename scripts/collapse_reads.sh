#!/usr/bin/env bash

##collapse same reads after removing adaptors - PCR duplicates removal
usage() {
    echo "bash $0 <INPUT Fastq> <Output Fastq> <Threads> <length to trim in 5end:15> <minimum length for output:25>"
    exit;
}

INPUT=$1
OUTPUT=$2
THREADS=$3
trim_length=$4
min_length=$5
##
if [ $# -ne 5 ]; then
    usage
fi


zcat -f $INPUT  |  paste - - - - |\
    sort -t$'\t' -k2,2 -S 8G -u |\
    mawk 'BEGIN{FS="\t";OFS="\t"}{print $1"\n"$2"\n"$3"\n"$4;}' |\
    paste - - - - |\
    mawk -v trimlen=$trim_length -v minlen=$min_length '
        BEGIN{FS="\t";OFS="\n"}
        {
            if(trimlen>0){
                newseq=substr($2,trimlen+1);
                newqseq=substr($4,trimlen+1);
                if(length(newseq)>=minlen)
                    {print $1,newseq,$3,newqseq;
            } else {
                print $0
            }
        }}' |\
    pigz -c -p ${THREADS} > $OUTPUT
