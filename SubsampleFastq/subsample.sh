#!/bin/bash -r

DIR=$(dirname $0)
VERSION="0.0.1"

USAGEMSG="usage: $(basename $0) -o outputDirectory -s FASTQSUFFIX -v FASTQ1 FASTQ2 RECORDS 

Subsample paired-end fastq files:

HiC assuming merged bam and properly sets flag
~/subsample.sh -o ~/tmp/ GMall_Ncol_R1.fastq.gz GMall_Ncol_R2.fastq.gz 10000

Author: Fabian Buske
inspired from https://www.biostars.org/p/76791/
Version: $VERSION

* FASTQ1 - read1 fastq file
* FASTQ2 - read2 fastq file
* RECORDS - number of records to subsample
* FASTQSUFFIX - file suffix e.g. fastq or fastq.gz
* o dir - where to put the data
"

[ $# -lt 2 ] && echo "$USAGEMSG" >&2 && exit 1
OUTDIR=
FASTQSUFFIX=

while getopts "o:s:v" opt;
do
	case ${opt} in
        o) OUTDIR="$OPTARG";;
	s) FASTQSUFFIX="$OPTARG";;
        v) VERBOSE="--verbose";;
        \?) print >&2 "$0: error - unrecognized option $1"
        exit 1;;
    esac
done

shift $(($OPTIND-2))
FASTQ1=$1
FASTQ2=$2
RECORDS=$3

#is ziped ?
CAT="cat"
if [[ ${FASTQ1##*.} == "gz" ]]; then 
     CAT="zcat";
elif [[ ${FASTQ1##*.} == "bz2" ]]; 
    then CAT="bzcat"; 
fi

if [ -z "$OUTDIR" ]; then
    OUTDIR=$(dirname $FASTQ1)
fi

paste <($CAT $FASTQ1) <($CAT $FASTQ2) |\ #merge the two fastqs
awk '{ printf("%s",$0); n++; if(n%4==0) { printf("\n");} else { printf("\t\t");} }' |\ #merge by group of 4 lines
shuf  |\ #shuffle
head -n $RECORDS |\
sed 's/\t\t/\n/g' |\ #restore the delimiters
awk -F"\t" '{print $1 > "${FASTQ1/%.gz/.tmp}"; print $2 > "${FASTQ2/%.gz/.tmp"}'
gzip -9 -c ${FASTQ1/%$FASTQSUFFIX/.tmp} > $OUTDIR/${FASTQ1/%$FASTQSUFFIX/_subsampled.fastq.gz}
gzip -9 -c ${FASTQ2/%$FASTQSUFFIX/.tmp} > $OUTDIR/${FASTQ2/%$FASTQSUFFIX/_subsampled.fastq.gz}

rm ${FASTQ1/%$FASTQSUFFIX/.tmp} ${FASTQ2/%$FASTQSUFFIX/.tmp}
