#!/bin/bash -e

DIR=$(dirname $0)
VERSION="0.0.1"

USAGEMSG="usage: $(basename $0) -s FASTQSUFFIX -v FASTQ1 FASTQ2 RECORDS OUTDIR

Will overwrite the 
Subsample paired-end fastq files:

HiC assuming merged bam and properly sets flag
~/subsample.sh GMall_Ncol_R1.fastq.gz GMall_Ncol_R2.fastq.gz 10000 ~/tmp/

Author: Fabian Buske
inspired from https://www.biostars.org/p/76791/
Version: $VERSION

* FASTQ1 - read1 fastq file
* FASTQ2 - read2 fastq file
* RECORDS - number of records to subsample
* FASTQSUFFIX - file suffix e.g. fastq or fastq.gz
* OUTDIR  - where to put the data
"

[ $# -lt 3 ] && echo "$USAGEMSG" >&2 && exit 1
OUTDIR=
FASTQSUFFIX=

while getopts "o:s:i:v" opt;
do
	case ${opt} in
	    s) FASTQSUFFIX="$OPTARG";;
        v) VERBOSE="--verbose";;
        \?) print >&2 "$0: error - unrecognized option $1"
        exit 1;;
    esac
done

shift $(($OPTIND-1))
FASTQ1=$1
FASTQ2=$2
RECORDS=$3
OUTDIR=$4

mkdir -p $OUTDIR

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

SAMPLE1=${FASTQ1##*/}
SAMPLE2=${FASTQ2##*/}
echo $OUTDIR/${SAMPLE1/%$FASTQSUFFIX/.tmp}


paste <($CAT $FASTQ1) <($CAT $FASTQ2) |\
awk '{ printf("%s",$0); n++; if(n%4==0) { printf("\n");} else { printf("\t\t");} }' |\
shuf -n $RECORDS |\
sed 's/\t\t/\n/g' |\
awk -F'\t' -v FN1=$OUTDIR/${SAMPLE1/%$FASTQSUFFIX/.tmp} -v FN2=$OUTDIR/${SAMPLE2/%$FASTQSUFFIX/.tmp} '{print $1 > FN1; print $2 > FN2}'
gzip -9 -c $OUTDIR/${SAMPLE1/%$FASTQSUFFIX/.tmp} > $OUTDIR/${SAMPLE1/%$FASTQSUFFIX/fastq.gz}
gzip -9 -c $OUTDIR/${SAMPLE2/%$FASTQSUFFIX/.tmp} > $OUTDIR/${SAMPLE2/%$FASTQSUFFIX/fastq.gz}

rm $OUTDIR/${SAMPLE1/%$FASTQSUFFIX/.tmp} $OUTDIR/${SAMPLE2/%$FASTQSUFFIX/.tmp}
