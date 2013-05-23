#!/bin/sh -e

USAGEMSG="usage: $(basename $0) IPstrength.output IPstrength.png ENCODE.output

writes a CHANCE ChIP-QC report in HTML

Author: Fabian Buske

* CHIP - the chip data (bam expected)
* CONTROL - the input/control data (bam expected)
* -v - print progress information (verbose).
"

DIR=$(dirname $0)
VERSION="0.0.1"

[ $# -lt 2 ] && echo "$USAGEMSG" >&2 && exit 1

CHIP=""
CONTROL=""
BUILD="hg19"
VERBOSE="--quiet"

while getopts "v" opt;
do
        case ${opt} in
                v) VERBOSE="--verbose";;
                \?) print >&2 "$0: error - unrecognized option $1"
                exit 1;;
        esac
done

shift $(($OPTIND-1))
IPSTRENGTH=$1
IPIMAGE=$2
ENCODE=$3

if [ ! -f $IPSTRENGTH ]; then
	echo "[ERROR] ChIP IPstrength output does not exist: $IPSTRENGTH"
	exit 1
fi
OUT=${IPSTRENGTH%.*}".html"

if [ ! -f $IPIMAGE ]; then
        echo "[WARN] ChIP IPstrength image does not exist: $IPIMAGE"
fi
IPIMAGE=${IPIMAGE##*/}

if [ ! -f $ENCODE ]; then
        echo "[WARN] ChIP Encode output does not exist: $ENCODE"
fi

CHIP=$(grep "^IP_sample_id," $IPSTRENGTH | awk -F\, '{print $2}')
INPUT=$(grep "^Input_sample_id," $IPSTRENGTH | awk -F\, '{print $2}')
OUTCOME=$(grep "^pass," $IPSTRENGTH | awk -F\, '{print $2}')
if [ "$OUTCOME" = 0 ]; then
	OUTCOME="No"
else
	OUTCOME="Yes"
fi
SCALING=$(grep "^input_scaling_factor," $IPSTRENGTH | awk -F\, '{print $2}')
ENRICHMENT=$(grep "^percent_genome_enriched," $IPSTRENGTH | awk -F\, '{print $2}')
DIFFENRICH=$(grep "^differential_percentage_enrichment," $IPSTRENGTH | awk -F\, '{print $2}')
FDR=$(grep "^fdr," $IPSTRENGTH | awk -F\, '{print $2}')
FDR_TFBS_NORMAL=$(grep "^tfbs_normal_fdr," $IPSTRENGTH | awk -F\, '{print $2}')
FDR_HIST_NORMAL=$(grep "^histone_normal_fdr," $IPSTRENGTH | awk -F\, '{print $2}')
FDR_TFBS_CANCER=$(grep "^tfbs_cancer_fdr," $IPSTRENGTH | awk -F\, '{print $2}')
FDR_HIST_CANCER=$(grep "^histone_cancer_fdr," $IPSTRENGTH | awk -F\, '{print $2}')


echo "<!DOCTYPE html PUBLIC '-//W3C//DTD XHTML 1.0 Transitional//EN' 'http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd'> \
<html xmlns='http://www.w3.org/1999/xhtml'> \
<html> \
<head><title>Chance output for $CHIP </title></head> \
<body> \
ChIP: $CHIP<br/> \
Input: $INPUT<br/> \
<h3>IPstrength</h3>\
<table border='1' cellspacing='0' cellpadding='2'>\
<tr>\
	<th>Significance</th>\
	<th>Input scaling</th>\
	<th>Genome enriched (%)</th>\
	<th>Diff. enrichment (%)</th>\
	<th>FDR (overall)</th>\
	<th>FDR TFBS (normal)</th>\
	<th>FDR Histone (normal)</th>\
        <th>FDR TFBS (cancer)</th>\
	<th>FDR Histone (normal)</th>\
</tr>\
<tr>\
	<td>$OUTCOME</td>\
	<td>$SCALING</td>\
	<td>$ENRICHMENT</td>\
	<td>$DIFFENRICH</td>\
	<td>$FDR</td>\
	<td>$FDR_TFBS_NORMAL</td>\
	<td>$FDR_HIST_NORMAL</td>\
	<td>$FDR_TFBS_CANCER</td>\
	<td>$FDR_HIST_CANCER</td>\
</tr>\
</table>\
<img src='$IPIMAGE'>" > $OUT

if [ -f $ENCODE ]; then
BUILD=$(grep "build," $ENCODE | awk -F\, '{print $2}')
SNR=$(grep "odds_ratio," $ENCODE | awk -F\, '{print $2}')
PROB=$(grep "probability," $ENCODE | awk -F\, '{print $2}')

echo "
<h3>Encode comparison</h3>\
<table border='1' cellspacing='0' cellpadding='2'>\
<tr>\
	<th> Build </th>\
	<th> Signal to noise ratio </th>\
	<th> Probability* </th>\
</tr>\
<tr>\
	<td>$BUILD</td>\
	<td>$SNR</td>\
	<td>$PROB</td>\
</tr>\
</table>\
*A small probability indicates your data differs greatly from ENCODE datasets" >> $OUT
fi

echo "</body>\
</html>" >> $OUT
