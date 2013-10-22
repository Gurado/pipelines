#!/bin/sh -e

# Script to obtain Repeat masked tracks from UCSC.
# author: Fabian Buske
# date: October 2013

function usage {
echo -e "usage: $(basename $0) [OPTIONS] GENOME SHORESIZE
* GENOME : genome assembly to use, default hg19 
-o OUTDIR : output location: default <GENOME>/
-g MAKEGTF : flag indicates that a gtf files should be created from the tracks
"

exit
}

if [ ! $# -ge 1 ]; then usage ; fi

GENOME=hg19
OUTDIR=
MAKEGTF=

#INPUTS                                                                                                           
while getopts "go:" opt;
do
	case ${opt} in
		o) OUTDIR="$OPTARG";;
		g) MAKEGTF="TRUE";;
		\?) print >&2 "$0: error - unrecognized option $1"
		exit 1;;
        esac
done

shift $(($OPTIND-1))
GENOME=$1

module load gi/ucsc_utils/283 gi/bedtools/2.17.0

if [ -z "$OUTDIR" ]; then
    OUTDIR=$GENOME
fi

mkdir -p $OUTDIR

mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select genoName, genoStart, genoEnd, repClass from $GENOME.rmsk" | tail -n+2 | awk '{OFS="\t";print $1,$2,$3,"Repeat_"NR,$4}' > ${OUTDIR}/Repeats.bed


for REPEAT in $(awk '{print $5}' ${OUTDIR}/Repeats.bed | grep -v '\?' | sort -u | tr '\n' ' '); do

	grep "$REPEAT" ${OUTDIR}/Repeats.bed | awk '{OFS="\t"; print $1,$2,$3,$4}' > ${OUTDIR}/$REPEAT.bed
	if [ "$MAKEGTF" = "TRUE" ]; then
    	bedToGenePred ${GENOME}/$REPEAT.bed stdout | genePredToGtf file stdin ${OUTDIR}/$REPEAT.gtf
	fi
done

mv $OUTDIR/Unknown.bed $OUTDIR/Unknown_repeat.bed

rm ${OUTDIR}/Repeats.bed
