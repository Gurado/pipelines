#!/bin/sh -e

# Script to obtain CpG islands track from UCSC.
# author: Fabian Buske
# date: October 2013

function usage {
echo -e "usage: $(basename $0) -o OUTDIR GENOME SHORESIZE
* GENOME : genome assembly to use, default hg19 
* SHORESIZE : size of the shores flanking CpG islands, default 2000bp 
-o OUTDIR : output location: default <GENOME>/
-g MAKEGTF : flag indicates that a gtf file should be created from the track
"

exit
}

if [ ! $# -ge 2 ]; then usage ; fi

GENOME=hg19
SHORESIZE=2000
OUTDIR=
MAKEGTF=

#INPUTS                                                                                                           
while getopts "o:" opt;
do
	case ${opt} in
		o) OUTDIR="$OPTARG";;
		\?) print >&2 "$0: error - unrecognized option $1"
		exit 1;;
        esac
done

shift $(($OPTIND-1))
GENOME=$1
SHORESIZE=$2

module load gi/ucsc_utils/283 gi/bedtools/2.17.0

if [ -z "$OUTDIR" ]; then
    OUTDIR=$GENOME
fi

mkdir -p $OUTDIR

mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, chromStart, chromEnd from $GENOME.cpgIslandExt" | tail -n+2| awk '{OFS="\t"; print $1,$2,$3,"CpGisland_"NR}' > ${OUTDIR}/CpGislands.bed

mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from $GENOME.chromInfo" > ${OUTDIR}/genome
	
bedtools slop -b $SHORESIZE -g ${OUTDIR}/genome -i ${OUTDIR}/CpGislands.bed > ${OUTDIR}/CpGislands$SHORESIZE.bed

bedtools subtract -a ${OUTDIR}/CpGislands$SHORESIZE.bed -b ${OUTDIR}/CpGislands.bed | bedtools sort | bedtools merge | awk '{OFS="\t";print $1,$2,$3,"CpGshore_"NR}'> ${OUTDIR}/CpGshores.bed

rm ${OUTDIR}/CpGislands$SHORESIZE.bed ${OUTDIR}/genome

if [ "$MAKEGTF" = "TRUE" ]; then
    bedToGenePred ${OUTDIR}/CpGislands.bed stdout | genePredToGtf file stdin ${OUTDIR}/CpGislands.gtf
    bedToGenePred ${OUTDIR}/CpGshores.bed stdout | genePredToGtf file stdin ${OUTDIR}/CpGshores.gtf
fi
