#!/bin/sh -e

# indicates what to do (set variable to 1)
DO_TRANSFERDATA=""
DO_CONVERTBAM2BEDGRAPH=""
DO_GENERATEARCHIVE=""
DO_TRAINSEGWAY=""
DO_PREDICTSEGWAY=""
DO_EVALUATE="1"

# General runtime parameter
OVERWRITEALL="1" # set to 1 to overwrite existing results

# target folder where all operation/results should operate on/in
TARGETDIR=/home/fabbus/research/HMEC/

# folder containing the original data dolder 
FILES_SOURCE=/Cancer-Epigenetics/Data/ClarkLab/Seq/ChIP-Seq/hg19/

# data files (tracks) 
FILES="USC20120622_HMEC_H3K27ac USC20120622_HMEC_H3K27me3 USC20120622_HMEC_H3K4me1 USC20120622_HMEC_H3K4me3 USC20120622_HMEC_RNApolIIp USC20120622_HMEC_Input"

# experimental identifier
EXPERIMENT="HMEC"

# data folder (to put or source the data from)
# Important: add final /
SEGWAY_DATA=${TARGETDIR}data/
SEGWAY_BIN=${TARGETDIR}bin/
SEGWAY_QOUT=${TARGETDIR}qout/
SEGWAY_RESULT=${TARGETDIR}result/
SEGWAY_TRAIN=${SEGWAY_RESULT}train-4-M-${EXPERIMENT}/
SEGWAY_PREDICT=${SEGWAY_RESULT}predict-4-M-${EXPERIMENT}/

# genome info
GENOME="hg19"
GENOMESEQ=/share/ClusterShare/biodata/galaxy_main/hg19/seq/${GENOME}.fa

# annotation data
ANNOTATION=gencode_release_15/gencode.v15.annotation.level_1_2.gtf
ANNOTATION_SOURCE=ftp://ftp.sanger.ac.uk/pub/gencode/release_15/gencode.v15.annotation.gtf.gz

## Segway parameters
LABELS=20
#LABELS=5
REGIONS="encodePilotRegions.hg19.bed"
SPECIAL="--clobber"
INSTANCES=3
TRAIN_REGIONS="${SEGWAY_DATA}encodePilotRegions.hg19.bed"

# some housekeeping
mkdir -p $SEGWAY_DATA $SEGWAY_BIN $SEGWAY_QOUT $SEGWAY_RESULT

# get encode pilot regions
if [ ! -f ${SEGWAY_DATA}/encodePilotRegions.hg19.bed ]; then
	wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/referenceSequences/encodePilotRegions.hg19.bed  
	mv encodePilotRegions.hg19.bed ${SEGWAY_DATA}
fi

if [ ! -f ${SEGWAY_DATA}${ANNOTATION} ]; then
	wget ${ANNOTATION_SOURCE}
	b=$(basename $f)
	mv ${b} ${SEGWAY_DATA}
	tar -xf ${SEGWAY_DATA}${b} 
fi

