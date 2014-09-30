#!/bin/sh -e

# target folder where all operation/results should operate on/in
TARGETDIR=/share/ClusterScratch/fabbus/seqway_prec/

# folder containing the original data dolder 
FILES_SOURCE=/Cancer-Epigenetics-Disco/ClarkLab/

# Temp dir
export TMPDIR=$TARGETDIR/tmp
# data files 
# tab separated file containing the following columns
# Mark	Cellline	Folder	Replicates	FILESUFFIX Fragmentsizes	Wiggler_smooting
EXPERIMENTS="experiments.txt"

# experimental identifier
EXPERIMENT="Prostate"

# Uncomment to use all track data defined in experiment, otherwise only use tracks matching $EXPERIMENT in the second column
# USE_ALL_TRACK_DATA="1"

# genome info
GENOME="hg19"
#GENOMESEQ=/share/ClusterShare/biodata/galaxy_main/hg19/seq/${GENOME}.fa
CHROMSIZES=/share/ClusterShare/biodata/contrib/ENCODE/encodeDCC/referenceSequences/male.hg19.chrom.sizes
SEQDIR=/share/ClusterShare/biodata/contrib/ENCODE/encodeDCC/referenceSequences/maleByChrom/
FASTASUFFIX=".fa"

# annotation data
ANNOTATION=/share/ClusterShare/biodata/contrib/GENCODE/release_19/gencode.v19.annotation.reduced.gtf

## Segway parameters
LABELS=15
#REGIONSOURCE=http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/referenceSequences/encodePilotRegions.hg19.bed
REGIONS=/share/ClusterShare/software/contrib/Cancer-Epigenetics/Annotation/hg19/Encode/encodePilotRegions.hg19.bed

#http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDacMapabilityConsensusExcludable.bed.gz
EXCLUDABLE=/share/ClusterShare/software/contrib/Cancer-Epigenetics/Annotation/hg19/Encode/wgEncodeDacMapabilityConsensusExcludable.bed

# additional options for segway 
# e.g. "--cluster-opt='-pe smp'"
SEGWAY_OPTIONS="--cluster-opt='-pe smp'"

TRAIN_EXPERIMENT=PrEC
SEGWAY_TRAIN_ADDPARAM=""
INSTANCES=3

# wiggler
# maps for 20-54 bp reads
WIGGLER_UMAPDIR_lt100=/share/ClusterShare/biodata/contrib/fabbus/umap/hg19_male/globalmap_k20tok54/
# maps for 100 bp reads
WIGGLER_UMAPDIR_ge100=/share/ClusterShare/biodata/contrib/fabbus/umap/hg19_male/globalmap_k101tok101/


# which set in the PredictOn column defined in EXPERIMENTS to work on
PREDICTON="1"
