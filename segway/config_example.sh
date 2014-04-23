#!/bin/sh -e

# target folder where all operation/results should operate on/in
TARGETDIR=/home/fabbus/research/tmp/

# folder containing the original data dolder 
FILES_SOURCE=/Cancer-Epigenetics/Data/ClarkLab/Seq/ChIP-Seq/hg19/

# data files 
# tab separated file containing the following columns
# Mark	Cellline	Folder	Replicates	Fragmentsizes	Wiggler_smooting    TrainOn PredictOn
EXPERIMENTS="experiments.txt"

# experimental identifier
EXPERIMENT="Prostate"

# genome info
GENOME="hg19"
#GENOMESEQ=/share/ClusterShare/biodata/galaxy_main/hg19/seq/${GENOME}.fa
CHROMSIZES=/share/ClusterShare/biodata/contrib/ENCODE/encodeDCC/referenceSequences/male.hg19.chrom.sizes
SEQDIR=/share/ClusterShare/biodata/contrib/ENCODE/encodeDCC/referenceSequences/maleByChrom/
FASTASUFFIX=".fa"

# annotation data
ANNOTATION=/share/ClusterShare/biodata/contrib/GENCODE/release_19/gencode.v14.annotation.gtf

## Segway parameters
LABELS=20
#REGIONSOURCE=http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/referenceSequences/encodePilotRegions.hg19.bed
REGIONS=/share/ClusterShare/software/contrib/Cancer-Epigenetics/Annotation/hg19/Encode/encodePilotRegions.hg19.bed
TRAIN_EXPERIMENT=PrEC
SEGWAY_TRAIN_ADDPARAM=""
INSTANCES=3

## Wiggler
# Wiggler smoothing for histone marks
WIGGLER_UMAPDIR=/share/ClusterShare/biodata/contrib/fabbus/umap/hg19_male/globalmap_k20tok54/
WIGGLER_SMOOTHING=300
