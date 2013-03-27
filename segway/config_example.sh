#!/bin/sh -e

# target folder where all operation/results should operate on/in
TARGETDIR=/home/fabbus/research/tmp/

# folder containing the original data dolder 
FILES_SOURCE=/Cancer-Epigenetics/Data/ClarkLab/Seq/ChIP-Seq/hg19/

# data files (tracks) - separate replicas by comma
FILES="USC20120622_HMEC_H3K27ac.bam USC20120622_HMEC_H3K27me3.bam USC20120622_HMEC_H3K4me1.bam USC20120622_HMEC_H3K4me3.bam USC20120622_HMEC_RNApolIIp.bam USC20120622_HMEC_Input.bam USC20130125_LNCaP_Input_1.bam,USC20130125_LNCaP_Input_2.bam"

# experimental identifier
EXPERIMENT="HMEC"

# genome info
GENOME="hg19"
#GENOMESEQ=/share/ClusterShare/biodata/galaxy_main/hg19/seq/${GENOME}.fa
CHROMSIZES=/share/ClusterShare/biodata/contrib/fabbus/encodeDCC/male.hg19.chrom.sizes
SEQDIR=/share/ClusterShare/biodata/contrib/fabbus/encodeDCC/maleByChrom/

# annotation data
ANNOTATION=/share/ClusterShare/software/contrib/Cancer-Epigenetics/Annotation/hg19/transcript/gencode.v14/gencode.v14.annotation.gtf
#ANNOTATION_SOURCE=ftp://ftp.sanger.ac.uk/pub/gencode/release_15/gencode.v15.annotation.gtf.gz

## Segway parameters
LABELS=20
#REGIONSOURCE=http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/referenceSequences/encodePilotRegions.hg19.bed
REGIONS=/share/ClusterShare/software/contrib/Cancer-Epigenetics/Annotation/hg19/Encode/encodePilotRegions.hg19.bed
SPECIAL=""
INSTANCES=3

## Wiggler
# Wiggler smoothing for histone marks
WIGGLER_UMAPDIR=/share/ClusterShare/biodata/contrib/fabbus/umap/hg19_male/globalmap_k20tok54/
WIGGLER_SMOOTHING=300
