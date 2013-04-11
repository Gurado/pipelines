#!/bin/sh -e

# target folder where all operation/results should operate on/in
TARGETDIR=/home/fabbus/research/HiC/

# folder containing the original data dolder 
RAW_FILES_SOURCE=/Cancer-Epigenetics-Data/HiC_RAW/
MAPPED_FILES_SOURCE=/Cancer-Epigenetics/Data/ClarkLab/HiC/Seq/HiC

# data files (tracks) - separate replicas by comma
FASTQ="TKCC20130321_HiC_LNCaP_1 TKCC20130321_HiC_LNCaP_2"
BAMFILES=""

# experimental identifier
EXPERIMENT="LnCAP"

# genome info
GENOME="hg19"
#GENOMESEQ=/share/ClusterShare/biodata/galaxy_main/hg19/seq/${GENOME}.fa
CHROMSIZES=/share/ClusterShare/biodata/contrib/fabbus/encodeDCC/male.hg19.chrom.sizes
SEQDIR=/share/ClusterShare/biodata/contrib/fabbus/encodeDCC/maleByChrom/


