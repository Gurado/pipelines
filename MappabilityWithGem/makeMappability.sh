#!/bin/bash
#
# SGE 
#$ -cwd
#$ -N GEM_indexer
#$ -l h_vmem=8G
#$ -b y
#$ -j y
#$ -V
#$ -pe smp 8

. ~/.profile
module load fabbus/gem gi/ucsc_utils
REFERENCE="/share/ClusterShare/biodata/contrib/genomeIndices_garvan/iGenomes/Mus_musculus/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa"
#gem-indexer -i $REFERENCE -o genome

TAGSIZE=75
if [ ! -s genome_$TAGSIZE.mappability ]; then
	gem-mappability -I genome.gem -o genome_$TAGSIZE -l $TAGSIZE -T 8
fi
if [ ! -s genome_$TAGSIZE.wig ]; then
	gem-2-wig -I genome.gem -i genome_$TAGSIZE.mappability -o genome_$TAGSIZE
	wigToBigWig genome_$TAGSIZE.wig genome_$TAGSIZE.sizes genome_$TAGSIZE.bw
fi

TAGSIZE=50
if [ ! -s genome_$TAGSIZE.mappability ]; then
        gem-mappability -I genome.gem -o genome_$TAGSIZE -l $TAGSIZE -T 8
fi
if [ ! -s genome_$TAGSIZE.wig ]; then
        gem-2-wig -I genome.gem -i genome_$TAGSIZE.mappability -o genome_$TAGSIZE
	wigToBigWig genome_$TAGSIZE.wig genome_$TAGSIZE.sizes genome_$TAGSIZE.bw
fi
