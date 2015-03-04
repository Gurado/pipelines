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

export MODULEPATH=/share/ClusterShare/Modules/modulefiles/contrib:/share/ClusterShare/Modules/modulefiles/centos6.2_x86_64:/share/ClusterShare/Modules/modulefiles/noarch:$MODULEPATH

module load fabbus/gem gi/ucsc_utils
REFERENCE="/share/ClusterShare/biodata/contrib/genomeIndices_garvan/iGenomes/Mus_musculus/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa"
#gem-indexer -i $REFERENCE -o genome

TAGSIZE=75
echo "make $TAGSIZE"
if [ ! -s genome_$TAGSIZE.mappability ]; then
	gem-mappability -I genome.gem -o genome_$TAGSIZE -l $TAGSIZE -T 8
fi
if [ ! -s genome_$TAGSIZE.wig ]; then
	gem-2-wig -I genome.gem -i genome_$TAGSIZE.mappability -o genome_$TAGSIZE
	wigToBigWig genome_$TAGSIZE.wig genome_$TAGSIZE.sizes genome_$TAGSIZE.bw
fi

TAGSIZE=50
echo "make $TAGSIZE"
if [ ! -s genome_$TAGSIZE.mappability ]; then
        gem-mappability -I genome.gem -o genome_$TAGSIZE -l $TAGSIZE -T 8
fi
if [ ! -s genome_$TAGSIZE.wig ]; then
        gem-2-wig -I genome.gem -i genome_$TAGSIZE.mappability -o genome_$TAGSIZE
	wigToBigWig genome_$TAGSIZE.wig genome_$TAGSIZE.sizes genome_$TAGSIZE.bw
fi

TAGSIZE=36
echo "make $TAGSIZE"
if [ ! -s genome_$TAGSIZE.mappability ]; then
        gem-mappability -I genome.gem -o genome_$TAGSIZE -l $TAGSIZE -T 8
fi
if [ ! -s genome_$TAGSIZE.wig ]; then
        gem-2-wig -I genome.gem -i genome_$TAGSIZE.mappability -o genome_$TAGSIZE
        wigToBigWig genome_$TAGSIZE.wig genome_$TAGSIZE.sizes genome_$TAGSIZE.bw
fi
