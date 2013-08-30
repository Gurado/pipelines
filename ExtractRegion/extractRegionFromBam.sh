#!/bin/sh -e

#$ -l walltime=1:30:00
#$ -l vmem=8gb   
#$ -l nodes=1:ppn=1
#$ -j oe
#$ -N toy
#$ -o toy.out
#$ -cwd
#$ -V

DIR=$(dirname $0)
VERSION="0.0.1"

USAGEMSG="usage: $(basename $0) -o outputDirectory -f outfilename -p -v BAM LOCATION 

Generate toy example fastq files by taking a mapped file (BAM) and returning the FASTQ data for a specific region of interest
examples:
./extractRegionFromBam.sh -o ~/research/Sandbox_ngd1/fastq/ChIPseq -n ChIPseq_CTCF_chr16 ~/research/integration/TFs/H1esc/bowtie/wgEncodeBroadHistoneH1hescCtcfStdRawDataRep1.asd.bam chr16:27184646-27472388
./extractRegionFromBam.sh -o ~/research/Sandbox_ngd1/fastq/ChIPseq -n ChIPseq_H3k9me3_chr16 ~/research/integration/TFs/H1esc/bowtie/wgEncodeBroadHistoneH1hescH3k09me3StdRawDataRep1.asd.bam chr16:27184646-27472388
./extractRegionFromBam.sh -o ~/research/Sandbox_ngd1/fastq/ChIPseq_input -n ChIPseq_Input_chr16 ~/research/integration/TFs/H1esc_control/bowtie/wgEncodeBroadHistoneH1hescControlStdRawData.asd.bam chr16:27184646-27472388

Author: Fabian Buske
Version: $VERSION

* BAM - the bam file
* LOCATION - chr location in the form chrx:start-end
* p - indicated that libary is paired end
* o dir - where to put the data
* f file - change filename to this prefix (suffix will be _Rx.fq.gz, with x in {1,2})
"

[ $# -lt 2 ] && echo "$USAGEMSG" >&2 && exit 1
ISPAIRED="FALSE"
OUTDIR=
FILENAME=

while getopts "o:n:pv" opt;
do
	case ${opt} in
		o) OUTDIR="$OPTARG";;
		n) FILENAME="$OPTARG";;
		p) ISPAIRED="TRUE";;
		v) VERBOSE="--verbose";;
		\?) print >&2 "$0: error - unrecognized option $1"
		exit 1;;
        esac
done

shift $(($OPTIND-1))
BAM=$1
LOCATION=$2

[ ! -f $BAM ] && echo "[ERROR] Bam file does not exist: $BAM"

[ -z "$OUTDIR" ]  && OUTDIR=$(dirname $BAM)
mkdir -p $OUTDIR 

if [ -z "$FILENAME" ]; then
    n=${BAM##*/}
    n=${n/.bam/_subset}
else
    n=$FILENAME
fi
echo $n

module load kevyin/java/1.7.0_25 gi/samtools/0.1.19 gi/picard-tools/1.91

if [ "$ISPAIRED" = "TRUE" ]; then
	samtools view -b -f 2 $BAM $LOCATION > $OUTDIR/$n.bam
        samtools index $OUTDIR/$n.bam
	java -Xmx5g -jar $(which SamToFastq.jar) INPUT=$OUTDIR/$n.bam FASTQ=$OUTDIR/${n}_R1.fq SECOND_END_FASTQ=$OUTDIR/${n}_R2.fq
	gzip -9 $OUTDIR/${n}_R1.fq $OUTDIR/${n}_R2.fq

else
        samtools view -b $BAM $LOCATION > $OUTDIR/$n.bam
        samtools index $OUTDIR/$n.bam
        java -Xmx5g -jar $(which SamToFastq.jar) INPUT=$OUTDIR/$n.bam FASTQ=$OUTDIR/${n}_R1.fq
        gzip -9 $OUTDIR/${n}_R1.fq

fi
