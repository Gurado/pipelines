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
./extractRegionFromBam.sh -o ~/research/Sandbox_ngd1/fastq/ChIPseq_input -n ChIPseq_Input_chr16 ~/research/integration/TFs/H1esc_control/bowtie/wgEncodeBroadHistoneH1hescControlStdRawData.asd.bam chr16:20000000-50000000

HiC assuming merged bam and properly sets flag
~/extractRegionFromBam.sh -1 R1 -2 R2 -p -o ~/tmp/ -n GMall_Ncol GMall_Ncol_uniques.bam chr16:20000000-30000000

Author: Fabian Buske
Version: $VERSION

* BAM - the bam file
* LOCATION - chr location in the form chrx:start-end
* 1 - read one identifier
* 2 - read 2 identifier
* p - indicates that libary is paired end
* i - include read pairs where the mate mapped outside the region of interest (otherwise the readpair is removed
* o dir - where to put the data
* n file - change filename to this prefix (suffix will be _Rx.fastq.gz, with x in {1,2})
"

[ $# -lt 2 ] && echo "$USAGEMSG" >&2 && exit 1
ISPAIRED="FALSE"
OUTDIR=
FILENAME=
LIBRARY="SINGLE"
OUTSITEMAPPINGMATES="EXCLUDE"

while getopts "1:2:o:n:ipv" opt;
do
	case ${opt} in
        1) READ1="$OPTARG";;
        2) READ2="$OPTARG";;
        o) OUTDIR="$OPTARG";;
        n) FILENAME="$OPTARG";;
        p) LIBRARY="PAIRED";;
        i) OUTSITEMAPPINGMATES="INCLUDE";;
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

module load gi/java/jdk1.7.0_45 gi/samtools/0.1.19 gi/picard-tools/1.103 gi/pigz/2.3

if [ "$LIBRARY" = "PAIRED" ]; then
    # sort by coordinates
    samtools sort $BAM $OUTDIR/$n.coord_sorted
    samtools index $OUTDIR/$n.coord_sorted.bam

    if [ "$OUTSITEMAPPINGMATES" == "INCLUDE" ]; then
        # get read ids for region of interest
        samtools view -f 3 $OUTDIR/$n.coord_sorted.bam $LOCATION | cut -f 1 | sort -u | awk -F'.' '{OFS="\t";print $1,$2}' | sort -k1,1 -k2,2g | awk '{OFS=".";print $1,$2}' | less > $OUTDIR/fastqIDs.txt
        
        # resort with read name
        java -Xmx30g -jar $(which SortSam.jar) INPUT=$BAM.bam OUTPUT=$BAM.name_sorted.bam SORT_ORDER=queryname
    
        # filter reads    
        java -Xmx30g -jar $(which FilterSamReads.jar) INPUT=$BAM.name_sorted.bam FILTER=includeReadList READ_LIST_FILE=$OUTDIR/fastqIDs.txt OUTPUT=$OUTDIR/$BAM.filtered.bam
    
        # get fastq
        java -Xmx30g -jar $(which SamToFastq.jar) INPUT=$OUTDIR/$BAM.filtered.bam FASTQ=${n}_R1.fastq SECOND_END_FASTQ=${n}_R2.fastq
    
        #cleanup
        rm $BAM.name_sorted.bam

    else
        samtools view -b -f 3 $OUTDIR/$n.coord_sorted.bam $LOCATION > $OUTDIR/$n.bam
        samtools index $OUTDIR/$n.bam
        java -Xmx30g -jar $(which SamToFastq.jar) VALIDATION_STRINGENCY=LENIENT INPUT=$OUTDIR/$n.bam FASTQ=$OUTDIR/${n}_R1.fastq SECOND_END_FASTQ=$OUTDIR/${n}_R2.fastq
	rm $OUTDIR/$n.bam $OUTDIR/$n.bam.bai
    fi    

    # zip
    pigz -11 $OUTDIR/${n}_R1.fastq $OUTDIR/${n}_R2.fastq

    rm $OUTDIR/$n.coord_sorted.bam $OUTDIR/$n.coord_sorted.bam.bai 

else
    samtools view -b $BAM $LOCATION > $OUTDIR/$n.bam
    samtools index $OUTDIR/$n.bam
    java -Xmx30g -jar $(which SamToFastq.jar) INPUT=$OUTDIR/$n.bam FASTQ=$OUTDIR/${n}_R1.fastq
    gzip -9 $OUTDIR/${n}_R1.fastq

fi
