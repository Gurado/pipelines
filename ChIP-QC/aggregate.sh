#!/bin/sh -e

USAGEMSG="usage: $(basename $0) INPUTFOLDER OUTPUTFOLDER

Aggregates results from NCIS and CHANCE QC calculations

Author: Fabian Buske

* INPUTFOLDER - contains the results of NCIS (*.txt) and CHANCE (*.IPstrength)
* OUTPUTFOLDER - specifies where the results are written to
"

DIR=$(dirname $0)
VERSION="0.0.1"

[ $# -lt 2 ] && echo "$USAGEMSG" >&2 && exit 1

CHIP=""
CONTROL=""
GAGRI=/Cancer-Epigenetics/Data/ClarkLab/Seq/ChIP-Seq/hg19/
NCIS=/home/fabbus/pipelines/NCIS/NCIS.R

while getopts "" opt;
do
        case ${opt} in
        \?) print >&2 "$0: error - unrecognized option $1"
                exit 1;;
        esac
done

shift $(($OPTIND-1))
SOURCE=$1
OUTPUT=$2

[ ! -d ${SOURCE} ] && echo "INPUT directory does not exist: ${INPUT}" && exit 1
mkdir -p ${OUTPUT}

$echo "ChIP	Control	NCIS_BackgroundRatio	NCIS_Binsize	NCIS_SeqDepth	NCIS_Normalization	CHANCE_Enriched	CHANCE_Cumulative	CHANCE_Scaling	CHANCE_FDR" > ${OUTPUT}/ChIP-QC_aggregate.txt


for line in $(ls -la ${SOURCE} | tail -n+4 | awk '{print $NF}' | awk -F\. '{print $1}' | sort -u ); do
	b=$(basename $line)
        chip=$(echo ${b} | awk -F\- '{print $1}')
        input=$(echo ${b} | awk -F\- '{print $2}')
	echo  "$chip $input"

	NCIS_BACKGROUNDRATIO=""
	NCIS_BINSIZE=""
	NCIS_NORMALIZATION=""
	NCIS_SEQDEPTH=""

	CHANCE_FDR=""
	CHANCE_SCALING=""
	CHANCE_ENRICHED=""
	CHANCE_CUMMULATIVE=""

	if [ -f ${line}.txt ]; then
		f=${line}.txt
		NCIS_BACKGROUNDRATIO=$(tail -n 1 $f | awk '{print $5}')
		NCIS_BINSIZE=$(tail -n 1 $f | awk '{print $3}')
		NCIS_NORMALIZATION=$(tail -n 1 $f | awk '{print $2}')
		NCIS_SEQDEPTH=$(tail -n 1 $f | awk '{print $4}')
	fi


	if [ -f ${line}.IPstrength ]; then
		f=${line}.IPstrength
	        CHANCE_FDR=$(cat $f | grep "^fdr," | awk -F\, '{print $2}')
		CHANCE_FDR_TFBS_NORMAL=$(cat $f | grep "^tfbs_normal_fdr" | awk -F\, '{print $2}')
                CHANCE_FDR_HISTONE_NORMAL=$(cat $f | grep "^tfbs_normal_fdr" | awk -F\, '{print $2}')
                CHANCE_FDR_TFBS_CANCER=$(cat $f | grep "^histone_cancer_fdr" | awk -F\, '{print $2}')
                CHANCE_FDR_HISTONE_CANCER=$(cat $f | grep "^histone_cancer_fdr" | awk -F\, '{print $2}')
        	CHANCE_SCALING=$(cat $f | grep "input_scaling_factor" | awk -F\, '{print $2}')
	        CHANCE_ENRICHED=$(cat $f | grep "percent_genome_enriched" | awk -F\, '{print $2}')
        	CHANCE_CUMMULATIVE=$(cat $f | grep "differential_percentage_enrichment" | awk -F\, '{print $2}')
	fi

	echo "${chip}	${input}	${NCIS_BACKGROUNDRATIO}	${NCIS_BINSIZE}	${NCIS_SEQDEPTH}	${NCIS_NORMALIZATION}	${CHANCE_ENRICHED}	${CHANCE_CUMMULATIVE}	${CHANCE_SCALING}	${CHANCE_FDR}	${CHANCE_FDR_TFBS_NORMAL}	${CHANCE_FDR_HISTONE_NORMAL}	${CHANCE_FDR_TFBS_CANCER}	${CHANCE_FDR_HISTONE_CANCER}" >> ${OUTPUT}/ChIP-QC_aggregate.txt
done

exit 1

RSCRIPT=${OUTPUT}/ChIP-QC_aggregate.R

echo "library(gridExtra)" > ${RSCRIPT}
echo "library(ggplot2)" >> ${RSCRIPT}
echo "library(scales)" >> ${RSCRIPT}
echo "data <- read.delim('${OUTPUT}/ChIP-QC_aggregate.txt', head=TRUE)" >> ${RSCRIPT}
echo "data <- cbind(data, 'CHANCE_BackgroundRatio'=0)" >> ${RSCRIPT}
echo "data['CHANCE_BackgroundRatio'] <- (100 - data['CHANCE_Cumulative'])/100" >> ${RSCRIPT}
echo "ChIP_cell <- sapply(strsplit(as.character(data[,1]), '_'), '[[', 2)" >> ${RSCRIPT}
echo "Control_cell <- sapply(strsplit(as.character(data[,2]), '_'), '[[', 2)" >> ${RSCRIPT}
echo "ChIP_mark <- sapply(strsplit(as.character(data[,1]), '_'), '[[', 3)" >> ${RSCRIPT}
echo "Control_mark <- sapply(strsplit(as.character(data[,2]), '_'), '[[', 3)" >> ${RSCRIPT}
echo "data <- cbind(data, ChIP_cell, Control_cell, ChIP_mark, Control_mark)" >> ${RSCRIPT}
echo "data <- data[with(data, order(ChIP_cell, ChIP_mark, Control_cell, Control_mark)),]"  >> ${RSCRIPT}
echo "data[,1] <- factor(data[,1], levels=unique(data[,1]))"  >> ${RSCRIPT}
echo "data[,2] <- factor(data[,2], levels=unique(data[,2]))"  >> ${RSCRIPT}
echo "data['CHANCE_FDR'] <- sapply(data['CHANCE_FDR'], log10)" >> ${RSCRIPT}
echo "data['CHANCE_FDR'] <- -data['CHANCE_FDR'] "  >>  ${RSCRIPT}
echo "fmt <- '%.2f'" >> ${RSCRIPT}
echo "p1 <- ggplot(data, aes(Control,ChIP)) + geom_tile(aes(fill = BackgroundRatio), colour = 'white') + scale_fill_gradient2(limits=c(0.25, 1.),midpoint=0.63, high=muted('red'), mid='yellow', low=muted('steelblue'), na.value='steelblue') + theme(axis.text.x = element_text(angle = 45, hjust=1, size=8), axis.text.y=element_text(size=8)) + theme(legend.direction = 'horizontal', legend.position = 'top')+geom_text(aes(label=sprintf(fmt, BackgroundRatio)),size=2)"  >> ${RSCRIPT}


echo "pdf(file='${OUTPUT}/ChIP-QC_aggregate.pdf')" >> ${RSCRIPT}
echo "p1" >> ${RSCRIPT}
echo "dev.off()" >> CHANCE_aggregate.R

/share/ClusterShare/software/contrib/Cancer-Epigenetics/tools/bin/Rscript ${OUTPUT}//ChIP-QC_aggregate.R
