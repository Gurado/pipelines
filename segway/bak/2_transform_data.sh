#!/bin/sh

DIR=`dirname $0`
source ${DIR}/0_config.sh

##
## transform the bam data into bedgraph
##

if [ -n $DO_CONVERTBAM2BEDGRAPH ]; then

        echo "get chromosome sizes for ${GENOME}"
        [ ! -f ${GENOME}.chrom.sizes ] && fetchChromSizes ${GENOME} > ${SEGWAY_DATA}/${GENOME}.chrom.sizes

	for F in `ls ${SEGWAY_DATA}/*.bam`; do

		FN=${F##*/}
		FB=${FN%.*}

		if [ ! -f ${SEGWAY_DATA}/${FB}.bedgraph.gz ]; then
			[ -f ${SEGWAY_QOUT}/td4${FB}.out ] && rm ${SEGWAY_QOUT}/td4${FB}.out

        	        echo 'echo job_id $JOB_ID startdata $(date)' > ${SEGWAY_BIN}/tdata${FB}.sh
	                echo 'echo convert ${FB}.bam to bedGraph' >>  ${SEGWAY_BIN}/tdata${FB}.sh
			echo "genomeCoverageBed -split -bg -ibam ${F} -g ${SEGWAY_DATA}/${GENOME}.chrom.sizes | gzip > ${SEGWAY_DATA}/${FB}.bedgraph.gz" >> ${SEGWAY_BIN}/tdata${FB}.sh	
	
			echo 'echo job_id $JOB_ID ending $(date)' >> ${SEGWAY_BIN}/gdata$i.sh
			chmod 777 ${SEGWAY_BIN}/tdata${FB}.sh
			# submit
			qsub -V -cwd -b y -j y -o ${SEGWAY_QOUT}/td4${FB}.out -N td4${FB} ${SEGWAY_BIN}/tdata${FB}.sh
		fi
	done
fi

