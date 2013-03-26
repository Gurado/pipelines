#!/bin/sh -e

DIR=`dirname $0`
source ${DIR}/0_config.sh

## load module

module load fabbus/segway/1.1.0

##
## generate the separate genome archives for the tissues.
## The archieve is then annotated with the bedgraph data  
##

if [ -n "$DO_GENERATEARCHIVE" ]; then
	[ -f ${SEGWAY_QOUT}/GnDt4M${EXPERIMENT}.out ] && rm ${SEGWAY_QOUT}/GnDt4M${EXPERIMENT}.out

	echo 'echo job_id $JOB_ID startdata $(date)' > ${SEGWAY_BIN}/gdata${EXPERIMENT}.sh
	echo 'echo get chromosome sizes for ${GENOME}' >> ${SEGWAY_BIN}/gdata${EXPERIMENT}.sh
        [ ! -f ${GENOME}.chrom.sizes ] && echo "fetchChromSizes ${GENOME} > ${SEGWAY_DATA}/${GENOME}.chrom.sizes" >> ${SEGWAY_BIN}/gdata${EXPERIMENT}.sh

	# genomedata-load call
	echo "genomedata-load --sizes -s ${SEGWAY_DATA}/${GENOME}.chrom.sizes \\" >> ${SEGWAY_BIN}/gdata${EXPERIMENT}.sh
	# add the -t <ID>=<FILE> sections for all tracks

	for f in $(ls $SEGWAY_DATA/*.gz ); do
	        b=$(basename $f)
		arrIN=(${b//./ })
		echo "-t "${arrIN[0]}=$f" \\" >> ${SEGWAY_BIN}/gdata${EXPERIMENT}.sh
	done
	echo "${SEGWAY_DATA}/${EXPERIMENT}.genomedata" >> ${SEGWAY_BIN}/gdata${EXPERIMENT}.sh
	echo 'echo job_id $JOB_ID ending $(date)' >> ${SEGWAY_BIN}/gdata${EXPERIMENT}.sh
	chmod 777 ${SEGWAY_BIN}/gdata${EXPERIMENT}.sh
	#submit
	qsub -V -cwd -b y -j y -o ${SEGWAY_QOUT}/GnDt4M${EXPERIMENT}.out -N GnDt4M${EXPERIMENT} ${SEGWAY_BIN}/gdata${EXPERIMENT}.sh

fi


