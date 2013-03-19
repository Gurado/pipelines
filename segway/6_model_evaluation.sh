#!/bin/sh -e

DIR=`dirname $0`
source ${DIR}/0_config.sh

## load module

module load fabbus/segway/1.1.0

##
## Evaluate
##

if [ -n "$DO_EVALUATE" ]; then

	echo 'echo job_id $JOB_ID startdata $(date)' > ${SEGWAY_BIN}segeval${EXPERIMENT}.sh     
	#preprocess file
	if [ -n $OVERWRITEALL ] || [ ! -f ${SEGWAY_PREDICT}segway.bed.gz.pkl.gz ]; then
		echo "echo '----------preprocess'" >> ${SEGWAY_BIN}segeval${EXPERIMENT}.sh     
		echo "segtools-preprocess ${SEGWAY_PREDICT}segway.bed.gz" >> ${SEGWAY_BIN}segeval${EXPERIMENT}.sh 
	fi
	
	echo "echo '----------lengthdist'" >> ${SEGWAY_BIN}segeval${EXPERIMENT}.sh
	echo "segtools-length-distribution ${SEGWAY_PREDICT}segway.bed.gz.pkl.gz --outdir=${SEGWAY_RESULT}length-dist/ --clobber" >> ${SEGWAY_BIN}segeval${EXPERIMENT}.sh 
	
	echo "echo '----------geneagg'" >> ${SEGWAY_BIN}segeval${EXPERIMENT}.sh
	echo "segtools-aggregation ${SEGWAY_PREDICT}segway.bed.gz.pkl.gz ${SEGWAY_DATA}${ANNOTATION} --normalize --mode=gene --outdir=${SEGWAY_RESULT}gencode-agg/ --clobber" >> ${SEGWAY_BIN}segeval${EXPERIMENT}.sh 
	
	echo "echo '----------gmtkparam'" >> ${SEGWAY_BIN}segeval${EXPERIMENT}.sh
	echo "segtools-gmtk-parameters ${SEGWAY_TRAIN}params/params.params --outdir=${SEGWAY_RESULT}gtmk-param/ --clobber" >> ${SEGWAY_BIN}segeval${EXPERIMENT}.sh 
	
	echo "echo '----------html'" >> ${SEGWAY_BIN}segeval${EXPERIMENT}.sh
	echo "cd ${SEGWAY_RESULT}" >> ${SEGWAY_BIN}segeval${EXPERIMENT}.sh
	echo "segtools-html-report -o segtools.html ${SEGWAY_PREDICT}segway.bed.gz.pkl.gz --clobber" >> ${SEGWAY_BIN}segeval${EXPERIMENT}.sh 
	echo "sed 's|${SEGWAY_RESULT}||g' segtools.html > segtools.html2" >> ${SEGWAY_BIN}segeval${EXPERIMENT}.sh 
	echo "mv segtools.html2 segtools.html" >> ${SEGWAY_BIN}segeval${EXPERIMENT}.sh 
	echo "cd $(pwd)" >> ${SEGWAY_BIN}segeval${EXPERIMENT}.sh 

	chmod 777 ${SEGWAY_BIN}segeval${EXPERIMENT}.sh

	echo "qsub -l mem_requested=16G -V -cwd -b y -j y -o ${SEGWAY_QOUT}sgevl4M${EXPERIMENT}.out -N sgevl4M${EXPERIMENT} ${SEGWAY_BIN}segeval${EXPERIMENT}.sh"
fi

