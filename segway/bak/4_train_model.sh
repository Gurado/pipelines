#!/bin/sh -e

DIR=`dirname $0`
source ${DIR}/0_config.sh

## load module

module load fabbus/segway/1.1.0

##
## train Seqway models
##

if [ -n "$DO_TRAINSEGWAY" ]; then
	[ -f ${SEGWAY_QOUT}sgtrn4M${EXPERIMENT}.out ] && rm ${SEGWAY_QOUT}sgtrn4M${EXPERIMENT}.out
	[ -d ${SEGWAY_TRAIN} ] && rm -r ${SEGWAY_TRAIN}

	echo 'echo job_id $JOB_ID startdata $(date)' > ${SEGWAY_BIN}segtrain${EXPERIMENT}.sh
	OPTIONS="--include-coords=${TRAIN_REGIONS} --num-labels=${LABELS} --num-instances=${INSTANCES} ${SPECIAL}"
	echo "segway $OPTIONS \\">> ${SEGWAY_BIN}segtrain${EXPERIMENT}.sh 

	# add the --track <ID> sections
	for f in $(ls $SEGWAY_DATA/*.bedgraph.gz ); do
		b=$(basename $f)
		arrIN=(${b//./ })
	echo "--track "${arrIN[0]}" \\" >> ${SEGWAY_BIN}segtrain${EXPERIMENT}.sh
	done
	echo "train ${SEGWAY_DATA}${EXPERIMENT}.genomedata ${SEGWAY_TRAIN}" >> ${SEGWAY_BIN}segtrain${EXPERIMENT}.sh
	echo 'echo job_id $JOB_ID ending $(date)' >> ${SEGWAY_BIN}segtrain${EXPERIMENT}.sh
	chmod 777 ${SEGWAY_BIN}segtrain${EXPERIMENT}.sh
	# submit
	#qsub -l mem_requested=16G -V -cwd -b y -j y -o ${SEGWAY_QOUT}/sgtrn4M${EXPERIMENT}.out -N sgtrn4M${EXPERIMENT} ${SEGWAY_BIN}/segtrain${EXPERIMENT}.sh
	${SEGWAY_BIN}segtrain${EXPERIMENT}.sh
fi
