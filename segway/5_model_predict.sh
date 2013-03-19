#!/bin/sh -e

DIR=`dirname $0`
source ${DIR}/0_config.sh

## load module

module load fabbus/segway/1.1.0

##
## predict usign a trained Seqway model
##

if [ -n "$DO_PREDICTSEGWAY" ]; then
	[ -f ${SEGWAY_QOUT}/sgprd4M${EXPERIMENT}.out ] && rm ${SEGWAY_QOUT}sgprd4M${EXPERIMENT}.out
	[ -d ${SEGWAY_PREDICT} ] && rm -r ${SEGWAY_PREDICT}

	echo 'echo job_id $JOB_ID startdata $(date)' > ${SEGWAY_BIN}segpredict${EXPERIMENT}.sh     

	echo 'export TMP=/tmp/' >> ${SEGWAY_BIN}segpredict${EXPERIMENT}.sh
	echo 'export TEMP=/tmp/' >> ${SEGWAY_BIN}segpredict${EXPERIMENT}.sh
	echo 'export TMPDIR=/tmp/' >> ${SEGWAY_BIN}segpredict${EXPERIMENT}.sh
        
	# segway call    
	echo "segway --num-labels=${LABELS} \\">> ${SEGWAY_BIN}segpredict${EXPERIMENT}.sh
        # add the --track <ID> sections
        for f in $(ls ${SEGWAY_DATA}/*.bedgraph.gz ); do
            b=$(basename $f)
            arrIN=(${b//./ })
            echo "--track "${arrIN[0]}" \\" >> ${SEGWAY_BIN}segpredict${EXPERIMENT}.sh
        done
        echo "identify ${SEGWAY_DATA}${EXPERIMENT}.genomedata ${SEGWAY_TRAIN} ${SEGWAY_PREDICT}" >> ${SEGWAY_BIN}segpredict${EXPERIMENT}.sh
        echo 'echo job_id $JOB_ID ending $(date)' >> ${SEGWAY_BIN}segpredict${EXPERIMENT}.sh
        chmod 777 ${SEGWAY_BIN}segpredict${EXPERIMENT}.sh
        # submit
#       echo "qsub -l mem_requested=16G -V -cwd -b y -j y -o ${SEGWAY_QOUT}sgprd4M${EXPERIMENT}.out -N sgprd4M${EXPERIMENT} ${SEGWAY_BIN}segpredict${EXPERIMENT}.sh"
	${SEGWAY_BIN}segpredict${EXPERIMENT}.sh
fi
