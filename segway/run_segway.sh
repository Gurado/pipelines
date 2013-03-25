#!/bin/sh -e

USAGEMSG="usage: $(basename $0) CONFIGFILE

Starts the segway pipeline given the config file

Author: Fabian Buske

Requirements (modules):
	module fabbus/segway/1.1.0 
	BEDtools (genomeCoverageBed, fetchChromSizes) in path

* CONFIGFILE - the config file describing the joba
* -1 step 1 - collect the bam data from gagri
* -2 step 2 - groom the bam data into bedGraph format
* -3 step 3 - put the data as tracks into a genomedata archive
* -4 step 4 - train the segway model
* -5 step 5 - predict using the segway model
* -6 step 6 - evaluate the output using segtools
* -f force - overwrite existing results
* -a armed - trigger the jobs after writing the scripts
* -v - print progress information (verbose).
"


PIPELINEDIR=`dirname $0`
VERSION="0.0.1"

DO_TRANSFERDATA=""
DO_CONVERTBAM2BEDGRAPH=""
DO_GENERATEARCHIVE=""
DO_TRAINSEGWAY=""
DO_PREDICTSEGWAY=""
DO_EVALUATE=""
CLOBBER=""

ARMED="FALSE"
OVERWRITEALL="FALSE"

[ $# -lt 1 ] && echo "$USAGEMSG" >&2 && exit 1

while getopts "123456afv" opt;
do
        case ${opt} in
	1) DO_TRANSFERDATA="TRUE";;
	2) DO_CONVERTBAM2BEDGRAPH="TRUE";;
	3) DO_GENERATEARCHIVE="TRUE";;
	4) DO_TRAINSEGWAY="TRUE";;
	5) DO_PREDICTSEGWAY="TRUE";;
	6) DO_EVALUATE="TRUE";;
        a) ARMED="TRUE";;
        f) OVERWRITEALL="TRUE";;
        v) VERBOSE="--verbose";;
        \?) print >&2 "$0: error - unrecognized option $1"
            exit 1;;
        esac
done

shift $(($OPTIND-1))
CONFIG=$1

source ${CONFIG}

if [ "$OVERWRITEALL" = "TRUE" ];then
	CLOBBER="--clobber "
fi 

if [ -d ${TARGETDIR} ]; then
	echo "[WARN] target directory exists already"
else
	mkdir -p ${TARGETDIR}
fi

#
# data folder (to put or source the data from)
# Important: add final /
#
SEGWAY_DATA=${TARGETDIR}data/
SEGWAY_BIN=${TARGETDIR}bin/
SEGWAY_QOUT=${TARGETDIR}qout/
SEGWAY_RESULT=${TARGETDIR}result/
SEGWAY_TRAIN=${SEGWAY_RESULT}train-${EXPERIMENT}
SEGWAY_PREDICT=${SEGWAY_RESULT}predict-${EXPERIMENT}

# some housekeeping
mkdir -p $SEGWAY_DATA $SEGWAY_BIN $SEGWAY_QOUT $SEGWAY_RESULT
cp ${REGIONS} ${SEGWAY_DATA}
TRAIN_REGIONS=${SEGWAY_DATA}$(basename ${REGIONS})

##
## collect the data (tracks) from gagri
##
if [ -n "$DO_TRANSFERDATA" ];  then
	# get all files
	echo "echo  'collect data tracks from gagri'" > ${SEGWAY_BIN}1_cdata.sh

	for F in ${FILES}; do
	        echo "echo 'datafile ${F}'" >> ${SEGWAY_BIN}1_cdata.sh
        	echo "smbclient \\\\\\\\gagri\\\\GRIW -A ~/.smbclient -c 'cd ${FILES_SOURCE}/${F}; get ${F}.bam' && mv ${F}.bam ${SEGWAY_DATA}" >> ${SEGWAY_BIN}1_cdata.sh
	        echo "smbclient \\\\\\\\gagri\\\\GRIW -A ~/.smbclient -c 'cd ${FILES_SOURCE}/${F}; get ${F}.bam.bai' && mv ${F}.bam.bai ${SEGWAY_DATA}" >> ${SEGWAY_BIN}1_cdata.sh

	done

	chmod 777 ${SEGWAY_BIN}1_cdata.sh
	if [ $ARMED = "TRUE" ]; then
	        ${SEGWAY_BIN}1_cdata.sh
	fi
fi

##      
## transform the bam data into bedgraph
##
if [ -n "$DO_CONVERTBAM2BEDGRAPH" ]; then
	echo "echo 'get chromosome sizes for ${GENOME}'" > ${SEGWAY_BIN}2_tdata.sh
        echo "fetchChromSizes ${GENOME} > ${SEGWAY_DATA}/${GENOME}.chrom.sizes" >>  ${SEGWAY_BIN}2_tdata.sh

	for F in `ls ${SEGWAY_DATA}/*.bam`; do
        	FN=${F##*/}
	        FB=${FN%.*}

        	if [ ! -f ${SEGWAY_DATA}${FB}.bedgraph.gz ]; then
                	[ -f ${SEGWAY_QOUT}td4${FB}.out ] && rm ${SEGWAY_QOUT}td4${FB}.out

	                echo 'echo job_id $JOB_ID startdata $(date)' > ${SEGWAY_BIN}tdata${FB}.sh
                	echo 'echo convert ${FB}.bam to bedGraph' >>  ${SEGWAY_BIN}tdata${FB}.sh
        	        echo "genomeCoverageBed -split -bg -ibam ${F} -g ${SEGWAY_DATA}/${GENOME}.chrom.sizes | gzip > ${SEGWAY_DATA}/${FB}.bedgraph.gz" >> ${SEGWAY_BIN}tdata${FB}.sh

	                echo 'echo job_id $JOB_ID ending $(date)' >> ${SEGWAY_BIN}tdata${FB}.sh
			chmod 777 ${SEGWAY_BIN}tdata${FB}.sh

        	        # submit
                	echo "qsub -V -cwd -b y -j y -o ${SEGWAY_QOUT}/td4${FB}.out -N td4${FB} ${SEGWAY_BIN}/tdata${FB}.sh" >>  ${SEGWAY_BIN}2_tdata.sh
	        fi
	done

	chmod 777 ${SEGWAY_BIN}2_tdata.sh
	if [ $ARMED = "TRUE" ]; then
        	${SEGWAY_BIN}2_tdata.sh
	fi
fi

##
## generate the separate genome archives for the tissues.
## The archieve is then annotated with the bedgraph data  
##
if [ -n "$DO_GENERATEARCHIVE" ]; then
	## load module
	echo "module load fabbus/segway/1.1.0" > ${SEGWAY_BIN}3_gdata.sh
	echo "[ -f ${SEGWAY_QOUT}/GnDt4M${EXPERIMENT}.out ] && rm ${SEGWAY_QOUT}/GnDt4M${EXPERIMENT}.out" >> ${SEGWAY_BIN}3_gdata.sh


	echo 'echo job_id $JOB_ID startdata $(date)' > ${SEGWAY_BIN}gdata${EXPERIMENT}.sh
	echo 'echo get chromosome sizes for ${GENOME}' >> ${SEGWAY_BIN}gdata${EXPERIMENT}.sh
	[ ! -f ${GENOME}.chrom.sizes ] && echo "fetchChromSizes ${GENOME} > ${SEGWAY_DATA}/${GENOME}.chrom.sizes" >> ${SEGWAY_BIN}gdata${EXPERIMENT}.sh

	# genomedata-load call
	echo "echo '*** create genomedata archive;"  >> ${SEGWAY_BIN}gdata${EXPERIMENT}.sh
	echo "genomedata-load --sizes -s ${SEGWAY_DATA}/${GENOME}.chrom.sizes \\" >> ${SEGWAY_BIN}gdata${EXPERIMENT}.sh
	# add the -t <ID>=<FILE> sections for all tracks

	for f in $(ls $SEGWAY_DATA/*.gz ); do
	        b=$(basename $f)
        	arrIN=(${b//./ })
        echo "-t "${arrIN[0]}=$f" \\" >> ${SEGWAY_BIN}gdata${EXPERIMENT}.sh
	done
	echo "${SEGWAY_DATA}/${EXPERIMENT}.genomedata" >> ${SEGWAY_BIN}gdata${EXPERIMENT}.sh
	echo 'echo job_id $JOB_ID ending $(date)' >> ${SEGWAY_BIN}gdata${EXPERIMENT}.sh
	chmod 777 ${SEGWAY_BIN}gdata${EXPERIMENT}.sh
	#submit
	echo "qsub -pe smp 2 -V -cwd -b y -j y -o ${SEGWAY_QOUT}/GnDt4M${EXPERIMENT}.out -N GnDt4M${EXPERIMENT} ${SEGWAY_BIN}/gdata${EXPERIMENT}.sh" >> ${SEGWAY_BIN}3_gdata.sh
	chmod 777  ${SEGWAY_BIN}3_gdata.sh

	if [ $ARMED = "TRUE" ]; then
	        ${SEGWAY_BIN}3_gdata.sh
	fi
fi

##
## train Seqway models
##
if [ -n "$DO_TRAINSEGWAY" ]; then
	echo "module load fabbus/segway/1.1.0" > ${SEGWAY_BIN}4_train.sh
        echo "[ -f ${SEGWAY_QOUT}sgtrn4M${EXPERIMENT}.out ] && rm ${SEGWAY_QOUT}sgtrn4M${EXPERIMENT}.out" >> ${SEGWAY_BIN}4_train.sh
        echo "[ -d ${SEGWAY_TRAIN} ] && rm -r ${SEGWAY_TRAIN}" >> ${SEGWAY_BIN}4_train.sh

        OPTIONS="--include-coords=${TRAIN_REGIONS} --num-labels=${LABELS} --num-instances=${INSTANCES} ${CLOBBER} ${SPECIAL}"
	echo "echo '*** train segway'" >> ${SEGWAY_BIN}segtrain${EXPERIMENT}.sh
        echo "segway $OPTIONS \\"> ${SEGWAY_BIN}segtrain${EXPERIMENT}.sh 

        # add the --track <ID> sections
        for f in $(ls $SEGWAY_DATA/*.bedgraph.gz ); do
                b=$(basename $f)
                arrIN=(${b//./ })
        echo "--track "${arrIN[0]}" \\" >> ${SEGWAY_BIN}segtrain${EXPERIMENT}.sh
        done
        echo "train ${SEGWAY_DATA}${EXPERIMENT}.genomedata ${SEGWAY_TRAIN}" >> ${SEGWAY_BIN}segtrain${EXPERIMENT}.sh

        chmod 777 ${SEGWAY_BIN}segtrain${EXPERIMENT}.sh
        #echo "qsub -l mem_requested=16G -V -cwd -b y -j y -o ${SEGWAY_QOUT}/sgtrn4M${EXPERIMENT}.out -N sgtrn4M${EXPERIMENT} ${SEGWAY_BIN}/segtrain${EXPERIMENT}.sh"  >> ${SEGWAY_BIN}4_train.sh
        echo "${SEGWAY_BIN}segtrain${EXPERIMENT}.sh"  >> ${SEGWAY_BIN}4_train.sh
	chmod 777  ${SEGWAY_BIN}4_train.sh

        if [ $ARMED = "TRUE" ]; then
                ${SEGWAY_BIN}4_train.sh
        fi
fi

##
## predict using a trained Seqway model
##
if [ -n "$DO_PREDICTSEGWAY" ]; then
	echo "module load fabbus/segway/1.1.0" > ${SEGWAY_BIN}5_predict.sh
        echo "[ -f ${SEGWAY_QOUT}/sgprd4M${EXPERIMENT}.out ] && rm ${SEGWAY_QOUT}sgprd4M${EXPERIMENT}.out" >> ${SEGWAY_BIN}5_predict.sh
        echo "[ -d ${SEGWAY_PREDICT} ] && rm -r ${SEGWAY_PREDICT}" >> ${SEGWAY_BIN}5_predict.sh

        echo 'echo job_id $JOB_ID startdata $(date)' > ${SEGWAY_BIN}segpredict${EXPERIMENT}.sh

        echo 'export TMP=/tmp/' >> ${SEGWAY_BIN}segpredict${EXPERIMENT}.sh
        echo 'export TEMP=/tmp/' >> ${SEGWAY_BIN}segpredict${EXPERIMENT}.sh
        echo 'export TMPDIR=/tmp/' >> ${SEGWAY_BIN}segpredict${EXPERIMENT}.sh

        # segway calli
	echo "echo '*** predict segmentation'" >>  ${SEGWAY_BIN}segpredict${EXPERIMENT}.sh
        echo "segway --num-labels=${LABELS} ${CLOBBER} \\">> ${SEGWAY_BIN}segpredict${EXPERIMENT}.sh
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
        echo "${SEGWAY_BIN}segpredict${EXPERIMENT}.sh" >> ${SEGWAY_BIN}5_predict.sh
	chmod 777 ${SEGWAY_BIN}5_predict.sh

	if [ $ARMED = "TRUE" ]; then
                ${SEGWAY_BIN}5_predict.sh
        fi
fi

##
## Evaluate prediction
##
if [ -n "$DO_EVALUATE" ]; then
	echo "module load fabbus/segway/1.1.0" > ${SEGWAY_BIN}6_evaluate.sh
        echo 'echo job_id $JOB_ID startdata $(date)' > ${SEGWAY_BIN}segeval${EXPERIMENT}.sh     
        #preprocess file
        if [ -n $OVERWRITEALL ] || [ ! -f ${SEGWAY_PREDICT}/segway.bed.gz.pkl.gz ]; then
                echo "echo '*** preprocess'" >> ${SEGWAY_BIN}segeval${EXPERIMENT}.sh     
                echo "segtools-preprocess ${SEGWAY_PREDICT}/segway.bed.gz" >> ${SEGWAY_BIN}segeval${EXPERIMENT}.sh 
        fi
        
        echo "echo '*** lengthdist'" >> ${SEGWAY_BIN}segeval${EXPERIMENT}.sh
        echo "segtools-length-distribution ${SEGWAY_PREDICT}/segway.bed.gz.pkl.gz --outdir=${SEGWAY_RESULT}length-dist/ ${CLOBBER}" >> ${SEGWAY_BIN}segeval${EXPERIMENT}.sh 
        
        echo "echo '*** geneagg'" >> ${SEGWAY_BIN}segeval${EXPERIMENT}.sh
        echo "segtools-aggregation ${SEGWAY_PREDICT}/segway.bed.gz.pkl.gz ${ANNOTATION} --normalize --mode=gene --outdir=${SEGWAY_RESULT}gencode-agg/ ${CLOBBER}" >> ${SEGWAY_BIN}segeval${EXPERIMENT}.sh 
        
        echo "echo '----------gmtkparam'" >> ${SEGWAY_BIN}segeval${EXPERIMENT}.sh
        echo "segtools-gmtk-parameters ${SEGWAY_TRAIN}/params/params.params --outdir=${SEGWAY_RESULT}gtmk-param/ ${CLOBBER}" >> ${SEGWAY_BIN}segeval${EXPERIMENT}.sh 
        
        echo "echo '----------html'" >> ${SEGWAY_BIN}segeval${EXPERIMENT}.sh
        echo "cd ${SEGWAY_RESULT}" >> ${SEGWAY_BIN}segeval${EXPERIMENT}.sh
        echo "segtools-html-report -L ${SEGWAY_PREDICT}/segway.layered.bed.gz --results-dir=${SEGWAY_RESULT} -o segtools.html ${SEGWAY_PREDICT}/segway.bed.gz.pkl.gz ${CLOBBER}" >> ${SEGWAY_BIN}segeval${EXPERIMENT}.sh 
        echo "sed 's|${SEGWAY_RESULT}||g' segtools.html > segtools.html2" >> ${SEGWAY_BIN}segeval${EXPERIMENT}.sh 
        echo "mv segtools.html2 segtools.html" >> ${SEGWAY_BIN}segeval${EXPERIMENT}.sh 
        echo "cd $(pwd)" >> ${SEGWAY_BIN}segeval${EXPERIMENT}.sh 

        chmod 777 ${SEGWAY_BIN}segeval${EXPERIMENT}.sh

        echo "qsub -l mem_requested=16G -V -cwd -b y -j y -o ${SEGWAY_QOUT}sgevl4M${EXPERIMENT}.out -N sgevl4M${EXPERIMENT} ${SEGWAY_BIN}segeval${EXPERIMENT}.sh" >> ${SEGWAY_BIN}6_evaluate.sh

	chmod 777 ${SEGWAY_BIN}6_evaluate.sh
        if [ $ARMED = "TRUE" ]; then
                ${SEGWAY_BIN}6_evaluate.sh
        fi
fi


