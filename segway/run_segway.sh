#!/bin/sh -e

USAGEMSG="usage: $(basename $0) CONFIGFILE

Starts the segway pipeline given the config file

Author: Fabian Buske

Requirements (modules):
	gi/gcc/4.8.2
	fabbus/segway_gbr/1.2.0 
	BEDtools (genomeCoverageBed) in path
	Samtools gi/samtools

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
VERSION="0.0.2"

DO_TRANSFERDATA=""
DO_CONVERTDATA2BEDGRAPH=""
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
        2) DO_CONVERTDATA2BEDGRAPH="TRUE";;
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

if [ -d ${TARGETDIR}/ ]; then
	echo "[WARN] target directory exists already"
else
	mkdir -p ${TARGETDIR}/
fi

#
# data folder (to put or source the data from)
# Important: add final /
#
SEGWAY_DATA=${TARGETDIR}/data/
SEGWAY_BIN=${TARGETDIR}/bin/
SEGWAY_QOUT=${TARGETDIR}/qout/
SEGWAY_RESULT=${TARGETDIR}/result_${LABELS}/
SEGWAY_TRAIN=${SEGWAY_RESULT}/train-${EXPERIMENT}
SEGWAY_PREDICT=${SEGWAY_RESULT}/predict-${EXPERIMENT}${PREDICTION}

# some housekeeping
mkdir -p $SEGWAY_DATA $SEGWAY_BIN $SEGWAY_QOUT $SEGWAY_RESULT
cp ${REGIONS} ${SEGWAY_DATA}
TRAIN_REGIONS=${SEGWAY_DATA}$(basename ${REGIONS})

module load gi/ucsc_utils fabbus/segway_gbr/ gi/samtools
##
## collect the data (tracks) from gagri
##
if [ -n "$DO_TRANSFERDATA" ];  then
	# get all files
	echo "echo  'collect data tracks from Gagri'" > ${SEGWAY_BIN}/1_cdata.sh
	echo "module load gi/samtools" > ${SEGWAY_BIN}/1_cdata.sh
	while IFS=$'\t' read -r -a DATA; do
        FILES_SOURCE=${DATA[2]}
        for F in $(echo ${DATA[3]} | tr ',' '\n'); do
            echo "echo 'datafile ${F}${DATA[4]}'" >> ${SEGWAY_BIN}/1_cdata.sh
            if [[ ! -f ${SEGWAY_DATA}/${DATA[4]} ]] || [ $OVERWRITEALL == 'TRUE' ]; then 
                echo "[ ! -f ${SEGWAY_DATA}/$F${DATA[4]} ] && smbclient //Gagri/diverse -A ~/.smbclient -c 'cd ${FILES_SOURCE}; get ${F}${DATA[4]}' && mv ${F}${DATA[4]} ${SEGWAY_DATA}" >> ${SEGWAY_BIN}/1_cdata.sh
                if [[ "${DATA[4]}" =~ ".bam" ]] && [[ ! -f ${SEGWAY_DATA}/$F${DATA[4]}.bai ]]; then
    	           echo "[ ! -f ${SEGWAY_DATA}/$F${DATA[4]}.bai ] && smbclient //Gagri/diverse -A ~/.smbclient -c 'cd ${FILES_SOURCE}; get ${F}${DATA[4]}.bai' && mv ${F}${DATA[4]}.bai ${SEGWAY_DATA}" >> ${SEGWAY_BIN}/1_cdata.sh
    	           echo "[ ! -f ${SEGWAY_DATA}/$F${DATA[4]}.bai ] && samtools index ${SEGWAY_DATA}/$F${DATA[4]}" >> ${SEGWAY_BIN}/1_cdata.sh
    	        fi
	        fi
        done
    done < $EXPERIMENTS

	chmod 777 ${SEGWAY_BIN}/1_cdata.sh
	if [ $ARMED = "TRUE" ]; then
        ${SEGWAY_BIN}/1_cdata.sh
	fi
fi

##      
## transform the bam data into bedgraph
##
if [ -n "$DO_CONVERTDATA2BEDGRAPH" ]; then
	if [ ! -d ${SEQDIR} ] || [ ! -d ${WIGGLER_UMAPDIR} ]; then
		echo "sequence dir or umap dir is missing/invalid"
		echo "seqDir: ${SEQDIR}"
		echo "umapDIR:${WIGGLER_UMAPDIR}"
		exit 1
	fi
	
	module load gi/bedtools gi/samtools fabbus/wiggler/2.0 gi/ucsc_utils gi/pigz
	echo "module load gi/bedtools gi/samtools fabbus/wiggler/2.0 gi/ucsc_utils gi/pigz" > ${SEGWAY_BIN}/2_tdata.sh
#	echo "echo 'get chromosome sizes for ${GENOME}'" >> ${SEGWAY_BIN}/2_tdata.sh

	while IFS=$'\t' read -r -a DATA; do
        if [[ "${DATA[4]}" =~ ".bam" ]]; then
            FRAGMENTS=""
            INPUTS=""
            REPNAMES=""
            for REPLICA in $(echo ${DATA[3]} | tr ',' '\n'); do
    	        [ ! -f ${SEGWAY_DATA}${REPLICA}${DATA[4]} ] && echo "[ERROR] file not found: ${SEGWAY_DATA}${REPLICA}${DATA[4]}" && exit 1
    		INPUTS="${INPUTS} -i=${SEGWAY_DATA}${REPLICA}${DATA[4]}"
                if [[ $(samtools view ${SEGWAY_DATA}${REPLICA}${DATA[4]} |  awk '{print length($10)}' | head -1000 | sort -u | head -n 1) -ge 100 ]]; then	
                    echo "using $WIGGLER_UMAPDIR_ge100"
                    WIGGLER_UMAP=$WIGGLER_UMAPDIR_ge100
                else
                    echo "using $WIGGLER_UMAPDIR_lt100"
                    WIGGLER_UMAP=$WIGGLER_UMAPDIR_lt100
                fi
            done
            for FRAGSIZE in $(echo ${DATA[5]} | tr ',' '\n'); do
                FRAGMENTS="${FRAGMENTS} -l=${FRAGSIZE}"
            done
    
    		REPNAMES=${DATA[0]}"_"${DATA[1]}
            echo $REPNAMES
            
            if [ ! -f ${SEGWAY_DATA}${REPNAMES}.bg.gz ] || [ "$OVERWRITEALL" = "TRUE" ]; then
                [ -f ${SEGWAY_QOUT}TrDa-${REPNAMES}.out ] && rm ${SEGWAY_QOUT}TrDa-${REPNAMES}.out
                [ -f ${SEGWAY_DATA}${REPNAMES}.bg.gz ] && rm ${SEGWAY_DATA}${REPNAMES}.bg.gz
                echo '#!/bin/bash' > ${SEGWAY_BIN}/tdata${REPNAMES}.sh
                echo 'echo job_id $JOB_ID startdata $(date)' >> ${SEGWAY_BIN}/tdata${REPNAMES}.sh
                echo "echo convert ${REPNAMES} to bedGraph using wiggler" >>  ${SEGWAY_BIN}/tdata${REPNAMES}.sh
                # wiggler
                
                echo "align2rawsignal -of=bg ${INPUTS} ${FRAGMENTS} -s=${SEQDIR} -u=${WIGGLER_UMAP} -n=5 -v=${SEGWAY_QOUT}wiggler-${REPNAMES}.log -k=tukey -w=${DATA[6]} -o=${SEGWAY_DATA}${REPNAMES}.bg" >> ${SEGWAY_BIN}/tdata${REPNAMES}.sh 
##		 apply asinh transformation to all signal values (done via (--distribution=asinh_norm) by default
#		echo "[ -f -c ${SEGWAY_DATA}${REPNAMES}.bg.gz ] && rm -c ${SEGWAY_DATA}${REPNAMES}.bg.gz" >> ${SEGWAY_BIN}/tdata${REPNAMES}.sh
#		echo "cat ${SEGWAY_DATA}${REPNAMES}.bg | awk 'function asinh(x){return log(x + sqrt(x*x +1))}{OFS=\"\\t\";print \$1,\$2,\$3,asinh(\$4)}' | pigz -9 -c > ${SEGWAY_DATA}${REPNAMES}.bg.gz " >> ${SEGWAY_BIN}/tdata${REPNAMES}.sh
#		echo "rm ${SEGWAY_DATA}${REPNAMES}.bg" >> ${SEGWAY_BIN}/tdata${REPNAMES}.sh
		echo "pigz -9 ${SEGWAY_DATA}${REPNAMES}.bg" >> ${SEGWAY_BIN}/tdata${REPNAMES}.sh
                echo 'echo job_id $JOB_ID ending $(date)' >> ${SEGWAY_BIN}/tdata${REPNAMES}.sh
                chmod 777 ${SEGWAY_BIN}/tdata${REPNAMES}.sh
                
                # submit
                echo "qsub -pe smp 4 -V -cwd -l h_rt=24:00:00 -j y -m e -S /bin/bash -o ${SEGWAY_QOUT}TrDa-${REPNAMES}.out -N TrDa-${REPNAMES} ${SEGWAY_BIN}/tdata${REPNAMES}.sh" >>  ${SEGWAY_BIN}/2_tdata.sh
            fi
        elif [[ "${DATA[4]}" =~ ".bw" ]]; then
            if [ "$(echo ${DATA[3]} | tr ',' '\n' | wc -l | cut -f 1)" -gt 1 ]; then
                echo "[ERROR] replicates for bigwigs not supported"
                exit 1
            fi

    	    REPNAMES=${DATA[0]}"_"${DATA[1]}
            echo $REPNAMES
        	if [ ! -f ${SEGWAY_DATA}${REPNAMES}.bg.gz ] || [ "$OVERWRITEALL" = "TRUE" ]; then
                [ -f ${SEGWAY_QOUT}TrDa-${REPNAMES}.out ] && rm ${SEGWAY_QOUT}TrDa-${REPNAMES}.out
                echo '#!/bin/bash' > ${SEGWAY_BIN}/tdata${REPNAMES}.sh
                echo 'echo job_id $JOB_ID startdata $(date)' >> ${SEGWAY_BIN}/tdata${REPNAMES}.sh
                echo "echo convert ${REPNAMES} bigwig to bedGraph " >>  ${SEGWAY_BIN}/tdata${REPNAMES}.sh
                # wiggler
                echo "bigWigToBedGraph ${SEGWAY_DATA}${DATA[3]}${DATA[4]} ${SEGWAY_DATA}${REPNAMES}.bg" >> ${SEGWAY_BIN}/tdata${REPNAMES}.sh  
                echo "pigz -9 ${SEGWAY_DATA}${REPNAMES}.bg" >> ${SEGWAY_BIN}/tdata${REPNAMES}.sh  
                echo 'echo job_id $JOB_ID ending $(date)' >> ${SEGWAY_BIN}/tdata${REPNAMES}.sh
                chmod 777 ${SEGWAY_BIN}/tdata${REPNAMES}.sh
                
                # submit
                echo "qsub -pe smp 4 -V -cwd -l h_rt=24:00:00 -j y -m e -S /bin/bash -o ${SEGWAY_QOUT}TrDa-${REPNAMES}.out -N TrDa-${REPNAMES} ${SEGWAY_BIN}/tdata${REPNAMES}.sh" >>  ${SEGWAY_BIN}/2_tdata.sh
            fi

        elif [[ "${DATA[4]}" =~ ".bg" ]]; then
            if [ "$(echo ${DATA[3]} | tr ',' '\n' | wc -l | cut -f 1)" -gt 1 ]; then
                echo "[ERROR] replicates for begraphs not supported"
                exit 1
            fi
    		REPNAMES=${DATA[0]}"_"${DATA[1]}
            echo $REPNAMES
        	if [ ! -f ${SEGWAY_DATA}${REPNAMES}.bg.gz ] || [ "$OVERWRITEALL" = "TRUE" ]; then
                [ -f ${SEGWAY_QOUT}TrDa-${REPNAMES}.out ] && rm ${SEGWAY_QOUT}TrDa-${REPNAMES}.out
                echo '#!/bin/bash' > ${SEGWAY_BIN}/tdata${REPNAMES}.sh
                echo 'echo job_id $JOB_ID startdata $(date)' >> ${SEGWAY_BIN}/tdata${REPNAMES}.sh
                echo "echo convert ${REPNAMES} bigwig to bedGraph " >>  ${SEGWAY_BIN}/tdata${REPNAMES}.sh
                # wiggler
                echo "pigz -9 -c ${SEGWAY_DATA}${DATA[3]}${DATA[4]} > ${SEGWAY_DATA}${REPNAMES}.bg.gz" >> ${SEGWAY_BIN}/tdata${REPNAMES}.sh  
                echo 'echo job_id $JOB_ID ending $(date)' >> ${SEGWAY_BIN}/tdata${REPNAMES}.sh
                chmod 777 ${SEGWAY_BIN}/tdata${REPNAMES}.sh
                
                # submit
                echo "qsub -pe smp 4 -V -cwd -l h_rt=24:00:00 -j y -m e -S /bin/bash -o ${SEGWAY_QOUT}TrDa-${REPNAMES}.out -N TrDa-${REPNAMES} ${SEGWAY_BIN}/tdata${REPNAMES}.sh" >>  ${SEGWAY_BIN}/2_tdata.sh
            fi
        fi
    done < $EXPERIMENTS

	chmod 777 ${SEGWAY_BIN}/2_tdata.sh
	if [ $ARMED = "TRUE" ]; then
    	${SEGWAY_BIN}/2_tdata.sh
    fi
fi

##
## generate the separate genome archives for the tissues.
## The archieve is then annotated with the bedgraph data  
##
if [ -n "$DO_GENERATEARCHIVE" ]; then
        mkdir -p ${SEGWAY_DATA}${EXPERIMENT}.genomedata
	if [ ! -n "$CHROMSIZES" ];then
		echo "[ERROR] Chromosome sizes not given: $CHROMSIZES"
		exit 1
	fi
	## load module
	module load fabbus/segway_gbr
	echo "module load fabbus/segway_gbr" > ${SEGWAY_BIN}/3_gdata.sh
	echo "[ -f ${SEGWAY_QOUT}/GnDt-${EXPERIMENT}.out ] && rm ${SEGWAY_QOUT}/GnDt-${EXPERIMENT}.out" >> ${SEGWAY_BIN}/3_gdata.sh

    for CHR in $(cut -f 1 ${CHROMSIZES}); do
        if [ ! -f $SEQDIR/$CHR$FASTASUFFIX ]; then 
            echo "[ERROR] chr file not found: $SEQDIR/$CHR$FASTASUFFIX"
            exit 1
        fi
        echo $SEQDIR/$CHR$FASTASUFFIX
        
    	echo 'echo job_id $JOB_ID startdata $(date)' > ${SEGWAY_BIN}/gdata${EXPERIMENT}${CHR}.sh
    	# genomedata-load call
    	echo "echo '*** create genomedata archive'"  >> ${SEGWAY_BIN}/gdata${EXPERIMENT}${CHR}.sh
    	echo "genomedata-load -s $SEQDIR/$CHR$FASTASUFFIX \\" >> ${SEGWAY_BIN}/gdata${EXPERIMENT}${CHR}.sh
    	# add the -t <ID>=<FILE> sections for all tracks

	while IFS=$'\t' read -r -a DATA; do
        REPNAMES=${DATA[0]}"_"${DATA[1]}
        if [[ -n "$USE_ALL_TRACK_DATA" ]] || [[ "${DATA[1]}" == "$EXPERIMENT" ]]; then
            echo "[USE ] $REPNAMES"
#                b=$(basename $f)
#            	arrIN=(${REPNAMES//./ })
                echo "-t" $REPNAMES=${SEGWAY_DATA}${REPNAMES}.bg.gz" \\" >> ${SEGWAY_BIN}/gdata${EXPERIMENT}${CHR}.sh
            else
       	        echo "[SKIP] $REPNAMES"
    	    fi
        done < $EXPERIMENTS
        # add dinucleotide
#        echo "-t dinucleotide \\" >> ${SEGWAY_BIN}/gdata${EXPERIMENT}${CHR}.sh
        
    	echo "${SEGWAY_DATA}${EXPERIMENT}.genomedata" >> ${SEGWAY_BIN}/gdata${EXPERIMENT}${CHR}.sh
    	echo 'echo job_id $JOB_ID ending $(date)' >> ${SEGWAY_BIN}/gdata${EXPERIMENT}${CHR}.sh
    	chmod 777 ${SEGWAY_BIN}/gdata${EXPERIMENT}${CHR}.sh
        sed -i 's|//*|/|g' ${SEGWAY_BIN}/gdata${EXPERIMENT}${CHR}.sh
    	#submit
    	echo "qsub -pe smp 1 -V -cwd -b y -j y -o ${SEGWAY_QOUT}GnDt-${EXPERIMENT}-${CHR}.out -N GnDt-${EXPERIMENT}-${CHR} ${SEGWAY_BIN}/gdata${EXPERIMENT}${CHR}.sh" >> ${SEGWAY_BIN}/3_gdata.sh
    done    
	chmod 777  ${SEGWAY_BIN}/3_gdata.sh

	if [ $ARMED = "TRUE" ]; then
             ${SEGWAY_BIN}/3_gdata.sh
	fi
fi

##
## train Seqway models
##
if [ -n "$DO_TRAINSEGWAY" ]; then

    module load gi/bedtools
    echo "module load fabbus/segway_gbr" > ${SEGWAY_BIN}/4_train.sh
    echo "[ -f ${SEGWAY_QOUT}SgTrn-${EXPERIMENT}.out ] && rm ${SEGWAY_QOUT}SgTrn-${EXPERIMENT}.out" >> ${SEGWAY_BIN}/4_train.sh
    echo "[ -d ${SEGWAY_TRAIN} ] && rm -r ${SEGWAY_TRAIN}" >> ${SEGWAY_BIN}/4_train.sh

    mkdir -p $SEGWAY_DATA/tmp/
    if [ -n "$EXCLUDABLE" ]; then
        EXLCUDECOORDS="--exclude-coords=$EXCLUDABLE"
    fi
    
    OPTIONS="$MODEL_ADDPARAM --include-coords=$TRAIN_REGIONS --num-labels=${LABELS} $EXLCUDECOORDS $CLUSTEROPT --num-instances=${INSTANCES} ${CLOBBER} ${SEGWAY_TRAIN_ADDPARAM}"
    echo "echo '*** train segway'" >> ${SEGWAY_BIN}/segtrain${EXPERIMENT}.sh
    echo "segway $OPTIONS \\"> ${SEGWAY_BIN}/segtrain${EXPERIMENT}.sh 
    if [[ -n "$USE_ALL_TRACK_DATA" ]]; then
        for e in $(cut -f 1 $EXPERIMENTS | sort -u); do
            REPNAMES=$(fgrep -w "$e" $EXPERIMENTS | cut -f1,2 | sort -k1,2 | tr '\t' '_' | tr '\n' ',' | sed 's/,*$//g')
            echo "[USE ] $REPNAMES"
            echo "--track=$REPNAMES \\" >> ${SEGWAY_BIN}/segtrain${EXPERIMENT}.sh
        done
    
    else
        while IFS=$'\t' read -r -a DATA; do
   	    REPNAMES=${DATA[0]}"_"${DATA[1]}
            if [[ "${DATA[1]}" == "$TRAIN_EXPERIMENT" ]]; then
               echo "[USE ] $REPNAMES"
               echo "--track=$REPNAMES \\" >> ${SEGWAY_BIN}/segtrain${EXPERIMENT}.sh
       	    else
       	       echo "[SKIP] $REPNAMES"
            fi
         done < $EXPERIMENTS
    fi
    ## add dinucleotide
#    echo "--track=dinucleotide \\" >> ${SEGWAY_BIN}/segtrain${EXPERIMENT}.sh
    echo "train ${SEGWAY_DATA}${EXPERIMENT}.genomedata ${SEGWAY_TRAIN}" >> ${SEGWAY_BIN}/segtrain${EXPERIMENT}.sh

    chmod 777 ${SEGWAY_BIN}/segtrain${EXPERIMENT}.sh
    #echo "qsub -l mem_requested=16G -V -cwd -b y -j y -o ${SEGWAY_QOUT}SgTrn0${EXPERIMENT}.out -N SgTrn-${EXPERIMENT} ${SEGWAY_BIN}//segtrain${EXPERIMENT}.sh"  >> ${SEGWAY_BIN}/4_train.sh
    echo "${SEGWAY_BIN}/segtrain${EXPERIMENT}.sh" >> ${SEGWAY_BIN}/4_train.sh
    # make sure there is no douple // in any path as segway doesn't like that
    sed -i 's|//*|/|g' ${SEGWAY_BIN}/segtrain${EXPERIMENT}.sh
    
    chmod 777  ${SEGWAY_BIN}/4_train.sh

    if [ $ARMED = "TRUE" ]; then
        ${SEGWAY_BIN}/4_train.sh
    fi
fi

##
## predict using a trained Seqway model
##
if [ -n "$DO_PREDICTSEGWAY" ]; then
    echo "module load fabbus/segway_gbr" > ${SEGWAY_BIN}/5_predict.sh
    echo "[ -f ${SEGWAY_QOUT}SgPrd-${EXPERIMENT}.out ] && rm ${SEGWAY_QOUT}SgPrd-${EXPERIMENT}.out" >> ${SEGWAY_BIN}/5_predict.sh
    echo "[ -d ${SEGWAY_PREDICT} ] && rm -r ${SEGWAY_PREDICT}" >> ${SEGWAY_BIN}/5_predict.sh

    echo 'echo job_id $JOB_ID startdata $(date)' > ${SEGWAY_BIN}/segpredict${EXPERIMENT}${PREDICTION}.sh

    echo "echo '*** predict segmentation'" >>  ${SEGWAY_BIN}/segpredict${EXPERIMENT}${PREDICTION}.sh
    if [ -n "$EXCLUDABLE" ]; then
        EXLCUDECOORDS="--exclude-coords=$EXCLUDABLE"
    fi
    echo "segway $MODEL_ADDPARAM --num-labels=${LABELS} $CLUSTEROPT $EXLCUDECOORDS ${CLOBBER} \\">> ${SEGWAY_BIN}/segpredict${EXPERIMENT}${PREDICTION}.sh
    # add the --track <ID> sections
    if [[ -n "$USE_ALL_TRACK_DATA" ]]; then
        for e in $(cut -f 1 $EXPERIMENTS | sort -u); do
            REPNAMES=$(fgrep -w "$e" $EXPERIMENTS | cut -f1,2 | sort -k1,2 | tr '\t' '_' | tr '\n' ',' | sed 's/,*$//g')
            echo "[USE ] $REPNAMES"
            echo "--track=$REPNAMES \\" >> ${SEGWAY_BIN}/segtrain${EXPERIMENT}.sh
        done

    else
    while IFS=$'\t' read -r -a DATA; do
	REPNAMES=${DATA[0]}"_"${DATA[1]}
        if [[ "${DATA[1]}" == "$PREDICTION" ]]; then
            echo "[USE ] $REPNAMES"
            echo "--track=$REPNAMES \\" >> ${SEGWAY_BIN}/segpredict${EXPERIMENT}${PREDICTION}.sh
    	else
    	   echo "[SKIP] $REPNAMES"
	fi
    done < $EXPERIMENTS
    fi
    ## add dinucleotide
#    echo "--track=dinucleotide \\" >> ${SEGWAY_BIN}/segpredict${EXPERIMENT}${PREDICTION}.sh
    echo "identify ${SEGWAY_DATA}${EXPERIMENT}.genomedata ${SEGWAY_TRAIN} ${SEGWAY_PREDICT}" >> ${SEGWAY_BIN}/segpredict${EXPERIMENT}${PREDICTION}.sh
    echo 'echo job_id $JOB_ID ending $(date)' >> ${SEGWAY_BIN}/segpredict${EXPERIMENT}${PREDICTION}.sh
    chmod 777 ${SEGWAY_BIN}/segpredict${EXPERIMENT}${PREDICTION}.sh
    # submit
#       echo "qsub -l mem_requested=16G -V -cwd -b y -j y -o ${SEGWAY_QOUT}SgPrd-${EXPERIMENT}.out -N SgPrd-${EXPERIMENT} ${SEGWAY_BIN}/segpredict${EXPERIMENT}${PREDICTION}.sh"   
    echo "${SEGWAY_BIN}/segpredict${EXPERIMENT}${PREDICTION}.sh" >> ${SEGWAY_BIN}/5_predict.sh
    # make sure there is no douple // in any path as segway doesn't like that
    sed -i 's|//*|/|g' ${SEGWAY_BIN}/segpredict${EXPERIMENT}${PREDICTION}.sh

    chmod 777 ${SEGWAY_BIN}/5_predict.sh
    
    if [ $ARMED = "TRUE" ]; then
        ${SEGWAY_BIN}/5_predict.sh
    fi
fi

##
## Evaluate prediction
##
if [ -n "$DO_EVALUATE" ]; then

    echo "#!/bin/bash -e" > ${SEGWAY_BIN}/6_evaluate.sh
    echo "unset module"
    echo "module load fabbus/segway_gbr gi/ucsc_utils/283 gi/bedtools" >> ${SEGWAY_BIN}/6_evaluate.sh

    for rf in ${SEGWAY_PREDICT}/segway.[0-9].bed.gz; do 
        COUNTER=$(echo $rf | sed 's/.*segway.\([0-9]\).bed.gz/\1/g')
        echo $COUNTER

    echo "#!/bin/bash -e" > ${SEGWAY_BIN}/segeval${EXPERIMENT}$COUNTER.sh
    echo "unset module" >> ${SEGWAY_BIN}/segeval${EXPERIMENT}$COUNTER.sh
    echo 'echo job_id $JOB_ID startdata $(date)' >> ${SEGWAY_BIN}/segeval${EXPERIMENT}$COUNTER.sh
    #preprocess file
    if [ -n $OVERWRITEALL ] || [ ! -f ${SEGWAY_PREDICT}/segway.$COUNTER.bed.gz.pkl.gz ]; then
        echo "echo '*** preprocess'" >> ${SEGWAY_BIN}/segeval${EXPERIMENT}$COUNTER.sh     
        echo "segtools-preprocess ${CLOBBER} ${SEGWAY_PREDICT}/segway.$COUNTER.bed.gz" >> ${SEGWAY_BIN}/segeval${EXPERIMENT}$COUNTER.sh 
    fi
    
    echo "echo '*** length disttribution analysis'" >> ${SEGWAY_BIN}/segeval${EXPERIMENT}$COUNTER.sh
    echo "segtools-length-distribution ${SEGWAY_PREDICT}/segway.$COUNTER.bed.gz.pkl.gz --outdir=${SEGWAY_RESULT}/length-dist$COUNTER/ ${CLOBBER}" >> ${SEGWAY_BIN}/segeval${EXPERIMENT}$COUNTER.sh 
    
    echo "echo '*** gene aggregation analysis'" >> ${SEGWAY_BIN}/segeval${EXPERIMENT}$COUNTER.sh
    echo "segtools-aggregation ${SEGWAY_PREDICT}/segway.$COUNTER.bed.gz.pkl.gz ${ANNOTATION} --normalize --mode=gene --outdir=${SEGWAY_RESULT}/gencode-agg$COUNTER/ ${CLOBBER}" >> ${SEGWAY_BIN}/segeval${EXPERIMENT}$COUNTER.sh 
    
    echo "echo '*** gmtk parameter generation'" >> ${SEGWAY_BIN}/segeval${EXPERIMENT}$COUNTER.sh
    echo "segtools-gmtk-parameters ${SEGWAY_TRAIN}/params/params.params --outdir=${SEGWAY_RESULT}/gtmk-param$COUNTER/ ${CLOBBER}" >> ${SEGWAY_BIN}/segeval${EXPERIMENT}$COUNTER.sh 
    
    echo "echo '*** html report generation'" >> ${SEGWAY_BIN}/segeval${EXPERIMENT}$COUNTER.sh
    echo "cd ${SEGWAY_RESULT}/" >> ${SEGWAY_BIN}/segeval${EXPERIMENT}$COUNTER.sh
    echo "segtools-html-report -L ${SEGWAY_PREDICT}/segway.$COUNTER.layered.bed.gz --results-dir=${SEGWAY_RESULT}/ -o segtools$COUNTER.html ${SEGWAY_PREDICT}/segway.$COUNTER.bed.gz.pkl.gz ${CLOBBER}" >> ${SEGWAY_BIN}/segeval${EXPERIMENT}$COUNTER.sh 
    echo "sed 's|${SEGWAY_RESULT}/||g' segtools$COUNTER.html > segtools$COUNTER.html2" >> ${SEGWAY_BIN}/segeval${EXPERIMENT}$COUNTER.sh 
    echo "mv segtools$COUNTER.html2 segtools$COUNTER.html" >> ${SEGWAY_BIN}/segeval${EXPERIMENT}$COUNTER.sh 
    echo "cd $(pwd)" >> ${SEGWAY_BIN}/segeval${EXPERIMENT}$COUNTER.sh
    echo "for LABEL in \$(seq 0 $(( $LABELS - 1 )) ); do zcat ${SEGWAY_PREDICT}/segway.$COUNTER.bed.gz | tail -n+2 | awk -v label=\$LABEL '{if (\$4 == label){OFS=\"\t\"; print \$1,\$2,\$3,\$4}}' | bedtools sort >  ${SEGWAY_PREDICT}/segway.$COUNTER.\$LABEL.bed" >> ${SEGWAY_BIN}/segeval${EXPERIMENT}$COUNTER.sh 
    echo "bedToBigBed -type=bed4 ${SEGWAY_PREDICT}/segway.$COUNTER.\$LABEL.bed $CHROMSIZES ${SEGWAY_PREDICT}/$EXPERIMENT.$COUNTER.\$LABEL.bb" >> ${SEGWAY_BIN}/segeval${EXPERIMENT}$COUNTER.sh 
    echo "rm  ${SEGWAY_PREDICT}/segway.$COUNTER.\$LABEL.bed" >> ${SEGWAY_BIN}/segeval${EXPERIMENT}$COUNTER.sh 
    echo "done" >> ${SEGWAY_BIN}/segeval${EXPERIMENT}$COUNTER.sh 
    # make sure there is no douple // in any path as segway doesn't like that
    sed -i 's|//*|/|g' ${SEGWAY_BIN}/segeval${EXPERIMENT}$COUNTER.sh
 
    chmod 777 ${SEGWAY_BIN}/segeval${EXPERIMENT}$COUNTER.sh
    
    echo "qsub -l mem_requested=16G -S /bin/bash -V -cwd -b y -j y -o ${SEGWAY_QOUT}SgEva-${EXPERIMENT}$COUNTER.out -N SgEva-${EXPERIMENT}$COUNTER ${SEGWAY_BIN}/segeval${EXPERIMENT}$COUNTER.sh" >> ${SEGWAY_BIN}/6_evaluate.sh
    
    done

    chmod 777 ${SEGWAY_BIN}/6_evaluate.sh
    if [ $ARMED = "TRUE" ]; then
        ${SEGWAY_BIN}/6_evaluate.sh
    fi
fi


