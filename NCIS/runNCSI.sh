#!/bin/sh -e

USAGEMSG="usage: $(basename $0) -g gagriChipSeqDir -w workingDir -N jobname -p 'gsub parameters' -d -f -v CHIP INPUT

Starts the Normalization for ChIP-Seq (NCBI) pipeline.

Author: Fabian Buske

Requirements (in PATH environment or specified):
        R with ShortRead library available

* CHIP - the chip data (bam expected)
* CONTROL - the input/control data (bam expected)
* -g data source directory - directory on gagri that points to the ChIP-seq data folders
* -w working directory - directory to put all the data, scripts, results
* -n NCIS source - source of the NCIS R script
* -N job name - give this job a name of your choicea
* -p parameters - cluster job parameters to use when submitting via qsub
* -f force - overwrite existing results
* -d dry - write scripts but dont trigger R job
* -v - print progress information (verbose).
"


DIR=`dirname $0`
VERSION="0.0.1"

[ $# -lt 2 ] && echo "$USAGEMSG" >&2 && exit 1

CHIP=""
CONTROL=""
GAGRI=/Cancer-Epigenetics/Data/ClarkLab/Seq/ChIP-Seq/hg19/
NCIS=http://pages.cs.wisc.edu/~kliang/NCIS/NCIS.R
WORKINGDIR=$PWD
DRYRUN="FALSE"
FORCE="FALSE"
NCORES=1
VERBOSE="--quiet"
JOBNAME=""
JOBPARAMS="-l h_vmem=25G,virtual_free=15G"

while getopts "g:w:n:N:p:dfv" opt;
do
        case ${opt} in
        g) GAGRI="$OPTARG";;
        n) NCIS="$OPTARG";;
	N) JOBNAME="$OPTARG";;
	p) JOBPARAMS="$OPTARG";;
        w) WORKINGDIR="$OPTARG";;
	d) DRYRUN="TRUE";;
	f) FORCE="TRUE";;
        v) VERBOSE="--verbose";;
        \?) print >&2 "$0: error - unrecognized option $1"
                exit 1;;
        esac
done

shift $(($OPTIND-1))
CHIP=$1
CONTROL=$2

if [ -z ${JOBNAME} ]; then
   JOBNAME="NCIS-${CHIP}-${CONTROL}"
fi

if [ ! -d ${WORKINGDIR} ]; then
   mkdir -p ${OUTPUT}
else
   echo "[WARN] working directory already exists, content will be overwritten" 
fi

DATA=${WORKINGDIR}/data/
RESULT=${WORKINGDIR}/results/
BIN=${WORKINGDIR}/bin/
LOG=${WORKINGDIR}/log/

mkdir -p ${DATA} ${BIN} ${RESULT} ${LOG}

##
## Check if results exists already or existing results are to be overwritten
##

if [ ${FORCE} = "FALSE" ] && [ -f ${RESULT}${CHIP}-${CONTROL}.txt ]; then
   echo "[NOTE] Results already exist: ${RESULT}${CHIP}_${CONTROL}.RData" >> ${LOG}/${JOBNAME}.log
   [ ${VERBOSE} = "--verbose" ] && tail -n 1 ${LOG}/${JOBNAME}.log
   exit 0
fi

##
## make log
##

echo "ChIPseq QC : v${VERSION}" > ${LOG}/${JOBNAME}.log
echo "Jobname    : $JOBNAME"    >> ${LOG}/${JOBNAME}.log
echo "chip       : $CHIP"       >> ${LOG}/${JOBNAME}.log
echo "control    : $CONTROL"    >> ${LOG}/${JOBNAME}.log
echo "gagri      : $GAGRI"      >> ${LOG}/${JOBNAME}.log

echo "working dir: $WORKINGDIR" >> ${LOG}/${JOBNAME}.log
echo "data       : $DATA"       >> ${LOG}/${JOBNAME}.log
echo "scripts    : $BIN"        >> ${LOG}/${JOBNAME}.log
echo "result     : $RESULT"     >> ${LOG}/${JOBNAME}.log
echo "logs       : $LOG"        >> ${LOG}/${JOBNAME}.log
echo "dry-run    : $DRYRUN"     >> ${LOG}/${JOBNAME}.log
echo "force      : $FORCE"      >> ${LOG}/${JOBNAME}.log

if [ ${VERBOSE} = "--verbose" ]; then
   cat ${LOG}/${JOBNAME}.log
fi

##
## Check if data already existing
##

if [ ${FORCE} = "TRUE" ] || [ ! -f ${DATA}${CHIP}.bam ]; then
   echo "** get ${CHIP} data from gagri" >> ${LOG}/${JOBNAME}.log
   if [ -f  ~/.smbclient ]; then
      smbclient \\\\gagri\\GRIW -A ~/.smbclient -c "cd ${GAGRI}/${CHIP}; get ${CHIP}.bam" && mv ${CHIP}.bam ${DATA}
   else
      smbclient \\\\gagri\\GRIW -U `whoami` -c "cd ${GAGRI}/${CHIP}; get ${CHIP}.bam" && mv ${CHIP}.bam ${DATA}
   fi
fi

if [ ${FORCE} = "TRUE" ] || [ ! -f ${DATA}${CONTROL}.bam ]; then
   echo "** get ${CONTROL} data from gagri" >> ${LOG}/${JOBNAME}.log
   if [ -f  ~/.smbclient ]; then
      smbclient \\\\gagri\\GRIW -A ~/.smbclient -c "cd ${GAGRI}/${CONTROL}; get ${CONTROL}.bam" && mv ${CONTROL}.bam ${DATA}
   else
      smbclient \\\\gagri\\GRIW -U `whoami` -c "cd ${GAGRI}/${CONTROL}; get ${CONTROL}.bam" && mv ${CONTROL}.bam ${DATA}
   fi
fi

##
## write R script
##

echo "** write R script" >> ${LOG}/${JOBNAME}.log
echo "source('http://pages.cs.wisc.edu/~kliang/NCIS/NCIS.R')" > ${BIN}/${CHIP}-${CONTROL}.R
echo "library(ShortRead)" >> ${BIN}/${CHIP}-${CONTROL}.R
echo "input <- readAligned(dirPath = '${DATA}', pattern ='${CONTROL}', type = 'BAM')" >> ${BIN}/${CHIP}-${CONTROL}.R
echo "chip <- readAligned(dirPath = '${DATA}', pattern = '${CHIP}', type = 'BAM')" >> ${BIN}/${CHIP}-${CONTROL}.R
echo "res <- NCIS(chip, input, data.type='AlignedRead')" >> ${BIN}/${CHIP}-${CONTROL}.R
echo "write.table(res, file='${RESULT}${CHIP}-${CONTROL}.txt', sep='\t')" >> ${BIN}/${CHIP}-${CONTROL}.R

##
## submit script to cluster
##

if [ ${DRYRUN} = "TRUE" ]; then

   echo "qsub -pe smp $NCORES -m e -o ${LOG} ${JOBPARAMS} -e ${LOG} -N ${JOBNAME} -M `whoami`@garvan.unsw.edu.au -wd ${WORKINGDIR} -b y /share/ClusterShare/software/contrib/Cancer-Epigenetics/tools/bin/Rscript --quiet '${BIN}${CHIP}-${CONTROL}.R'" >> ${LOG}/${JOBNAME}.log
   tail -n 1 ${LOG}/${JOBNAME}.log

else
   echo "** submit job" >> ${LOG}/${JOBNAME}.log
   qsub -pe smp $NCORES -m e -o ${LOG} ${JOBPARAMS}  -e ${LOG} -N ${JOBNAME} -M `whoami`@garvan.unsw.edu.au -wd ${WORKINGDIR} -b y /share/ClusterShare/software/contrib/Cancer-Epigenetics/tools/bin/Rscript --quiet "${BIN}${CHIP}-${CONTROL}.R"
fi

