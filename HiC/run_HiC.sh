#!/bin/sh -e

USAGEMSG="usage: $(basename $0) CONFIGFILE

Starts the segway pipeline given the config file

Author: Fabian Buske

Requirements (modules):
	module fabbus/segway/1.1.0 
	BEDtools (genomeCoverageBed, fetchChromSizes) in path

* CONFIGFILE - the config file describing the joba
* -1 step 1 - collect the bam data from gagri
* -2 step 2 - run FastQC on data
* -3 step 3 - map data using bowtie
* -4 step 4 - 
* -5 step 5 - 
* -6 step 6 - 
* -f force - overwrite existing results
* -a armed - trigger the jobs after writing the scripts
* -v - print progress information (verbose).
"


PIPELINEDIR=`dirname $0`
VERSION="0.0.1"

DO_TRANSFERDATA=""
DO_FASTQC=""
DO_MAPPING=""
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
	2) DO_FASTQC="TRUE";;
	3) DO_MAPPING="TRUE";;
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
HIC_DATA=${TARGETDIR}data/
HIC_BIN=${TARGETDIR}bin/
HIC_QOUT=${TARGETDIR}qout/
HIC_RESULT=${TARGETDIR}result/

# some housekeeping
mkdir -p $HIC_DATA $HIC_BIN $HIC_QOUT $HIC_RESULT

##
## collect the data (tracks) from gagri
##
if [ -n "$DO_TRANSFERDATA" ];  then
	# get all files
	echo "echo  'collect data tracks from gagri'" > ${HIC_BIN}1_cdata.sh
	FILEARR=$(echo ${FASTQ} | tr " ," "\n")
	for F in ${FILEARR}; do
		FN=${F##*/}
                FB=${FN%.*}
	        echo "echo 'datafile ${F}'" >> ${HIC_BIN}1_cdata.sh
        	echo "smbclient \\\\\\\\gagri\\\\GRIW -A ~/.smbclient -c 'prompt; recurse; cd ${RAW_FILES_SOURCE}${FB}; mget ${F}*.gz' && mv ${F}*.gz ${HIC_DATA}" >> ${HIC_BIN}1_cdata.sh
	done

	chmod 777 ${HIC_BIN}1_cdata.sh
	if [ $ARMED = "TRUE" ]; then
	        ${HIC_BIN}1_cdata.sh
	fi
fi

##      
## transform the bam data into bedgraph
##
if [ -n "$DO_FASTQC" ]; then
	echo "echo  'run fastqc'" > ${HIC_BIN}fastqc.sh
	echo "mkdir -p ${HIC_RESULT}fastqc" >> ${HIC_BIN}fastqc.sh
	echo "fastqc -f fastq -t 6 -o ${HIC_RESULT}fastqc ${HIC_DATA}*.gz " >> ${HIC_BIN}fastqc.sh

	chmod 777 ${HIC_BIN}fastqc.sh	
	echo "qsub -V -cwd -l h_rt=01:00:00 -pe smp 6 -j y -m e -M `whoami`@garvan.unsw.edu.au -S /bin/bash -o ${HIC_QOUT}FastQC.out -N FASTQC ${HIC_BIN}fastqc.sh" >  ${HIC_BIN}2_fastqc.sh
	chmod 777 ${HIC_BIN}2_fastqc.sh

        if [ $ARMED = "TRUE" ]; then
		${HIC_BIN}2_fastqc.sh
        fi

fi
