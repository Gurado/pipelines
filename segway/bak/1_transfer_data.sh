#!/bin/sh -e

DIR=`dirname $0`
source ${DIR}/0_config.sh


##
## transfer the data (tracks) from gagri
##

if [ -n "$DO_TRANSFERDATA" ]; then

	# get all files
	for F in ${FILES}; do
		echo "get datafile ${F} from gagri"
		[ ! -f ${F}.bam ] && smbclient \\\\gagri\\GRIW -A ~/.smbclient -c "cd ${FILES_SOURCE}/${F}; get ${F}.bam" && mv ${F}.bam ${SEGWAY_DATA}
		[ ! -f ${F}.bam.bai ] && smbclient \\\\gagri\\GRIW -A ~/.smbclient -c "cd ${FILES_SOURCE}/${F}; get ${F}.bam.bai" && mv ${F}.bam.bai ${SEGWAY_DATA}

	done

fi
