#!/bin/sh -e

USAGEMSG="usage: $(basename $0) CONFIGFILE

Rezips a file with pigs -11 option. Using the cluster with 64 CPUs.

Author: Fabian Buske

Requirements (modules):
        gi/pigz/2.3
"
module load gi/pigz/2.3

[ $# -lt 1 ] && echo "$USAGEMSG" >&2 && exit 1

while getopts "v" opt;
do
        case ${opt} in
        v) VERBOSE="--verbose";;
        \?) print >&2 "$0: error - unrecognized option $1"
            exit 1;;
        esac
done


for arg
do
	F="${arg[$i]}"

	if [[ ${F##*.} != "gz" ]]; then 
		echo "file not zipped? $F skipped" 
		continue
	fi

	echo "#!/bin/sh" > $F.pigztmp.sh
	echo "ls -la $F" >> $F.pigztmp.sh
	echo "unpigz $F" >> $F.pigztmp.sh
	echo "pigz -11 ${F/.gz/}" >> $F.pigztmp.sh
	echo "ls -la $F" >> $F.pigztmp.sh
	chmod 777 $F.pigztmp.sh

	qsub -V -S /bin/bash -j y -o $F.qout -cwd -pe smp 64 -l h_vmem=40G -N pigz_$F -l h_rt=4:00:00 $F.pigztmp.sh
done
