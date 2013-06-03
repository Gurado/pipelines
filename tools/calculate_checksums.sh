#!/bin/sh -e


USAGEMSG="usage: $(basename $0) FOLDER 

calculates checksum for all files in FOLDER

Author: Fabian Buske

* FOLDER - target folder containing relevant files
* -c FILE - verifies the files in the folder given the provided checksum FILE
"

DIR=$(dirname $0)
VERSION="0.0.1"

[ $# -lt 1 ] && echo "$USAGEMSG" >&2 && exit 1

CHECKSUMS=""

while getopts "c:" opt;
do
    case ${opt} in
    c) CHECKSUMS="$OPTARG";;
    \?) print >&2 "$0: error - unrecognized option $1"
        exit 1;;
    esac
done

shift $(($OPTIND-1))
FOLDER=$1

if [ -n "$CHECKSUMS" ]; then
    echo "checking files"
    while read MD5 FILE
    do
        MD5new=$(md5 -r $FILE | awk '{print $1}')
        if [[ "$MD5" != "$MD5new" ]]; then
            echo "[ERROR] $FILE md5 mismatch: $MD5new should have been $MD5"
        else
            echo "[OK] $FILE verified"
        fi
        
    done < $CHECKSUMS

else
    for f in $(find $FOLDER -type f | grep -v ".DS_Store" | grep -v ".md5"); do
        MD5=$(md5 -r $f)
        echo $MD5
    done
fi