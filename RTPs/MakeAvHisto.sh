#!/bin/bash


## $1 is the file
## $2 is the starting time
## $3 is Lx
## $4 is Ly

file=$1
output=$1"-av"
time=$2
Lx=$3
Ly=$4

if [[ $# -ne 4 ]] ; then
    echo 'usage: file time Lx Ly'
    exit 0
fi

gawk -v time=$2 -v Lx=$3 -v Ly=$4 -f MakeavHisto.awk < $file > $output

