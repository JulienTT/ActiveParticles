#!/bin/bash
# this file generates a pdf of a histogram

if [[ $# -ne 7 ]] ; then
    echo 'usage: file time dt output rho v Dr'
    exit 0
fi

file=$1
time=$2
dt=$3
fileout=$4
rho=$5
v=$6
Dr=$7

./MakeAverageHisto.sh $file $time $dt $fileout

rhomax=$(awk '$2>rhomax {rhomax=$2} END {print rhomax}' $fileout)

echo "Peak of distribution: $rhomax"

xmax=$(awk -v rhomax=$rhomax '$2>(rhomax/100.) {xmax=$1} END{print xmax}' $fileout)

echo "largest relevant rho is $xmax"

./plot-histo-title.py $fileout $xmax $rho $v $Dr
