#!/bin/bash

#This code compute the pressure starting from a time T up to the end
#of the system, for T between 'first' and the final time, every
#'increment'

if [[ $# -ne 3 ]] ;
then
    echo 'usage: ./ComputePressure.sh file first-time increment'
    exit 1;
fi

file=$1
first=$2
increment=$3

Last=$(awk 'NF>1 {time=$1} END{print time-'"$increment"'}' $file)

rm  $file-av;
for T in `seq $first $increment $Last`;
do
    #echo "Time $T"
    awk -v time=$T '$1>time {p1+=$2;p2+=$3;rho+=$4;n++} END{print time,p1/n,p2/n,rho/n}' $file  >> $file-av
done
