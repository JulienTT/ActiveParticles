#!/bin/bash

if [[ $# -ne 3 ]] ; then
    echo 'usage: file xmin xmax'
    exit 0
fi

file=$1
x1=$2
x2=$3

awk -v x1=$x1 -v x2=$x2 '(($1>x1)&&($1<x2)&&($2>max)) {max=$2;xmax=$1} END{print xmax,max}' $file

