#!/bin/bash

Dr=$1

echo "#rhomax1 p(rhomax1) rhomax2 p(rhomax2)"> DiagPhase-Dr-"$Dr".txt

for i in *Dr-"$Dr"*histo-plot
do
    echo $i | sed -e 's/test-//; s/histo-plot//; s/^\(.\)/\#\1/';
    v=$(echo $i | sed -e 's/[^v]*v-//; s/-.*//')
    awk -v vitesse=$v -f FindExtrema.awk < $i;
done >> DiagPhase-Dr-"$Dr".txt

