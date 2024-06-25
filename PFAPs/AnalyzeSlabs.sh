#!/bin/bash
rm  DiagPhase-Dr-2-slab.txt
for v in 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.01 1.02 1.03 1.04 1.05 1.06 1.07 1.08 1.09 1.1 1.11;
do
    ./MakeProfile.sh Data-slab-v-"$v"-Dr-2-Lx-120-Ly-40-rhoG-0.24-rhoL-5.32-Tf-4000-pos 150 1 120 40
    
    var=$(awk '$1<20||($1>100&&$1<118) {rhom+=$2;min++} $1>45&&$1<75 {rhoM+=$2;max++}END{print rhom/min*.37*.37,rhoM/max*.37*.37}' Data-slab-v-"$v"-Dr-2-Lx-120-Ly-40-rhoG-0.24-rhoL-5.32-Tf-4000-pos-Av)
    set -- $var
    
    echo $v $1 $2  >> DiagPhase-Dr-2-slab.txt
done;
