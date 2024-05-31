#!/bin/bash

t1=1000;
t2=800;
v=10;
Dr=1;
Dt=0;

for dt in 1e-3 5e-3 1e-2 5e-2 1e-1; do
    ./ComputePressure.sh DataPressure/test-pressure-non-interacting-rho-0.5-Lx-100-Ly-50-dt-"$dt"-v-"$v"-Dr-"$Dr"-Dt-"$Dt"-mu-1-pressure 100 100;
    
    var=`awk '$1~/'"$t1"'/' DataPressure/test-pressure-non-interacting-rho-0.5-Lx-100-Ly-50-dt-"$dt"-v-10-Dr-1-Dt-0-mu-1-pressure-av`
    set -- $var
    p1=$2
    p2=$3
    realrho=$4
    p=$(echo $p1 $p2 | awk '{print .5*($1+$2)}')
    p_theo=$(echo "" | awk '{print '"$realrho"'*('"$Dt"'+'"$v"'*'"$v"'/2./'"$Dr"')}')
    sigma=$( echo $p1 $p2 | awk '{print ($1>$2)?($1-$2):($2-$1)}' )
    
    echo "$dt $v $realrho $Dr $Dt $p_theo $p $sigma"
done > StatPressureNonInteracting;

    # Old results with translational diffusivity

# echo "#v rho Dr Dt p_theo p erreur" > StatPressureNI
# for v in 1 2 3; do
#     for rho in .5 1 1.5; do
# 	for Dr in 1 2 3; do
# 	    for Dt in 1 2 3; do
# 		if [ -f test-pressure-v-"$v"-rho-"$rho"-Dr-"$Dr"-Dt-"$Dt"-pressure ]; then
# 		    echo "Analyse"
# 		    ./ComputePressure.sh test-pressure-v-"$v"-rho-"$rho"-Dr-"$Dr"-Dt-"$Dt"-pressure 50 1;
# 		    var=`awk '$1~/'"$t1"'/' test-pressure-v-"$v"-rho-"$rho"-Dr-"$Dr"-Dt-"$Dt"-pressure-av`
# 		    set -- $var
# 		    p1=$2
# 		    p2=$3
# 		    realrho=$4
# 		    p=$(echo $p1 $p2 | awk '{print .5*($1+$2)}')
# 		    p_theo=$(echo "" | awk '{print '"$realrho"'*('"$Dt"'+'"$v"'*'"$v"'/2./'"$Dr"')}')
# 		    sigma=$( echo $p1 $p2 | awk '{print ($1>$2)?($1-$2):($2-$1)}' )

# 		    echo "$v $realrho $Dr $Dt $p_theo $p $sigma" >>  StatPressureNI
# 		fi
# 	    done;
# 	done;
#     done;
# done;

	 
