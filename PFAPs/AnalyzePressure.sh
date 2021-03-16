#!/bin/bash

t1=100;
t2=800;

echo "#v rho Dr Dt p_theo p erreur" > StatPressureNI
for v in 1 2 3; do
    for rho in .5 1 1.5; do
	for Dr in 1 2 3; do
	    for Dt in 1 2 3; do
		if [ -f test-pressure-v-"$v"-rho-"$rho"-Dr-"$Dr"-Dt-"$Dt"-pressure ]; then
		    echo "Analyse"
		    ./ComputePressure.sh test-pressure-v-"$v"-rho-"$rho"-Dr-"$Dr"-Dt-"$Dt"-pressure 50 1;
		    var=`awk '$1~/'"$t1"'/' test-pressure-v-"$v"-rho-"$rho"-Dr-"$Dr"-Dt-"$Dt"-pressure-av`
		    set -- $var
		    p1=$2
		    p2=$3
		    realrho=$4
		    p=$(echo $p1 $p2 | awk '{print .5*($1+$2)}')
		    p_theo=$(echo "" | awk '{print '"$realrho"'*('"$Dt"'+'"$v"'*'"$v"'/2./'"$Dr"')}')
		    sigma=$( echo $p1 $p2 | awk '{print ($1>$2)?($1-$2):($2-$1)}' )

		    echo "$v $realrho $Dr $Dt $p_theo $p $sigma" >>  StatPressureNI
		fi
	    done;
	done;
    done;
done;

	 
