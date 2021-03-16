#!/bin/bash

t1=100;

echo "# dt rho P sigma" > StatPressureWCA-rho-0.45
for dt in 1e-6 2.5e-6 5e-6 1e-5 2e-5 4e-5 8e-5 8e-4 1e-3 2e-3 4e-3 8e-3 1e-2 2e-2 4e-2 8e-2 1e-1 2e-1;
do
    ./ComputePressure.sh test-pressure-WCA-rho-0.45-Lx-50-Ly-22-dt-"$dt"-pressure 50 1;
    var=`awk '$1~/'"$t1"'/' test-pressure-WCA-rho-0.45-Lx-50-Ly-22-dt-"$dt"-pressure-av`
    set -- $var
    p1=$2
    p2=$3
    realrho=$4
    p=$(echo $p1 $p2 | awk '{print .5*($1+$2)}')
    sigma=$( echo $p1 $p2 | awk '{print ($1>$2)?($1-$2):($2-$1)}' )
    echo "$dt $realrho $p $sigma" >>  StatPressureWCA-rho-0.45
done;

echo "# dt rho P sigma" > StatPressureWCA-rho-0.80
for dt in 2e-5 4e-5 8e-5 2e-4 4e-4 8e-4 1e-3;
do
    ./ComputePressure.sh test-pressure-WCA-rho-0.8-Lx-50-Ly-22-dt-"$dt"-pressure 50 1;
    var=`awk '$1~/'"$t1"'/' test-pressure-WCA-rho-0.8-Lx-50-Ly-22-dt-"$dt"-pressure-av`
    set -- $var
    p1=$2
    p2=$3
    realrho=$4
    p=$(echo $p1 $p2 | awk '{print .5*($1+$2)}')
    sigma=$( echo $p1 $p2 | awk '{print ($1>$2)?($1-$2):($2-$1)}' )
    echo "$dt $realrho $p $sigma" >>  StatPressureWCA-rho-0.80
done;

echo "# dt rho P sigma" > StatPressureWCA-rho-0.80-Lx-90
for dt in 5e-6 1e-5 2e-5 4e-5 8e-5 1e-4 2e-4 4e-4 8e-4 1e-3 1e-1;
do
    ./ComputePressure.sh test-pressure-WCA-rho-0.8-Lx-90-Ly-22-dt-"$dt"-pressure 50 1;
    var=`awk '$1~/'"$t1"'/' test-pressure-WCA-rho-0.8-Lx-90-Ly-22-dt-"$dt"-pressure-av`
    set -- $var
    p1=$2
    p2=$3
    realrho=$4
    p=$(echo $p1 $p2 | awk '{print .5*($1+$2)}')
    sigma=$( echo $p1 $p2 | awk '{print ($1>$2)?($1-$2):($2-$1)}' )
    echo "$dt $realrho $p $sigma" >>  StatPressureWCA-rho-0.80
done;
