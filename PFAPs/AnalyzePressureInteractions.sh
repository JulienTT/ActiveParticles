#!/bin/bash

t1=1000;
v=10
Dr=1

echo "#realrho rho p pNI sigma" >  StatPressureWCA-v-10-Dr-1-vs-rho;
for rho in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8;
do
    ./ComputePressure.sh DataPressure/test-pressure-WCA-rho-"$rho"-Lx-100-Ly-50-dt-5e-5-v-10-Dr-1-Dt-0-mu-1-pressure 100 100;
    var=`awk '$1~/^'"$t1"'$/' DataPressure/test-pressure-WCA-rho-"$rho"-Lx-100-Ly-50-dt-5e-5-v-10-Dr-1-Dt-0-mu-1-pressure-av`
    set -- $var
    p1=$2
    p2=$3
    realrho=$4
    p=$(echo $p1 $p2 | awk '{print .5*($1+$2)}')
    pNI=$(echo $v $Dr | awk -v v0=$v -v Dr=$Dr '{print '"$realrho"'*.5*v0*v0/Dr}')
    sigma=$( echo $p1 $p2 | awk '{print ($1>$2)?($1-$2):($2-$1)}' )
    echo "$realrho $rho $p $pNI $sigma" 
done >> StatPressureWCA-v-10-Dr-1-vs-rho;

rho=0.8
t1=200;
echo "#dt realrho rho p pNI sigma" >   StatPressureWCA-v-10-Dr-1-rho-0.8-vs-dt;
for dt in 2.5e-6 5e-6 1e-5 2.5e-5 5e-5 1e-4 2.5e-4 5e-4; 
do
    ./ComputePressure.sh DataPressure/test-pressure-WCA-rho-"$rho"-Lx-100-Ly-50-dt-"$dt"-v-10-Dr-1-Dt-0-mu-1-pressure 100 100;
    var=`awk '$1~/^'"$t1"'$/' DataPressure/test-pressure-WCA-rho-"$rho"-Lx-100-Ly-50-dt-"$dt"-v-10-Dr-1-Dt-0-mu-1-pressure-av`
    set -- $var
    p1=$2
    p2=$3
    realrho=$4
    p=$(echo $p1 $p2 | awk '{print .5*($1+$2)}')
    pNI=$(echo $v $Dr | awk -v v0=$v -v Dr=$Dr '{print '"$realrho"'*.5*v0*v0/Dr}')
    sigma=$( echo $p1 $p2 | awk '{print ($1>$2)?($1-$2):($2-$1)}' )
    echo "$dt $realrho $rho $p $pNI $sigma" 
done >> StatPressureWCA-v-10-Dr-1-rho-0.8-vs-dt;














## This was to test the convergence as dt->0

# echo "# dt rho P sigma" > StatPressureWCA-rho-0.45
# for dt in 1e-6 2.5e-6 5e-6 1e-5 2e-5 4e-5 8e-5 8e-4 1e-3 2e-3 4e-3 8e-3 1e-2 2e-2 4e-2 8e-2 1e-1 2e-1;
# do
#     ./ComputePressure.sh test-pressure-WCA-rho-0.45-Lx-50-Ly-22-dt-"$dt"-pressure 50 1;
#     var=`awk '$1~/'"$t1"'/' test-pressure-WCA-rho-0.45-Lx-50-Ly-22-dt-"$dt"-pressure-av`
#     set -- $var
#     p1=$2
#     p2=$3
#     realrho=$4
#     p=$(echo $p1 $p2 | awk '{print .5*($1+$2)}')
#     sigma=$( echo $p1 $p2 | awk '{print ($1>$2)?($1-$2):($2-$1)}' )
#     echo "$dt $realrho $p $sigma" >>  StatPressureWCA-rho-0.45
# done;

# echo "# dt rho P sigma" > StatPressureWCA-rho-0.80
# for dt in 2e-5 4e-5 8e-5 2e-4 4e-4 8e-4 1e-3;
# do
#     ./ComputePressure.sh test-pressure-WCA-rho-0.8-Lx-50-Ly-22-dt-"$dt"-pressure 50 1;
#     var=`awk '$1~/'"$t1"'/' test-pressure-WCA-rho-0.8-Lx-50-Ly-22-dt-"$dt"-pressure-av`
#     set -- $var
#     p1=$2
#     p2=$3
#     realrho=$4
#     p=$(echo $p1 $p2 | awk '{print .5*($1+$2)}')
#     sigma=$( echo $p1 $p2 | awk '{print ($1>$2)?($1-$2):($2-$1)}' )
#     echo "$dt $realrho $p $sigma" >>  StatPressureWCA-rho-0.80
# done;

# echo "# dt rho P sigma" > StatPressureWCA-rho-0.80-Lx-90
# for dt in 5e-6 1e-5 2e-5 4e-5 8e-5 1e-4 2e-4 4e-4 8e-4 1e-3 1e-1;
# do
#     ./ComputePressure.sh test-pressure-WCA-rho-0.8-Lx-90-Ly-22-dt-"$dt"-pressure 50 1;
#     var=`awk '$1~/'"$t1"'/' test-pressure-WCA-rho-0.8-Lx-90-Ly-22-dt-"$dt"-pressure-av`
#     set -- $var
#     p1=$2
#     p2=$3
#     realrho=$4
#     p=$(echo $p1 $p2 | awk '{print .5*($1+$2)}')
#     sigma=$( echo $p1 $p2 | awk '{print ($1>$2)?($1-$2):($2-$1)}' )
#     echo "$dt $realrho $p $sigma" >>  StatPressureWCA-rho-0.80
# done;
