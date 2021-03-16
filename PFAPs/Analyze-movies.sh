#!/bin/bash

################ MOVIES ############################
dt=2

Dr=2
rho=2.5

for v in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9;
do
    echo "";
    echo "rho = $rho      Dr = $Dr      v = $v ";
    echo "";
    ./Make_movie_points_colored_by_density-JT.sh test-rho-"$rho"-v-"$v"-Dr-"$Dr"-L-80-pos 80 80 .185 $dt
done


#### Rho 4 #######

rho=4
L=72

Dr=2

for v in 1 2 12.25 12.5 12.75 13 13.25 13.5 13.75 14.5 15;
do
    echo "";
    echo "rho = $rho      Dr = $Dr      v = $v ";
    echo "";
    ./Make_movie_points_colored_by_density-JT.sh test-rho-"$rho"-v-"$v"-Dr-"$Dr"-L-72-pos $L $L .185 $dt
done

Dr=2

#create histograms
for v in 0 2 3 4 5 6 7 8 9 10 12 14 16 18 20 22 24 26 35 40 45 50 60;
do
    echo "";
    echo "rho = $rho      Dr = $Dr      v = $v ";
    echo "";
    ./Make_movie_points_colored_by_density-JT.sh test-rho-"$rho"-v-"$v"-Dr-"$Dr"-pos $L $L .185 $dt
done
