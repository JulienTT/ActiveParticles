#!/bin/bash


Dr=2

for v in 1 2 3 4 5 6 7 8 9 10 12 14 16 18 20 22 24 26 28 30;
do
    screen -d -m ./abp-lj test-rho-4-v-"$v"-Dr-"$Dr" 4 72 72 2.5e-5 $v $Dr .42 .37 1 $RANDOM 2000 0 1 3 30 .5 10
done

Dr=3

for v in 1 2 3 4 5 6 7 8 9 10 12 14 16 18 20 22 24 26 28 30;
do
    screen -d -m ./abp-lj test-rho-4-v-"$v"-Dr-"$Dr" 4 72 72 2.5e-5 $v $Dr .42 .37 1 $RANDOM 2000 0 1 3 30 .5 10
done

Dr=2

for v in 35 40 45 50 60 70 80 90 100;
do
    screen -d -m ./abp-lj test-rho-4-v-"$v"-Dr-"$Dr" 4 72 72 1e-5 $v $Dr .42 .37 1 $RANDOM 2000 0 1 3 30 .5 10
done

Dr=3

for v in 35 40 45 50 60 70 80 90 100;
do
    screen -d -m ./abp-lj test-rho-4-v-"$v"-Dr-"$Dr" 4 72 72 1e-5 $v $Dr .42 .37 1 $RANDOM 2000 0 1 3 30 .5 10
done
