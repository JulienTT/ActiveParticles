#!/bin/bash


Dr=2

for v in 2 4 6 8 10 15 20;
do
    screen -d -m ./abp-lj test-rho-2.5-v-"$v"-Dr-"$Dr" 2.5 51 51 2.5e-5 $v $Dr .42 .37 1 $RANDOM 2000 0 1 3 30 .5 10
done

Dr=3

for v in 2 5 10 15 20 25 30;
do
    screen -d -m ./abp-lj test-rho-2.5-v-"$v"-Dr-"$Dr" 2.5 51 51 2.5e-5 $v $Dr .42 .37 1 $RANDOM 2000 0 1 3 30 .5 10
done
