#!/bin/bash

Dr=2

for v in 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9;
do
    screen -d -m ./abp-lj-input test-rho-2.5-v-"$v"-Dr-"$Dr"-L-80 2.5 80 80 2.5e-5 $v $Dr .42 .37 1 $RANDOM 2000 0 4 5 25 .5 10 test-rho-2.5-v-0-Dr-1-L80-input 16000
done

for v in 12.5 13 13.5 15;
do
    screen -d -m ./abp-lj-input test-rho-4-v-"$v"-Dr-"$Dr"-L-72 4 72 72 2.5e-5 $v $Dr .42 .37 1 $RANDOM 500 0 4 4 30 .5 10 test-rho-4-v-20-Dr-2-input 20736
done

