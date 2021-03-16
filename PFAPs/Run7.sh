#!/bin/bash

Dr=2

for v in 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9;
do
    screen -d -m ./abp-lj-input test-rho-2.5-v-"$v"-Dr-"$Dr"-L-80 2.5 80 80 2.5e-5 $v $Dr .42 .37 1 $RANDOM 2000 0 4 5 25 .5 10 test-rho-2.5-v-0-Dr-1-L80-input 16000
done

