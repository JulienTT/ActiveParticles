#!/bin/bash

Dr=1

for v in 0.3 0.4 0.6 0.7 0.8 0.9 1.25;
do
    screen -d -m ./abp-lj-input test-rho-2.5-v-"$v"-Dr-"$Dr" 2.5 51 51 2.5e-5 $v $Dr .42 .37 1 $RANDOM 2000 0 1 3 25 .5 10 test-rho-2.5-v-0-Dr-2-input 6502
done


Dr=2

for v in 0.3 0.4 0.6 0.7 0.8 0.9 1.25;
do
    screen -d -m ./abp-lj-input test-rho-2.5-v-"$v"-Dr-"$Dr" 2.5 51 51 2.5e-5 $v $Dr .42 .37 1 $RANDOM 2000 0 1 3 25 .5 10 test-rho-2.5-v-0-Dr-2-input 6502
done

