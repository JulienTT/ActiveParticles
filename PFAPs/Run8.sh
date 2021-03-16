#!/bin/bash

Dr=2

for v in 12.25 12.75 13.25 13.75 14.5;
do
    screen -d -m ./abp-lj-input test-rho-4-v-"$v"-Dr-"$Dr"-L-72 4 72 72 2.5e-5 $v $Dr .42 .37 1 $RANDOM 500 0 4 4 30 .5 10 test-rho-4-v-20-Dr-2-input 20736
done

