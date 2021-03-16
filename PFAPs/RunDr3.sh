#!/bin/bash

for v in 40 80 100;
do
    ./abp-lj test-rho-3-v-"$v"-Dr-"$Dr" 3 80 80 1e-5 $v $Dr .42 .37 1 $RANDOM 1000 0 1 4 30 .5 10
done
