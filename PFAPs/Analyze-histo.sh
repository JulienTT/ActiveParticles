#!/bin/bash

################ HISTOGRAMS ############################

time_offset=41

### rho 3 ###

Dr=1
rho=3

## create each histogram
for v in 0 1 5 10 20 40 80 100;
do
    echo "";
    echo "rho = $rho      Dr = $Dr      v = $v ";
    echo "";
    time=$(awk 'NF>1 {time=$1} END{print time-'"$time_offset"'}' test-rho-"$rho"-v-"$v"-Dr-"$Dr"-pos)
    echo time "$time"
    ./MakeHistogram.sh test-rho-"$rho"-v-"$v"-Dr-"$Dr"-histo $time 1 test-rho-"$rho"-v-"$v"-Dr-"$Dr"-histo-plot $rho $v $Dr
done

## bundle them all in one pdf
cp test-rho-"$rho"-v-0-Dr-"$Dr"-histo-plot.pdf test-rho-"$rho"-Dr-"$Dr"-histo-plot.pdf
for v in 1 5 10 20 40 80 100;
do
    pdftk test-rho-"$rho"-Dr-"$Dr"-histo-plot.pdf test-rho-"$rho"-v-"$v"-Dr-"$Dr"-histo-plot.pdf cat output temp.pdf
    mv temp.pdf test-rho-"$rho"-Dr-"$Dr"-histo-plot.pdf 
done


Dr=2
#Create histograms
for v in 1 5 10 20 40 80 100;
do
    echo "";
    echo "rho = $rho      Dr = $Dr      v = $v ";
    echo "";
    echo time "$time"
    time=$(awk 'NF>1 {time=$1} END{print time-'"$time_offset"'}' test-rho-"$rho"-v-"$v"-Dr-"$Dr"-pos)
    ./MakeHistogram.sh test-rho-"$rho"-v-"$v"-Dr-"$Dr"-histo $time 1 test-rho-"$rho"-v-"$v"-Dr-"$Dr"-histo-plot $rho $v $Dr
done

## bundle them all in one pdf
cp test-rho-"$rho"-v-0-Dr-1-histo-plot.pdf test-rho-"$rho"-Dr-"$Dr"-histo-plot.pdf
for v in 1 5 10 20 40 80 100;
do
    pdftk test-rho-"$rho"-Dr-"$Dr"-histo-plot.pdf test-rho-"$rho"-v-"$v"-Dr-"$Dr"-histo-plot.pdf cat output temp.pdf
    mv temp.pdf test-rho-"$rho"-Dr-"$Dr"-histo-plot.pdf 
done

Dr=3

#create histograms
for v in 1 5 10 20 40 80 100;
do
    echo "";
    echo "rho = $rho      Dr = $Dr      v = $v ";
    echo "";
    echo time "$time"
    time=$(awk 'NF>1 {time=$1} END{print time-'"$time_offset"'}' test-rho-"$rho"-v-"$v"-Dr-"$Dr"-pos)
    ./MakeHistogram.sh test-rho-"$rho"-v-"$v"-Dr-"$Dr"-histo $time 1 test-rho-"$rho"-v-"$v"-Dr-"$Dr"-histo-plot $rho $v $Dr
done

## bundle them all in one pdf
cp test-rho-"$rho"-v-0-Dr-1-histo-plot.pdf test-rho-"$rho"-Dr-"$Dr"-histo-plot.pdf
for v in 1 5 10 20 40 80 100;
do
    pdftk test-rho-"$rho"-Dr-"$Dr"-histo-plot.pdf test-rho-"$rho"-v-"$v"-Dr-"$Dr"-histo-plot.pdf cat output temp.pdf
    mv temp.pdf test-rho-"$rho"-Dr-"$Dr"-histo-plot.pdf 
done

### rho 4 ###
rho=4

Dr=1

#create histograms
for v in 0.1 0.2 0.5 1 2 4;
do
    echo "";
    echo "rho = $rho      Dr = $Dr      v = $v ";
    echo "";
    echo time "$time"
    time=$(awk 'NF>1 {time=$1} END{print time-'"$time_offset"'}' test-rho-"$rho"-v-"$v"-Dr-"$Dr"-pos)
    ./MakeHistogram.sh test-rho-"$rho"-v-"$v"-Dr-"$Dr"-histo $time 1 test-rho-"$rho"-v-"$v"-Dr-"$Dr"-histo-plot $rho $v $Dr
done

#bundle them all in one pdf
cp test-rho-"$rho"-v-0-Dr-2-histo-plot.pdf test-rho-"$rho"-Dr-"$Dr"-histo-plot.pdf
for v in 0.1 0.2 0.5 1 2 4;
do
    pdftk test-rho-"$rho"-Dr-"$Dr"-histo-plot.pdf test-rho-"$rho"-v-"$v"-Dr-"$Dr"-histo-plot.pdf cat output temp.pdf
    mv temp.pdf test-rho-"$rho"-Dr-"$Dr"-histo-plot.pdf 
done

Dr=2

#create histograms
for v in 0 1 2 3 4 5 6 7 8 9 10 12 14 16 18 20 22 24 26 28 30 35 40 45 50 60;
do
    echo "";
    echo "rho = $rho      Dr = $Dr      v = $v ";
    echo "";
    echo time "$time"
    time=$(awk 'NF>1 {time=$1} END{print time-'"$time_offset"'}' test-rho-"$rho"-v-"$v"-Dr-"$Dr"-pos)
    ./MakeHistogram.sh test-rho-"$rho"-v-"$v"-Dr-"$Dr"-histo $time 1 test-rho-"$rho"-v-"$v"-Dr-"$Dr"-histo-plot $rho $v $Dr
done

#bundle them all in one pdf
cp test-rho-"$rho"-v-0-Dr-2-histo-plot.pdf test-rho-"$rho"-Dr-"$Dr"-histo-plot.pdf
for v in 1 2 3 4 5 6 7 8 9 10 12 14 16 18 20 22 24 26 28 30 35 40 45 50 60 70 80 90 100;
do
    pdftk test-rho-"$rho"-Dr-"$Dr"-histo-plot.pdf test-rho-"$rho"-v-"$v"-Dr-"$Dr"-histo-plot.pdf cat output temp.pdf
    mv temp.pdf test-rho-"$rho"-Dr-"$Dr"-histo-plot.pdf 
done

Dr=3

#create histograms
for v in 1 2 3 4 5 6 7 8 9 10 12 14 16 18 20 22 24 26 28 30 35 40 45 50 60 70;
do
    echo "";
    echo "rho = $rho      Dr = $Dr      v = $v ";
    echo "";
    echo time "$time"
    time=$(awk 'NF>1 {time=$1} END{print time-'"$time_offset"'}' test-rho-"$rho"-v-"$v"-Dr-"$Dr"-pos)
    ./MakeHistogram.sh test-rho-"$rho"-v-"$v"-Dr-"$Dr"-histo $time 1 test-rho-"$rho"-v-"$v"-Dr-"$Dr"-histo-plot $rho $v $Dr
done

#bundle them all in one pdf
cp test-rho-"$rho"-v-0-Dr-2-histo-plot.pdf test-rho-"$rho"-Dr-"$Dr"-histo-plot.pdf
for v in 2 3 4 5 6 7 8 9 10 12 14 16 18 20 22 24 26 28 30 35 40 45 50 60 70;
do
    pdftk test-rho-"$rho"-Dr-"$Dr"-histo-plot.pdf test-rho-"$rho"-v-"$v"-Dr-"$Dr"-histo-plot.pdf cat output temp.pdf
    mv temp.pdf test-rho-"$rho"-Dr-"$Dr"-histo-plot.pdf 
done

##### rho 2.5 #######

Dr=2
rho=2.5

#create histograms
for v in 0 0.1 0.2 0.5 1 2 4 6 8 10 15 20;
do
    echo "";
    echo "rho = $rho      Dr = $Dr      v = $v ";
    echo "";
    time=$(awk 'NF>1 {time=$1} END{print time-'"$time_offset"'}' test-rho-"$rho"-v-"$v"-Dr-"$Dr"-pos)
    echo time "$time"
    ./MakeHistogram.sh test-rho-"$rho"-v-"$v"-Dr-"$Dr"-histo $time 1 test-rho-"$rho"-v-"$v"-Dr-"$Dr"-histo-plot $rho $v $Dr
done

#bundle them all in one pdf
cp test-rho-"$rho"-v-0-Dr-2-histo-plot.pdf test-rho-"$rho"-Dr-"$Dr"-histo-plot.pdf
for v in 0.1 0.2 0.5 1 2 4 6 8 10 15 20;
do
    pdftk test-rho-"$rho"-Dr-"$Dr"-histo-plot.pdf test-rho-"$rho"-v-"$v"-Dr-"$Dr"-histo-plot.pdf cat output temp.pdf
    mv temp.pdf test-rho-"$rho"-Dr-"$Dr"-histo-plot.pdf 
done

Dr=3

#create histograms
for v in 2 5 10 15 20 25 30;
do
    echo "";
    echo "rho = $rho      Dr = $Dr      v = $v ";
    echo "";
    time=$(awk 'NF>1 {time=$1} END{print time-'"$time_offset"'}' test-rho-"$rho"-v-"$v"-Dr-"$Dr"-pos)
    echo time "$time"
    ./MakeHistogram.sh test-rho-"$rho"-v-"$v"-Dr-"$Dr"-histo $time 1 test-rho-"$rho"-v-"$v"-Dr-"$Dr"-histo-plot $rho $v $Dr
done

#bundle them all in one pdf
cp test-rho-"$rho"-v-0-Dr-2-histo-plot.pdf test-rho-"$rho"-Dr-"$Dr"-histo-plot.pdf
for v in 5 10 15 20 25 30;
do
    pdftk test-rho-"$rho"-Dr-"$Dr"-histo-plot.pdf test-rho-"$rho"-v-"$v"-Dr-"$Dr"-histo-plot.pdf cat output temp.pdf
    mv temp.pdf test-rho-"$rho"-Dr-"$Dr"-histo-plot.pdf 
done

##########################T 0.40####################################
time_offset=41
rho=2.5
Dr=1
T=0.40
#create histograms
for v in 0 0.5;
do
    echo "";
    echo "T = $T      rho = $rho      Dr = $Dr      v = $v ";
    echo "";
    time=$(awk 'NF>1 {time=$1} END{print time-'"$time_offset"'}' Data/data-rho-"$rho"-v-"$v"-Dr-"$Dr"-T-"$T"-pos)
    echo time "$time"
    ./MakeHistogram.sh Data/data-rho-"$rho"-v-"$v"-Dr-"$Dr"-T-"$T"-histo $time 1 Data/data-rho-"$rho"-v-"$v"-Dr-"$Dr"-T-"$T"-histo-plot $rho $v $Dr
done

#bundle them all in one pdf
cp Data/data-rho-"$rho"-v-0-Dr-"$Dr"-T-"$T"-histo-plot.pdf data-rho-"$rho"-Dr-"$Dr"-T-"$T"-histo-plot.pdf
for v in 0.5;
do
    pdftk data-rho-"$rho"-Dr-"$Dr"-T-"$T"-histo-plot.pdf Data/data-rho-"$rho"-v-"$v"-Dr-"$Dr"-T-"$T"-histo-plot.pdf cat output temp.pdf
    mv temp.pdf data-rho-"$rho"-Dr-"$Dr"-T-"$T"-histo-plot.pdf 
done

##########################T 0.41####################################
time_offset=41
rho=2.5
Dr=1
T=0.40
#create histograms
for v in 0 0.5;
do
    echo "";
    echo "T = $T      rho = $rho      Dr = $Dr      v = $v ";
    echo "";
    time=$(awk 'NF>1 {time=$1} END{print time-'"$time_offset"'}' Data/data-rho-"$rho"-v-"$v"-Dr-"$Dr"-T-"$T"-pos)
    echo time "$time"
    ./MakeHistogram.sh Data/data-rho-"$rho"-v-"$v"-Dr-"$Dr"-T-"$T"-histo $time 1 Data/data-rho-"$rho"-v-"$v"-Dr-"$Dr"-T-"$T"-histo-plot $rho $v $Dr
done

#bundle them all in one pdf
cp Data/data-rho-"$rho"-v-0-Dr-"$Dr"-T-"$T"-histo-plot.pdf data-rho-"$rho"-Dr-"$Dr"-T-"$T"-histo-plot.pdf
for v in 0.5;
do
    pdftk data-rho-"$rho"-Dr-"$Dr"-T-"$T"-histo-plot.pdf Data/data-rho-"$rho"-v-"$v"-Dr-"$Dr"-T-"$T"-histo-plot.pdf cat output temp.pdf
    mv temp.pdf data-rho-"$rho"-Dr-"$Dr"-T-"$T"-histo-plot.pdf 
done
