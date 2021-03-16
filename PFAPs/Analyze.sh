#!/bin/bash


################ HISTOGRAMS ############################

time_offset=300

### rho 3 ###

Dr=1
rho=3

## create each histogram
for v in 100;
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
for v in 1 5 10 20 40 100;
do
    pdftk test-rho-"$rho"-Dr-"$Dr"-histo-plot.pdf test-rho-"$rho"-v-"$v"-Dr-"$Dr"-histo-plot.pdf cat output temp.pdf
    mv temp.pdf test-rho-"$rho"-Dr-"$Dr"-histo-plot.pdf 
done


Dr=2
#Create histograms
for v in 100;
do
    echo "";
    echo "rho = $rho      Dr = $Dr      v = $v ";
    echo "";
    time=$(awk 'NF>1 {time=$1} END{print time-'"$time_offset"'}' test-rho-"$rho"-v-"$v"-Dr-"$Dr"-pos)
    echo time "$time"
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

### rho 4 ###
rho=4

Dr=1

Dr=2

#create histograms
for v in 0 1 2 3 4 5 6 7 8 9 10 12 14 16 18 20 22 24 26 35 40 45 50 60;
do
    echo "";
    echo "rho = $rho      Dr = $Dr      v = $v ";
    echo "";
    time=$(awk 'NF>1 {time=$1} END{print time-'"$time_offset"'}' test-rho-"$rho"-v-"$v"-Dr-"$Dr"-pos)
    echo time "$time"
    ./MakeHistogram.sh test-rho-"$rho"-v-"$v"-Dr-"$Dr"-histo $time 1 test-rho-"$rho"-v-"$v"-Dr-"$Dr"-histo-plot $rho $v $Dr
done

for v in 12.25 12.5 12.75 13 13.25 13.5 13.75 14.5 15;
do
    echo "";
    echo "rho = $rho      Dr = $Dr      v = $v ";
    echo "";
    time=$(awk 'NF>1 {time=$1} END{print time-'"$time_offset"'}' test-rho-"$rho"-v-"$v"-Dr-"$Dr"-L-72-pos)
    echo time "$time"
    ./MakeHistogram.sh test-rho-"$rho"-v-"$v"-Dr-"$Dr"-histo $time 1 test-rho-"$rho"-v-"$v"-Dr-"$Dr"-L-72-histo-plot $rho $v $Dr
done

#bundle them all in one pdf
cp test-rho-"$rho"-v-0-Dr-2-histo-plot.pdf test-rho-"$rho"-Dr-"$Dr"-histo-plot.pdf
for v in 1 2 3 4 5 6 7 8 9 10 12;
do
    pdftk test-rho-"$rho"-Dr-"$Dr"-histo-plot.pdf test-rho-"$rho"-v-"$v"-Dr-"$Dr"-histo-plot.pdf cat output temp.pdf
    mv temp.pdf test-rho-"$rho"-Dr-"$Dr"-histo-plot.pdf 
done

for v in 12.25 12.5 12.75 13 13.25 13.5 13.75;
do
    pdftk test-rho-"$rho"-Dr-"$Dr"-histo-plot.pdf test-rho-"$rho"-v-"$v"-Dr-"$Dr"-L-72-histo-plot.pdf cat output temp.pdf
    mv temp.pdf test-rho-"$rho"-Dr-"$Dr"-histo-plot.pdf 
done

for v in 14;
do
    pdftk test-rho-"$rho"-Dr-"$Dr"-histo-plot.pdf test-rho-"$rho"-v-"$v"-Dr-"$Dr"-histo-plot.pdf cat output temp.pdf
    mv temp.pdf test-rho-"$rho"-Dr-"$Dr"-histo-plot.pdf 
done

for v in 14.5 15;
do
    pdftk test-rho-"$rho"-Dr-"$Dr"-histo-plot.pdf test-rho-"$rho"-v-"$v"-Dr-"$Dr"-L-72-histo-plot.pdf cat output temp.pdf
    mv temp.pdf test-rho-"$rho"-Dr-"$Dr"-histo-plot.pdf 
done

for v in 16 18 20 22 24 26 28 30 35 40 45 50 60;
do
    pdftk test-rho-"$rho"-Dr-"$Dr"-histo-plot.pdf test-rho-"$rho"-v-"$v"-Dr-"$Dr"-histo-plot.pdf cat output temp.pdf
    mv temp.pdf test-rho-"$rho"-Dr-"$Dr"-histo-plot.pdf 
done

Dr=3


##### rho 2.5 #######

rho=2.5
Dr=1
#create histograms
for v in 0.3 0.4 0.6 0.7 0.8 0.9 1.25;
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
for v in 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.25;
do
    pdftk test-rho-"$rho"-Dr-"$Dr"-histo-plot.pdf test-rho-"$rho"-v-"$v"-Dr-"$Dr"-histo-plot.pdf cat output temp.pdf
    mv temp.pdf test-rho-"$rho"-Dr-"$Dr"-histo-plot.pdf 
done

Dr=2
#create histograms
for v in 0.3 0.4 0.6 0.7 0.8 0.9 1.25 ;
do
    echo "";
    echo "rho = $rho      Dr = $Dr      v = $v ";
    echo "";
    time=$(awk 'NF>1 {time=$1} END{print time-'"$time_offset"'}' test-rho-"$rho"-v-"$v"-Dr-"$Dr"-pos)
    echo time "$time"
    ./MakeHistogram.sh test-rho-"$rho"-v-"$v"-Dr-"$Dr"-histo $time 1 test-rho-"$rho"-v-"$v"-Dr-"$Dr"-histo-plot $rho $v $Dr
done

#create histograms
for v in 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2;
do
    echo "";
    echo "rho = $rho      Dr = $Dr      v = $v ";
    echo "";
    time=$(awk 'NF>1 {time=$1} END{print time-'"$time_offset"'}' test-rho-"$rho"-v-"$v"-Dr-"$Dr"-L-80-pos)
    echo time "$time"
    ./MakeHistogram.sh test-rho-"$rho"-v-"$v"-Dr-"$Dr"-L-80-histo $time 1 test-rho-"$rho"-v-"$v"-Dr-"$Dr"-L-80-histo-plot $rho $v $Dr
done

#bundle them all in one pdf
cp test-rho-"$rho"-v-0-Dr-2-histo-plot.pdf test-rho-"$rho"-Dr-"$Dr"-histo-plot.pdf
for v in  0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.25 2 4 6 8 10 15 20;
do
    pdftk test-rho-"$rho"-Dr-"$Dr"-histo-plot.pdf test-rho-"$rho"-v-"$v"-Dr-"$Dr"-histo-plot.pdf cat output temp.pdf
    mv temp.pdf test-rho-"$rho"-Dr-"$Dr"-histo-plot.pdf 
done

cp test-rho-"$rho"-v-0-Dr-2-L-80-histo-plot.pdf test-rho-"$rho"-Dr-"$Dr"-L-80-histo-plot.pdf
for v in 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2;
do
    pdftk test-rho-"$rho"-Dr-"$Dr"-L-80-histo-plot.pdf test-rho-"$rho"-v-"$v"-Dr-"$Dr"-L-80-histo-plot.pdf cat output temp.pdf
    mv temp.pdf test-rho-"$rho"-Dr-"$Dr"-L-80-histo-plot.pdf 
done

Dr=3

#bundle them all in one pdf
cp test-rho-"$rho"-v-0-Dr-2-histo-plot.pdf test-rho-"$rho"-Dr-"$Dr"-histo-plot.pdf
for v in 2 5 10 15 20 25 30;
do
    pdftk test-rho-"$rho"-Dr-"$Dr"-histo-plot.pdf test-rho-"$rho"-v-"$v"-Dr-"$Dr"-histo-plot.pdf cat output temp.pdf
    mv temp.pdf test-rho-"$rho"-Dr-"$Dr"-histo-plot.pdf 
done

################ MOVIES ############################

Dr=1
dt=2
#### Rho 3 #######

rho=3
for v in 100;
do
    echo "";
    echo "rho = $rho      Dr = $Dr      v = $v ";
    echo "";
    
    ./Make_movie_points_colored_by_density-JT.sh test-rho-"$rho"-v-"$v"-Dr-"$Dr"-pos 80 80 .185 $dt
done

Dr=2

for v in 100;
do
    echo "";
    echo "rho = $rho      Dr = $Dr      v = $v ";
    echo "";
    ./Make_movie_points_colored_by_density-JT.sh test-rho-"$rho"-v-"$v"-Dr-"$Dr"-pos 80 80 .185 $dt
done

Dr=3


#### Rho 4 #######

rho=4

Dr=2
L=72
for v in 0 1 2 3 4 5 6 7 8 9 10 12 14 16 18 20 22 24 26 35 40 45 50 60;
do
    echo "";
    echo "rho = $rho      Dr = $Dr      v = $v ";
    echo "";
    ./Make_movie_points_colored_by_density-JT.sh test-rho-"$rho"-v-"$v"-Dr-"$Dr"-pos $L $L .185 $dt
done

for v in 12.25 12.5 12.75 13 13.25 13.5 13.75 14.5 15;
do
    echo "";
    echo "rho = $rho      Dr = $Dr      v = $v ";
    echo "";
    ./Make_movie_points_colored_by_density-JT.sh test-rho-"$rho"-v-"$v"-Dr-"$Dr"-L-72-pos $L $L .185 $dt
done
	 
Dr=3

#### Rho 2.5 #######

rho=2.5

Dr=1
L=51
#create histograms
for v in 0.3 0.4 0.6 0.7 0.8 0.9 1.25;
do
    echo "";
    echo "rho = $rho      Dr = $Dr      v = $v ";
    echo "";
    ./Make_movie_points_colored_by_density-JT.sh test-rho-"$rho"-v-"$v"-Dr-"$Dr"-pos $L $L .185 $dt
done;

Dr=2
L=51
for v in 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.25  4 6 8 10 15 20;
do
    echo "";
    echo "rho = $rho      Dr = $Dr      v = $v ";
    echo "";
    ./Make_movie_points_colored_by_density-JT.sh test-rho-"$rho"-v-"$v"-Dr-"$Dr"-pos $L $L .185 $dt
done

L=80
for v in 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2;
do
    echo "";
    echo "rho = $rho      Dr = $Dr      v = $v ";
    echo "";
    ./Make_movie_points_colored_by_density-JT.sh test-rho-"$rho"-v-"$v"-Dr-"$Dr"-L-80-pos $L $L .185 $dt
done



Dr=3
L=51
for v in 2 5 10 15 20 25 30;
do
    echo "";
    echo "rho = $rho      Dr = $Dr      v = $v ";
    echo "";
    ./Make_movie_points_colored_by_density-JT.sh test-rho-"$rho"-v-"$v"-Dr-"$Dr"-pos $L $L .185 $dt
done


