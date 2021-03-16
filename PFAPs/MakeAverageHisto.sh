#!/bin/bash
# This script generates a file that can be plotted to get the histogram

rm temp data*

file=$1
time=$2
dt=$3
fileout=$4

awk -v time=$time '$1>time {m=1} m==1 {print $0}' $file > temp

# i is the marker of the file

# t is a time of the next block to read. eps is a small increment used to "compare" the times.

# iread is turned off when a new file is made and turned to 1 when data are stored in this file

# NF>2 {if($1>t-eps) {iread=1;print $2,$3 >> file}}
# -> if it is a data line (NF>2), then if $1 is larger than the current time ($1>t-eps) then you should record

# NF<2 {if(iread==1) {i+=1;iread=0;file=sprintf("data%0'"$pad"'d",i)}}
# if you are between two blocks and you have just read the last block
# (iread==1) then you should increment i and open a new file, and
# increment t up to the next time you want to record

awk 'BEGIN{iread=0;i=0;file=sprintf("data%0'"$pad"'d",i)}
NF<2  {if(iread==1) {i+=1;iread=0;file=sprintf("data%0'"$pad"'d",i)}}
NF>2 {iread=1;print $2,$3  >> file}
' temp

awk '{print $1}' data0 > temp2

for i in data*;
do
    awk '{print $2}' $i > temp;
    paste temp2 temp > temp3;
    mv temp3 temp2
done

awk '{m=0;for(i=2;i<=NF;i++) {m+=$i}; print $1, m/(NF-1)}' temp2 > $fileout

rm temp data*


    
