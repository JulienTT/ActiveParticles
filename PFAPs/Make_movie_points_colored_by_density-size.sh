#!/bin/bash

## $1 is the filename
## $2 is Lx
## $3 is Ly
## $4 is size
## $5 is increment

sigma=$4
dt=$5

file=$1

rm data*

N=$(awk 'NF<2 {m++} END{print m}' $file)
echo "N=$N"
pad=$(echo ${#N} | awk '{print $1+1}')
echo "pad $pad"

#Calcul of max density
rho=$(awk 'BEGIN{rho=0} $6>rho {rho=$6} END{print rho}' $file)

# i is the marker of the file

# t is a time of the next block to read. eps is a small increment used to "compare" the times.

# iread is turned off when a new file is made and turned to 1 when data are stored in this file

# NF>2 {if($1>t-eps) {iread=1;print $1,$2,$3,$4,$6 >> file}}
# -> if it is a data line (NF>2), then if $1 is larger than the current time ($1>t-eps) then you should record

# NF<2 {if(iread==1) {i+=1;iread=0;file=sprintf("data%0'"$pad"'d",i)}}
# if you are between two blocks and you have just read the last block
# (iread==1) then you should increment i and open a new file, and
# increment t up to the next time you want to record

awk 'BEGIN{iread=0;i=0;t=0;t_increment='"$dt"';eps=0.000001;file=sprintf("data%0'"$pad"'d",i)}
NF<2  {if(iread==1) {i+=1;iread=0;t+=t_increment;file=sprintf("data%0'"$pad"'d",i)}}
NF>2 {if($1>t-eps) {iread=1;print $1,$2,$3,$4,$6 >> file}}
' $file

for i in data*;
do
    echo "file $i";
    Time=$(head -n 2 $i | tail -n 1 | awk '{printf("%.1f",$1)}')
    echo "time $Time"
    ./draw_points_colored_by_density-size.py $i $rho $Time $2 $3 $sigma; 
done

ffmpeg -r 10 -i data%0"$pad"d.png -b:a 16M -vcodec libx264 "$file"-python.mp4 

# ffmepg makes nicer movies

#rm data0*

