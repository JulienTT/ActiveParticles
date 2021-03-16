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

echo "rho max is $rho"

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
NF>2 {if($1>t-eps) {iread=1;print $1,$2,$3,$4,$5,$6 >> file}}
' $file

for i in data*;
do
    Time=$(head -n 2 $i | tail -n 1 | awk '{printf("%.1f",$1)}')
    echo "file $i $Time";
    gnuplot <<EOF
    set title 'Time $Time'
    set cbrange [0:$rho]
    set xr[0:$2]
    set yr[0:$3]
    set size square
    set terminal png size 1600,1600
    set output "$i.png"
    unset key
    size=$sigma
    set style fill solid
    pl "$i" us 3:4 w p pt 7 ps .2 lt -1,"$i" us 3:4:(size):(\$6*4/3.1415) w circles lc palette
    set output
EOF
done

ffmpeg -y -r 10 -framerate 10 -i data%0"$pad"d.png -b:a 16M -vcodec libx264 "$file"-gnuplot.mp4 

# ffmepg makes nicer movies

rm data*

