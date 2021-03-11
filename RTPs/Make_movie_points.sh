#!/bin/bash

## $1 is the filename
## $2 is Lx
## $3 is Ly
## $4 is increment

if [[ $# -ne 4 ]] ; then
    echo 'usage: file Lx Ly dt'
    exit 0
fi

file=$1
Lx=$2
Ly=$3
#sigma=$4
dt=$4
filepos=$file"-pos"
fileobs=$file"-obs"

rm data*

#Compute the number of time blocks by counting the number of empty lines

#This stores the result of the command `awk 'NF<2 {m++} END{print m}' $filepos' in the variable N
N=$(awk 'NF<2 {m++} END{print m}' $filepos)
echo "number of time frames in the data file $N"

pad=$(echo ${#N} | awk '{print $1+1}')
echo "pad $pad"

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
' $filepos

for i in data*;
do
    Time=$(head -n 2 $i | tail -n 1 | awk '{printf("%.1f",$1)}')
    echo "file $i $Time";
    gnuplot <<EOF
    set style rect fc lt -1 fs solid 0.01 noborder
    set object rect from -1,-1 to 0,$Ly+1 fs solid fc rgb "white" front	
    set object rect from -1,-1 to $Lx+1,0 fs solid fc rgb "white" front
    set object rect from -1,$Ly to $Lx+1,$Ly+1 fs solid fc rgb "white" front	 
    set object rect from $Lx,-1 to $Lx+1,$Ly+1 fs solid fc rgb "white" front
    #    set key box opaque height 2

    set title 'Time $Time'
    set xr[-1:$Lx+1]
    set yr[-1:$Ly+1]
    set size square
    set object rectangle from  0,0 to $Lx,$Ly front
    set terminal png size 1600,1600
    set output "$i.png"
    unset key
    #size=$sigma
    set style fill solid
    pl "$fileobs" us 1:2:3 w circles,\
    "$i" us 3:4 w p pt 7 ps 2 lt -1
    set output
EOF
done

gnuplot <<EOF
    set style rect fc lt -1 fs solid 0.01 noborder
    set object rect from -1,-1 to 0,$Ly+1 fs solid fc rgb "white" front	
    set object rect from -1,-1 to $Lx+1,0 fs solid fc rgb "white" front
    set object rect from -1,$Ly to $Lx+1,$Ly+1 fs solid fc rgb "white" front	 
    set object rect from $Lx,-1 to $Lx+1,$Ly+1 fs solid fc rgb "white" front

    set title 'Full trajectory'
    set xr[-1:$Lx+1]
    set yr[-1:$Ly+1]
    set object rectangle from  0,0 to $Lx,$Ly
    set size square
    set terminal png size 1600,1600
    set output "$file-full-traj.png"
    unset key
    #size=$sigma
    set style fill solid
    pl "$fileobs" us 1:2:3 w circles,\
    "$filepos" us 3:4 w p pt 7 ps 2 lt -1
    set output
EOF

ffmpeg -y -r 10 -framerate 10 -i data%0"$pad"d.png -b:a 16M -vcodec libx264 "$file"-gnuplot.mp4 

# ffmepg makes nicer movies
# use imagemagick to produce a gif instead ?

rm data*

