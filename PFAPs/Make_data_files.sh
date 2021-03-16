#!/bin/bash

## $1 is the filename
## $2 is Lx
## $3 is Ly
## $4 is size

sigma=$4

file=$1

#Calcul of max density
rho=$(awk 'BEGIN{rho=0} $6>rho {rho=$6} END{print rho}' $file)

echo "file $file";
Time=$(head -n 2 $file | tail -n 1 | awk '{printf("%.1f",$1)}')
echo "time $Time"
./draw_points_colored_by_density-size-inset.py $file $rho $Time $2 $3 $sigma; 


