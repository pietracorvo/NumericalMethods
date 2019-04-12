#!/usr/bin/gnuplot -persist

reset
set term gif animate
set output "propagation.gif"
set yrange [-1:1.5]
set arrow from 0,-1 to 0,1.5 nohead


#set label "y=x" at 1,2
do for [i=0:999] {plot 'data_psi.dat' u 1:2 index i w l, '' u 1:3 index i w l} 

set output
