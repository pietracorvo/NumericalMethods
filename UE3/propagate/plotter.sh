#!/usr/bin/gnuplot -persist

# führe mich in bash mit 'gnuplot plotter.sh' aus!

reset
set term gif animate size 1200, 800
set output "propagation.gif"

set yrange [-1.1:1]
set xlabel "x"
set arrow from 0,-1.1 to 0,1 nohead

do for [i=0:199] {	
					set title sprintf('t = %i s', i)
					plot 'data_psi.dat' u 1:2 index i w l title 'psi²(x)', \
					'' u 1:3 index i w l title 'V(x)' }
