#!/usr/bin/gnuplot -persist

reset
set term png size 1200, 800
set output "laplace.png"

set title 'Heatmap of Potetnial (500x500 gridpoints)'

set xrange[0:499]
set yrange[0:499]

set xlabel "x [cm]"
set ylabel "y [cm]"

plot 'laplace.dat' matrix with image
