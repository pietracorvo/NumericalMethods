#!/usr/bin/gnuplot -persist

reset
set term png size 1500,1000
set output "wannier_EV.png"

set title "Electron density over potential V (E=0)"
set xlabel "x [a.u.]"
set ylabel "|psi(x)|²"
set y2label "V(x)"
set yrange [-0.01:0.1]
set y2range [-1.01:2]
set ytics nomirror
set y2tics

#plot 'wannier_EV.dat' i 0 u 1:2 w l ax x1y2 title 'V', \
plot for [IND=3:22]{'wannier_EV.dat' i 0 u 1:IND w l ax x1y1 title sprintf('EV %i', IND-2)}
