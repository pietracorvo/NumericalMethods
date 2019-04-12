#!/usr/bin/gnuplot -persist

reset
set term gif size 1200, 800
set term gif animate delay 100
set output "wannier.gif"

set xlabel "x [a.u.]"
set ylabel "|psi(x)|²"
set y2label "V(x)"
set yrange [-0.01:0.05]
set y2range [-1.01:3]
set ytics nomirror
set y2tics

energies="0.0 0.001 0.002 0.004 0.006 0.008 0.01 0.012 0.014 0.016 0.018 0.02"


do for [IND=0:11]{
		set title sprintf('Electron density at E = %s', word(energies, IND+1));
		plot 'wannierEV.dat' i IND u 1:2 w l ax x1y2 title 'V', \
		'' i IND u 1:3 lw 3 w l ax x1y1 title 'EV1', \
		'' i IND u 1:4 w l ax x1y1 title 'EV2', \
		'' i IND u 1:5 w l ax x1y1 title 'EV3', \
		'' i IND u 1:6 w l ax x1y1 title 'EV4', \
		'' i IND u 1:7 w l ax x1y1 title 'EV5', \
		'' i IND u 1:8 w l ax x1y1 title 'EV6', \
		'' i IND u 1:9 w l ax x1y1 title 'EV7', \
		'' i IND u 1:10 w l ax x1y1 title 'EV8', \
		'' i IND u 1:11 w l ax x1y1 title 'EV9', \
		'' i IND u 1:12 w l ax x1y1 title 'EV10', \
		'' i IND u 1:13 w l ax x1y1 title 'EV11', \
		'' i IND u 1:14 w l ax x1y1 title 'EV12', \
		'' i IND u 1:15 w l ax x1y1 title 'EV13', \
		'' i IND u 1:16 w l ax x1y1 title 'EV14', \
		'' i IND u 1:17 w l ax x1y1 title 'EV15', \
		'' i IND u 1:18 w l ax x1y1 title 'EV16', \
		'' i IND u 1:19 w l ax x1y1 title 'EV17', \
		'' i IND u 1:20 w l ax x1y1 title 'EV18', \
		'' i IND u 1:21 w l ax x1y1 title 'EV19', \
		'' i IND u 1:22 w l ax x1y1 title 'EV20';
		}


