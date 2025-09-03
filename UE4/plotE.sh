#!/usr/bin/gnuplot -persist

reset
set term png
set output "energy.png"

set title 'Energy of system'

set xlabel "t [d]"
set ylabel "E []"

plot 'EarthOrbit_Energy_Euler_1.dat' u 1:2 title 'Euler (dt=1)' w l, 'EarthOrbit_Energy_Euler_001.dat' u 1:2 title 'Euler (dt=0.01)' w l, 'EarthOrbit_Energy_RK4_1.dat' u 1:2 title 'RK4 (dt=1)' w l, 'EarthOrbit_Energy_RK4_001.dat' u 1:2 title 'RK4 (dt=0.01)' w l
