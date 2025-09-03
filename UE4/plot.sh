#!/usr/bin/gnuplot -persist

# führe mich in bash mit 'gnuplot plotlap.sh' aus!

reset
set term png size 800, 800
set output "orbit.png"

set title 'Orbit with Euler and RK4 method'

set xrange [-1.5:1.5]
set yrange [-1.5:1.5]

set xlabel "x [AU]"
set ylabel "y [AU]"

plot 'EarthOrbit_Euler_1.dat' u 2:3 title 'Euler (dt=1)' w l, 'EarthOrbit_Euler_001.dat' u 2:3 title 'Euler (dt=0.01)' w l, 'EarthOrbit_RK4_1.dat' u 2:3 title 'RK4 (dt=1)' w l, 'EarthOrbit_RK4_001.dat' u 2:3 title 'RK4 (dt=0.01)' w l
