plot "soma_v.dat" u 1:2 w lines lw 3 lc 7 title 'soma v'
set ylabel ("soma voltage (mV)")
set xlabel ("time (ms)")
set xrange [0:1200]
set yrange [-80:60]
set border lw 3
set output 'soma_v.png'
set terminal pngcairo size 1024, 343 enhanced font 'Verdana, 15'
replot
