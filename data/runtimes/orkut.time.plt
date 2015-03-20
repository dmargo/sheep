#!/usr/bin/gnuplot
set terminal eps font 'Verdana,14'
set output 'orkut.time.eps'

set title 'Orkut partition scaling'
set ylabel "Time (seconds)"
set xlabel "#-Partitions"

set key center right
set yrange [0:]

plot 'orkut.time' u 1:2 t 'metis', '' u 1:3 t 'fennel', '' u 1:4 t 'sheep: 6 workers'
