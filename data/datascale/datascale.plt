#!/usr/bin/gnuplot
set terminal eps font 'Verdana,14'
set output 'datascale.eps'

set xtics nomirror
set x2tics
set autoscale xfix
set autoscale x2fix

set title 'Data scaling: time vs. input size'
set ylabel 'Time (seconds)'
set xlabel 'Vertices + edges (millions)'
set xrange [0:140]

set key left top

plot 'datascale.dat' using ($4/1000000):5:x2tic(1) with points pointtype 9 lc rgb "red" title '1 worker', \
                  '' using ($4/1000000):6:x2tic(1) with points pointtype 11 lc rgb "blue" title '6 workers'

