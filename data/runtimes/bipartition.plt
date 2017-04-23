#!/usr/bin/gnuplot
set terminal eps font 'Verdana,14'
set output 'bipartition.eps'

set xtics nomirror
set x2tics
set autoscale xfix
set autoscale x2fix

set title 'Comparative scaling: time vs. input size'
set ylabel 'Time (seconds)'
set xlabel 'Vertices + edges (millions)'
set xrange [0:140]

set key top left

plot 'bipartition.time' u ($5/1000000):4:x2tic(1) with points pt 5 title 'metis',\
                     '' u ($5/1000000):2:x2tic(1) with points pt 7 title 'fennel',\
                     '' u ($5/1000000):3:x2tic(1) with points pt 11 title 'sheep: 6 workers'

