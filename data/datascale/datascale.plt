#!/usr/bin/gnuplot
set terminal eps font 'Verdana,14'
set output 'datascale.eps'

set title 'Data scaling: time vs. input size'
set ylabel 'Time (seconds)'
set xlabel 'Vertices + edges (millions)'

set key left top

plot 'datascale.dat' using ($4/1000000):5 with lines title '1 worker', \
                  '' using ($4/1000000):5:1 with labels notitle, \
                  '' using ($4/1000000):6 with lines title '6 workers', \
                  '' using ($4/1000000):6:1 with labels notitle

