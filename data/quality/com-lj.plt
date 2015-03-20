#!/usr/bin/gnuplot
set terminal eps font 'Verdana,14'
set output 'com-lj-cost.eps'

set title 'com-LiveJournal partition costs'
set ylabel 'Edge partition CV'
set xlabel '#-Partitions'

set key top left reverse Left
set yrange [0:]

plot 'com-lj.cost' u 1:5 t 'fennel-VP',\
  '' u 1:7 t 'ordered fennel-VP' lt 8,\
  '' u 1:2 t 'fennel-EP',\
  '' u 1:3 t 'powergraph',\
  '' u 1:4 t 'sheep' lt 7,\
  '' u 1:6 t 'metis' lt 5
