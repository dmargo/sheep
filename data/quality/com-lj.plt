#!/usr/bin/gnuplot
set terminal eps font 'Verdana,14'
set output 'com-lj-cost.eps'

set title 'com-LiveJournal partition costs'
set ylabel 'Edge partition CV (millions)'
set xlabel '#-Partitions'

set key top left reverse Left
set yrange [0:]

plot 'com-lj.cost' u 1:($5/1000000) t 'fennel-VP',\
  '' u 1:($7/1000000) t 'ordered fennel-VP' lt 8,\
  '' u 1:($2/1000000) t 'fennel-EP',\
  '' u 1:($3/1000000) t 'powergraph',\
  '' u 1:($4/1000000) t 'sheep' lt 7,\
  '' u 1:($6/1000000) t 'metis' lt 5
