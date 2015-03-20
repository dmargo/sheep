#!/usr/bin/gnuplot
set terminal eps font 'Verdana,14'
set output 'hep.cost.eps'

set title 'HEP partition costs'
set ylabel 'Edge partition CV'
set xlabel '#-partitions'

set key top left reverse Left
set yrange [0:]

plot 'hep.cost' u 1:2 t 'sheep-degree' lt 7, '' u 1:3 t 'sheep-BC' lt 9, '' u 1:4 t 'metis' lt 5
