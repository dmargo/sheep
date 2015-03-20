#!/usr/bin/gnuplot
set terminal eps font 'Verdana,14'
set output 'slurm-twitter-cost.eps'

set title 'Twitter partition costs'
set ylabel 'Edge partition CV'
set xlabel '#-Partitions'

set key top left reverse Left
set yrange [0:]

plot 'slurm-twitter.cost' u 1:2 t 'ordered fennel-VP' lt 8, '' u 1:3 t 'sheep' lt 7
