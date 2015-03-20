#!/usr/bin/gnuplot
set terminal eps font 'Verdana,14'
set output "slurm-twitter-times.eps"

set style data histograms
set style histogram rowstacked
set style fill solid 1.0 border -1
set boxwidth 1 relative

set title "twitter parallel scaling"
set ylabel "Time (seconds)"
set xlabel "f#: fennel with # partitions, s#: sheep with # cores"

plot "slurm-25.avg" using 2 title "Load", \
  '' using 3 title "Sort", \
  '' using 4 title "Map", \
  '' using 5 title "Reduce", \
  '' using 6:xticlabel(1) title "Partition"

