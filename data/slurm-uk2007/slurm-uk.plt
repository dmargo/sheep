#!/usr/bin/gnuplot
set terminal eps font 'Verdana,14'
set output 'slurm-uk.eps'

set style data histograms
set style histogram rowstacked
set style fill solid 1.0 border -1
set boxwidth 1 relative

set title 'Sheep distributed scaling on uk_2007'
set ylabel 'Time (seconds)
set xlabel '#-workers x #-nodes'

set yrange [0:]

plot "slurm-uk.summary.withtime" using 2 title "Load", \
  '' using 3 title "Sort", \
  '' using 4 title "Map", \
  '' using 5:xticlabel(1) title "Reduce",\
  '' using 6 title "Partition"

