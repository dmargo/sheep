#!/usr/bin/gnuplot
set terminal eps font 'Verdana,14'

set style data histograms
set style histogram rowstacked
set style fill solid 1.0 border -1
set boxwidth 1 relative

set output "twitter.fennel.time.eps"
set title 'Fennel partition scaling on Twitter'
set ylabel 'Time (seconds)'
set xlabel 'Fennel: #-Partitions'

set key top left
set yrange [0:1800]

plot 'slurm-fen.avg' using 2 title 'Load', \
  '' using 6:xticlabel(1) title 'Partition'  lt 5


set output 'sheep.fennel.time.eps'
set title 'Sheep parallel scaling on Twitter'
set xlabel 'Sheep: #-Workers'

set key top right
set yrange [0:600]

plot "slurm-sheep.avg" using 2 title "Load", \
  '' using 3 title "Sort", \
  '' using 4 title "Map", \
  '' using 5 title "Reduce", \
  '' using 6:xticlabel(1) title "Partition"

