#!/bin/bash

MAKE_DATA=false
PLOT_DATA=true

DDIR="${HOME}/data"
GRAPHS=()
GRAPHS+=("${DDIR}/amazon/com-amazon.ungraph.dat")
GRAPHS+=("${DDIR}/dblp/com-dblp.ungraph.dat")
GRAPHS+=("${DDIR}/youtube/com-youtube.ungraph.dat")
GRAPHS+=("${DDIR}/cit-patents/cit-Patents.dat")
GRAPHS+=("${DDIR}/com-lj/com-lj.ungraph.dat")
GRAPHS+=("${DDIR}/soc-lj/soc-LiveJournal1.dat")
GRAPHS+=("${DDIR}/orkut/com-orkut.ungraph.dat")

RDIR='data/quality'

if [ ! -d scripts ]; then
  exit 1
fi



if $MAKE_DATA; then
  mkdir -p $RDIR
  for G in ${GRAPHS[@]}; do
    DIRNAME=$(dirname $G)
    NAME=$(basename $DIRNAME)
    RAW="${RDIR}/${NAME}.raw"
    rm -f $RAW
    scripts/dist-partition.sh -w 1 $G $(seq 2 40) 2>&1 | tee -a "${RDIR}/${NAME}.raw"
  done
fi



if $PLOT_DATA; then
  RAW_DATA=( ${RDIR}/*.raw )
  for RAW in ${RAW_DATA[@]}; do
    NAME=$(basename $RAW .raw)
    grep 'ECV(down)' $RAW | egrep -o "[0][\.][[:digit:]]+" > "${RDIR}/${NAME}.dat"

gnuplot <<EOF
set output "${RDIR}/${NAME}.png"
set terminal png

set title "${NAME} partition quality"
set xlabel "#-Partitions"
set ylabel "Edge CV per edge"
set yrange [0:1]

plot "${RDIR}/${NAME}.dat" using (2+\$0):1 title "${NAME}"
EOF

  done
fi

