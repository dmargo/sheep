#!/bin/bash
if [ ! -d scripts ]; then
  exit 1
fi

TRUE=0
FALSE=1

MAKE_DATA=$FALSE
PLOT_DATA=$FALSE
RDIR='data/parallel'
TRIALS=6

while getopts "mpo:t:airc:" opt; do
  case $opt in
    m)
      MAKE_DATA=$TRUE;;
    p)
      PLOT_DATA=$TRUE;;
    o)
      RDIR=$OPTARG;;
    t)
      TRIALS=$OPTARG;;
    a)
      VERTICAL='-a';;
    i)
      MPI_SORT='-i';;
    r)
      MPI_REDUCE='-r';;
    c)
      CORES="-c $OPTARG";;
    :)
      echo "Option -$OPTARG requires an argument."
      exit 1;;
    \?)
      echo "Invalid option: -$OPTARG"
      exit 1;;
  esac
done

DDIR="${HOME}/data"
GRAPHS=()
GRAPHS+=("${DDIR}/amazon/com-amazon.ungraph.dat")
GRAPHS+=("${DDIR}/dblp/com-dblp.ungraph.dat")
GRAPHS+=("${DDIR}/youtube/com-youtube.ungraph.dat")
GRAPHS+=("${DDIR}/cit-patents/cit-Patents.dat")
GRAPHS+=("${DDIR}/com-lj/com-lj.ungraph.dat")
GRAPHS+=("${DDIR}/soc-lj/soc-LiveJournal1.dat")
GRAPHS+=("${DDIR}/orkut/com-orkut.ungraph.dat")
WORKER_LIST=( 1 2 4 6 8 10 12 )



if [ $MAKE_DATA -eq $TRUE ]; then
  mkdir -p $RDIR

  for G in ${GRAPHS[@]}; do
    DIRNAME=$(dirname $G)
    NAME=$(basename $DIRNAME)
    RAW="${RDIR}/${NAME}.raw"
    rm -f $RAW

    for WORKERS in ${WORKER_LIST[@]}; do
      for i in $(seq 1 $TRIALS); do
        echo "Starting with $WORKERS workers..." | tee -a $RAW
        scripts/dist-partition.sh $VERTICAL $MPI_SORT $MPI_REDUCE $CORES -w $WORKERS $G 0 | tee -a $RAW
        echo | tee -a $RAW
      done
    done
  done
fi



if [ $PLOT_DATA -eq $TRUE ]; then
  RAW_DATA=( ${RDIR}/*.raw )
  for RAW in ${RAW_DATA[@]}; do
    NAME=$(basename $RAW .raw)

    egrep "^Starting with[[:blank:]]" $RAW | egrep -o "[[:digit:]]+" > "/tmp/${NAME}.workers"
    egrep "^Loaded graph[[:blank:]]" $RAW | egrep -o "[[:digit:]]*\.[[:digit:]]+" > "/tmp/${NAME}.load"
    egrep "^Sorted[[:blank:]]" $RAW | egrep -o "[[:digit:]]*\.[[:digit:]]+" > "/tmp/${NAME}.sort"
    egrep "^Mapped[[:blank:]]" $RAW | egrep -o "[[:digit:]]*\.[[:digit:]]+" > "/tmp/${NAME}.map"
    egrep "^Reduced[[:blank:]]" $RAW | egrep -o "[[:digit:]]*\.[[:digit:]]+" > "/tmp/${NAME}.red"

    paste /tmp/${NAME}.workers /tmp/${NAME}.load /tmp/${NAME}.sort /tmp/${NAME}.map /tmp/${NAME}.red > ${RDIR}/${NAME}.dat
    rm /tmp/${NAME}.workers /tmp/${NAME}.load /tmp/${NAME}.sort /tmp/${NAME}.map /tmp/${NAME}.red

    rm -f "${RDIR}/${NAME}.avg"
    for W in $(awk '{print $1}' ${RDIR}/${NAME}.dat | sort -nu); do
      echo -n "$W " >> "${RDIR}/${NAME}.avg"
      egrep "^$W[[:blank:]]" ${RDIR}/${NAME}.dat | awk 'NR > 1' |
          awk '{ls += $2; ss += $3; ms += $4; rs += $5} END {print ls/NR" "ss/NR" "ms/NR" "rs/NR}' >> "${RDIR}/${NAME}.avg"
    done

gnuplot <<EOF
set terminal eps font 'Verdana,14'
set output "${RDIR}/${NAME}.eps"

set style data histograms
set style histogram rowstacked
set style fill solid 1.0 border -1
set boxwidth 1 relative

set title "${NAME} parallel scaling"
set ylabel "Time (seconds)"
set xlabel "#-Workers"

plot "${RDIR}/${NAME}.avg" using 2 title "Load", \
  '' using 3 title "Sort", \
  '' using 4 title "Map", \
  '' using 5:xticlabel(1) title "Reduce"
EOF

  done
fi

