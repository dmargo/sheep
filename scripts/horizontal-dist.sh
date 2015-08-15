#!/bin/bash

# SETUP
if [ $SEQ_FILE = '-' ]; then
  export SEQ_FILE="${PREFIX}.seq"
  if [ $USE_MPI_SORT -eq $FALSE ]; then
    source scripts/sort-worker.sh
  fi
fi



# MAP
FAST_PART=$( [ $USE_MPI_REDUCE -eq $TRUE ] && [ "$OUT_FILE" != '' ] && [ "$PARTS" != 0 ] && \
  echo $TRUE || echo $FALSE )

if [ $USE_MPI_SORT -eq $FALSE ] && [ $USE_MPI_REDUCE -eq $FALSE ]; then
  echo "Loaded in 0.0 seconds."
  BEG=$(date +%s%N)

  for ID_NUM in $( seq 0 $(( $WORKERS - 1 )) ); do
    $RUN scripts/map-worker.sh $ID_NUM &
    if [ $(( ($ID_NUM + 1) % $CORES )) -eq 0 ]; then wait; fi
  done
  wait

  END=$(date +%s%N)
  ELAPSED=$(echo "scale=8; ($END - $BEG) / 1000000000" | bc)
  echo "Mapped in $ELAPSED seconds."
else
  MPI_SORT=$( [ $USE_MPI_SORT -eq $TRUE ] && echo '-i' || echo '')
  MPI_REDUCE=$( [ $USE_MPI_REDUCE -eq $TRUE ] && echo '-r' || echo '')
  if [ $FAST_PART -eq $TRUE ]; then
    echo 'Using fast partition path...'
    mpiexec -n $WORKERS ./graph2tree $GRAPH -s $SEQ_FILE -o $OUT_FILE -p $PARTS $MPI_SORT $MPI_REDUCE $VERBOSE
  else
    mpiexec -n $WORKERS ./graph2tree $GRAPH -s $SEQ_FILE -o $PREFIX $MPI_SORT $MPI_REDUCE $VERBOSE
  fi
fi



# REDUCE
if [ $USE_MPI_REDUCE -eq $FALSE ]; then
  BEG=$(date +%s%N)

  export STEP=0
  export STEP_SIZE=$WORKERS
  export WORKERS=$(( ($WORKERS + $REDUCTION - 1) / $REDUCTION ))
  while [ $STEP_SIZE -ne 1 ]; do
    for ID_NUM in $( seq 0 $(( $WORKERS - 1 )) ); do
      $RUN scripts/reduce-worker.sh $ID_NUM &
      if [ $(( ($ID_NUM + 1) % $CORES )) -eq 0 ]; then wait; fi
    done
    wait

    export STEP=$(( $STEP + 1 ))
    export STEP_SIZE=$WORKERS
    export WORKERS=$(( ($WORKERS + $REDUCTION - 1) / $REDUCTION ))
  done
  
  END=$(date +%s%N)
  ELAPSED=$(echo "scale=8; ($END - $BEG) / 1000000000" | bc)
  echo "Reduced in $ELAPSED seconds."
  mv "${PREFIX}00r${STEP}.tre" "${PREFIX}.tre"
elif [ $FAST_PART -eq $FALSE ]; then
  mv $PREFIX "${PREFIX}.tre"
fi



# PARTITION
if [ $FAST_PART -eq $FALSE ]; then
  source scripts/part-worker.sh
fi
