#!/bin/bash

# XXX: MUST BE DEFINED
#USE_INOTIFY
#VERBOSE

#GRAPH
#DIR
#PREFIX
#PARTS

#REDUCTION
#WORKERS
#ID_NUM
ID_NUM=${ID_NUM:-$1}

if [ $ID_NUM -eq 0 ]; then
  BEG=$(date +%s%N)
fi



# MAP
source scripts/map-worker.sh



# REDUCE
STEP=0
STEP_SIZE=$WORKERS
WORKERS=$(( ($WORKERS + $REDUCTION - 1) / $REDUCTION ))
while [ $STEP_SIZE -ne 1 ] && [ $ID_NUM -lt $WORKERS ]; do

  source scripts/reduce-worker.sh

  STEP=$(( $STEP + 1 ))
  STEP_SIZE=$WORKERS
  WORKERS=$(( ($WORKERS + $REDUCTION - 1) / $REDUCTION ))
done

if [ $ID_NUM -eq 0 ]; then
  mv "${PREFIX}00r${STEP}.tre" "${PREFIX}.tre"

  END=$(date +%s%N)
  ELAPSED=$(echo "scale=8; ($END - $BEG) / 1000000000" | bc)
  echo "Mapped in $ELAPSED seconds."
  echo "Reduced in 0.0 seconds."

  # PARTITION
  source scripts/part-worker.sh
fi

