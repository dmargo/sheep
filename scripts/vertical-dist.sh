#!/bin/bash

# SETUP
if [ $SEQ_FILE = '-' ]; then
  export SEQ_FILE="${PREFIX}.seq"
  source scripts/sort-worker.sh
fi

# LAUNCH WORKERS
for ID_NUM in `seq 0 $(( $WORKERS - 1 ))`; do
  $RUN scripts/vertical-worker.sh $ID_NUM &
done
wait

