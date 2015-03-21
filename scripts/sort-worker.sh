#!/bin/bash

# XXX: MUST BE DEFINED:
#VERBOSE

#GRAPH
#PREFIX

#WORKERS

if [ "$VERBOSE" = "-v" ]; then
  echo "SPLIT: $(hostname)"
fi



# SORT & SPLIT
BEG=$(date +%s%N)

./degree_sequence $GRAPH "${SEQ_FILE}.tmp" > /dev/null

mv "${SEQ_FILE}.tmp" $SEQ_FILE

END=$(date +%s%N)
ELAPSED=$(echo "scale=8; ($END - $BEG) / 1000000000" | bc)
echo "Sorted in $ELAPSED seconds."
