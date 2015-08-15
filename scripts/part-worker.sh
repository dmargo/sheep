#!/bin/bash

# XXX: MUST BE DEFINED
#USE_INOTIFY
#VERBOSE

#GRAPH
#DIR
#PREFIX
#PARTS

if [ "$PARTS" != 0 ]; then
  if [ "$VERBOSE" = "-v" ]; then
    echo "PARTITION: $(hostname)"
  fi

  INPUT_TREE="${PREFIX}.tre"
  while [ ! -f $INPUT_TREE ]; do
    [ $USE_INOTIFY -eq 0 ] && inotifywait -qqt 1 -e create -e moved_to $DIR || sleep 1
  done

  BEG=$(date +%s%N)

  if [ "$OUT_FILE" = '' ]; then
    ./partition_tree -f -g $GRAPH $SEQ_FILE $INPUT_TREE $PARTS
  else
    ./partition_tree -f -g $GRAPH $SEQ_FILE $INPUT_TREE $PARTS -o $OUT_FILE
  fi

  END=$(date +%s%N)
  ELAPSED=$(echo "scale=8; ($END - $BEG) / 1000000000" | bc)
  echo "Partitioned in $ELAPSED seconds."
fi

