#!/bin/bash

# XXX: MUST BE DEFINED
#USE_INOTIFY
#VERBOSE

#DIR
#PREFIX

#STEP
#STEP_SIZE
#WORKERS
#ID_NUM
ID_NUM=${ID_NUM:-$1}
printf -v ID_STR '%02d' $ID_NUM

if [ "$VERBOSE" = "-v" ]; then
  echo "REDUCE: $(hostname)"
fi



# REDUCE
INPUT_LIST=$( seq -f "${PREFIX}%02gr${STEP}.tre" -s ' ' $ID_NUM $WORKERS $(( $STEP_SIZE - 1 )) )

INPUT_ARRAY=($INPUT_LIST)
for INPUT_FILE in ${INPUT_ARRAY[*]}; do
  while [ ! -f $INPUT_FILE ]; do
    [ $USE_INOTIFY -eq 0 ] && inotifywait -qqt 1 -e create -e moved_to $DIR || sleep 1
  done
done

OUTPUT_FILE="${PREFIX}${ID_STR}r$(( $STEP + 1 )).tre"

if [ ${#INPUT_ARRAY[@]} -eq 1 ]; then
  mv $INPUT_LIST $OUTPUT_FILE
else
  ./merge_trees $INPUT_LIST -o "${OUTPUT_FILE}.tmp" $VERBOSE
  mv "${OUTPUT_FILE}.tmp" $OUTPUT_FILE
fi
