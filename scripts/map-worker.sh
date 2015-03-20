#!/bin/bash

# XXX: MUST BE DEFINED
#USE_INOTIFY
#VERBOSE

#GRAPH
#DIR
#PREFIX

#WORKERS
#ID_NUM
ID_NUM=${ID_NUM:-$1}
printf -v ID_STR '%02d' $ID_NUM

if [ "$VERBOSE" = "-v" ]; then
  echo "MAP: $(hostname)"
fi



# MAP
#INPUT_FILE="${PREFIX}${ID_STR}.net"
#while [ ! -f $INPUT_FILE ]; do
#  [ $USE_INOTIFY -eq 0 ] && inotifywait -qqt 1 -e create -e moved_to $DIR || sleep 1
#done

while [ ! -f $SEQ_FILE ]; do
  [ $USE_INOTIFY -eq 0 ] && inotifywait -qqt 1 -e create -e moved_to $DIR || sleep 1
done

OUTPUT_FILE="${PREFIX}${ID_STR}"
./graph2tree $GRAPH -l "$(( $ID_NUM + 1 ))/$WORKERS" -s $SEQ_FILE -o $OUTPUT_FILE $VERBOSE
mv $OUTPUT_FILE "${OUTPUT_FILE}r0.tre"
