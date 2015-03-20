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

./degree_sequence $GRAPH -o "${SEQ_FILE}.tmp" > /dev/null
#split $GRAPH $PREFIX -d --additional-suffix '.net.tmp' -n "l/$WORKERS" \
#  --filter "egrep $'^[0-9]+[ \\t][0-9]+[\\r]?$' > \$FILE"

mv "${SEQ_FILE}.tmp" $SEQ_FILE
#find $DIR -regex "${PREFIX}..\.net\.tmp" | while read FILE; do mv $FILE ${FILE%.tmp}; done

END=$(date +%s%N)
ELAPSED=$(echo "scale=8; ($END - $BEG) / 1000000000" | bc)
echo "Sorted in $ELAPSED seconds."
