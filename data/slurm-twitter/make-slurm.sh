#!/bin/bash
OUT_DATA='./slurm-25.avg'
rm -f $OUT_DATA

RAW_DATA=( ./slurm-25-*.out )
for RAW in ${RAW_DATA[@]}; do
  NAME=$(basename $RAW .out)
  WORKERS=$(echo $NAME | egrep -o '[[:digit:]]+$' )

  egrep "^Loaded graph[[:blank:]]" $RAW | egrep -o "[[:digit:]]*\.[[:digit:]]+" > "/tmp/${NAME}.load"
  egrep "^Sorted[[:blank:]]" $RAW | egrep -o "[[:digit:]]*\.[[:digit:]]+" > "/tmp/${NAME}.sort"
  egrep "^Mapped[[:blank:]]" $RAW | egrep -o "[[:digit:]]*\.[[:digit:]]+" > "/tmp/${NAME}.map"
  egrep "^Reduced[[:blank:]]" $RAW | egrep -o "[[:digit:]]*\.[[:digit:]]+" > "/tmp/${NAME}.red"

  paste /tmp/${NAME}.load /tmp/${NAME}.sort /tmp/${NAME}.map /tmp/${NAME}.red > ${NAME}.dat
  rm    /tmp/${NAME}.load /tmp/${NAME}.sort /tmp/${NAME}.map /tmp/${NAME}.red

  echo -n "$WORKERS " >> $OUT_DATA
  awk '{ls += $1; ss += $2; ms += $3; rs += $4} END {print ls/NR" "ss/NR" "ms/NR" "rs/NR" "(ls+ss+ms+rs)/NR}' ${NAME}.dat >> $OUT_DATA
done

