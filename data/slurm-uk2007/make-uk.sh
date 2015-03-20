#!/bin/bash
OUTPUT=slurm-uk.summary
rm -f $OUTPUT

for FILE in $@; do
  PREFIX=$(echo $FILE | egrep -o '[[:digit:]]+(x[[:digit:]]+)?')
  echo -n "$PREFIX " >> $OUTPUT
  cat $FILE | awk '{s1+=$1;s2+=$2;s3+=$3;s4+=$4} END {s1=s1/NR;s2=s2/NR;s3=s3/NR;s4=s4/NR; print s1" "s2" "s3" "s4}' >> $OUTPUT
done
