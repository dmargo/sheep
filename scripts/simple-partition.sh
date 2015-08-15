#!/bin/bash

JTREE_HOME=${JTREE_HOME:-$(pwd)}
USE_INOTIFY=${USE_INOTIFY:-$(command -v inotifywait > /dev/null)$?}
VERBOSE=${VERBOSE:-''}

GRAPH=${GRAPH:-${1:-'data/hep-th.dat'}}
DIR=${DIR:-$(dirname $GRAPH)}
PREFIX=${PREFIX:-${GRAPH%.net}}

shift 1
PARTS=${PARTS:-${@:-2}}



cd $JTREE_HOME

USE_SEQ=$( [ $SEQ_FILE != '-' ] && echo "-s $SEQ_FILE" || echo '' )
if [ "$OUT_FILE" != '' ] && [ "$PARTS" != '0' ]; then
  echo 'Using fast partition path...'
  ./graph2tree $GRAPH $USE_SEQ -o $OUT_FILE -p $PARTS $VERBOSE
  echo "Reduced in 0.0 seconds."
else
  ./graph2tree $GRAPH $USE_SEQ -o "${PREFIX}.tre" $VERBOSE
  echo "Reduced in 0.0 seconds"
  source scripts/part-worker.sh
fi

