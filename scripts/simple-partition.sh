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

./graph2tree $GRAPH -o "${PREFIX}.tre" $VERBOSE

echo "Reduced in 0.0 seconds."

SEQ_FILE='-'
source scripts/part-worker.sh

