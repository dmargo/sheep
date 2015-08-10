#!/bin/bash

# INVARIANTS
TRUE=0
FALSE=1

# These are currently invariants, but may become options.
export USE_INOTIFY=$(command -v inotifywait > /dev/null)$?
export REDUCTION=2


# OPTIONS
USE_SLURM=$FALSE
JTREE_HOME=$(pwd)
TRIALS=1

USE_VERTICAL=$FALSE
USE_MPI_SORT=$FALSE
USE_MPI_REDUCE=$FALSE
KEEP_DATA=$FALSE

export VERBOSE=''
export SEQ_FILE='-'
export OUT_FILE=''
INITIAL_WORKERS=2

while getopts "lh:t:airkvs:o:w:c:" opt; do
  case $opt in
    l)
      USE_SLURM=$TRUE;;
    h)
      JTREE_HOME=$OPTARG;;
    t)
      TRIALS=$OPTARG;;
    a)
      USE_VERTICAL=$TRUE;;
    i)
      USE_MPI_SORT=$TRUE;;
    r)
      USE_MPI_REDUCE=$TRUE;;
    k)
      KEEP_DATA=$TRUE;;
    v)
      export VERBOSE='-v';;
    s)
      export SEQ_FILE=$OPTARG;;
    o)
      export OUT_FILE=$OPTARG;;
    w)
      INITIAL_WORKERS=$OPTARG;;
    c)
      CORES=$OPTARG;;
    :)
      echo "Option -$OPTARG requires an argument."
      exit 1;;
    \?)
      echo "Invalid option: -$OPTARG"
      exit 1;;
  esac
done

export CORES=${CORES:-$INITIAL_WORKERS}

if [ $USE_SLURM -eq $TRUE ]; then
  DEFAULT_GRAPH='/n/regal/seltzer_lab/data/as20graph/as20graph.dat'
  RUN='srun -n 1'
else
  DEFAULT_GRAPH='data/hep-th.dat'
  RUN=''
fi


# ARGUMENTS
shift $(( $OPTIND - 1 ))
export GRAPH=${1:-$DEFAULT_GRAPH}
shift 1
export PARTS=${@:-2}

if [ $USE_SLURM -eq $FALSE ] && [ ! -f $GRAPH ]; then
  echo "$GRAPH does not exist."
  exit 1
fi

echo "Starting dist-partition on $GRAPH with $INITIAL_WORKERS workers..."
echo "s:$USE_SLURM a:$USE_VERTICAL i:$USE_MPI_SORT r:$USE_MPI_REDUCE c:$CORES"
if [ $USE_SLURM -eq $TRUE ]; then
  echo "p:$SBATCH_PARTITION nodes:$SLURM_JOB_NUM_NODES $SLURM_TASKS_PER_NODE"
fi


#ENVIRONMENT SETUP
cd $JTREE_HOME

BASEDIR=$(dirname $GRAPH)

if [ $USE_SLURM -eq $TRUE ]; then
  if [ $SLURM_JOB_NUM_NODES -eq 1 ]; then
    SBCP='cp -f -v'
  else
    SBCP='sbcast -f -v'
  fi

  TMP_GRAPH="/scratch/$(basename $GRAPH)"
  $SBCP $GRAPH $TMP_GRAPH
  if [ -f "${GRAPH}.ini" ]; then
    $SBCP "${GRAPH}.ini" "${TMP_GRAPH}.ini"
  fi
  export GRAPH=$TMP_GRAPH
fi


# DIRTY DEEDS DONE DIRT CHEAP
for t in $(seq $TRIALS); do
  export DIR="$BASEDIR/$(date +%s%N)"
  export PREFIX="$DIR/$(basename $GRAPH .dat)"
  mkdir $DIR

  export WORKERS=$INITIAL_WORKERS
  if [ $WORKERS -eq 1 ]; then
    source scripts/simple-partition.sh
  elif [ $USE_VERTICAL -eq $TRUE ]; then
    source scripts/vertical-dist.sh
  else
    source scripts/horizontal-dist.sh
  fi

  # CLEANUP
  if [ $KEEP_DATA -eq $FALSE ]; then
    rm -rf $DIR
  fi
done
if [ $USE_SLURM -eq $TRUE ]; then
  rm -rf $TMP_GRAPH
fi
