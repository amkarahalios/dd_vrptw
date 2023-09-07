#!/bin/bash

export LD_LIBRARY_PATH="$HOME/dd_vrptw/Cliquer/src/lib/:$LD_LIBRARY_PATH"

if [ $# -ne 3 ]; then
  echo "Invalid arguments"
  echo "Argument 1: instances dir"
  echo "Argument 2: COL_GEN or COL_ELIM"
  echo "Argument 3: param file name"
  exit 2
fi

OUTPUT_DIR=$HOME/dd_vrptw/logs/$1_$3/
mkdir -p $OUTPUT_DIR

FILES=$HOME/dd_vrptw/instances/$1/*
PARAM_FILE=$HOME/dd_vrptw/param_files/$3.csv
for f in $FILES
do
  echo "Processing $f file..."

  OUTPUT_FILE=$OUTPUT_DIR/$(basename $f).log
  if [ -f "$OUTPUT_FILE" ]; then
    echo "$OUTPUT_FILE does exist."
  else
    echo "Running for $OUTPUT_FILE"
    timeout 3600 .././solver $f $2 $PARAM_FILE > $OUTPUT_FILE
  fi
done
