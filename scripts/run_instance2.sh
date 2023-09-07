#!/bin/bash

export LD_LIBRARY_PATH="$HOME/dd_vrptw/Cliquer/src/lib/:$LD_LIBRARY_PATH"

if [ $# -ne 9 ]; then
  echo "Invalid arguments"
  echo "Argument 1: instances dir"
  echo "Argument 2: COL_GEN or COL_ELIM"
  echo "Argument 3: DD/DP/LP/LAG"
  echo "Argument 4: Q/NG"
  echo "Argument 5: s/k value"
  echo "Argument 6: max s value"
  echo "Argument 7: use cuts"
  echo "Argument 8: timeout"
  echo "Argument 9: description"
  exit 2
fi

OUTPUT_DIR=$HOME/dd_vrptw/logs/$9/
META_FILE=$OUTPUT_DIR/meta.log
mkdir -p $OUTPUT_DIR
echo "Arguments $1, $2, $3, $4, $5, $6, $7, $8, $9" > $META_FILE

FILES=$HOME/dd_vrptw/$1/*
for f in $FILES
do
  echo "Processing $f file..."

  OUTPUT_FILE=$OUTPUT_DIR/$(basename $f).log
  if [ -f "$OUTPUT_FILE" ]; then
    echo "$OUTPUT_FILE does exist."
  else
    echo "Running for $OUTPUT_FILE"
    timeout $8 .././solver_nomussp $f $2 $3 $4 $5 $6 $7 $8 > $OUTPUT_FILE
  fi
done
