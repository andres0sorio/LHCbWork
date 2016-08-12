#!/bin/sh
# Andres Osorio - aooliver@ph.ed.ac.uk

if [ $# -lt 3 ]
    then
    echo "usage:: $0 <ni> <nf> <opt>"
    exit 1
fi

MIN=$1
MAX=$2
OPT=$3

INPUT=inputParameters'_'$OPT'_'1.dat

BASENAME="inputParameters_"$OPT

for i in `seq $MIN $MAX`;
  do
  
  target=$BASENAME'_'$i.dat
  cp $INPUT $target
  
done

echo "prepareGen> input files done!"

