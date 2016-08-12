#!/bin/sh
# prepareGen.sh generates n scripts (csh) )to be submitted to a batch queue
# The scripts will run the fit to MCToy data
# Andres Osorio - aooliver@ph.ed.ac.uk

if [ $# -lt 3 ]
    then
    echo "usage:: $0 <ni> <nf> <opt>"
    exit 1
fi

MIN=$1
MAX=$2
OPT=$3

INPUT=prototypeFit.job



for i in `seq $MIN $MAX`;
  do
  
  OUTPUT=./fit'_'$OPT'_'$i.job
  
  INDATA=results'_'$OPT'_'$i.tar
  
  sed -e "s/set INPUT=/set INPUT=$INDATA/
          s/set OPTION=/set OPTION=$OPT/
          s/set FITID=/set FITID=$i/" $INPUT > $OUTPUT 

  chmod +x $OUTPUT
  
done

echo "prepareGen> done!"

