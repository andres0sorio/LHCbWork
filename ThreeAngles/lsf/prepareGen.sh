#!/bin/sh
# prepareGen.sh generates n scripts (csh) )to be submitted to a batch queue
# The scripts will run the generation of MCToy data
# Andres Osorio - aooliver@ph.ed.ac.uk

if [ $# -lt 3 ]
    then
    echo "usage:: $0 <ni> <nf> <opt>"
    exit 1
fi

MIN=$1
MAX=$2
OPT=$3

INPUT=prototypeGen.job

BASENAME="inputParameters_"$OPT

for i in `seq $MIN $MAX`;
  do
  
  OUTPUT=./gen'_'$OPT'_'$i.job
  FILE=$BASENAME'_'$i.dat
  
  sed -e "s/set DATA=/set DATA=$FILE/
          s/set OPTION=/set OPTION=$OPT/
          s/set GENID=/set GENID=$i/" $INPUT > $OUTPUT 

  chmod +x $OUTPUT
  
done

echo "prepareGen> done!"

