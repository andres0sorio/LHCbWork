#!/bin/sh
# Andres Osorio - aooliver@ph.ed.ac.uk

if [ $# -lt 4 ]
    then
    echo "usage:: $0 <ni> <nf> <opt> <nevt>"
    exit 1
fi

MIN=$1
MAX=$2
OPT=$3
JUMP=$4
EVT='0'

INPUT=phis_study'_'$OPT.root

rm fitting_steps.log

touch fitting_steps.log

for i in `seq $MIN $MAX`;
  do

  EVT=$(($JUMP + $JUMP * ($i-1) ))

  ./stepOne $INPUT $EVT >> fitting_steps.log

done

echo "runAnalysis> done!"

