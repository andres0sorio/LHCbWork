#!/bin/sh
# prepareGen.sh generates n scripts (jdl_csh) )to be submitted to a batch
# The scripts will run the generation of MCToy data
# system - like the GRID
# Andres Osorio - aooliver@ph.ed.ac.uk

if [ $# -lt 4 ]
    then
    echo "usage:: $0 <ni> <nf> <opt> <nevt>"
    exit 1
fi

MIN=$1
MAX=$2
OPT=$3
EVT=$4

INPUT=phis_study'_'$OPT

INITIAL=initialSO'_'$OPT

for i in `seq $MIN $MAX`;
  do

  paramfile=$INITIAL'_'$i.dat

  if [ ! -f $paramfile ]
      then
      cp initialSO_xxT1_0.dat $paramfile
  fi
  
  ./stepOne $INPUT'_'$i.root $EVT
  
done

echo "runAnalysis> done!"

