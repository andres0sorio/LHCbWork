#!/bin/sh
# prepareRun.sh generates n scripts (jdl_csh) )to be submitted to a batch
# system - like the GRID
# Andres Osorio - aooliver@ph.ed.ac.uk

if [ $# -lt 4 ]
    then
    echo "usage:: $0 <ni> <nf> <opt> <exec>"
    exit 1
fi

MIN=$1
MAX=$2
OPT=$3
EXE=$4

INPUT=prototest.csh

BASENAME="initialSO_test_"$OPT

for i in `seq $MIN $MAX`;
  do
  
  OUTPUT=./test'_'$OPT'_'$i.csh
  FILE=$BASENAME'_'$i.dat
  
  sed -e "s/set DATA=/set DATA=$FILE/" $INPUT > $OUTPUT
  
  chmod +x $OUTPUT
  
done

###################
INPUT=prototype.jdl

for i in `seq $MIN $MAX`;
  do
  
  OUTPUT=./test'_'$OPT'_'$i.jdl

  script=test'_'$OPT'_'$i.csh
  data=phis_study_xxRt_0.root
  parfile=$BASENAME'_'$i.dat
  results=resultsStepOne.out

  sed -e "s/Executable =/Executable = \"$script\";/;
          s/Arguments =/Arguments = \"$EXE\";/;
          s/OutputSandbox =/OutputSandbox = \{\"std.out\",\"std.err\",\"$results\"\};/;
          s/InputSandbox =/InputSandbox = \{\"$script\",\"$data\",\"$EXE\",\"setROOT.csh\",\"findROOT.sh\",\"$parfile\",\"$results\"\};/" $INPUT > $OUTPUT

  chmod +x $OUTPUT
  
done





echo "prepareRun> done!"

