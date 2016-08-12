#!/bin/sh
# the script generates n scripts (jdl+csh)to be submitted to a batch
# system - like the GRID
# Andres Osorio - aooliver@ph.ed.ac.uk

if [ $# -lt 3 ]
    then
    echo "usage:: $0 <ni> <nf> <opt>"
    exit 1
fi

MIN=$1
MAX=$2
OPT=$3
EXE=stepOne

INPUT=prototype.csh

BASENAME="phis_study_"$OPT

for i in `seq $MIN $MAX`;
  do
  
  OUTPUT=./fit'_'$OPT'_'$i.csh
  FILE=$BASENAME'_'$i.root

  sed -e "s/set DATA=/set DATA=$FILE/" $INPUT > $OUTPUT
  
  chmod +x $OUTPUT

done

###################
INPUT=prototype.jdl

for i in `seq $MIN $MAX`;
  do
  
  OUTPUT=./fit'_'$OPT'_'$i.jdl

  script=fit'_'$OPT'_'$i.csh
  data=$BASENAME'_'$i.root
  results=resultsStepOne.out

  sed -e "s/Executable =/Executable = \"$script\";/;
          s/Arguments =/Arguments = \"$EXE\";/;
          s/OutputSandbox =/OutputSandbox = \{\"std.out\",\"std.err\",\"$results\"\};/;
          s/InputSandbox =/InputSandbox = \{\"$script\",\"$data\",\"$EXE\",\"setROOT.csh\",\"findROOT.sh\"};/" $INPUT > $OUTPUT

  chmod +x $OUTPUT
  
done





echo "prepareRun> done!"

