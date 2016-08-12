#!/bin/sh
# prepareGen.sh generates n scripts (jdl_csh) )to be submitted to a batch
# The scripts will run the generation of MCToy data
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

INPUT=prototypeGen.csh

BASENAME="inputParameters_"$OPT

for i in `seq $MIN $MAX`;
  do
  
  OUTPUT=./gen'_'$OPT'_'$i.csh
  FILE=$BASENAME'_'$i.dat
  
  sed -e "s/set DATA=/set DATA=$FILE/
          s/set OPTION=/set OPTION=$OPT/" $INPUT > $OUTPUT 
  chmod +x $OUTPUT
  
done

###################
INPUT=prototype.jdl

for i in `seq $MIN $MAX`;
  do
  
  OUTPUT=./gen'_'$OPT'_'$i.jdl

  script=gen'_'$OPT'_'$i.csh
  parfile=$BASENAME'_'$i.dat
  results=results'_'$OPT.tar

  sed -e "s/Executable =/Executable = \"$script\";/;
          s/Arguments =/Arguments = \"$EXE\";/;
          s/OutputSandbox =/OutputSandbox = \{\"std.out\",\"std.err\",\"$results\"\};/;
          s/InputSandbox =/InputSandbox = \{\"$script\",\"$data\",\"$EXE\",\"setROOT.csh\",\"findROOT\",\"$parfile\"\};/" $INPUT > $OUTPUT

  chmod +x $OUTPUT
  
done


echo "prepareGen> done!"

