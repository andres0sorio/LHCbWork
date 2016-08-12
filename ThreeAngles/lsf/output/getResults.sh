#!/bin/bash
## Get results from stepOne of the fit
## they are stored in a single file (unfortunately with constant name)
## Andres Osorio - aooliver@ph.ed.ac.uk

if [ $# -lt 1 ]
    then
    echo "usage:: $0 <log file>"
    exit 1
fi

INFILE=$1
LIM=`wc -l < $INFILE`

var='1'

while [ "$var" -le "$LIM" ]
  do
  LINE=`sed -e $var!d $INFILE`
  NJOB=`echo $LINE`
  dir=${NJOB%%res*out}

  cd $dir
  rootfile=`ls *.root`
  echo $rootfile
  noext=${rootfile%.*root}
  opt=${noext:11}
  cp resultsStepOne.out ../results/resultsSO'_'$opt.out
  cd ../

  var=$(($var+1))

done


