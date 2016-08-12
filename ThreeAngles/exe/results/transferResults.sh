#!/bin/sh
# Andres Osorio - aooliver@ph.ed.ac.uk

if [ $# -lt 4 ]
    then
    echo "usage:: $0 <ni> <nf> <opt1> <opt2>"
    exit 1
fi

MIN=$1
MAX=$2
OPT=$3
OPT2=$4

DIR=$OPT'_'$OPT2

BASENAME="results_"$OPT

for i in `seq $MIN $MAX`;
  do

   if [ ! -d $DIR ]
    then
    mkdir $DIR
   fi
  
   file=$BASENAME'_'$i.tar
   echo $file
   rfcp $CASTOR_HOME/ToyMCfits/$file ./$DIR/
  
done

echo "transferResults>  done!"

