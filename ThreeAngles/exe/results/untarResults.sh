#!/bin/sh
# Andres Osorio - aooliver@ph.ed.ac.uk

if [ $# -lt 4 ]
    then
    echo "usage:: $0 <ni> <nf> <opt> <dir>"
    exit 1
fi

MIN=$1
MAX=$2
OPT=$3
DIR=$4

BASENAME="results_"$OPT

cd $DIR

for i in `seq $MIN $MAX`;
  do

  file=$BASENAME'_'$i.tar
  echo $file
  tar xvf $file
  
done

rm initialSO*

cd ../

echo "untarResults>  done!"

