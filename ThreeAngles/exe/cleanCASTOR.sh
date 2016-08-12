#!/bin/sh
# Andres Osorio - aooliver@ph.ed.ac.uk

if [ $# -lt 5 ]
    then
    echo "usage:: $0 <ni> <nf> <file> <ext> <dir>"
    exit 1
fi

MIN=$1
MAX=$2
BASE=$3
EXT=$4
CASTORDIR=$5

for i in `seq $MIN $MAX`;
  do
  
  FILE=$BASE'_'$i.$EXT
  
  rfrm $CASTOR_HOME/$CASTORDIR/$FILE
  
done

rfdir $CASTOR_HOME/$CASTORDIR

echo "Removed files from CASTOR: done!"
