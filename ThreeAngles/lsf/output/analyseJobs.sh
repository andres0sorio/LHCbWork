#!/bin/bash

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
  echo $dir
 
  var=$(($var+1))
  

done


