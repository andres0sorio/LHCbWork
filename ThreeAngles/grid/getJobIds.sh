#!/bin/bash
## getJobIds.sh: extract and list from the main log file the job IDs
## Andres Osorio - aooliver@ph.ed.ac.uk

if [ $# -lt 1 ]
    then
    echo "usage:: $0 <log file>"
    exit 1
fi

INFILE=$1
TMPFILE=temp.out
grep "https" $INFILE > $TMPFILE
LIMIT=`wc -l < $TMPFILE`

var='1'

while [ "$var" -le "$LIMIT" ]
  do

  LINE=`sed -e $var!d $TMPFILE`
  NLINE=`echo $LINE`
  ENVAR=${NLINE:2}
  echo $ENVAR
  var=$(($var+1))

done

rm $TMPFILE
