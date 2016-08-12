#!/bin/bash
## Match job script (jdl) with the jobID
## Andres Osorio - aooliver@ph.ed.ac.uk
##

if [ $# -lt 1 ]
    then
    echo "usage:: $0 <log file>"
    echo "log file: full log file" 
    exit 1
fi

INFILE="$1"

######################################################

jobid=(`grep "https" $INFILE | awk '{print $2}'`)
jdls=(`grep "jdl" $INFILE | awk '{print $2}'`)

len1=`echo ${#jobid[*]}`
len2=`echo ${#jdls[*]}`

var='0'

if [ "$len1" -ne "$len2" ]
    then 
    echo "matchJobId> incompatible sizes!"
    echo $len1 $len2
    exit 1
fi

while [ "$var" -lt "$len1" ]
  do
  
  echo ${jdls[$var]} ${jobid[$var]} 
  
  var=$(($var+1))
  
done

