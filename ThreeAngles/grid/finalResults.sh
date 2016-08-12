#!/bin/bash
## Match job script (jdl) with the jobID and Job Status
## Andres Osorio - aooliver@ph.ed.ac.uk
## using matchJobId.sh and makeTable.sh

if [ $# -lt 2 ]
    then
    echo "usage:: $0 <log file> <opt>"
    echo "log file: full log file" 
    echo "opt: gdrb01 or rb107"
    exit 1
fi

INFILE=$1
OPT=$2
######################################################

statfile=jobStatus'_'$OPT.txt
jdls=(`./matchJobId.sh $INFILE | grep "$OPT" | awk '{print $1}'`)
jobid=(`./matchJobId.sh $INFILE | grep "$OPT" | awk '{print $2}'`)
status=(`./makeTable.sh $statfile | awk '{print $2}'`)

len1=`echo ${#jobid[*]}`
len2=`echo ${#jdls[*]}`
len3=`echo ${#status[*]}`

var='0'

if [ "$len1" -ne "$len2" -o "$len1" -ne "$len3" ]
    then 
    echo "matchJobId> incompatible sizes!"
    echo $len1 $len2
    exit 1
fi

while [ "$var" -lt "$len1" ]
  do
  
  echo ${jdls[$var]} ${jobid[$var]} ${status[$var]}
  
  var=$(($var+1))
  
done

