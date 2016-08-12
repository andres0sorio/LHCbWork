#!/bin/bash
## Make a nice table with job id, status, CE
## Andres Osorio - aooliver@ph.ed.ac.uk
##

if [ $# -lt 1 ]
    then
    echo "usage:: $0 <log file>"
    echo "log file> jobStatus_[gdrb01,rb107].txt"
    exit 1
fi

INFILE="$1"

######################################################
# gdrb01x

jobid=(`grep "Status info" $INFILE | awk '{print $7}'`)
status=(`grep "Current Status" $INFILE | awk '{print $3}'`)
place=(`grep "Destination" $INFILE | awk '{print $2}'`)

len1=`echo ${#jobid[*]}`
len2=`echo ${#status[*]}`
len3=`echo ${#place[*]}`

var='0'

if [ "$len1" -ne "$len2" -o "$len1" -ne "$len3" ]
   then 
   exit 1
fi

while [ "$var" -lt "$len1" ]
do

echo ${jobid[$var]} ${status[$var]} ${place[$var]}

var=$(($var+1))
done

