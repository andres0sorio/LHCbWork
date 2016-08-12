#!/bin/bash
# script to retrieve jobs from the grid
# Andres Osorio - aooliver@ph.ed.ac.uk
# provide a log file - after using analyseJobs.sh

if [ $# -lt 1 ]
    then
    echo "usage:: $0 <log file>"
    exit 1
fi

if [ ! -d output ]
   then
   echo "$0> output directory is missing. Cannot continue."
   exit 1
fi

INFILE="$1"

######################################################
jobid=(`grep "Status info" $INFILE | awk '{print $7}'`)
status=(`grep "Current Status" $INFILE | awk '{print $3}'`)

len1=`echo ${#jobid[*]}`
len2=`echo ${#status[*]}`

var='0'

if [ "$len1" -ne "$len2" ]
   then 
   echo "retrieveJobs> there are problems in your logfile"
   exit 1
fi

while [ "$var" -lt "$len1" ]
do

if [ "${status[$var]}" == "Done" ]
  then
  echo "retrieveJobs> ${jobid[$var]}"
  edg-job-get-output --dir output ${jobid[$var]}
fi


var=$(($var+1))
done

