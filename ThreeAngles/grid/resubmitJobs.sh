#!/bin/bash
# script to resubmit CPStudies jobs to the grid
# for instance when jobs are aboerted due to a expired proxy 
# Andres Osorio - aooliver@ph.ed.ac.uk

if [ $# -lt 2 ]
    then
    echo "usage: $0 <jdl list> <destination list>"
    exit 1	
fi

jdllist=$1
destlist=$2

st=`date '+%j%H%M'`
jobs=(`cat $jdllist`)
len=`echo ${#jobs[*]}` 

dests=(`cat $destlist`)
max=`echo ${#dests[*]}`

logfile=Jobsub'_'$st.log
echo "jobs submitted on the:" > $logfile
date >> $logfile

var='0'
max=$(($max-1))
len=$(($len-1))

for i in `seq 0 $len`;
  do

  echo "submitJobs> ${jobs[$i]}" >> $logfile
  echo edg-job-submit --vo lhcb -r ${dests[$var]} ${jobs[$i]} >> $logfile
#  sleep 1

  var=$(($var+1))

  if [ "$var" -gt "$max" ]
	then
	var='0'
  fi
  
done

###################

echo "resubmitJobs.sh> done!"

