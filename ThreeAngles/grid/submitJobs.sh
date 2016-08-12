#!/bin/bash
# script to submit CPStudies jobs to the grid
# Andres Osorio - aooliver@ph.ed.ac.uk
# the destinations are read from a file

if [ $# -lt 5 ]
    then
    echo "usage: $0 <ni> <nf> <opt 1> <opt 2> <destination list>"
    exit 1	
fi

MIN=$1
MAX=$2
OPT1=$3
OPT2=$4
DESTLIST=$5

st=`date '+%j%H%M'`
logfile=Jobsub'_'$st.log
echo "jobs submitted on the:" > $logfile
date >> $logfile

dests=(`cat $DESTLIST`)
max=`echo ${#dests[*]}`
var='0'
max=$(($max-1))

for i in `seq $MIN $MAX`;
  do
  
  job=$OPT1'_'$OPT2'_'$i.jdl
  echo "submitJobs> $job" >> $logfile
  edg-job-submit --vo lhcb -r ${dests[$var]} $job >> $logfile
  sleep 1
  
  var=$(($var+1))
  if [ "$var" -gt "$max" ]
      then
      var='0'
  fi
  
done

###################

echo "submitJobs.sh> done!"

