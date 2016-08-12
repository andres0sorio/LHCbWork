#!/bin/bash
#

if [ $# -lt 4 ]
    then
    echo "usage:: $0 <ni> <nf> <name> <queue>"
    echo "name is the base name of the job e.g. gen_xxxx"
    exit 1	
fi

MIN=$1
MAX=$2
TYPE=$3
QUEUE=$4

for i in `seq $MIN $MAX`;
  do
  
  JOB=$TYPE'_'$i.job
  echo $JOB
  bsub -q $QUEUE $JOB

  sleep 5
  
done

###################

echo "submitJobs.sh> done!"

