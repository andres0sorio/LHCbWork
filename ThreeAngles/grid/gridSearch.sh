#!/bin/bash
#

if [ $# -lt 1 ]
    then
    echo "usage:: $0 <jdl-script>"
    exit 1
fi

SCRIPT=$1

OUTPUT=job_submission.log

date > $OUTPUT

for i in `seq 1 100`;
  do

  edg-job-submit --vo lhcb $SCRIPT >> $OUTPUT
  
done

###################

echo "gridSeach> done!"

