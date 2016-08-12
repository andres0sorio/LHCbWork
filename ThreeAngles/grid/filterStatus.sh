#!/bin/bash
## Match job script (jdl) with a Job Status
## Andres Osorio - aooliver@ph.ed.ac.uk
## uses finalResults.sh 

if [ $# -lt 2 ]
    then
    echo "usage: $0 <run summary> <status>"
    echo "run summary: cat all the log files (in asc. order) from submission"
    echo "-> file contains jdls and jobid|"
    echo "status: one of possible status (Ready, Cleared, ...)"
    exit 1
fi

summary=$1
status=$2
######################################################

./finalResults.sh $summary gdrb01 | grep "$status" | awk '{print $1}'
./finalResults.sh $summary rb107 | grep "$status" | awk '{print $1}'
