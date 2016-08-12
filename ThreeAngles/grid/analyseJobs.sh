#!/bin/bash

if [ $# -lt 1 ]
    then
    echo "usage:: $0 <log file>"
    echo "where <log file> contains only job ids (use getJobIds.sh)"
    exit 1
fi

INFILE=$1
TMPF1=temp1.out
grep "gdrb01" $INFILE > $TMPF1
LIM1=`wc -l < $TMPF1`

TMPF2=temp2.out
grep "rb107" $INFILE > $TMPF2
LIM2=`wc -l < $TMPF2`

echo "Jobs sent to gdrb01" > jobStatus_gdrb01.txt 
var='1'

while [ "$var" -le "$LIM1" ]
  do

  LINE=`sed -e $var!d $TMPF1`
  NJOB=`echo $LINE`
  edg-job-status $NJOB >> jobStatus_gdrb01.txt
  var=$(($var+1))

done

echo "Jobs sent to rb107" > jobStatus_rb107.txt
var='1'

##############################

while [ "$var" -le "$LIM2" ]
  do

  LINE=`sed -e $var!d $TMPF2`
  NJOB=`echo $LINE`
  edg-job-status $NJOB >> jobStatus_rb107.txt
  var=$(($var+1))

done

rm $TMPF1 $TMPF2

##check for status

codes[1]="W"
codes[2]="r"
codes[3]="s"
codes[4]="S"
codes[5]="R"
codes[6]="D"
codes[7]="A"
codes[8]="C"

col1[1]=`grep "Waiting" jobStatus_gdrb01.txt | wc -l`
col1[2]=`grep "Ready" jobStatus_gdrb01.txt | wc -l`
col1[3]=`grep "Scheduled" jobStatus_gdrb01.txt | wc -l`
col1[4]=`grep "Submitted" jobStatus_gdrb01.txt | wc -l`
col1[5]=`grep "Running" jobStatus_gdrb01.txt | wc -l`
col1[6]=`grep "Done" jobStatus_gdrb01.txt | wc -l`
col1[7]=`grep "Aborted" jobStatus_gdrb01.txt | wc -l`
col1[8]=`grep "Cleared" jobStatus_gdrb01.txt | wc -l`

col2[1]=`grep "Waiting" jobStatus_rb107.txt | wc -l`
col2[2]=`grep "Ready" jobStatus_rb107.txt | wc -l`
col2[3]=`grep "Scheduled" jobStatus_rb107.txt | wc -l`
col2[4]=`grep "Submitted" jobStatus_rb107.txt | wc -l`
col2[5]=`grep "Running" jobStatus_rb107.txt | wc -l`
col2[6]=`grep "Done" jobStatus_rb107.txt | wc -l`
col2[7]=`grep "Aborted" jobStatus_rb107.txt | wc -l`
col2[8]=`grep "Cleared" jobStatus_rb107.txt | wc -l`

for i in `seq 1 8`;
  do
  echo "${codes[$i]}     ${col1[$i]}     ${col2[$i]}"
done
