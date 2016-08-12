#!/bin/bash
## Get results from stepOne of the fit
## they are stored in a single file (unfortunately all of them have the same name)
## requiring special attention
## Andres Osorio - aooliver@ph.ed.ac.uk

if [ $# -lt 1 ]
    then
    echo "usage:: $0 <opt>"
    exit 1
fi

tempf=temp.out
target=resultsStepOne.out

search=fit'_'$1
grep -r "$search" * | awk '{ print $1 }' > $tempf

LIM=`wc -l < $tempf`

var='1'

while [ "$var" -le "$LIM" ]
  do
  LINE=`sed -e $var!d $tempf`
  NJOB=`echo $LINE`
  dir=${NJOB%%/*csh:}
  dataset=${NJOB##aos*/}
  dataset=${dataset:4}
  dataset=${dataset%.*:}

  echo $dir
  ##echo $dataset  
  
  file=$dir/$target
  
  if [ -f $file ];
      then
      newfile=resultsSO'_'$dataset.out
      if [ -d results ];
	  then
	  cp $file results/$newfile
      else
	  mkdir results
	  cp $file results/$newfile
      fi
  else
      echo "getResults> $target not found"      
  fi
  
  var=$(($var+1))
  
done

rm $tempf