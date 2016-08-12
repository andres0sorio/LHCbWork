#!/bin/bash

for iexp in `seq 0 199`;
do

 if [ $iexp -lt "10" ] 
  then 
  datafile=main2d'_000'$iexp.dat
  else 
   if [ $iexp -lt "100" ]
     then
     datafile=main2d'_00'$iexp.dat
     else   
       if [ $iexp -lt "1000" ]
       then
       datafile=main2d'_0'$iexp.dat
       fi
    fi
  fi


./fitData $datafile initialfitparams.dat 150000

done

echo "fit to data completed"


