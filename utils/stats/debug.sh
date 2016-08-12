#!/bin/bash

make clean

make

./statsReader cutSummary_NoCuts_Signal.out option

gv StudyOnCutsoption1.eps 

echo "Done." 
