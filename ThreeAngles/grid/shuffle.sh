#!/bin/bash
# shuffle function - takes a list as argument and outputs the list at random order (shuffle)
# Andres Osorio - aooliver@ph.ed.ac.uk

if [ $# -lt 1 ];
    then
    echo "usage: $0 <infile>"
    exit 0
fi

init=(`cat $1`)
output=( '' )
val='0'
max=`echo ${#init[*]}`

shuffle ()
{
    local temp
    temp=( `echo "$1"` )
    #pass an array as argument use its length as max range
    range=`echo ${#temp[*]}`
    number=$RANDOM
    let "number %= $range"
    elem=`echo $number`
    output[$val]=${temp[$elem]}
    unset temp[$elem]
    nsize=`echo ${#output[*]}`
    if [ "$nsize" ==  "$max" ]
	then
	val=$(($val+1))
	output[$val]=${temp[0]}
	return 1
    else
	val=$(($val+1))
	argument=`echo ${temp[@]}`
	shuffle "$argument"
    fi
    
}

#echo "initial = ${init[@]}"
argument=`echo ${init[@]}`
shuffle "$argument"

for i in `seq 0 $(($max-1))`;
  do
  echo ${output[$i]}
done
