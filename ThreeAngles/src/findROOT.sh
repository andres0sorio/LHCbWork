#!/bin/sh
# this script tries to find the ROOT installation
# on the CE where the job is submitted
# Andres Osorio - aooliver@ph.ed.ac.uk
# returns ROOTSYS

#if lhcb env then look inside this directory for a ROOT installation

if [ $# -lt 1 ]
    then
    echo "usage:: $0 <version>"
    exit 1
fi

version=$1
config=$CMTCONFIG

#test at CERN
if [ "$CERN" == "/cern" ]
    then 
    VO_LHCB_SW_DIR=/afs/cern.ch/sw
    LCG_EXT=$VO_LHCB_SW_DIR/lcg/external
else
    ##on the GRID
    LCG_EXT=$VO_LHCB_SW_DIR/lib/lcg/external
fi


rootdir=`ls $LCG_EXT/root`

if [ "$rootdir" == "" ];
    then
    echo "root dir not found"
    exit 0
else
    for i in $rootdir
      do
      if [ "$i" == "$version" ];
	  then
	  found='1'
      fi
    done	
fi

if [ $found ];
    then
    echo $LCG_EXT/root/$version/$config/root
fi
