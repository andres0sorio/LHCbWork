# script setting up ROOT
# on the CE where the job is submitted
# Andres Osorio - aooliver@ph.ed.ac.uk


LCGEXT=/scratch/lhcb-soft/lcg/external

export ROOTSYS=${LCGEXT}/root/5.10.00c/slc3_ia32_gcc323/root
export PATH=${ROOTSYS}/bin:${PATH}

ISLDPATH=`/bin/env | /bin/grep "LD_LIBRARY_PATH"`

if [ "$ISLDPATH" == "" ];
   then
   export LD_LIBRARY_PATH=${ROOTSYS}/lib
   else
   export LD_LIBRARY_PATH=${ROOTSYS}/lib:${LD_LIBRARY_PATH}
fi
