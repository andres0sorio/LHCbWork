# script setting up ROOT
# on the CE where the job is submitted
# Andres Osorio - aooliver@ph.ed.ac.uk

##set LCGEXT=/opt/exp_soft/lhcb/lib/lcg/external
##set LCGEXT=/gridstorage/exptsw/lhcb/lib/lcg/external
##setenv ROOTSYS ${LCGEXT}/root/5.11.02/slc3_ia32_gcc323/root

echo "searching for: 5.12.00"
set rootsys=`./findROOT 5.12.00`
echo ${rootsys}

setenv ROOTSYS ${rootsys}

setenv PATH ${ROOTSYS}/bin:${PATH}

set ISLDPATH = `/bin/env | /bin/grep "LD_LIBRARY_PATH"`

if ( "$ISLDPATH" == "" ) then
   setenv LD_LIBRARY_PATH ${ROOTSYS}/lib
else
   setenv LD_LIBRARY_PATH ${ROOTSYS}/lib:${LD_LIBRARY_PATH}
endif

echo "defined variables:"
echo "ROOTSYS: ${ROOTSYS}"
echo "PATH: ${PATH}"
echo "LD_PATH: ${LD_LIBRARY_PATH}"
echo "setROOT.csh> finished successfuly"
