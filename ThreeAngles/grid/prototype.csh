#!/bin/csh
#######################################################
# ML fit studies
# Andres Osorio - aooliver@ph.ed.ac.uk


if ($#argv == 0 ) then
  echo "Usage: $0 <executable> <option>"
  exit 0
else
  set EXEC = $1

endif

#Important: set ROOT
chmod +x findROOT
source setROOT.csh

#####################################
set NEVT='130000'
set DATA=

# == run the executable ================================
chmod +x ${EXEC}

./${EXEC} ${DATA} ${NEVT}

echo "${0}: Done"

exit 0

