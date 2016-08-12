#!/bin/bash

rm test.out

valgrind --run-libc-freeres=no --skin=memcheck --leak-check=yes --show-reachable=yes ./simdetAnal pythia_eemmmm_E.sim.gz 1000 > test.out 2>&1



echo "Done." 
