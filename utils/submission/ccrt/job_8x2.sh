#!/bin/bash
#MSUB -r dumses_test_run_8x2         # Job name
#MSUB -n 8                           # Total number of proc. used
#MSUB -c 2                           # Number of core per proc. to use
#MSUB -T 600                         # Elapsed time limit in s.
#MSUB -o dumses_test_run_8x2.%I.o    # Standard output
#MSUB -o dumses_test_run_8x2.%I.e    # Error output
#MSUB -q standard                    # Queue name

##MSUB -A paxxx                      # Project ID

set -x

cd ${BRIDGE_MSUB_PWD}
mkdir test_run_8x2.${BRIDGE_MSUB_JOBID}
cd test_run_8x2.${BRIDGE_MSUB_JOBID}
cp $CCCWORKDIR/tests_dumses/magnetic_loop/input_8 input
cp $CCCWORKDIR/tests_dumses/magnetic_loop/dumses .

export OMP_NUM_THREADS=2

ccc_mprun ./dumses
