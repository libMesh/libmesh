#!/bin/sh

# This shell runs some parallel tests.

# $Id: run_par_tests.sh,v 1.2 2007/12/20 16:25:26 ed Exp $

# Even for successful runs, mpiexec seems to set a non-zero return
# code!  
#set -e
echo ""
echo "Testing parallel I/O with HDF5..."

mpiexec -n 1 ./tst_h_par
mpiexec -n 2 ./tst_h_par
mpiexec -n 4 ./tst_h_par
echo "SUCCESS!!!"

exit 0
