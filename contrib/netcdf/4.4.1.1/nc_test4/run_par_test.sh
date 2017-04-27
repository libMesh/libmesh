#!/bin/sh

# This shell runs some parallel I/O tests.

set -e
echo
echo "Testing MPI parallel I/O with various other mode flags..."
mpiexec -n 1 ./tst_mode
echo
echo "Testing MPI parallel I/O without netCDF..."
mpiexec -n 4 ./tst_mpi_parallel
echo
echo "Testing very simple parallel I/O with 4 processors..."
mpiexec -n 4 ./tst_parallel
echo
echo "Testing simple parallel I/O with 16 processors..."
mpiexec -n 16 ./tst_parallel3
echo
echo "num_proc   time(s)  write_rate(B/s)"
mpiexec -n 1 ./tst_parallel4
mpiexec -n 2 ./tst_parallel4
mpiexec -n 4 ./tst_parallel4
mpiexec -n 8 ./tst_parallel4

#mpiexec -n 16 ./tst_parallel4
#mpiexec -n 32 ./tst_parallel4
#mpiexec -n 64 ./tst_parallel4
echo
echo "Parallel Performance Test for NASA"
mpiexec -n 4 ./tst_nc4perf

echo
echo "Parallel I/O test for Collective I/O, contributed by HDF Group."
mpiexec -n 1 ./tst_simplerw_coll_r
mpiexec -n 2 ./tst_simplerw_coll_r
