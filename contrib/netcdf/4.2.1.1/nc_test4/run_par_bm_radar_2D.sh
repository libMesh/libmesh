#!/bin/sh

# This shell file runs benchmarks on the 2D radar data on parallel platforms.

# $Id: run_par_bm_radar_2D.sh,v 1.1 2007/12/18 01:16:25 ed Exp $

set -e
echo ""
echo "Getting radar 2D data file from Unidata FTP site..."
file=20070803-2300_tile1-2d.nc3
if ! test -f $file; then
    wget ftp://ftp.unidata.ucar.edu/pub/netcdf/sample_data/$file
fi

echo "*** Running bm_file for parallel access on $file..."
header="-h"
chunksizes="1501:2001"
for numproc in 1 4 16
  do
  mpiexec -n $numproc ./bm_file -p -d ${header} -s 16 -f 4 -o tst_r2d.nc -c 0:-1:0:1501:2001 $file
  header=
done
echo '*** SUCCESS!!!'

exit 0