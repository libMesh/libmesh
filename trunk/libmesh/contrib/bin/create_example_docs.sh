#!/bin/bash
#
# This shell script is a wrapper for the script ex2html.sh
# also located in this directory.  This script calls ex2html.sh
# for all of the example programs.

if [ "x$LIBMESH_DIR" = x ]; then
  export LIBMESH_DIR=../..
fi

cd $LIBMESH_DIR/examples/
exampleslist=`ls -d ex*/ */ex*/ | sed 's#/$##g'`
echo $exampleslist
cd ../contrib/bin/

for exdir in $exampleslist; do
    mybasename=`basename $exdir`
    echo "Processing example $exdir";
    ./ex2html.sh $LIBMESH_DIR/examples/$exdir/$mybasename.C "make --no-print-directory -C $LIBMESH_DIR/examples/$exdir run"
done
