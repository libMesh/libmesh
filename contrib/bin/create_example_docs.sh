#!/bin/sh

# This shell script is a wrapper for the script ex2html.sh
# also located in this directory.  This script calls ex2html.sh
# for all of the example programs.

if test "$LIBMESH_DIR" = ""; then
  export LIBMESH_DIR=../..
fi

cd "$LIBMESH_DIR"/examples/ || exit
exampleslist=$(find . -type d -name "*ex*")
echo "$exampleslist"
cd ../contrib/bin/ || exit

for exdir in $exampleslist; do
    mybasename=$(basename $exdir)
    echo "Processing example $exdir, file $LIBMESH_DIR/examples/$exdir/$mybasename.C";
    ./ex2html.sh $LIBMESH_DIR/examples/$exdir/$mybasename.C "make --no-print-directory -C $LIBMESH_DIR/examples/$exdir run"
done
