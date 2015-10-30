#!/bin/bash

# Note: script was /bin/sh but I don't think -n, -e are POSIX
# arguments to echo

if [ "x$LIBMESH_DIR" = x ]; then
  export LIBMESH_DIR=../..
fi

headers_to_test=`ls $LIBMESH_DIR/include/*/*.h`

if test $# -ge 1; then
    headers_to_test=$*
fi

# Terminal commands to goto specific columns
rescol=65;
gotocolumn="\e["$rescol"G";
#testnamecolumn=15;
#gotoname="\e["$testnamecolumn"G";

# Terminal commands for setting the color
white="\e[01;37m";
green="\e[01;32m";
red="\e[01;31m";
#grey="\e[00;37m";

# Terminal command to reset to terminal default
colorreset="\e[m";

# Errors during the tests will be printed here
errlog=test_headers.log

for i in $headers_to_test; do
    header_name=`basename $i`
    source_name=TestHeader_$header_name.C # Use .C here, take advantage of our make rule
    source_file=$LIBMESH_DIR/contrib/bin/$source_name
    app_name=TestHeader_$header_name.o
    app_file=$LIBMESH_DIR/contrib/bin/$app_name

    rm -f $source_file $app_file
    
    echo "#include \"$header_name\"" >> $source_file
    echo "int main () { return 0; }" >> $source_file

    echo -n "Testing Header File $header_name ... "

    # Use make -s for silent operation
    if make -s -C $LIBMESH_DIR contrib/bin/$app_name 2>> $errlog; then
	echo -e $gotocolumn $white"["$green"   OK   "$white"]";
    else
	echo -e $gotocolumn $white"["$red" FAILED "$white"]";
    fi
    
    echo -e -n $colorreset;    
    rm -f $source_file $app_file
done
