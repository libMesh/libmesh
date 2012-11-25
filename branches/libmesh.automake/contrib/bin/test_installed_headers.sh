#!/bin/bash
set -e
#env

# Terminal commands to goto specific columns
rescol=65;

# Terminal commands for setting the color
gotocolumn=;
white=;
green=;
red=;
grey=;
colorreset=;
if (test "X$TERM" != Xdumb); then
  gotocolumn="\033["$rescol"G";
  white="\033[01;37m";
  green="\033[01;32m";
  red="\033[01;31m";
  grey="\033[00;37m";
  colorreset="\033[m"; # Terminal command to reset to terminal default
fi

#echo "CXX=$CXX"

if (test "x$PKG_CONFIG" != "xno"); then
    installed_CXXFLAGS=`pkg-config libmesh --cflags`

elif (test -x $LIBMESH_CONFIG_PATH/libmesh-config); then
    installed_CXXFLAGS=`$LIBMESH_CONFIG_PATH/libmesh-config --cppflags --cxxflags --include`

else
    echo "Cannot query package installation!!"
    exit 1
fi



returnval=0
for header_to_test in $HEADERS_TO_TEST ; do
    echo -n "Testing Header $header_to_test ... "

    header_name=`basename $header_to_test`
    app_file=`mktemp -t $header_name.XXXXXXXXXX`
    source_file=$app_file.cxx
    object_file=$app_file.o
    errlog=$app_file.log
    
    echo "#include \"libmesh/$header_name\"" >> $source_file
    echo "int foo () { return 0; }" >> $source_file

    #echo $CXX $installed_CXXFLAGS $source_file -o $app_file
    if $CXX $installed_CXXFLAGS $source_file -c -o $object_file >$errlog 2>&1 ; then
 	echo -e $gotocolumn $white"["$green"   OK   "$white"]";
	echo -e -n $colorreset;    
    else
 	echo -e $gotocolumn $white"["$red" FAILED "$white"]";
	echo -e -n $colorreset;
	echo "Source file:"
	cat $source_file
	echo ""
	echo "Command line:"
	echo $CXX $installed_CXXFLAGS $source_file -c -o $object_file
	echo ""
	echo "Output:"
	cat $errlog
	echo ""
	returnval=1
    fi

    #echo -e -n $colorreset;    
    #cat $source_file
    rm -f $source_file $app_file $object_file $errlog
done

exit $returnval

