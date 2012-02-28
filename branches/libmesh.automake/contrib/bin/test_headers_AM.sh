#!/bin/bash
set -e

# Terminal commands to goto specific columns
rescol=65;


# [JWP] My Mac did not recognize the \e commands below, but it did
# recognize \033 commands, and they seem to work on Linux as well.  

# Terminal commands for setting the color
# gotocolumn="\e["$rescol"G";
# white="\e[01;37m";
# green="\e[01;32m";
# red="\e[01;31m";
# #grey="\e[00;37m";
# colorreset="\e[m"; # Terminal command to reset to terminal default


gotocolumn="\033["$rescol"G";
white="\033[01;37m";
green="\033[01;32m";
red="\033[01;31m";
#grey="\033[00;37m";
colorreset="\033[m";

echo "CXX=$CXX"
echo "INCLUDEDIR_BASE=$INCLUDEDIR"

include_I_path=""

if (test "x$TBB_DIR" != x); then
    if (test -d $TBB_DIR/include); then
	include_I_path="-I$TBB_DIR/include $include_I_path"
    fi
fi


for dir in $INCLUDEDIR/* ; do
   echo $dir
   if (test -d $INCLUDEDIR); then
     include_I_path="-I$dir $include_I_path"
   fi
done

echo "include_I_path=$include_I_path"
errlog=test_headers.log
returnval=0
for header_to_test in $HEADERS_TO_TEST ; do
    echo -n "Testing Header $header_to_test ... "

    header_name=`basename $header_to_test`
    app_file=`mktemp -t $header_name.XXXXXXXXXX`
    source_file=$app_file.cxx
    object_file=$app_file.o

    echo "#include \"$header_name\"" >> $source_file
    echo "int foo () { return 0; }" >> $source_file

    #echo $CXX $include_I_path $source_file -o $app_file
    if $CXX $include_I_path $source_file -c -o $object_file >$errlog 2>&1 ; then
 	echo -e $gotocolumn $white"["$green"   OK   "$white"]";
	echo -e -n $colorreset;    
    else
 	echo -e $gotocolumn $white"["$red" FAILED "$white"]";
	echo -e -n $colorreset;    
	cat $errlog
	returnval=1
    fi

    #echo -e -n $colorreset;    
    #cat $source_file
    rm -f $source_file $app_file $object_file
done

exit $returnval

# # Note: script was /bin/sh but I don't think -n, -e are POSIX
# # arguments to echo

# if [ "x$LIBMESH_DIR" = x ]; then
#   export LIBMESH_DIR=../..
# fi

# headers_to_test=`ls $LIBMESH_DIR/include/*/*.h`

# if test $# -ge 1; then
#     headers_to_test=$*
# fi


# # Errors during the tests will be printed here
# errlog=test_headers.log

# returnval=0

# for i in $headers_to_test; do
#     header_name=`basename $i`
#     source_name=TestHeader_$header_name.C # Use .C here, take advantage of our make rule
#     source_file=$LIBMESH_DIR/contrib/bin/$source_name
#     app_name=TestHeader_$header_name.o
#     app_file=$LIBMESH_DIR/contrib/bin/$app_name

#     rm -f $source_file $app_file
    
#     echo "#include \"$header_name\"" >> $source_file
#     echo "int main () { return 0; }" >> $source_file

#     echo -n "Testing Header File $header_name ... "

#     # Use make -s for silent operation
#     if make -s -C $LIBMESH_DIR contrib/bin/$app_name 2>> $errlog; then
# 	echo -e $gotocolumn $white"["$green"   OK   "$white"]";
#     else
# 	echo -e $gotocolumn $white"["$red" FAILED "$white"]";
#         returnval=1
#     fi
    
#     echo -e -n $colorreset;    
#     rm -f $source_file $app_file
# done

# exit $returnval
