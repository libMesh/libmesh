#!/bin/bash
#set -e

# This script is designed to be run as part of "make installcheck"
# but you can also run it manually from the source tree if you have a
# libMesh installed in $LIBMESH_DIR.

# You can run the script on a single header file by doing:
# test_CXXFLAGS="`$LIBMESH_DIR/bin/libmesh-config --cppflags --cxxflags --include`" HEADERS_TO_TEST=exact_solution.h ./contrib/bin/test_installed_headers.sh

# To run this script on *every* header file in an installed libMesh:
# test_CXXFLAGS="`$LIBMESH_DIR/bin/libmesh-config --cppflags --cxxflags --include`" HEADERS_TO_TEST="`find $LIBMESH_DIR/include/libmesh -name "*.h" -type f -exec basename {} \;`" ./contrib/bin/test_installed_headers.sh

# Respect the JOBS environment variable, if it is set
if [ -n "$JOBS" ]; then
    n_concurrent=$JOBS
else
    n_concurrent=20
fi

#echo MAKEFLAGS=$MAKEFLAGS

# Terminal commands to goto specific columns
rescol=65;

# Terminal commands for setting the color
gotocolumn=;
white=;
green=;
red=;
grey=;
colorreset=;
if (test "X$TERM" != Xdumb && { test -t 1; } 2>/dev/null); then
  gotocolumn="\033["$rescol"G";
  white="\033[01;37m";
  green="\033[01;32m";
  red="\033[01;31m";
  grey="\033[00;37m";
  colorreset="\033[m"; # Terminal command to reset to terminal default
fi


#echo "CXX=$CXX"

testing_installed_tree="no"

if (test "x$test_CXXFLAGS" = "x"); then

    testing_installed_tree="yes"

    if (test "x$PKG_CONFIG" != "xno"); then
        test_CXXFLAGS=`pkg-config libmesh --cflags`

    elif (test -x $LIBMESH_CONFIG_PATH/libmesh-config); then
        test_CXXFLAGS=`$LIBMESH_CONFIG_PATH/libmesh-config --cppflags --cxxflags --include`

    else
        echo "Cannot query package installation!!"
        exit 1
    fi
fi


# this function handles the I/O and compiling of a particular header file
# by encapsulating this in a function we can fork it and run multiple builds
# simultaneously
function test_header()
{
    myreturn=0
    header_to_test=$1
    header_name=`basename $header_to_test`
    app_file=`mktemp -t $header_name.XXXXXXXXXX`
    source_file=$app_file.cxx
    object_file=$app_file.o
    errlog=$app_file.log
    stdout=$app_file.stdout

    printf '%s' "Testing Header $header_to_test ... " > $stdout
    echo "#include \"libmesh/$header_name\"" >> $source_file
    echo "int foo () { return 0; }" >> $source_file

    #echo $CXX $test_CXXFLAGS $source_file -o $app_file
    if $CXX $test_CXXFLAGS $source_file -c -o $object_file >$errlog 2>&1 ; then
        # See color codes above.  We:
        # .) skip to column 65
        # .) print [ in white
        # .) print OK in green
        # .) print ] in white
        # .) reset the terminal color
        # .) print a newline
        printf '\e[65G\e[1;37m[\e[1;32m%s\e[1;37m]\e[m\e[m\n' "   OK   " >> $stdout
    else
        # See comment above for OK status
        printf '\e[65G\e[1;37m[\e[1;31m%s\e[1;37m]\e[m\e[m\n' " FAILED " >> $stdout
        echo "Source file:" >> $stdout
        cat $source_file  >> $stdout
        echo ""  >> $stdout
        echo "Command line:" >> $stdout
        echo $CXX $test_CXXFLAGS $source_file -c -o $object_file  >> $stdout
        echo ""  >> $stdout
        echo "Output:" >> $stdout
        cat $errlog >> $stdout
        echo "" >> $stdout
        myreturn=1
    fi

    cat $stdout
    rm -f $source_file $app_file $object_file $errlog $stdout

    return $myreturn
}


if [ "x$HEADERS_TO_TEST" = "x" ]; then
    HEADERS_TO_TEST=$DEFAULT_HEADERS_TO_TEST
fi


# loop over each header and fork tests
returnval=0
nrunning=0
runninglist=""
for header_to_test in $HEADERS_TO_TEST ; do

    # skip the files that live in contrib that we are installing when testing
    # from the source tree - paths will not be correct
    if (test "x$testing_installed_tree" = "xno"); then
        if (test "x`dirname $header_to_test`" = "xcontrib"); then
            continue
        fi
    fi

    if [ $nrunning -lt $n_concurrent ]; then
        test_header $header_to_test &
        runninglist="$runninglist $!"
        nrunning=$(($nrunning+1))
    else
        for pid in $runninglist ; do
            wait $pid
            # accumulate the number of failed tests
            returnval=$(($returnval+$?))
        done
        nrunning=0
        runninglist=""
    fi
done

wait

exit $returnval
