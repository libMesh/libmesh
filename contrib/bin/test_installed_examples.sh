#!/bin/bash
#set -e
#env

# Respect the JOBS environment variable, if it is set
if [ -n "$JOBS" ]; then
    n_concurrent=$JOBS
else
    n_concurrent=5
fi

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


echo "Testing examples in $examples_install_path"

if (test "x$PKG_CONFIG" != "xno"); then
    installed_CXXFLAGS=`pkg-config libmesh --cflags`
    installed_LIBS=`pkg-config libmesh --libs`

elif (test -x $LIBMESH_CONFIG_PATH/libmesh-config); then
    installed_CXXFLAGS=`$LIBMESH_CONFIG_PATH/libmesh-config --cppflags --cxxflags --include`
    installed_LIBS=`$LIBMESH_CONFIG_PATH/libmesh-config --libs`

else
    echo "Cannot query package installation!!"
    exit 1
fi



#echo "installed_CXXFLAGS=$installed_CXXFLAGS"
#echo "installed_LIBS=$installed_LIBS"

# this function handles the I/O and compiling of a particular example.
# by encapsulating this in a function we can fork it and run multiple builds
# simultaneously
function test_example()
{
    myreturn=0
    app_to_link="ex_app"
    errlog=$app_to_link.log

    stdout=`mktemp -t stdout.XXXXXXXXXX`

    echo -n "Testing example installation $exdir ... " > $stdout
    cd $examples_install_path/$exdir

    if $CXX $installed_CXXFLAGS *.C -o $app_to_link $installed_LIBS > $errlog 2>&1 ; then

 	echo -e $gotocolumn $white"["$green"   OK   "$white"]" >> $stdout
	echo -e -n $colorreset >> $stdout

    else
 	echo -e $gotocolumn $white"["$red" FAILED "$white"]" >> $stdout
	echo -e -n $colorreset >> $stdout
	echo "Compile line:"  >> $stdout
	echo $CXX $installed_CXXFLAGS *.C -o $app_to_link $installed_LIBS >> $stdout
	echo ""  >> $stdout
	echo "Output:"  >> $stdout
	cat $errlog  >> $stdout
	echo ""  >> $stdout
	myreturn=1
    fi

    cat $stdout
    rm -f $app_to_link $errlog $stdout

    return $myreturn
}


# loop over each example and fork tests

cd $examples_install_path

returnval=0
nrunning=0
runninglist=""

for exdir in */* ; do

    if [ $nrunning -lt $n_concurrent ]; then
        test_example $exdir &
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
