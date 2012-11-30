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

returnval=0

app_to_link="ex_app"
errlog=$app_to_link.log

#echo " "
#echo "$app_to_link"
#echo "$errlog"
#echo " "

cd $examples_install_path
for exdir in */* ; do

    echo -n "Testing example installation $exdir ... "
    cd $examples_install_path/$exdir
    
    if $CXX $installed_CXXFLAGS *.C -o $app_to_link $installed_LIBS > $errlog 2>&1 ; then
	 
 	echo -e $gotocolumn $white"["$green"   OK   "$white"]";
	echo -e -n $colorreset;
	
    else
 	echo -e $gotocolumn $white"["$red" FAILED "$white"]";
	echo -e -n $colorreset;
	echo "Compile line:"
	echo $CXX $installed_CXXFLAGS *.C -o $app_to_link $installed_LIBS
	echo ""
	echo "Output:"
	cat $errlog
	echo ""
	returnval=1	
    fi
    
    rm -f $app_to_link $errlog    
done

exit $returnval
