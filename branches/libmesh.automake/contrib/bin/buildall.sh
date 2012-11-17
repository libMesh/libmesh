#!/bin/bash

#./bootstrap

top_dir=`pwd`


#
# Figure out the installation prefix
#
if (test "x$AEROLAB_SYSTEM_CLASS" = "x"); then
    install_dir=$HOME/codes/install
    AEROLAB_SYSTEM_CLASS="local"
else
    install_dir=/lustre/work/benkirk/codes/install/$AEROLAB_SYSTEM_CLASS
fi

if (test "x$MPI_ID_STRING" != "x"); then
    install_dir="$install_dir-$MPI_ID_STRING"
fi
if (test "x$COMPILER_ID_STRING" != "x"); then
    install_dir="$install_dir-$COMPILER_ID_STRING"
fi
echo "installing in $install_dir"

#rm -rf $install_dir



#
# honor METHODS, but set a default
#  e.g. could use 
#  $ METHODS="opt devel" ./buildall.sh 
#
if (test "x$METHODS" = "x"); then
    METHODS="opt devel dbg"
fi
echo "building methods \"$METHODS\""


for METHOD in $METHODS ; do

    cd $top_dir
    builddir=$METHOD-$AEROLAB_SYSTEM_CLASS
    rm -rf $builddir && mkdir $builddir && cd $builddir

    echo " "
    echo " Configuring METHOD=$METHOD via ../configure --with-cxx=`which mpicxx` --with-cc=`which mpicc` --with-f77=`which mpif77` --with-fc=`which mpif90` --prefix=$install_dir --disable-dependency-tracking $@"
    echo " "

    ../configure --with-cxx=`which mpicxx` --with-cc=`which mpicc` --with-f77=`which mpif77` --with-fc=`which mpif90` --prefix=$install_dir --disable-dependency-tracking $@
    make --no-print-directory -j6 && make --no-print-directory -j6 install || exit 1
done
