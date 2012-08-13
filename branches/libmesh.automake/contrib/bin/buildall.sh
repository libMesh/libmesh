#!/bin/bash

#./bootstrap


if (test "x$AEROLAB_SYSTEM_CLASS" = "x"); then
    echo "ERROR:  expecting a AEROLAB_SYSTEM_CLASS env var!"
    exit 1
fi

if (test "x$METHODS" = "x"); then
    METHODS="opt devel dbg"
fi
echo "building methods \"$METHODS\""

top_dir=`pwd`
install_dir=/lustre/work/benkirk/codes/install/$AEROLAB_SYSTEM_CLASS

if (test "x$MPI_ID_STRING" != "x"); then
    install_dir="$install_dir-$MPI_ID_STRING"
fi
if (test "x$COMPILER_ID_STRING" != "x"); then
    install_dir="$install_dir-$COMPILER_ID_STRING"
fi
echo "installing in $install_dir"

#rm -rf $install_dir

for METHOD in $METHODS ; do

    cd $top_dir
    builddir=$METHOD-$AEROLAB_SYSTEM_CLASS
    rm -rf $builddir && mkdir $builddir && cd $builddir

    echo " "
    echo " Configuring METHOD=$METHOD via ../configure --with-cxx=`which mpicxx` --with-cc=`which mpicc` --with-f77=`which mpif77` --with-fc=`which mpif90` --prefix=$install_dir --disable-dependency-tracking $@"
    echo " "

    ../configure --with-cxx=`which mpicxx` --with-cc=`which mpicc` --with-f77=`which mpif77` --with-fc=`which mpif90` --prefix=$install_dir --disable-dependency-tracking $@
    make --no-print-directory -j 8 && make --no-print-directory -j 8 install
done
    