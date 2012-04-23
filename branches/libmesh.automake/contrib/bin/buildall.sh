#!/bin/bash

./bootstrap

top_dir=`pwd`
install_dir=/lustre/work/benkirk/codes/install/aerolab_workstations

rm -rf $install_dir

for METHOD in dbg devel opt ; do

    cd $top_dir
    rm -rf $METHOD && mkdir $METHOD && cd $METHOD

    echo " "
    echo " Configuring via ../configure --prefix=$install_dir $@"
    echo " "

    ../configure --prefix=$install_dir $@
    make -j 24 install
done
    