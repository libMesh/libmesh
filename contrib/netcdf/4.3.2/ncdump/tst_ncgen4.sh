#!/bin/sh
verbose=1
set -e

# To add a new test,
# 1. put the .cdl file in the 'cdl4' directory
# 2. put the result of running ncgen then ncdump
#    into the directory 'expected4' as .dmp
# 3. Modify the file tst_ncgen4_shared.sh to add
#    the test to the end of the TESTS4 variable
# 4. Add the new files into cdl4/Makefile.am
#    and expected4/Makefile.am 

if test "x$builddir" = "x"; then builddir=`pwd`; fi
if test "x$srcdir" = "x"; then srcdir=`dirname $0`; fi

# Make buildir absolute
cd $builddir
builddir=`pwd`

# Make srcdir be absolute
cd $srcdir
srcdir=`pwd`
cd $builddir

export verbose
export srcdir
export builddir

KFLAG=1 ; export KFLAG
echo "*** Performing diff tests: k=1"
sh ${srcdir}/tst_ncgen4_diff.sh
echo "*** Performing cycle tests: k=1"
sh  ${srcdir}/tst_ncgen4_cycle.sh
KFLAG=3 ; export KFLAG
echo "*** Performing diff tests: k=3"
sh  ${srcdir}/tst_ncgen4_diff.sh
echo "*** Performing cycle tests: k=3"
sh  ${srcdir}/tst_ncgen4_cycle.sh
KFLAG=4 ; export KFLAG
echo "*** Performing diff tests: k=4"
sh  ${srcdir}/tst_ncgen4_diff.sh
echo "*** Performing cycle tests: k=4"
sh  ${srcdir}/tst_ncgen4_cycle.sh
exit

