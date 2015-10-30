#!/bin/bash
# The purpose of this script is to automatically test .C files
# for unnecessary includes by attempting to compile the source
# with specific headers removed.  The expected behavior
# is for the compilation to fail, so if it succeeds that may
# mean you can remove a header from the .C to speed up compilation
# times.  Be sure to double-check or at least know what you are
# doing before just automatically removing header files.
# author: John W. Peterson, 2005

print_usage () {
printf '%s\n' "Usage: $0 filename.C"
}

# Test input args
if test $# -ne 1; then
    print_usage;
    exit;
fi

# Script should be run from the contrib/bin directory so
# we know where libmesh-config is.
if [ ! `echo $PWD | grep "contrib/bin"` ]; then
  echo "Please run this script from the contrib/bin directory."
  exit;
fi

########################################################
# Terminal commands to goto specific columns
rescol=65;
gotocolumn="\033["$rescol"G";
white="\033[01;37m";
green="\033[01;32m";
red="\033[01;31m";
colorreset="\033[m";
########################################################

# Informational message re: which file we are testing
printf '%s\n' "Testing header removal in $1";

# The base filename
fn=`basename $1`;

# The base filename without the .C extension
bn=`echo $fn | cut -d"." -f1`;

# Find line numbers with local include statements.
# Ignore lines for libmesh_config.h and libmesh_common.h,
# assume they are definitely needed if they are present.
inc_line_nums=`cat $1 | grep -n "^#include \"" | grep -v libmesh_config | grep -v libmesh_common | cut -d":" -f1`;

# This script requires that you have libmesh installed in $LIBMESH_DIR
# with a valid libmesh-config script.
LIBMESH_CONFIG=${LIBMESH_DIR}/bin/libmesh-config

# Save the compiler name determined by libmesh-config
compiler_name=`$LIBMESH_CONFIG --cxx`

# Save the compile flags determined by libmesh-config
compile_flags=`$LIBMESH_CONFIG --cppflags --cxxflags --include`

for i in $inc_line_nums; do
  # Grab the i'th line from the file
  line_i=`head -$i $1 | tail -1`;

  printf '%s' "Testing with ($line_i) removed ... ";
  
  # Name for temporary file to test line deletion
  temp_file=`printf '%s' "line_${i}_$bn"`;
  
  # Create a temporary file with line $i deleted
  cat $1 | sed "${i}d" > $temp_file.C;

  # Attempt to compile $temp_file
  $compiler_name $compile_flags -c $temp_file.C &> /dev/null

  # If an object file is successfully created, the compilation succeeds!
  if [ -f $temp_file.o ]; then
      printf '\e[65G\e[1;37m[\e[1;32m%s\e[1;37m]\e[m\e[m\n' "   OK   "
      rm $temp_file.o
  else
      printf '\e[65G\e[1;37m[\e[1;31m%s\e[1;37m]\e[m\e[m\n' " FAILED "
  fi

  # Clean up temporary file(s)
  rm $temp_file.C
done;

