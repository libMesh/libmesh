#!/bin/bash

# This script detabs (and reindents) a file using emacs-style indentation.
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 file.C"
    exit 1
fi

# You can specify what Emacs to use whe you run this script with e.g.
#
# EMACS=/Applications/Emacs.app/Contents/MacOS/Emacs detab.sh foo.C
#
# otherwise /usr/bin/emacs will be used.  The /usr/bin/emacs on OSX is
# very old (circa 2007 22.1.1 on Yosemite in 2015) and there are some
# bugs in the way that it indents C++ code, so we definitely recommend
# using something newer if possible.
if [ -z $EMACS ]; then
  EMACS=/usr/bin/emacs
fi

# Use perl to get a list of line numbers that have <TAB> characters
tabbed_lines=(`perl -ne "print \"$. \" if /\t/" $1`)

# Now actually use perl to remove those tab characters
perl -pli -e 's/\t//g' $1

# How much work is there to be done?
n=${#tabbed_lines[@]}

if [ $n -gt 0 ]; then
  # Print status message
  echo "Reindenting $n lines in $1."

  # Initialize the list of eval commands with the command which causes
  # Emacs to use spaces for tabs.
  eval_commands="--eval=\"(setq-default indent-tabs-mode nil)\" --eval=\"(c++-mode)\""

  # Make a long string of --eval commands
  for i in "${tabbed_lines[@]}"
  do
      eval_commands+=" --eval \"(goto-line $i)\" --eval \"(indent-for-tab-command)\""
  done

  # Construct the command you want to run as a string, and then use
  # 'eval' to run it.  Note that this is *not* the same thing as
  # running, say, `$cmd`.
  cmd="$EMACS -batch $1 $eval_commands -f save-buffer &> /dev/null"
  eval $cmd
else
  echo "Nothing to be done in $1."
fi
