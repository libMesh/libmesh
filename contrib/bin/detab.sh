#!/bin/bash

# This script detabs (and reindents) a file using emacs-style indentation.
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 file.C"
    exit 1
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
  cmd="/usr/bin/emacs -batch $1 $eval_commands -f save-buffer &> /dev/null"
  eval $cmd
else
  echo "Nothing to be done in $1."
fi
