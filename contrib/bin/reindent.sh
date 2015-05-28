#!/bin/bash

# This script reindents a file using emacs-style indentation,
# with no indentation in namespaces.
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 file.C"
    exit 1
fi

# You can specify what Emacs to use whe you run this script with e.g.
#
# EMACS=/Applications/Emacs.app/Contents/MacOS/Emacs reindent.sh foo.C
#
# otherwise /usr/bin/emacs will be used.  The /usr/bin/emacs on OSX is
# very old (circa 2007 22.1.1 on Yosemite in 2015) and there are some
# bugs in the way that it indents C++ code, so we definitely recommend
# using something newer if possible.
if [ -z $EMACS ]; then
  EMACS=/usr/bin/emacs
fi

# Print file we are working on
echo "Indenting $1"

# The following command:
# .) Uses only spaces for tabs
# .) Forces emacs to use C++ mode
# .) Sets the indentation level within namespaces to 0
# .) Runs indent-region on the entire file
# .) Saves the buffer
$EMACS -batch $1  \
  --eval="(setq-default indent-tabs-mode nil)" \
  --eval="(c++-mode)" \
  --eval="(c-set-offset 'innamespace 0)" \
  --eval="(indent-region (point-min) (point-max) nil)" \
  -f save-buffer &> /dev/null
