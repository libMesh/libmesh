#!/bin/bash

# This script reindents a file using emacs-style indentation,
# with no indentation in namespaces.
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 file.C"
    exit 1
fi

# The following command:
# .) Uses only spaces for tabs
# .) Forces emacs to use C++ mode
# .) Sets the indentation level within namespaces to 0
# .) Runs indent-region on the entire file
# .) Saves the buffer
/usr/bin/emacs -batch $1  \
  --eval="(setq-default indent-tabs-mode nil)" \
  --eval="(c++-mode)" \
  --eval="(c-set-offset 'innamespace 0)" \
  --eval="(indent-region (point-min) (point-max) nil)" \
  -f save-buffer
