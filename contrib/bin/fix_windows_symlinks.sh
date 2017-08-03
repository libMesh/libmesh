#!/bin/sh


LIBMESH_ROOT=$(shell git rev-parse --show-toplevel)
SYMLINKS=$(git ls-files -s "$LIBMESH_ROOT" | awk '/120000/{print $4}')

for sl in $SYMLINKS
do
    TARGET=$(dirname $sl)/$(cat $sl)
    echo -n "Replacing symlink $sl by a copy of $TARGET... "

    # Do not track these changes in git
    git update-index --assume-unchanged "$sl"

    # Remove symlink and copy
    rm "$sl"
    cp -r "$TARGET" $sl

    echo "done"
done
