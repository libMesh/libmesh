#!/bin/sh

# LIBMESH_ROOT should be set to the location of the libmesh source tree
if test "$LIBMESH_ROOT" = ""; then
  export LIBMESH_ROOT=../..
fi

our_headers=$(find $LIBMESH_ROOT/include -name "*.h" -type f | sort | uniq)

for file in "$@" ; do
    echo "Processing file $file"
    for header_name_with_path in $our_headers ; do

        header_name=$(basename $header_name_with_path)

        # This grep command is not correct because it doesn't handle
        # the fact that "." is a regular expression matching any
        # single character. So, for example, when we search for
        # "mesh.h", it matches the following line which is just a
        # comment that contains the substring "mesh h":
        # // the elements in the underlying elements in the mesh have changed,
        count=$(grep $header_name $file | grep -v libmesh/)

        # Debugging
        # if test "$count" != ""; then
        #     echo "count = $count"
        # fi

        if test "$count" != ""; then
            string_to_replace="#include \"$header_name\""
            replacement_string="#include \"libmesh/$header_name\""

            search_regex="include \"$header_name\""
            replacement_regex="include \"libmesh\/$header_name\""

            echo "        Replacing "
            echo "          "$string_to_replace" with"
            echo "          "$replacement_string


           sed "s/$search_regex/$replacement_regex/" $file > $file.new
           mv $file.new $file

           string_to_replace="#include <$header_name>"
           replacement_string="#include <libmesh/$header_name>"

           search_regex="include <$header_name>"
           replacement_regex="include <libmesh\/$header_name>"

           echo "        Replacing "
           echo "          "$string_to_replace" with"
           echo "          "$replacement_string

           sed "s/$search_regex/$replacement_regex/" $file > $file.new
           mv $file.new $file
        fi
    done    
done
