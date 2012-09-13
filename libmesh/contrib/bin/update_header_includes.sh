#!/bin/bash
set -e

# Note: script was /bin/sh but I don't think -n, -e are POSIX
# arguments to echo

if [ "x$LIBMESH_DIR" = x ]; then
  export LIBMESH_DIR=../..
fi


echo " "


for file in $@ ; do
    echo "Processing file $file"
    for header_name_with_path in `find $LIBMESH_ROOT/include -name "*.h" | sort | uniq` ; do
	
	header_name=`basename $header_name_with_path`
	
	string_to_replace="#include \"$header_name\""
	replacement_string="#include \"libmesh/$header_name\""
	replacement_regex="#include \"libmesh\/$header_name\""
	
	echo "  Replacing "
	echo "          "$string_to_replace" with"
	echo "          "$replacement_string
	
	
	#echo $file
	#echo "sed \"s/$string_to_replace/$replacement_regex/\" $file > $file.new"
	sed "s/$string_to_replace/$replacement_regex/" $file > $file.new
	mv $file.new $file
    done    
done

