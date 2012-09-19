#!/bin/bash

headers=`find . -name "*.h" -type f | sort`
#echo $headers

include_headers="include_HEADERS = "'\n'

echo -n "include_HEADERS = "
for header_with_path in $headers ; do
    
    header=`basename $header_with_path`
    source=`echo $header_with_path | gsed 's/.\///' -`

    echo " \\"
    echo -n "        "$source

done


echo " "

