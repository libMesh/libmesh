#!/bin/sh

headers_to_test=`ls ../../include/*.h`

if test $# -ge 1; then
    headers_to_test=$*
fi

for i in $headers_to_test; do
    header_name=`basename $i`
    source_file=TestHeader_$header_name.cc
    app_file=TestHeader_$header_name.o

    rm -f $source_file $app_file
    
    echo "#include \"$header_name\"" >> $source_file
    echo "int main () { return 0; }" >> $source_file

    #cat ../../src/apps/TestHeader_$header_name.cc

    echo " "
    echo "-------------------------------------------------"
    echo "Testing Header File $header_name"
        
    make -C ../.. contrib/bin/$app_file > /dev/null
    rm -f $source_file $app_file
done
