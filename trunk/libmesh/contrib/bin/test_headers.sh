#!/bin/sh

for i in `ls ../../include/*.h`; do
    header_name=`basename $i`
    source_file=../../src/apps/TestHeader_$header_name.cc
    app_file=bin/TestHeader_$header_name

    rm -f $source_file ../../$app_file
    
    echo "#include \"$header_name\"" >> $source_file
    echo "int main () { return 0; }" >> $source_file

    #cat ../../src/apps/TestHeader_$header_name.cc

    echo " "
    echo "-------------------------------------------------"
    echo "Testing Header File $header_name"
        
    make -C ../.. $app_file > /dev/null
    rm -f $source_file ../../$app_file
done
