#!/bin/sh

built_sources=""

headers=`find .. -name "*.h" -type f | sort`

for header_with_path in $headers ; do
    
    #echo $header_with_path
    header=`basename $header_with_path`
    #echo $header
    built_sources="$built_sources $header"
done

specializations=`find .. -name "*specializations" -type f | sort`

for specialization_with_path in $specializations ; do
    
    #echo $specialization_with_path
    specialization=`basename $specialization_with_path`
    #echo $specialization
    built_sources="$built_sources $specialization"
done


cat <<EOF > Makefile.am
# special handholding for prefix_config.m4 generated files
# so that 'make clean ; make' works as does 'make distcheck'
# libmesh_config.h is made by ./configure, so it should get 
# cleaned by 'make distclean'
DISTCLEANFILES = libmesh_config.h

#
# special things to do when running 'make dist'
dist-hook:
	rm -rf \$(distdir)/libmesh_config.h

#
# include the magic script!
EXTRA_DIST = rebuild_makefile.sh

BUILT_SOURCES = $built_sources

DISTCLEANFILES += \$(BUILT_SOURCES)

EOF


for header_with_path in $headers $specializations ; do  
    header=`basename $header_with_path`
    source=`echo $header_with_path | sed 's/../$(top_srcdir)\/include/' -`
    #echo $source
    cat <<EOF >> Makefile.am
$header: $source
	\$(AM_V_GEN)\$(LN_S) $source $header

EOF
done
#cat Makefile.am