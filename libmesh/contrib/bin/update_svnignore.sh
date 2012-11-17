#!/bin/bash

dirs=`find . -type d | grep -v .svn`

for dir in $dirs ; do
    if (test -d $dir/.svn); then 
	cd $dir 
	echo "----------------------------------------------------------------"
	echo $dir
	echo "----------------------------------------------------------------"
	svn propset svn:ignore . --file /home/benkirk/codes/libmesh/svnignore.txt
# 	    svn status
	for unversioned in `svn status | grep -v '/'| grep ? | cut -d'?' -f2`; do
	    echo $unversioned >> .svnignore
	done
	if (test -f .svnignore); then
	    cat .svnignore | sort | uniq > .svnignore2
	    mv .svnignore2 .svnignore
	    svn propget svn:ignore . >> .svnignore
	    cat .svnignore
	    
	    svn propset svn:ignore . --file .svnignore
	    rm -f .svnignore
	fi
	
# 	if (test -f Makefile.); then
# 	    svn propset svn:ignore . --file /home/benkirk/codes/libmesh/svnignore.txt
# 	fi
	    cd -
    fi
done
