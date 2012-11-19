#!/bin/bash 

# Write pretty status message
message_running() {
    echo "***************************************************************"
    echo "* Running Example " $LIBMESH_RUN $1 $2 $LIBMESH_OPTIONS
    echo "***************************************************************"
    echo " "
}

# Write pretty status message
message_done_running() {
    echo " "
    echo "***************************************************************"
    echo "* Done Running Example " $LIBMESH_RUN $1 $2 $LIBMESH_OPTIONS
    echo "***************************************************************"
}

run_example() {
    $LIBMESH_RUN ./$1 $2 $LIBMESH_OPTIONS || exit 1
}

# link_if_needed() {
#     if (test -f `basename $1`); then
# 	exit 0
#     fi
#     echo "--> linking $1 for setup"
#     ln -s $1 .
# }

# discard_link() {
#     file=`basename $1`
#     if (test -L $file); then
# 	echo "--> unlinking $file for cleanup"
# 	rm $file
#     fi
# }
