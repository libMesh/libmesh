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
    $LIBMESH_RUN ./$1 $2 || exit 1
}