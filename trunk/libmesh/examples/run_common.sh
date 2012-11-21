#!/bin/bash 

# Write pretty status message
message_running() {

    executable=$1
    if (test ! -x $executable -a -x example-$METHOD); then
	executable=example-$METHOD
    fi

    echo "***************************************************************"
    echo "* Running Example " $LIBMESH_RUN $executable $2 $LIBMESH_OPTIONS
    echo "***************************************************************"
    echo " "
}

# Write pretty status message
message_done_running() {

    executable=$1
    if (test ! -x $executable -a -x example-$METHOD); then
	executable=example-$METHOD
    fi

    echo " "
    echo "***************************************************************"
    echo "* Done Running Example " $LIBMESH_RUN $executable $2 $LIBMESH_OPTIONS
    echo "***************************************************************"
}

run_example() {

    executable=$1
    if (test ! -x $executable -a -x example-$METHOD); then
	executable=example-$METHOD
    fi

    $LIBMESH_RUN ./$executable $2 $LIBMESH_OPTIONS || exit 1
}
