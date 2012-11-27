#!/bin/bash 

# Write pretty status message
message_running() {
    example_name=$1
    shift
    executable=$1
    shift
    options=$@

    echo "***************************************************************"
    echo "* Running Example $example_name:" 
    echo "*  $LIBMESH_RUN $executable $options $LIBMESH_OPTIONS"
    echo "***************************************************************"
    echo " "
}

# Write pretty status message
message_done_running() {

    example_name=$1
    shift
    executable=$1
    shift
    options=$@

    echo " "
    echo "***************************************************************"
    echo "* Done Running Example $example_name:" 
    echo "*  $LIBMESH_RUN $executable $options $LIBMESH_OPTIONS"
    echo "***************************************************************"
}

run_example() {

    example_name=$1
    shift
    options=$@

    for method in ${METHODS}; do
	
	case "${method}" in
	    opt)   executable=example-opt ;;
	    dbg)   executable=example-dbg ;;
	    devel) executable=example-devel ;;
	    *) echo "ERROR: unknown method: ${method}!" ; exit 1 ;;
	esac

	if (test ! -x ${executable}); then
	    echo "ERROR: cannot find ${executable}!"
	    exit 1
	fi
	
	message_running $example_name $executable $options

	$LIBMESH_RUN ./$executable $options $LIBMESH_OPTIONS || exit 1
	
	message_done_running $example_name $executable $options
    done
}
