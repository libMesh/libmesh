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
    echo "*  $LIBMESH_RUN ./$executable $options $LIBMESH_OPTIONS"
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
    echo "*  $LIBMESH_RUN ./$executable $options $LIBMESH_OPTIONS"
    echo "***************************************************************"
}

run_example() {

    example_name=$1
    shift
    options=$@

    # when run outside of the automake envionment make sure we get METHODS set
    # to something useful
    if (test "x${METHODS}" = "x"); then
	if (test "x${METHOD}" = "x"); then
	    METHODS=opt
	else
	    METHODS="$METHOD"
	fi
    fi             
    
    for method in ${METHODS}; do
	
	case "${method}" in
	    optimized|opt)      executable=example-opt   ;;
	    debug|dbg)          executable=example-dbg   ;;
	    devel)              executable=example-devel ;;
	    profiling|pro|prof) executable=example-prof  ;;
	    oprofile|oprof)     executable=example-oprof ;;
	    *) echo "ERROR: unknown method: ${method}!" ; exit 1 ;;
	esac

	if (test ! -x ${executable}); then
	    echo "ERROR: cannot find ${executable}!"
	    exit 1
	fi
	
	message_running $example_name $executable $options

	$LIBMESH_RUN ./$executable $options $LIBMESH_OPTIONS || exit $?
	
	message_done_running $example_name $executable $options
    done
}

run_example_no_extra_options() {
  LIBMESH_OPTIONS='' run_example $@
}
