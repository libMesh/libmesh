#!/bin/sh

#set -x

. "$LIBMESH_DIR"/examples/run_common.sh

example_name=adaptivity_ex3
example_dir=examples/adaptivity/$example_name

# Test with a bunch of different refinement types; while we're at it
# we can get a little test coverage of our nodal partitioning options.
#
# --load-balanced-nodes-bfs doesn't appear to be safe to use on a
# distributed mesh!  I see dbg assertion failures from a sweep over as
# few as 15 processors!

run_example "$example_name" refinement_type=h
# Some solvers still give us trouble with too much p
run_example "$example_name" refinement_type=p max_r_steps=4
run_example "$example_name" refinement_type=hp max_r_steps=6 --load-balanced-nodes-linear
if [ "x$petscmajor" == "x" ]; then
run_example "$example_name" refinement_type=matchedhp max_r_steps=4 --load-balanced-nodes-petscpartitioner
else
run_example "$example_name" refinement_type=matchedhp max_r_steps=4
fi

# Let's get some 1D coverage too; that's cheap.
# Load balancing nodes should be pointless here but we can make sure
# it doesn't regress in this special case.
run_example "$example_name" dimension=1 refinement_type=h --load-balanced-nodes-linear
# After enough p refinement we start breaking non-pivoting ILU preconditioners
run_example "$example_name" dimension=1 refinement_type=p max_r_steps=8
run_example "$example_name" dimension=1 refinement_type=hp max_r_steps=6
run_example "$example_name" dimension=1 refinement_type=matchedhp max_r_steps=4

# And try out another element type
if [ "x$petscmajor" == "x" ]; then
run_example "$example_name" dimension=1 approx_type=HERMITE approx_order=3 refinement_type=h --load-balanced-nodes-petscpartitioner
else
run_example "$example_name" dimension=1 approx_type=HERMITE approx_order=3 refinement_type=h
fi
run_example "$example_name" dimension=1 approx_type=HERMITE approx_order=3 refinement_type=p max_r_steps=6
run_example "$example_name" dimension=1 approx_type=HERMITE approx_order=3 refinement_type=hp max_r_steps=8 --load-balanced-nodes-linear
run_example "$example_name" dimension=1 approx_type=HERMITE approx_order=3 refinement_type=matchedhp max_r_steps=4

# Examples to use for test coverage
run_example "$example_name" dimension=2 element_type=simplex approx_type=LAGRANGE approx_order=2 refinement_type=h max_r_steps=3
run_example "$example_name" dimension=2 element_type=simplex complete=true approx_type=LAGRANGE approx_order=3 refinement_type=h max_r_steps=3
run_example "$example_name" dimension=2 element_type=simplex approx_type=HIERARCHIC approx_order=3 refinement_type=h max_r_steps=3
run_example "$example_name" dimension=2 element_type=simplex approx_type=CLOUGH singularity=false approx_order=3 refinement_type=h max_r_steps=3
run_example "$example_name" dimension=3 element_type=simplex extrusion=true approx_type=LAGRANGE approx_order=1 refinement_type=h max_r_steps=3
run_example "$example_name" dimension=3 element_type=simplex extrusion=true approx_type=LAGRANGE approx_order=2 refinement_type=h max_r_steps=3
run_example "$example_name" dimension=3 element_type=simplex extrusion=true complete=true approx_type=LAGRANGE approx_order=2 refinement_type=h max_r_steps=3
run_example "$example_name" dimension=3 element_type=simplex extrusion=true complete=true approx_type=LAGRANGE approx_order=3 refinement_type=h max_r_steps=3

# Examples to use for benchmarking
benchmark_example 1 "$example_name" dimension=3 approx_type=LAGRANGE approx_order=2 refinement_type=h max_r_steps=12
benchmark_example 1 "$example_name" dimension=3 approx_type=HIERARCHIC approx_order=3 refinement_type=h
