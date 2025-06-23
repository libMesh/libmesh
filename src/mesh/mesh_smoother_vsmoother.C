// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#include "libmesh/libmesh_config.h"
#if defined(LIBMESH_ENABLE_VSMOOTHER)

// Local includes
#include "libmesh/mesh_smoother_vsmoother.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/elem.h"
#include "libmesh/unstructured_mesh.h"
#include "libmesh/utility.h"
#include "libmesh/boundary_info.h"
#include "libmesh/equation_systems.h"
#include "libmesh/distributed_mesh.h"
#include "libmesh/steady_solver.h"
#include "libmesh/diff_solver.h"
#include "libmesh/variational_smoother_constraint.h"
#include "libmesh/parallel_ghost_sync.h"

// C++ includes
#include <time.h> // for clock_t, clock()
#include <cstdlib> // *must* precede <cmath> for proper std:abs() on PGI, Sun Studio CC
#include <cmath>
#include <iomanip>
#include <limits>

namespace libMesh
{

// Optimization at -O2 or greater seem to break Intel's icc. So if we are
// being compiled with icc let's dumb-down the optimizations for this file
#ifdef __INTEL_COMPILER
#  pragma optimize ( "", off )
#endif

// Member functions for the Variational Smoother
VariationalMeshSmoother::VariationalMeshSmoother(UnstructuredMesh & mesh,
                                                 Real dilation_weight,
                                                 const bool preserve_subdomain_boundaries) :
  MeshSmoother(mesh),
  _dilation_weight(dilation_weight),
  _preserve_subdomain_boundaries(preserve_subdomain_boundaries)
{}


void VariationalMeshSmoother::smooth(unsigned int)
{
  // Check for multiple dimensions
  if (_mesh.elem_dimensions().size() > 1)
    libmesh_not_implemented_msg("Meshes containing elements of differing dimension are not yet supported.");

  // Create a new mesh, EquationSystems, and System
  DistributedMesh mesh_copy(_mesh);
  EquationSystems es(mesh_copy);
  VariationalSmootherSystem & sys = es.add_system<VariationalSmootherSystem>("variational_smoother_system");

  // Set this to something > 0 to add more quadrature points than the default
  // rule that integrates order 2 * nq_points + 1 polynomials exactly.
  // Using higher quadrature orders has not had a significant effect on observed solutions.
  //sys.extra_quadrature_order = 0;

  // Uncomment these to debug
  //sys.print_element_solutions=true;
  //sys.print_element_residuals=true;
  //sys.print_element_jacobians=true;

  // Add boundary node and hanging node constraints
  VariationalSmootherConstraint constraint(sys, _preserve_subdomain_boundaries);
  sys.attach_constraint_object(constraint);

  // Set system parameters
  sys.get_dilation_weight() = _dilation_weight;

  // Set up solver
  sys.time_solver =
    std::make_unique<SteadySolver>(sys);

  // Uncomment this line and use -snes_test_jacobian and -snes_test_jacobian_view
  // flags to compare the hand-coded jacobian in VariationalSmootherSystem
  // to finite difference jacobians.
  //sys.time_solver->diff_solver() = std::make_unique<PetscDiffSolver>(sys);

  es.init();

  // More debugging options
  DiffSolver & solver = *(sys.time_solver->diff_solver().get());
  //solver.quiet = false;
  solver.verbose = true;
  sys.time_solver->diff_solver()->relative_residual_tolerance = TOLERANCE*TOLERANCE;

  sys.solve();

  // Update _mesh from mesh_copy
  for (auto * node_copy : mesh_copy.local_node_ptr_range())
  {
    auto & node = _mesh.node_ref(node_copy->id());
    for (const auto d : make_range(mesh_copy.mesh_dimension()))
      node(d) = (*node_copy)(d);
  }

  SyncNodalPositions sync_object(_mesh);
  Parallel::sync_dofobject_data_by_id (_mesh.comm(), _mesh.nodes_begin(), _mesh.nodes_end(), sync_object);
}

} // namespace libMesh

#endif // defined(LIBMESH_ENABLE_VSMOOTHER)
