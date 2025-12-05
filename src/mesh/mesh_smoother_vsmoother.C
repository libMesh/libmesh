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
VariationalMeshSmoother::VariationalMeshSmoother(
    UnstructuredMesh &mesh, Real dilation_weight,
    const bool preserve_subdomain_boundaries,
    const Real relative_residual_tolerance,
    const Real absolute_residual_tolerance, const unsigned int verbosity)
    : MeshSmoother(mesh), _verbosity(verbosity),
      _dilation_weight(dilation_weight),
      _preserve_subdomain_boundaries(preserve_subdomain_boundaries),
      _setup_called(false),
      _relative_residual_tolerance(relative_residual_tolerance),
      _absolute_residual_tolerance(absolute_residual_tolerance) {}

void VariationalMeshSmoother::setup()
{
  // Check for multiple dimensions
  if (_mesh.elem_dimensions().size() > 1)
    libmesh_not_implemented_msg("Meshes containing elements of differing dimension are not yet supported.");

  /*
  Ideally we'd want to update _mesh directly even if it already has an
  EquationSystems attached and doing work on it ... and that's thwarted by the
  problems:

  1. We can't easily tell if there's already an EquationSystems attached.
  2. We can't attach a second EquationSystems safely (if it already has a system 0)
     because the DoF indexing will need to be overwritten by our system.
  3. The destructor of es won't even clean up after itself (we generally expect
     a mesh to go unused after its EquationSystems is destroyed), much less know
     how to restore anything from a previous EquationSystems.

  To avoid these issues, we'll just construct a new DistributedMesh _mesh_copy
  from _mesh, then do the solve on _mesh_copy, then copy its node locations back
  to _mesh after the solve is done. That'll be slightly less memory-efficient
  and somewhat more CPU-efficient in the case where _mesh is serial, though
  it'll be significantly less memory-efficient when _mesh is already distributed,
  but either way the robustness is probably worth it.
  */

  // Create a new mesh, EquationSystems, and System
  _mesh_copy = std::make_unique<DistributedMesh>(_mesh);

  // If the _mesh wasn't prepared, that's fine (we'll just be moving
  // its nodes), but we do need the copy to be prepared before our
  // solve does things like looking at neighbors.  We'll disable
  // repartitioning and renumbering first to make sure that we can
  // transfer our geometry changes back to the original mesh.
  _mesh_copy->allow_renumbering(false);
  _mesh_copy->skip_partitioning(true);
  _mesh_copy->complete_preparation();

  _equation_systems = std::make_unique<EquationSystems>(*_mesh_copy);
  _system =
      &(_equation_systems->add_system<VariationalSmootherSystem>("variational_smoother_system"));

  // Set this to something > 0 to add more quadrature points than the default
  // rule that integrates order 2 * fe_order + 1 polynomials exactly.
  // Using higher quadrature orders has not had a significant effect on observed solutions.
  //system()->extra_quadrature_order = 0;

  // Uncomment these to debug
  //system()->print_element_solutions=true;
  //system()->print_element_residuals=true;
  //system()->print_element_jacobians=true;

  // Add boundary node and hanging node constraints
  _constraint = std::make_unique<VariationalSmootherConstraint>(
      *_system, _preserve_subdomain_boundaries, _verbosity);
  system()->attach_constraint_object(*_constraint);

  // Set system parameters
  system()->set_verbosity(_verbosity);
  system()->get_dilation_weight() = _dilation_weight;

  // Set up solver
  system()->time_solver = std::make_unique<SteadySolver>(*_system);

  // Uncomment this line and use -snes_test_jacobian and -snes_test_jacobian_view
  // flags to compare the hand-coded jacobian in VariationalSmootherSystem
  // to finite difference jacobians.
  //system()->time_solver->diff_solver() = std::make_unique<PetscDiffSolver>(*_system);

  _equation_systems->init();

  // Solver verbosity
  if (_verbosity > 15) {
    DiffSolver &solver = *(system()->time_solver->diff_solver().get());
    solver.quiet = false;
    solver.verbose = true;
  }

  // Solver convergence tolerances
  system()->time_solver->diff_solver()->relative_residual_tolerance = _relative_residual_tolerance;
  system()->time_solver->diff_solver()->absolute_residual_tolerance = _absolute_residual_tolerance;

  _setup_called = true;
}

void VariationalMeshSmoother::smooth()
{
  if (!_setup_called)
    setup();

  libmesh_try
    {
      system()->solve();
    }

  libmesh_catch (ConvergenceFailure& e)
    {
#ifdef LIBMESH_ENABLE_EXCEPTIONS
      throw ConvergenceFailure("The VariationalMeshSmoother solve failed to converge.");
#endif
    }

  // Update _mesh from _mesh_copy
  for (auto * node_copy : _mesh_copy->local_node_ptr_range())
    {
      auto & node = _mesh.node_ref(node_copy->id());
      for (const auto d : make_range(_mesh_copy->mesh_dimension()))
        node(d) = (*node_copy)(d);
    }

  SyncNodalPositions sync_object(_mesh);
  Parallel::sync_dofobject_data_by_id (_mesh.comm(), _mesh.nodes_begin(), _mesh.nodes_end(), sync_object);

  // Release memory occupied by _mesh_copy
  // Destruct this before _mesh_copy because it references _mesh_copy
  _equation_systems.reset();
  _mesh_copy.reset();
  // We'll need to call setup again since we'll have to reconstruct _mesh_copy
  _setup_called = false;
}

const MeshQualityInfo & VariationalMeshSmoother::get_mesh_info() const
{
  libmesh_error_msg_if(!_setup_called, "Need to first call the setup() method of "
      << "this VariationalMeshSmoother object, then call get_mesh_info() prior "
      << "to calling smooth().");
  return _system->get_mesh_info();
}

} // namespace libMesh

#endif // defined(LIBMESH_ENABLE_VSMOOTHER)
