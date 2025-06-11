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
#if defined(LIBMESH_ENABLE_VSMOOTHER) && LIBMESH_DIM > 1

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
  _percent_to_move(1),
  _dim(mesh.mesh_dimension()),
  _dilation_weight(dilation_weight),
  _n_nodes(0),
  _n_cells(0),
  _preserve_subdomain_boundaries(preserve_subdomain_boundaries)
{}


Real VariationalMeshSmoother::smooth(unsigned int)
{
  // Update the mesh dimension, since the mesh may have changed since initialization
  _dim = _mesh.mesh_dimension();

  // Check for multiple dimensions
  if (_mesh.elem_dimensions().size() > 1)
    libmesh_not_implemented_msg("Meshes containing elements of differing dimension are not yet supported.");

  // Records the relative "distance moved"
  Real dist_norm = 0.;

  Real dilation_weight = _dilation_weight;

  // Initialize the _n_nodes and _n_cells member variables
  this->_n_nodes = _mesh.n_nodes();
  this->_n_cells = _mesh.n_active_elem();

  Array2D<Real> R(_n_nodes, _dim);

  // initial grid
  readgr(R);

  // grid optimization
  full_smooth(R, dilation_weight);

  // save result
  dist_norm = writegr(R);

  return dist_norm;
}



// save grid
Real VariationalMeshSmoother::writegr(const Array2D<Real> & R)
{
  libMesh::out << "Starting writegr" << std::endl;

  Real dist_norm = 0.;
  // Adjust nodal coordinates to new positions
  {
    int i = 0;
    for (auto & node : _mesh.node_ptr_range())
      {
        Real total_dist = 0.;

        // Get a reference to the node
        Node & node_ref = *node;

        // For each node set its X Y [Z] coordinates
        for (unsigned int j=0; j<_dim; j++)
          {
            Real distance = R[i][j] - node_ref(j);

            // Save the squares of the distance
            total_dist += Utility::pow<2>(distance);

            node_ref(j) += distance * _percent_to_move;
          }

        libmesh_assert_greater_equal (total_dist, 0.);

        // Add the distance this node moved to the global distance
        dist_norm += total_dist;

        i++;
      }

    // Relative "error"
    dist_norm = std::sqrt(dist_norm/_mesh.n_nodes());
  }

  libMesh::out << "Finished writegr" << std::endl;
  return dist_norm;
}



// reading grid from input file
int VariationalMeshSmoother::readgr(Array2D<Real> & R)
{
  libMesh::out << "Starting readgr" << std::endl;

  int i = 0;
  for (auto & node : _mesh.node_ptr_range())
  {
    // Get a reference to the node
    Node & node_ref = *node;

    // For each node grab its X [Y] [Z] coordinates
    for (unsigned int j=0; j<_dim; j++)
      R[i][j] = node_ref(j);
  }

  return 0;
}


void VariationalMeshSmoother::readNodesIntoArray(Array2D<Real> & R, const VariationalSmootherSystem & system) const
{
  unsigned int i = 0;
  for (auto * node : system.get_mesh().node_ptr_range())
  {
    // For each node grab its X Y [Z] coordinates
    for (unsigned int j=0; j < _dim; j++)
      R[i][j] = (*node)(j);

    ++i;
  }
}


// Preprocess mesh data and control smoothing/untangling iterations
void VariationalMeshSmoother::full_smooth(Array2D<Real> & R, Real dilation_weight)
{
  Real vol = 1.;
 
  // Create a new mesh, EquationSystems, and System
  DistributedMesh mesh(_mesh);
  EquationSystems es(mesh);
  VariationalSmootherSystem & sys = es.add_system<VariationalSmootherSystem>("variational_smoother_system");
  //sys.print_element_solutions=true;
  //sys.print_element_residuals=true;
  //sys.print_element_jacobians=true;
  sys.extra_quadrature_order = 0;

  // Add boundary node and hanging node constraints
  VariationalSmootherConstraint constraint(sys, _preserve_subdomain_boundaries);
  sys.attach_constraint_object(constraint);

  // Set system parameters
  sys.get_ref_vol() = vol;
  sys.get_dilation_weight() = dilation_weight;

  // Set up solver
  sys.time_solver =
    std::make_unique<SteadySolver>(sys);
  // Uncomment this line and use -snes_test_jacobian and -snes_test_jacobian_view
  // flags to compare the hand-coded jacobian in VariationalSmootherSystem
  // to finite difference jacobians.
  //sys.time_solver->diff_solver() = std::make_unique<PetscDiffSolver>(sys);

  es.init();

  DiffSolver & solver = *(sys.time_solver->diff_solver().get());
  sys.time_solver->diff_solver()->relative_residual_tolerance = 1e-10;
  solver.quiet = false;
  solver.verbose = true;
  solver.relative_step_tolerance = 1e-10;

  sys.solve();

  // Update _mesh with solution
  readNodesIntoArray(R, sys);
}

} // namespace libMesh

#endif // defined(LIBMESH_ENABLE_VSMOOTHER) && LIBMESH_DIM > 1
