// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



// <h1>Adaptivity Example 1 - Solving 1D PDE Using Adaptive Mesh Refinement</h1>
// \author Benjamin S. Kirk
// \date 2003
//
// This example demonstrates how to solve a simple 1D problem
// using adaptive mesh refinement. The PDE that is solved is:
// -epsilon*u''(x) + u(x) = 1, on the domain [0,1] with boundary conditions
// u(0) = u(1) = 0 and where epsilon << 1.
//
// The approach used to solve 1D problems in libMesh is virtually identical to
// solving 2D or 3D problems, so in this sense this example represents a good
// starting point for new users. Note that many concepts are used in this
// example which are explained more fully in subsequent examples.

// Libmesh includes
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/edge_edge3.h"
#include "libmesh/gnuplot_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/fe.h"
#include "libmesh/getpot.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dof_map.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/error_vector.h"
#include "libmesh/kelly_error_estimator.h"
#include "libmesh/mesh_refinement.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

void assemble_1D(EquationSystems & es,
                 const std::string & system_name);

int main(int argc, char ** argv)
{
  // Initialize the library.  This is necessary because the library
  // may depend on a number of other libraries (i.e. MPI and PETSc)
  // that require initialization before use.  When the LibMeshInit
  // object goes out of scope, other libraries and resources are
  // finalized.
  LibMeshInit init (argc, argv);

  // Skip adaptive examples on a non-adaptive libMesh build
#ifndef LIBMESH_ENABLE_AMR
  libmesh_example_requires(false, "--enable-amr");
#else

  // Create a mesh, with dimension to be overridden later, on the
  // default MPI communicator.
  Mesh mesh(init.comm());

  GetPot command_line (argc, argv);

  int n = 4;
  if (command_line.search(1, "-n"))
    n = command_line.next(n);

  // Build a 1D mesh with 4 elements from x=0 to x=1, using
  // EDGE3 (i.e. quadratic) 1D elements. They are called EDGE3 elements
  // because a quadratic element contains 3 nodes.
  MeshTools::Generation::build_line(mesh, n, 0., 1., EDGE3);

  // Define the equation systems object and the system we are going
  // to solve. See Introduction Example 2 for more details.
  EquationSystems equation_systems(mesh);
  LinearImplicitSystem & system = equation_systems.add_system
    <LinearImplicitSystem>("1D");

  // Add a variable "u" to the system, using second-order approximation
  system.add_variable("u", SECOND);

  // Give the system a pointer to the matrix assembly function. This
  // will be called when needed by the library.
  system.attach_assemble_function(assemble_1D);

  // Define the mesh refinement object that takes care of adaptively
  // refining the mesh.
  MeshRefinement mesh_refinement(mesh);

  // These parameters determine the proportion of elements that will
  // be refined and coarsened. Any element within 30% of the maximum
  // error on any element will be refined, and any element within 30%
  // of the minimum error on any element might be coarsened
  mesh_refinement.refine_fraction()  = 0.7;
  mesh_refinement.coarsen_fraction() = 0.3;
  // We won't refine any element more than 5 times in total
  mesh_refinement.max_h_level()      = 5;

  // Initialize the data structures for the equation system.
  equation_systems.init();

  // Refinement parameters
  const unsigned int max_r_steps = 5; // Refine the mesh 5 times

  // Define the refinement loop
  for (unsigned int r_step=0; r_step<=max_r_steps; r_step++)
    {
      // Solve the equation system
      equation_systems.get_system("1D").solve();

      // We need to ensure that the mesh is not refined on the last iteration
      // of this loop, since we do not want to refine the mesh unless we are
      // going to solve the equation system for that refined mesh.
      if (r_step != max_r_steps)
        {
          // Error estimation objects, see Adaptivity Example 2 for details
          ErrorVector error;
          KellyErrorEstimator error_estimator;

          // Compute the error for each active element
          error_estimator.estimate_error(system, error);

          // Output error estimate magnitude
          libMesh::out << "Error estimate\nl2 norm = "
                       << error.l2_norm()
                       << "\nmaximum = "
                       << error.maximum()
                       << std::endl;

          // Flag elements to be refined and coarsened
          mesh_refinement.flag_elements_by_error_fraction (error);

          // Perform refinement and coarsening
          mesh_refinement.refine_and_coarsen_elements();

          // Reinitialize the equation_systems object for the newly refined
          // mesh. One of the steps in this is project the solution onto the
          // new mesh
          equation_systems.reinit();
        }
    }

  // Construct gnuplot plotting object, pass in mesh, title of plot
  // and boolean to indicate use of grid in plot. The grid is used to
  // show the edges of each element in the mesh.
  GnuPlotIO plot(mesh, "Adaptivity Example 1", GnuPlotIO::GRID_ON);

  // Write out script to be called from within gnuplot:
  // Load gnuplot, then type "call 'gnuplot_script'" from gnuplot prompt
  plot.write_equation_systems("gnuplot_script", equation_systems);
#endif // #ifndef LIBMESH_ENABLE_AMR

  // All done.  libMesh objects are destroyed here.  Because the
  // LibMeshInit object was created first, its destruction occurs
  // last, and it's destructor finalizes any external libraries and
  // checks for leaked memory.
  return 0;
}




// Define the matrix assembly function for the 1D PDE we are solving
void assemble_1D(EquationSystems & es,
                 const std::string & libmesh_dbg_var(system_name))
{

#ifdef LIBMESH_ENABLE_AMR

  // It is a good idea to check we are solving the correct system
  libmesh_assert_equal_to (system_name, "1D");

  // Get a reference to the mesh object
  const MeshBase & mesh = es.get_mesh();

  // The dimension we are using, i.e. dim==1
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to the system we are solving
  LinearImplicitSystem & system = es.get_system<LinearImplicitSystem>("1D");

  // Get a reference to the DofMap object for this system. The DofMap object
  // handles the index translation from node and element numbers to degree of
  // freedom numbers. DofMap's are discussed in more detail in future examples.
  const DofMap & dof_map = system.get_dof_map();

  // Get a constant reference to the Finite Element type for the first
  // (and only) variable in the system.
  FEType fe_type = dof_map.variable_type(0);

  // Build a finite element object of the specified type. The build
  // function dynamically allocates memory so we use a UniquePtr in this case.
  // A UniquePtr is a pointer that cleans up after itself. See examples 3 and 4
  // for more details on UniquePtr.
  UniquePtr<FEBase> fe(FEBase::build(dim, fe_type));

  // Tell the finite element object to use fifth order Gaussian quadrature
  QGauss qrule(dim, FIFTH);
  fe->attach_quadrature_rule(&qrule);

  // Here we define some references to cell-specific data that will be used to
  // assemble the linear system.

  // The element Jacobian * quadrature weight at each integration point.
  const std::vector<Real> & JxW = fe->get_JxW();

  // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector<Real> > & phi = fe->get_phi();

  // The element shape function gradients evaluated at the quadrature points.
  const std::vector<std::vector<RealGradient> > & dphi = fe->get_dphi();

  // Declare a dense matrix and dense vector to hold the element matrix
  // and right-hand-side contribution
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  // This vector will hold the degree of freedom indices for the element.
  // These define where in the global system the element degrees of freedom
  // get mapped.
  std::vector<dof_id_type> dof_indices;

  // We now loop over all the active elements in the mesh in order to calculate
  // the matrix and right-hand-side contribution from each element. Use a
  // const_element_iterator to loop over the elements. We make
  // el_end const as it is used only for the stopping condition of the loop.
  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();

  // Note that ++el is preferred to el++ when using loops with iterators
  for ( ; el != el_end; ++el)
    {
      // It is convenient to store a pointer to the current element
      const Elem * elem = *el;

      // Get the degree of freedom indices for the current element.
      // These define where in the global matrix and right-hand-side this
      // element will contribute to.
      dof_map.dof_indices(elem, dof_indices);

      // Compute the element-specific data for the current element. This
      // involves computing the location of the quadrature points (q_point)
      // and the shape functions (phi, dphi) for the current element.
      fe->reinit(elem);

      // Store the number of local degrees of freedom contained in this element
      const int n_dofs = dof_indices.size();

      // We resize and zero out Ke and Fe (resize() also clears the matrix and
      // vector). In this example, all elements in the mesh are EDGE3's, so
      // Ke will always be 3x3, and Fe will always be 3x1. If the mesh contained
      // different element types, then the size of Ke and Fe would change.
      Ke.resize(n_dofs, n_dofs);
      Fe.resize(n_dofs);


      // Now loop over quadrature points to handle numerical integration
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
          // Now build the element matrix and right-hand-side using loops to
          // integrate the test functions (i) against the trial functions (j).
          for (std::size_t i=0; i<phi.size(); i++)
            {
              Fe(i) += JxW[qp]*phi[i][qp];

              for (std::size_t j=0; j<phi.size(); j++)
                {
                  Ke(i,j) += JxW[qp]*(1.e-3*dphi[i][qp]*dphi[j][qp] +
                                      phi[i][qp]*phi[j][qp]);
                }
            }
        }


      // At this point we have completed the matrix and RHS summation. The
      // final step is to apply boundary conditions, which in this case are
      // simple Dirichlet conditions with u(0) = u(1) = 0.

      // Define the penalty parameter used to enforce the BC's
      double penalty = 1.e10;

      // Loop over the sides of this element. For a 1D element, the "sides"
      // are defined as the nodes on each edge of the element, i.e. 1D elements
      // have 2 sides.
      for (unsigned int s=0; s<elem->n_sides(); s++)
        {
          // If this element has a NULL neighbor, then it is on the edge of the
          // mesh and we need to enforce a boundary condition using the penalty
          // method.
          if (elem->neighbor_ptr(s) == libmesh_nullptr)
            {
              Ke(s,s) += penalty;
              Fe(s)   += 0*penalty;
            }
        }

      // This is a function call that is necessary when using adaptive
      // mesh refinement. See Adaptivity Example 2 for more details.
      dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

      // Add Ke and Fe to the global matrix and right-hand-side.
      system.matrix->add_matrix(Ke, dof_indices);
      system.rhs->add_vector(Fe, dof_indices);
    }
#endif // #ifdef LIBMESH_ENABLE_AMR
}
