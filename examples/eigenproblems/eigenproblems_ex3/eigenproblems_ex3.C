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



// <h1>Eigenproblems Example 3 - Can you "hear the shape" of a drum?</h1>
// \author David Knezevic
// \date 2012
//
// The sound that a drum makes is determined by it's resonant frequencies,
// which are given by the eigenvalues of the Laplacian. "Can One Hear the
// Shape of a Drum?" was the title of an article by Mark Kac
// in the American Mathematical Monthly in 1966, where he raised the question:
// If we know all the eigenvalues of a drum, can we uniquely determine it's shape?
// This question was resolved in 1992, when Gordon, Webb, and Wolpert constructed
// a pair of regions in 2D that have different shapes but identical eigenvalues.
// So the answer to Kac's question is no: the spectrum of the Laplacian does
// not uniquely determine the shape of the domain.

// In this example, we compute the first few eigenvalues of the two domains proposed
// by Gordon, Webb and Wolpert. This amounts to solving a generalized eigenvalue
// problem in each case. The computed eigenvalues are stored in drum1_evals.txt and
// drum2_evals.txt. We can compare these to the (highly accurate) values reported in:
// T.A. Driscoll, "Eigenmodes of Isospectral Drums", SIAM Review, Vol. 39, No. 1, pp. 1-17, 1997.
//
// The first five eigenvalues from Driscoll are listed below (the author states that
// "all digits shown are believed to be correct"):
//  2.53794399980
//  3.65550971352
//  5.17555935622
//  6.53755744376
//  7.24807786256


// C++ include files
#include <fstream>

// libMesh include files.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/condensed_eigen_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dof_map.h"
#include "libmesh/fe_interface.h"
#include "libmesh/getpot.h"
#include "libmesh/elem.h"
#include "libmesh/zero_function.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/slepc_macro.h"
#include "libmesh/enum_eigen_solver_type.h"

#define BOUNDARY_ID 100

// Bring in everything from the libMesh namespace
using namespace libMesh;


// Function prototype.  This is the function that will assemble
// the eigen system. Here, we will simply assemble a mass matrix.
void assemble_matrices(EquationSystems & es,
                       const std::string & system_name);

// We store the Dirichlet dofs in a set in order to impose the boundary conditions
void get_dirichlet_dofs(EquationSystems & es,
                        const std::string & system_name,
                        std::set<unsigned int> & global_dirichlet_dofs_set);


int main (int argc, char ** argv)
{
  // Initialize libMesh and the dependent libraries.
  LibMeshInit init (argc, argv);

  // This example uses an ExodusII input file
#ifndef LIBMESH_HAVE_EXODUS_API
  libmesh_example_requires(false, "--enable-exodus");
#endif

  // This example is designed for the SLEPc eigen solver interface.
#ifndef LIBMESH_HAVE_SLEPC
  if (init.comm().rank() == 0)
    libMesh::err << "ERROR: This example requires libMesh to be\n"
                 << "compiled with SLEPc eigen solvers support!"
                 << std::endl;

  return 0;
#else

#ifdef LIBMESH_DEFAULT_SINGLE_PRECISION
  // SLEPc currently gives us a nasty crash with Real==float
  libmesh_example_requires(false, "--disable-singleprecision");
#endif

#if defined(LIBMESH_USE_COMPLEX_NUMBERS) && SLEPC_VERSION_LESS_THAN(3,6,2)
  // SLEPc used to give us an "inner product not well defined" with
  // Number==complex; but this problem seems to be solved in newer versions.
  libmesh_example_requires(false, "--disable-complex or use SLEPc>=3.6.2");
#endif

  // We use Dirichlet boundary conditions here
#ifndef LIBMESH_ENABLE_DIRICHLET
  libmesh_example_requires(false, "--enable-dirichlet");
#endif

  // Tell the user what we are doing.
  {
    libMesh::out << "Running " << argv[0];

    for (int i=1; i<argc; i++)
      libMesh::out << " " << argv[i];

    libMesh::out << std::endl << std::endl;
  }

  // Skip this 2D example if libMesh was compiled as 1D-only.
  libmesh_example_requires(2 <= LIBMESH_DIM, "2D support");

  // Use GetPot to parse the command line arguments
  GetPot command_line (argc, argv);

  // Read the mesh name from the command line
  const std::string mesh_name =
    libMesh::command_line_next("-mesh_name", std::string());

  // Also, read in the index of the eigenvector that we should plot
  // (zero-based indexing, as usual!)
  const unsigned int plotting_index =
    libMesh::command_line_next("-plotting_index", 0u);

  // Finally, read in the number of eigenpairs we want to compute!
  const unsigned int n_evals =
    libMesh::command_line_next("-n_evals", 0u);

  // Append the .e to mesh_name
  std::ostringstream mesh_name_exodus;
  mesh_name_exodus << mesh_name << "_mesh.e";

  // Create a mesh, with dimension to be overridden by the file, on
  // the default MPI communicator.
  Mesh mesh(init.comm());

  mesh.read(mesh_name_exodus.str());

  // Add boundary IDs to this mesh so that we can use DirichletBoundary
  // Each processor should know about each boundary condition it can
  // see, so we loop over all elements, not just local elements.
  for (const auto & elem : mesh.element_ptr_range())
    for (auto side : elem->side_index_range())
      if (elem->neighbor_ptr (side) == nullptr)
        mesh.get_boundary_info().add_side(elem, side, BOUNDARY_ID);

  mesh.get_boundary_info().regenerate_id_sets();

  // Print information about the mesh to the screen.
  mesh.print_info();

  // Create an equation systems object.
  EquationSystems equation_systems (mesh);

  // Create a CondensedEigenSystem named "Eigensystem" and (for convenience)
  // use a reference to the system we create.
  CondensedEigenSystem & eigen_system =
    equation_systems.add_system<CondensedEigenSystem> ("Eigensystem");

  // Declare the system variables.
  // Adds the variable "p" to "Eigensystem".   "p"
  // will be approximated using second-order approximation.
  eigen_system.add_variable("p", SECOND);

  // Give the system a pointer to the matrix assembly
  // function defined below.
  eigen_system.attach_assemble_function (assemble_matrices);

  // Set the number of requested eigenpairs n_evals and the number
  // of basis vectors used in the solution algorithm.
  equation_systems.parameters.set<unsigned int>("eigenpairs")    = n_evals;
  equation_systems.parameters.set<unsigned int>("basis vectors") = n_evals*3;

  // Set the solver tolerance and the maximum number of iterations.
  equation_systems.parameters.set<Real>("linear solver tolerance") = pow(TOLERANCE, 5./3.);
  equation_systems.parameters.set<unsigned int>
    ("linear solver maximum iterations") = 1000;

  // Set the type of the problem, here we deal with
  // a generalized Hermitian problem.
  eigen_system.set_eigenproblem_type(GHEP);

  // Set the target eigenvalue
  eigen_system.get_eigen_solver().set_position_of_spectrum(0., TARGET_REAL);

  {
    ZeroFunction<> zf;

#ifdef LIBMESH_ENABLE_DIRICHLET
    // Most DirichletBoundary users will want to supply a "locally
    // indexed" functor
    DirichletBoundary dirichlet_bc({BOUNDARY_ID}, {0}, zf,
                                   LOCAL_VARIABLE_ORDER);

    eigen_system.get_dof_map().add_dirichlet_boundary(dirichlet_bc);
#endif
  }

  // Initialize the data structures for the equation system.
  equation_systems.init();

  // Prints information about the system to the screen.
  equation_systems.print_info();

  eigen_system.initialize_condensed_dofs();

  // Solve the system "Eigensystem".
  eigen_system.solve();

  // Get the number of converged eigen pairs.
  unsigned int nconv = eigen_system.get_n_converged();

  libMesh::out << "Number of converged eigenpairs: "
               << nconv
               << "\n"
               << std::endl;

  if (plotting_index > n_evals)
    {
      libMesh::out << "WARNING: Solver did not converge for the requested eigenvector!" << std::endl;
    }

  // write out all of the computed eigenvalues and plot the specified eigenvector
  std::ostringstream eigenvalue_output_name;
  eigenvalue_output_name << mesh_name << "_evals.txt";
  std::ofstream evals_file(eigenvalue_output_name.str().c_str());

  for (unsigned int i=0; i<nconv; i++)
    {
      std::pair<Real,Real> eval = eigen_system.get_eigenpair(i);

      // The eigenvalues should be real!
      libmesh_assert_less (eval.second, TOLERANCE);
      evals_file << eval.first << std::endl;

      // plot the specified eigenvector
      if (i == plotting_index)
        {
#ifdef LIBMESH_HAVE_EXODUS_API
          // Write the eigen vector to file.
          std::ostringstream eigenvector_output_name;
          eigenvector_output_name << mesh_name << "_evec.e";
          ExodusII_IO (mesh).write_equation_systems (eigenvector_output_name.str(), equation_systems);
#endif // #ifdef LIBMESH_HAVE_EXODUS_API
        }
    }

  evals_file.close();

#endif // LIBMESH_HAVE_SLEPC

  // All done.
  return 0;
}



void assemble_matrices(EquationSystems & es,
                       const std::string & libmesh_dbg_var(system_name))
{

  // It is a good idea to make sure we are assembling
  // the proper system.
  libmesh_assert_equal_to (system_name, "Eigensystem");

#ifdef LIBMESH_HAVE_SLEPC

  // Get a constant reference to the mesh object.
  const MeshBase & mesh = es.get_mesh();

  // The dimension that we are running.
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to our system.
  EigenSystem & eigen_system = es.get_system<EigenSystem> ("Eigensystem");

  // Get a constant reference to the Finite Element type
  // for the first (and only) variable in the system.
  FEType fe_type = eigen_system.get_dof_map().variable_type(0);

  // A reference to the two system matrices
  SparseMatrix<Number> & matrix_A = eigen_system.get_matrix_A();
  SparseMatrix<Number> & matrix_B = eigen_system.get_matrix_B();

  // Build a Finite Element object of the specified type.  Since the
  // FEBase::build() member dynamically creates memory we will
  // store the object as a std::unique_ptr<FEBase>.  This can be thought
  // of as a pointer that will clean up after itself.
  std::unique_ptr<FEBase> fe (FEBase::build(dim, fe_type));

  // A  Gauss quadrature rule for numerical integration.
  // Use the default quadrature order.
  QGauss qrule (dim, fe_type.default_quadrature_order());

  // Tell the finite element object to use our quadrature rule.
  fe->attach_quadrature_rule (&qrule);

  // The element Jacobian * quadrature weight at each integration point.
  const std::vector<Real> & JxW = fe->get_JxW();

  // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector<Real>> & phi = fe->get_phi();

  // The element shape function gradients evaluated at the quadrature
  // points.
  const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();

  // A reference to the DofMap object for this system.  The DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.
  const DofMap & dof_map = eigen_system.get_dof_map();

  // The element mass and stiffness matrices.
  DenseMatrix<Number> Me;
  DenseMatrix<Number> Ke;

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;


  // Now we will loop over all the elements in the mesh that
  // live on the local processor. We will compute the element
  // matrix and right-hand-side contribution.  In case users
  // later modify this program to include refinement, we will
  // be safe and will only consider the active elements;
  // hence we use a variant of the active_elem_iterator.
  for (const auto & elem : mesh.active_local_element_ptr_range())
    {
      // Get the degree of freedom indices for the
      // current element.  These define where in the global
      // matrix and right-hand-side this element will
      // contribute to.
      dof_map.dof_indices (elem, dof_indices);

      // Compute the element-specific data for the current
      // element.  This involves computing the location of the
      // quadrature points (q_point) and the shape functions
      // (phi, dphi) for the current element.
      fe->reinit (elem);

      // Zero the element matrices before
      // summing them.  We use the resize member here because
      // the number of degrees of freedom might have changed from
      // the last element.  Note that this will be the case if the
      // element type is different (i.e. the last element was a
      // triangle, now we are on a quadrilateral).
      const unsigned int n_dofs =
        cast_int<unsigned int>(dof_indices.size());
      Ke.resize (n_dofs, n_dofs);
      Me.resize (n_dofs, n_dofs);

      // Now loop over the quadrature points.  This handles
      // the numeric integration.
      //
      // We will build the element matrix.  This involves
      // a double loop to integrate the test functions (i) against
      // the trial functions (j).
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        for (unsigned int i=0; i<n_dofs; i++)
          for (unsigned int j=0; j<n_dofs; j++)
            {
              Me(i,j) += JxW[qp]*phi[i][qp]*phi[j][qp];
              Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
            }

      // The calls to constrain_element_matrix below have no effect in
      // the current example. However, if users modify this example to
      // include hanging nodes due to mesh refinement, for example,
      // then it is essential to call constrain_element_matrix.
      // As a result we include constrain_element_matrix here to
      // ensure this example is ready to be used with hanging nodes.
      // (Note that constrained rows/cols will be eliminated from
      // the eigenproblem by the CondensedEigenSystem.)
      dof_map.constrain_element_matrix(Ke, dof_indices, false);
      dof_map.constrain_element_matrix(Me, dof_indices, false);

      // Finally, simply add the element contribution to the
      // overall matrices A and B.
      matrix_A.add_matrix (Ke, dof_indices);
      matrix_B.add_matrix (Me, dof_indices);
    } // end of element loop


#else
  // Avoid compiler warnings
  libmesh_ignore(es);
#endif // LIBMESH_HAVE_SLEPC
}
