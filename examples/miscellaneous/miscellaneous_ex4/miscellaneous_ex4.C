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



// <h1>Miscellaneous Example 4 - Using a shell matrix</h1>
// \author Tim Kroger
// \date 2008
//
// This example solves the equation
//
// \f$-\Delta u+\int u = 1\f$
//
// with homogeneous Dirichlet boundary conditions.  This system has
// a full system matrix which can be written as the sum of of sparse
// matrix and a rank 1 matrix.  The shell matrix concept is used to
// solve this problem.
//
// The problem is solved in parallel on a non-uniform grid in order
// to demonstrate all the techniques that are required for this.
// The grid is fixed, however, i.e. no adaptive mesh refinement is
// used, so that the example remains simple.
//
// The example is 2d; extension to 3d is straight forward.

// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <cstdlib> // *must* precede <cmath> for proper std:abs() on PGI, Sun Studio CC
#include <cmath>

// Basic include file needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/vtk_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/sum_shell_matrix.h"
#include "libmesh/tensor_shell_matrix.h"
#include "libmesh/sparse_shell_matrix.h"
#include "libmesh/mesh_refinement.h"

#include "libmesh/getpot.h"

// This example will solve a linear transient system,
// so we need to include the TransientLinearImplicitSystem definition.
#include "libmesh/transient_system.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/vector_value.h"

// The definition of a geometric element
#include "libmesh/elem.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Function prototype.  This function will assemble the system matrix
// and right-hand-side.
void assemble (EquationSystems & es,
               const std::string & system_name);

// Begin the main program.  Note that the first
// statement in the program throws an error if
// you are in complex number mode, since this
// example is only intended to work with real
// numbers.
int main (int argc, char ** argv)
{
  // Initialize libMesh.
  LibMeshInit init (argc, argv);

#if !defined(LIBMESH_ENABLE_AMR)
  libmesh_example_requires(false, "--enable-amr");
#else
  libmesh_example_requires(libMesh::default_solver_package() == PETSC_SOLVERS, "--enable-petsc");

  // Brief message to the user regarding the program name
  // and command line arguments.

  libMesh::out << "Running: " << argv[0];

  for (int i=1; i<argc; i++)
    libMesh::out << " " << argv[i];

  libMesh::out << std::endl << std::endl;

  // Skip this 2D example if libMesh was compiled as 1D-only.
  libmesh_example_requires(2 <= LIBMESH_DIM, "2D support");

  // Create a mesh, with dimension to be overridden later, distributed
  // across the default MPI communicator.
  Mesh mesh(init.comm());

  // Create an equation systems object.
  EquationSystems equation_systems (mesh);

  MeshTools::Generation::build_square (mesh,
                                       16,
                                       16,
                                       -1., 1.,
                                       -1., 1.,
                                       QUAD4);

  LinearImplicitSystem & system =
    equation_systems.add_system<LinearImplicitSystem>
    ("System");

  // Adds the variable "u" to "System".  "u"
  // will be approximated using first-order approximation.
  system.add_variable ("u", FIRST);

  // Also, we need to add two vectors.  The tensor matrix v*w^T of
  // these two vectors will be part of the system matrix.
  system.add_vector("v", false);
  system.add_vector("w", false);

  // We need an additional matrix to be used for preconditioning since
  // a shell matrix is not suitable for that.
  system.add_matrix("Preconditioner");

  // Give the system a pointer to the matrix assembly function.
  system.attach_assemble_function (assemble);

  // Initialize the data structures for the equation system.
  equation_systems.init ();

  // Prints information about the system to the screen.
  equation_systems.print_info();

  equation_systems.parameters.set<unsigned int>
    ("linear solver maximum iterations") = 250;
  equation_systems.parameters.set<Real>
    ("linear solver tolerance") = TOLERANCE;

  // Refine arbitrarily some elements.
  for (unsigned int i=0; i<2; i++)
    {
      MeshRefinement mesh_refinement(mesh);
      MeshBase::element_iterator       elem_it  = mesh.elements_begin();
      const MeshBase::element_iterator elem_end = mesh.elements_end();
      for (; elem_it != elem_end; ++elem_it)
        {
          Elem * elem = *elem_it;
          if (elem->active())
            {
              if ((elem->id()%20)>8)
                {
                  elem->set_refinement_flag(Elem::REFINE);
                }
              else
                {
                  elem->set_refinement_flag(Elem::DO_NOTHING);
                }
            }
          else
            {
              elem->set_refinement_flag(Elem::INACTIVE);
            }
        }
      mesh_refinement.refine_elements();
      equation_systems.reinit();
    }

  // Prints information about the system to the screen.
  equation_systems.print_info();

  // Before the assemblation of the matrix, we have to clear the two
  // vectors that form the tensor matrix (since this is not performed
  // automatically).
  system.get_vector("v").init(system.n_dofs(), system.n_local_dofs());
  system.get_vector("w").init(system.n_dofs(), system.n_local_dofs());

  // We need a shell matrix to solve.  There is currently no way to
  // store the shell matrix in the system.  We just create it locally
  // here (a shell matrix does not occupy much memory).
  SumShellMatrix<Number> shellMatrix(system.comm());
  TensorShellMatrix<Number> shellMatrix0(system.get_vector("v"), system.get_vector("w"));
  shellMatrix.matrices.push_back(&shellMatrix0);
  SparseShellMatrix<Number> shellMatrix1(*system.matrix);
  shellMatrix.matrices.push_back(&shellMatrix1);

  // Attach that to the system.
  system.attach_shell_matrix(&shellMatrix);

  // Reset the preconditioning matrix to zero (for the system matrix,
  // the same thing is done automatically).
  system.get_matrix("Preconditioner").zero();

  // Assemble & solve the linear system
  system.solve();

  // Detach the shell matrix from the system since it will go out of
  // scope.  Nobody should solve the system outside this function.
  system.detach_shell_matrix();

  // Print a nice message.
  libMesh::out << "Solved linear system in "
               << system.n_linear_iterations()
               << " iterations, residual norm is "
               << system.final_linear_residual()
               << "."
               << std::endl;

#if defined(LIBMESH_HAVE_VTK) && !defined(LIBMESH_ENABLE_PARMESH)
  // Write result to file.
  VTKIO(mesh).write_equation_systems ("out.pvtu", equation_systems);
#endif // #ifdef LIBMESH_HAVE_VTK

#endif // #ifndef LIBMESH_ENABLE_AMR

  return 0;
}



// This function defines the assembly routine.  It is responsible for
// computing the proper matrix entries for the element stiffness
// matrices and right-hand sides.
void assemble (EquationSystems & es,
               const std::string & libmesh_dbg_var(system_name))
{
#ifdef LIBMESH_ENABLE_AMR
  // It is a good idea to make sure we are assembling
  // the proper system.
  libmesh_assert_equal_to (system_name, "System");

  // Get a constant reference to the mesh object.
  const MeshBase & mesh = es.get_mesh();

  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to the Convection-Diffusion system object.
  LinearImplicitSystem & system =
    es.get_system<LinearImplicitSystem> ("System");

  // Get the Finite Element type for the first (and only)
  // variable in the system.
  FEType fe_type = system.variable_type(0);

  // Build a Finite Element object of the specified type.  Since the
  // FEBase::build() member dynamically creates memory we will
  // store the object as a UniquePtr<FEBase>.  This can be thought
  // of as a pointer that will clean up after itself.
  UniquePtr<FEBase> fe      (FEBase::build(dim, fe_type));
  UniquePtr<FEBase> fe_face (FEBase::build(dim, fe_type));

  // A Gauss quadrature rule for numerical integration.
  // Let the FEType object decide what order rule is appropriate.
  QGauss qrule (dim,   fe_type.default_quadrature_order());
  QGauss qface (dim-1, fe_type.default_quadrature_order());

  // Tell the finite element object to use our quadrature rule.
  fe->attach_quadrature_rule      (&qrule);
  fe_face->attach_quadrature_rule (&qface);

  // Here we define some references to cell-specific data that
  // will be used to assemble the linear system.  We will start
  // with the element Jacobian * quadrature weight at each integration point.
  const std::vector<Real> & JxW      = fe->get_JxW();
  const std::vector<Real> & JxW_face = fe_face->get_JxW();

  // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector<Real> > & phi = fe->get_phi();
  const std::vector<std::vector<Real> > & psi = fe_face->get_phi();

  // The element shape function gradients evaluated at the quadrature
  // points.
  const std::vector<std::vector<RealGradient> > & dphi = fe->get_dphi();

  // The XY locations of the quadrature points used for face integration
  //const std::vector<Point>& qface_points = fe_face->get_xyz();

  // A reference to the DofMap object for this system.  The DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.  We will talk more about the DofMap
  // in future examples.
  const DofMap & dof_map = system.get_dof_map();

  // Define data structures to contain the element matrix
  // and right-hand-side vector contribution.  Following
  // basic finite element terminology we will denote these
  // "Ke" and "Fe".
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  // Analogous data structures for thw two vectors v and w that form
  // the tensor shell matrix.
  DenseVector<Number> Ve;
  DenseVector<Number> We;

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;

  // Now we will loop over all the elements in the mesh that
  // live on the local processor. We will compute the element
  // matrix and right-hand-side contribution.  Since the mesh
  // will be refined we want to only consider the ACTIVE elements,
  // hence we use a variant of the active_elem_iterator.
  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
    {
      // Store a pointer to the element we are currently
      // working on.  This allows for nicer syntax later.
      const Elem * elem = *el;

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

      // Zero the element matrix and right-hand side before
      // summing them.  We use the resize member here because
      // the number of degrees of freedom might have changed from
      // the last element.  Note that this will be the case if the
      // element type is different (i.e. the last element was a
      // triangle, now we are on a quadrilateral).
      Ke.resize (dof_indices.size(),
                 dof_indices.size());

      Fe.resize (dof_indices.size());
      Ve.resize (dof_indices.size());
      We.resize (dof_indices.size());

      // Now we will build the element matrix and right-hand-side.
      // Constructing the RHS requires the solution and its
      // gradient from the previous timestep.  This myst be
      // calculated at each quadrature point by summing the
      // solution degree-of-freedom values by the appropriate
      // weight functions.
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
          // Now compute the element matrix and RHS contributions.
          for (std::size_t i=0; i<phi.size(); i++)
            {
              // The RHS contribution
              Fe(i) += JxW[qp]*phi[i][qp];

              for (std::size_t j=0; j<phi.size(); j++)
                {
                  // The matrix contribution
                  Ke(i,j) += JxW[qp]*(
                                      // Stiffness matrix
                                      (dphi[i][qp]*dphi[j][qp])
                                      );
                }

              // V and W are the same for this example.
              Ve(i) += JxW[qp]*phi[i][qp];
              We(i) += JxW[qp]*phi[i][qp];
            }
        }

      // At this point the interior element integration has
      // been completed.  However, we have not yet addressed
      // boundary conditions.  For this example we will only
      // consider simple Dirichlet boundary conditions imposed
      // via the penalty method.
      //
      // The following loops over the sides of the element.
      // If the element has no neighbor on a side then that
      // side MUST live on a boundary of the domain.
      {
        // The penalty value.
        const Real penalty = 1.e10;

        // The following loops over the sides of the element.
        // If the element has no neighbor on a side then that
        // side MUST live on a boundary of the domain.
        for (unsigned int s=0; s<elem->n_sides(); s++)
          if (elem->neighbor_ptr(s) == libmesh_nullptr)
            {
              fe_face->reinit(elem, s);

              for (unsigned int qp=0; qp<qface.n_points(); qp++)
                {
                  // Matrix contribution
                  for (std::size_t i=0; i<psi.size(); i++)
                    for (std::size_t j=0; j<psi.size(); j++)
                      Ke(i,j) += penalty*JxW_face[qp]*psi[i][qp]*psi[j][qp];
                }
            }
      }


      // We have now built the element matrix and RHS vector in terms
      // of the element degrees of freedom.  However, it is possible
      // that some of the element DOFs are constrained to enforce
      // solution continuity, i.e. they are not really "free".  We need
      // to constrain those DOFs in terms of non-constrained DOFs to
      // ensure a continuous solution.  The
      // DofMap::constrain_element_matrix_and_vector() method does
      // just that.

      // However, constraining both the sparse matrix (and right hand
      // side) plus the rank 1 matrix is tricky.  The dof_indices
      // vector has to be backuped for that because the constraining
      // functions modify it.

      std::vector<dof_id_type> dof_indices_backup(dof_indices);
      dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);
      dof_indices = dof_indices_backup;
      dof_map.constrain_element_dyad_matrix(Ve, We, dof_indices);

      // The element matrix and right-hand-side are now built
      // for this element.  Add them to the global matrix and
      // right-hand-side vector.  The SparseMatrix::add_matrix()
      // and NumericVector::add_vector() members do this for us.
      system.matrix->add_matrix (Ke, dof_indices);
      system.get_matrix("Preconditioner").add_matrix (Ke, dof_indices);
      system.rhs->add_vector (Fe, dof_indices);
      system.get_vector("v").add_vector(Ve, dof_indices);
      system.get_vector("w").add_vector(We, dof_indices);
    }
  // Finished computing the sytem matrix and right-hand side.

  // Matrices and vectors must be closed manually.  This is necessary
  // because the matrix is not directly used as the system matrix (in
  // which case the solver closes it) but as a part of a shell matrix.
  system.matrix->close();
  system.get_matrix("Preconditioner").close();
  system.rhs->close();
  system.get_vector("v").close();
  system.get_vector("w").close();

#endif // #ifdef LIBMESH_ENABLE_AMR
}
