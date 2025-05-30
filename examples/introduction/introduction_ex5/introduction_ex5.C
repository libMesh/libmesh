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



// <h1>Introduction Example 5 - Run-Time Quadrature Rule Selection</h1>
// \author Benjamin S. Kirk
// \date 2003
//
// This is the fifth example program.  It builds on
// the previous two examples, and extends the use
// of the std::unique_ptr as a convenient build method to
// determine the quadrature rule at run time.


// C++ include files that we need
#include <iostream>
#include <sstream>
#include <algorithm>
#include <math.h>

// Basic include file needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/equation_systems.h"

// Define the Finite Element object.
#include "libmesh/fe.h"

// Define the base quadrature class, with which
// specialized quadrature rules will be built.
#include "libmesh/quadrature.h"

// Define useful datatypes for finite element
// matrix and vector components.
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"

// Define the DofMap, which handles degree of freedom
// indexing.
#include "libmesh/dof_map.h"

// To impose Dirichlet boundary conditions
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/analytic_function.h"

// The definition of a geometric element
#include "libmesh/elem.h"
#include "libmesh/enum_solver_package.h"

// Reading a quadrature rule from user arguments
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/string_to_enum.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;



// Function prototype, as before.
void assemble_poisson(EquationSystems & es,
                      const std::string & system_name);

// Exact solution function prototype, as before.
Real exact_solution (const Real x,
                     const Real y,
                     const Real z = 0.);

// Define a wrapper for exact_solution that will be needed below
void exact_solution_wrapper (DenseVector<Number> & output,
                             const Point & p,
                             const Real)
{
  output(0) = exact_solution(p(0), p(1), p(2));
}


// The quadrature type the user requests.
QuadratureType quad_type=INVALID_Q_RULE;



// Begin the main program.
int main (int argc, char ** argv)
{
  // Initialize libMesh and any dependent libraries, like in example 2.
  LibMeshInit init (argc, argv);

  // This example requires a linear solver package.
  libmesh_example_requires(libMesh::default_solver_package() != INVALID_SOLVER_PACKAGE,
                           "--enable-petsc, --enable-trilinos, or --enable-eigen");

  // Check for proper usage.  The quadrature rule
  // must be given at run time.
  libmesh_error_msg_if(argc < 3,
                       "Usage: " << argv[0] << " -q <rule>\n"
                       "  where <rule> is one of QGAUSS, QSIMPSON, or QTRAP.");

  // Tell the user what we are doing.
  libMesh::out << "Running " << argv[0];

  for (int i=1; i<argc; i++)
    libMesh::out << " " << argv[i];

  libMesh::out << std::endl << std::endl;

  // Set the quadrature rule type that the user wants
  quad_type = Utility::string_to_enum<QuadratureType>
    (libMesh::command_line_next("-q", std::string("QGAUSS")));

  // Skip this 3D example if libMesh was compiled as 1D-only.
  libmesh_example_requires(3 <= LIBMESH_DIM, "3D support");

  // We use Dirichlet boundary conditions here
#ifndef LIBMESH_ENABLE_DIRICHLET
  libmesh_example_requires(false, "--enable-dirichlet");
#endif

  // The following is identical to example 4, and therefore
  // not commented.  Differences are mentioned when present.
  Mesh mesh(init.comm());

  // We will use a linear approximation space in this example,
  // hence 8-noded hexahedral elements are sufficient.  This
  // is different than example 4 where we used 27-noded
  // hexahedral elements to support a second-order approximation
  // space.
  MeshTools::Generation::build_cube (mesh,
                                     16, 16, 16,
                                     -1., 1.,
                                     -1., 1.,
                                     -1., 1.,
                                     HEX8);

  mesh.print_info();

  EquationSystems equation_systems (mesh);

  equation_systems.add_system<LinearImplicitSystem> ("Poisson");

  unsigned int u_var = equation_systems.get_system("Poisson").add_variable("u", FIRST);

  equation_systems.get_system("Poisson").attach_assemble_function (assemble_poisson);

  // Construct a Dirichlet boundary condition object

  // Indicate which boundary IDs we impose the BC on
  // We either build a line, a square or a cube, and
  // here we indicate boundaries covering each case
  std::set<boundary_id_type> boundary_ids {0,1,2,3,4,5};

  // Create an AnalyticFunction object that we use to project the BC
  // This function just calls the function exact_solution via exact_solution_wrapper
  AnalyticFunction<> exact_solution_object(exact_solution_wrapper);

#ifdef LIBMESH_ENABLE_DIRICHLET
  // In general, when reusing a system-indexed exact solution, we want
  // to use the default system-ordering constructor for
  // DirichletBoundary, so we demonstrate that here.  In this case,
  // though, we have only one variable, so system- and local-
  // orderings are the same.
  DirichletBoundary dirichlet_bc
    (boundary_ids, {u_var}, exact_solution_object);

  // We must add the Dirichlet boundary condition _before_
  // we call equation_systems.init()
  equation_systems.get_system("Poisson").get_dof_map().add_dirichlet_boundary(dirichlet_bc);
#endif

  equation_systems.init();

  equation_systems.print_info();

  equation_systems.get_system("Poisson").solve();

  // "Personalize" the output, with the
  // number of the quadrature rule appended.
  std::ostringstream f_name;
  f_name << "out_" << quad_type << ".e";

#ifdef LIBMESH_HAVE_EXODUS_API
  ExodusII_IO(mesh).write_equation_systems (f_name.str(),
                                            equation_systems);
#endif // #ifdef LIBMESH_HAVE_EXODUS_API

  // All done.
  return 0;
}




void assemble_poisson(EquationSystems & es,
                      const std::string & libmesh_dbg_var(system_name))
{
  libmesh_assert_equal_to (system_name, "Poisson");

  const MeshBase & mesh = es.get_mesh();

  const unsigned int dim = mesh.mesh_dimension();

  LinearImplicitSystem & system = es.get_system<LinearImplicitSystem>("Poisson");

  const DofMap & dof_map = system.get_dof_map();

  FEType fe_type = dof_map.variable_type(0);

  // Build a Finite Element object of the specified type.  Since the
  // FEBase::build() member dynamically creates memory we will
  // store the object as a std::unique_ptr<FEBase>.  Below, the
  // functionality of std::unique_ptr's is described more detailed in
  // the context of building quadrature rules.
  std::unique_ptr<FEBase> fe (FEBase::build(dim, fe_type));

  // Now this deviates from example 4.  we create a
  // 5th order quadrature rule of user-specified type
  // for numerical integration.  Note that not all
  // quadrature rules support this order.
  std::unique_ptr<QBase> qrule(QBase::build(quad_type, dim, THIRD));

  // Tell the finite element object to use our
  // quadrature rule.  Note that a std::unique_ptr<QBase> returns
  // a QBase* pointer to the object it handles with get().
  // However, using get(), the std::unique_ptr<QBase> qrule is
  // still in charge of this pointer. I.e., when qrule goes
  // out of scope, it will safely delete the QBase object it
  // points to.  This behavior may be overridden using
  // std::unique_ptr<Xyz>::release(), but is currently not
  // recommended.
  fe->attach_quadrature_rule (qrule.get());

  // Declare a special finite element object for
  // boundary integration.
  std::unique_ptr<FEBase> fe_face (FEBase::build(dim, fe_type));

  // As already seen in example 3, boundary integration
  // requires a quadrature rule.  Here, however,
  // we use the more convenient way of building this
  // rule at run-time using quad_type.  Note that one
  // could also have initialized the face quadrature rules
  // with the type directly determined from qrule, namely
  // through:
  // \verbatim
  // std::unique_ptr<QBase>  qface (QBase::build(qrule->type(),
  // dim-1,
  // THIRD));
  // \endverbatim
  // And again: using the std::unique_ptr<QBase> relaxes
  // the need to delete the object afterward,
  // they clean up themselves.
  std::unique_ptr<QBase>  qface (QBase::build(quad_type,
                                              dim-1,
                                              THIRD));

  // Tell the finite element object to use our
  // quadrature rule.  Note that a std::unique_ptr<QBase> returns
  // a QBase* pointer to the object it handles with get().
  // However, using get(), the std::unique_ptr<QBase> qface is
  // still in charge of this pointer. I.e., when qface goes
  // out of scope, it will safely delete the QBase object it
  // points to.  This behavior may be overridden using
  // std::unique_ptr<Xyz>::release(), but is not recommended.
  fe_face->attach_quadrature_rule (qface.get());

  // This is again identical to example 4, and not commented.
  const std::vector<Real> & JxW = fe->get_JxW();

  const std::vector<Point> & q_point = fe->get_xyz();

  const std::vector<std::vector<Real>> & phi = fe->get_phi();

  const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();

  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;
  std::vector<dof_id_type> dof_indices;

  // The global system matrix
  SparseMatrix<Number> & matrix = system.get_system_matrix();

  // Now we will loop over all the elements in the mesh.
  // See example 3 for details.
  for (const auto & elem : mesh.active_local_element_ptr_range())
    {
      dof_map.dof_indices (elem, dof_indices);

      const unsigned int n_dofs =
        cast_int<unsigned int>(dof_indices.size());

      fe->reinit (elem);

      libmesh_assert_equal_to (n_dofs, phi.size());

      Ke.resize (n_dofs, n_dofs);

      Fe.resize (n_dofs);

      // Now loop over the quadrature points.  This handles
      // the numeric integration.  Note the slightly different
      // access to the QBase members!
      for (unsigned int qp=0; qp<qrule->n_points(); qp++)
        {
          // Add the matrix contribution
          for (unsigned int i=0; i != n_dofs; i++)
            for (unsigned int j=0; j != n_dofs; j++)
              Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);

          // fxy is the forcing function for the Poisson equation.
          // In this case we set fxy to be a finite difference
          // Laplacian approximation to the (known) exact solution.
          //
          // We will use the second-order accurate FD Laplacian
          // approximation, which in 2D on a structured grid is
          //
          // u_xx + u_yy = (u(i-1,j) + u(i+1,j) +
          //                u(i,j-1) + u(i,j+1) +
          //                -4*u(i,j))/h^2
          //
          // Since the value of the forcing function depends only
          // on the location of the quadrature point (q_point[qp])
          // we will compute it here, outside of the i-loop
          const Real x = q_point[qp](0);
          const Real y = q_point[qp](1);
          const Real z = q_point[qp](2);
          const Real eps = 1.e-3;

          const Real uxx = (exact_solution(x-eps, y, z) +
                            exact_solution(x+eps, y, z) +
                            -2.*exact_solution(x, y, z))/eps/eps;

          const Real uyy = (exact_solution(x, y-eps, z) +
                            exact_solution(x, y+eps, z) +
                            -2.*exact_solution(x, y, z))/eps/eps;

          const Real uzz = (exact_solution(x, y, z-eps) +
                            exact_solution(x, y, z+eps) +
                            -2.*exact_solution(x, y, z))/eps/eps;

          const Real fxy = - (uxx + uyy + ((dim==2) ? 0. : uzz));


          // Add the RHS contribution
          for (unsigned int i=0; i != n_dofs; i++)
            Fe(i) += JxW[qp]*fxy*phi[i][qp];
        }

      // If this assembly program were to be used on an adaptive mesh,
      // we would have to apply any hanging node constraint equations
      // Call heterogenously_constrain_element_matrix_and_vector to impose
      // non-homogeneous Dirichlet BCs
      dof_map.heterogenously_constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

      // The element matrix and right-hand-side are now built
      // for this element.  Add them to the global matrix and
      // right-hand-side vector.  The SparseMatrix::add_matrix()
      // and NumericVector::add_vector() members do this for us.
      matrix.add_matrix (Ke, dof_indices);
      system.rhs->add_vector    (Fe, dof_indices);

    } // end of element loop
}
