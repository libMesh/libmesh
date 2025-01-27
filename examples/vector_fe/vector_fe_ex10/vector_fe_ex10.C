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


// <h1>Vector Finite Elements Example 10 - Raviart-Thomas elements (grad-div)</h1>
// \author Nuno Nobre
// \date 2025
//
// This example uses Raviart-Thomas elements to solve a model grad-div problem
// in H(div) in both 2d and 3d: -\nabla (\nabla \cdot \vec{u}) + \vec{u} = f.
// Note that, unlike the problem in Vector Finite Elements Example 6, this is
// _not_ an elliptic problem.

// Basic utilities.
#include "libmesh/string_to_enum.h"

// The solver packages supported by libMesh.
#include "libmesh/enum_solver_package.h"

// The mesh object and mesh generation and modification utilities.
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_modification.h"

// Matrix and vector types.
#include "libmesh/dense_matrix.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/numeric_vector.h"

// The finite element object and the geometric element type.
#include "libmesh/fe.h"
#include "libmesh/elem.h"

// Gauss quadrature rules.
#include "libmesh/quadrature_gauss.h"

// The dof map, which handles degree of freedom indexing.
#include "libmesh/dof_map.h"

// The system of equations.
#include "libmesh/equation_systems.h"
#include "libmesh/linear_implicit_system.h"

// The exact solution and error computation.
#include "libmesh/exact_solution.h"
#include "libmesh/enum_norm_type.h"
#include "solution_function.h"

// I/O utilities.
#include "libmesh/getpot.h"
#include "libmesh/exodusII_io.h"


// Bring in everything from the libMesh namespace.
using namespace libMesh;

// Function prototype. This is the function that will assemble
// the linear system for our grad-div problem. Note that the
// function will take the EquationSystems object and the
// name of the system we are assembling as input. From the
// EquationSystems object we have access to the Mesh and
// other objects we might need.
void assemble_graddiv(EquationSystems & es,
                      const std::string & system_name);

int main (int argc, char ** argv)
{
  // Initialize libMesh.
  LibMeshInit init (argc, argv);

  // This example requires a linear solver package.
  libmesh_example_requires(libMesh::default_solver_package() != INVALID_SOLVER_PACKAGE,
                           "--enable-petsc, --enable-trilinos, or --enable-eigen");

  // Parse the input file.
  GetPot infile("vector_fe_ex10.in");

  // But allow the command line to override it.
  infile.parse_command_line(argc, argv);

  // Read in parameters from the command line and the input file.
  const unsigned int dimension = infile("dim", 2);
  const unsigned int grid_size = infile("grid_size", 15);

  // Skip higher-dimensional examples on a lower-dimensional libMesh build.
  libmesh_example_requires(dimension <= LIBMESH_DIM, dimension << "D support");

  // Create a mesh, with dimension to be overridden later, distributed
  // across the default MPI communicator.
  Mesh mesh(init.comm());

  // Use the MeshTools::Generation mesh generator to create a uniform
  // grid on the cube [-1,1]^D. To accomodate Raviart-Thomas elements, we must
  // use TRI6/7 or QUAD8/9 elements in 2d, or TET14 or HEX27 in 3d.
  const std::string elem_str = infile("element_type", std::string("TRI6"));

  libmesh_error_msg_if((dimension == 2 && elem_str != "TRI6" && elem_str != "TRI7" && elem_str != "QUAD8" && elem_str != "QUAD9") ||
                       (dimension == 3 && elem_str != "TET14" && elem_str != "HEX27"),
                       "You selected " << elem_str <<
                       " but this example must be run with TRI6, TRI7, QUAD8, or QUAD9 in 2d" <<
                       " or with TET14, or HEX27 in 3d.");

  if (dimension == 2)
    MeshTools::Generation::build_square (mesh,
                                         grid_size,
                                         grid_size,
                                         -1., 1.,
                                         -1., 1.,
                                         Utility::string_to_enum<ElemType>(elem_str));
  else if (dimension == 3)
    MeshTools::Generation::build_cube (mesh,
                                       grid_size,
                                       grid_size,
                                       grid_size,
                                       -1., 1.,
                                       -1., 1.,
                                       -1., 1.,
                                       Utility::string_to_enum<ElemType>(elem_str));

  // Make sure the code is robust against nodal reorderings.
  MeshTools::Modification::permute_elements(mesh);

  // Print information about the mesh to the screen.
  mesh.print_info();

  // Create an equation systems object.
  EquationSystems equation_systems (mesh);

  // Declare the system  "GradDiv" and its variable.
  LinearImplicitSystem & system = equation_systems.add_system<LinearImplicitSystem>("GradDiv");

  // Set the FE approximation order for the vector field variable.
  const Order vector_order = static_cast<Order>(infile("order", 1u));

  libmesh_error_msg_if(vector_order < FIRST || vector_order > ((dimension == 3) ? FIRST : FIFTH),
                       "You selected: " << vector_order <<
                       " but this example must be run with either 1 <= order <= 5 in 2d"
                       " or with order 1 in 3d.");

  // Adds the variable "u" to "GradDiv". "u" will be our vector field.
  system.add_variable("u", vector_order, RAVIART_THOMAS);

  // Give the system a pointer to the matrix assembly
  // function. This will be called when needed by the library.
  system.attach_assemble_function(assemble_graddiv);

  // Initialize the data structures for the equation system.
  equation_systems.init();

  // Prints information about the system to the screen.
  equation_systems.print_info();

  // Solve the system "GradDiv". Note that calling this
  // member will assemble the linear system and invoke
  // the default numerical solver.
  system.solve();

  ExactSolution exact_sol(equation_systems);

  SolutionFunction soln_func;
  SolutionGradient soln_grad;

  // Build FunctionBase* containers to attach to the ExactSolution object.
  std::vector<FunctionBase<Number> *> sols(1, &soln_func);
  std::vector<FunctionBase<Gradient> *> grads(1, &soln_grad);

  exact_sol.attach_exact_values(sols);
  exact_sol.attach_exact_derivs(grads);

  // Use higher quadrature order for more accurate error results.
  int extra_error_quadrature = infile("extra_error_quadrature", 2);
  exact_sol.extra_quadrature_order(extra_error_quadrature);

  // Compute the error.
  exact_sol.compute_error("GradDiv", "u");

  // Print out the error values.
  libMesh::out << "~~ Vector field (u) ~~"
               << std::endl;
  libMesh::out << "L2 error is: "
               << exact_sol.l2_error("GradDiv", "u")
               << std::endl;
  libMesh::out << "HDiv semi-norm error is: "
               << exact_sol.error_norm("GradDiv", "u", HDIV_SEMINORM)
               << std::endl;
  libMesh::out << "HDiv error is: "
               << exact_sol.hdiv_error("GradDiv", "u")
               << std::endl;

#ifdef LIBMESH_HAVE_EXODUS_API

  // We write the file in the ExodusII format.
  ExodusII_IO(mesh).write_equation_systems("out.e", equation_systems);

#endif // #ifdef LIBMESH_HAVE_EXODUS_API

  // All done.
  return 0;
}



// We now define the matrix assembly function for the
// grad-div system. We need to first compute element
// matrices and right-hand sides, and then take into
// account the boundary conditions, which will be handled
// via a penalty method.
void assemble_graddiv(EquationSystems & es,
                      const std::string & libmesh_dbg_var(system_name))
{

  // It is a good idea to make sure we are assembling
  // the proper system.
  libmesh_assert_equal_to (system_name, "GradDiv");

  // Get a constant reference to the mesh object.
  const MeshBase & mesh = es.get_mesh();

  // The dimension that we are running.
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to the LinearImplicitSystem we are solving.
  LinearImplicitSystem & system = es.get_system<LinearImplicitSystem>("GradDiv");

  // A reference to the  DofMap object for this system.  The  DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.
  const DofMap & dof_map = system.get_dof_map();

  // Get a constant reference to the Finite Element type
  // for the variable in the system.
  FEType vector_fe_type = dof_map.variable_type(system.variable_number("u"));

  // Build the Finite Element object. Since the
  // FEBase::build() member dynamically creates memory we will
  // store the object as a std::unique_ptr<FEBase>. This can be thought
  // of as a pointer that will clean up after itself. Introduction Example 4
  // describes some advantages of  std::unique_ptr's in the context of
  // quadrature rules.
  std::unique_ptr<FEVectorBase> vector_fe (FEVectorBase::build(dim, vector_fe_type));

  // A just-high-enough Gauss quadrature rule for numerical integration.
  QGauss qrule (dim, vector_fe_type.default_quadrature_order());

  // Tell the finite element object to use our quadrature rule.
  vector_fe->attach_quadrature_rule (&qrule);

  // Declare a special finite element object for boundary integration.
  std::unique_ptr<FEVectorBase> vector_fe_face (FEVectorBase::build(dim, vector_fe_type));

  // Boundary integration requires one quadrature rule with dimensionality one
  // less than the dimensionality of the element.
  QGauss qface(dim-1, vector_fe_type.default_quadrature_order());

  // Tell the finite element object to use our quadrature rule.
  vector_fe_face->attach_quadrature_rule (&qface);

  // Here we define some references to cell-specific data that
  // will be used to assemble the linear system.
  //
  // The element Jacobian * quadrature weight at each integration point.
  const std::vector<Real> & JxW = vector_fe->get_JxW();

  // The physical XY locations of the quadrature points on the element.
  // These might be useful for evaluating spatially varying material
  // properties at the quadrature points.
  const std::vector<Point> & q_point = vector_fe->get_xyz();

  // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector<RealGradient>> & vector_phi = vector_fe->get_phi();

  // The divergence of the element vector shape functions evaluated at the
  // quadrature points.
  const std::vector<std::vector<Real>> & div_vector_phi = vector_fe->get_div_phi();

  // Define data structures to contain the element matrix
  // and right-hand-side vector contribution.  Following
  // basic finite element terminology we will denote these
  // "Ke" and "Fe".  These datatypes are templated on
  //  Number, which allows the same code to work for real
  // or complex numbers.
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  // These vectors will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;

  // The global system matrix
  SparseMatrix<Number> & matrix = system.get_system_matrix();

  // Now we will loop over all the elements in the mesh.
  // We will compute the element matrix and right-hand-side
  // contribution.
  //
  // Element ranges are a nice way to iterate through all the
  // elements, or all the elements that have some property.  The
  // range will iterate from the first to the last element on
  // the local processor.
  // It is smart to make this one const so that we don't accidentally
  // mess it up!  In case users later modify this program to include
  // refinement, we will be safe and will only consider the active
  // elements; hence we use a variant of the
  // active_local_element_ptr_range.
  for (const auto & elem : mesh.active_local_element_ptr_range())
    {
      // Get the degree of freedom indices for the
      // current element.  These define where in the global
      // matrix and right-hand-side this element will
      // contribute to.
      dof_map.dof_indices (elem, dof_indices);

      // Cache the total number of degrees of freedom on this element,
      // for use as array and loop bounds later.
      // We use cast_int to explicitly convert from size() (which may be
      // 64-bit) to unsigned int (which may be 32-bit but which is definitely
      // enough to count *local* degrees of freedom.
      const unsigned int n_dofs =
        cast_int<unsigned int>(dof_indices.size());

      // Compute the element-specific data for the current
      // element.  This involves computing the location of the
      // quadrature points (q_point) and the shape functions
      // and their divergences for the current element.
      vector_fe->reinit (elem);

      // We should also have the same number of degrees of freedom as
      // shape functions for our variable.
      libmesh_assert_equal_to (n_dofs, vector_phi.size());

      // Zero the element matrix and right-hand side before
      // summing them.  We use the resize member here because
      // the number of degrees of freedom might have changed from
      // the last element.  Note that this will be the case if the
      // element type is different (i.e. the last element was a
      // triangle, now we are on a quadrilateral).

      // The  DenseMatrix::resize() and the  DenseVector::resize()
      // members will automatically zero out the matrix  and vector.
      Ke.resize (n_dofs, n_dofs);
      Fe.resize (n_dofs);

      // Now loop over the quadrature points.  This handles
      // the numeric integration.
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {

          // Now we will build the element matrix.
          // This a double loop to integrate the vector test functions (i)
          // against the vector trial functions (j) and their divergences
          for (unsigned int i = 0; i != n_dofs; i++)
            for (unsigned int j = 0; j != n_dofs; j++)
              {
                Ke(i, j) += JxW[qp]*(div_vector_phi[i][qp]*div_vector_phi[j][qp]+
                                     vector_phi[i][qp]*vector_phi[j][qp]);
              }

          // This is the end of the matrix summation loop
          // Now we build the element right-hand-side contribution.
          // This involves a single loop in which we integrate the "forcing
          // function" in the PDE against the vector test functions (k).
          {
            // The location of the current quadrature point.
            const Real x = q_point[qp](0);
            const Real y = q_point[qp](1);
            const Real z = q_point[qp](2);

            // "f" is the forcing function, given by the well-known
            // "method of manufactured solutions".
            RealGradient f = GradDivExactSolution().forcing(x, y, z);

            // Loop to integrate the vector test functions (k) against the
            // forcing function.
            for (unsigned int k = 0; k != n_dofs; k++)
              {
                Fe(k) += JxW[qp]*f*vector_phi[k][qp];
              }
          }
        }

      // We have now reached the end of the quadrature point loop, so
      // the interior element integration has been completed. However, we have
      // not yet addressed boundary conditions.
      {

        // The following loop is over the sides of the element.
        // If the element has no neighbor on a side then that
        // side MUST live on a boundary of the domain.
        for (auto side : elem->side_index_range())
          if (elem->neighbor_ptr(side) == nullptr)
            {
              // The value of the shape functions at the quadrature points.
              const std::vector<std::vector<RealGradient>> & vector_phi_face = vector_fe_face->get_phi();

              // The Jacobian * Quadrature Weight at the quadrature
              // points on the face.
              const std::vector<Real> & JxW_face = vector_fe_face->get_JxW();

              // The XYZ locations (in physical space) of, and the normals at,
              // the quadrature points on the face.  This is where
              // we will interpolate the boundary value function.
              const std::vector<Point> & qface_point = vector_fe_face->get_xyz();
              const std::vector<Point> & normals = vector_fe_face->get_normals();

              // Compute the vector shape function values on the element face.
              vector_fe_face->reinit(elem, side);

              // Some shape functions will be 0 on the face, but for ease of
              // indexing and generality of code we loop over them anyway.
              libmesh_assert_equal_to (n_dofs, vector_phi_face.size());

              // Loop over the face quadrature points for integration.
              for (unsigned int qp=0; qp<qface.n_points(); qp++)
                {
                  // The location on the boundary of the current
                  // face quadrature point.
                  const Real xf = qface_point[qp](0);
                  const Real yf = qface_point[qp](1);
                  const Real zf = qface_point[qp](2);

                  // The boundary value for the vector variable.
                  RealGradient vector_value = GradDivExactSolution()(xf, yf, zf);

                  // We use the penalty method to set the flux of the vector
                  // variable at the boundary, i.e. the RT vector boundary dof.
                  const Real penalty = 1.e10;

                  // A double loop to integrate the normal component of the
                  // vector test functions (i) against the normal component of
                  // the vector trial functions (j).
                  for (unsigned int i = 0; i != n_dofs; i++)
                    for (unsigned int j = 0; j != n_dofs; j++)
                      {
                        Ke(i, j) += JxW_face[qp]*penalty*vector_phi_face[i][qp]*
                                    normals[qp]*vector_phi_face[j][qp]*normals[qp];
                      }

                  // Loop to integrate the normal component of the vector test
                  // functions (i) against the normal component of the
                  // exact solution for the vector variable.
                  for (unsigned int i = 0; i != n_dofs; i++)
                    {
                      Fe(i) += JxW_face[qp]*penalty*vector_phi_face[i][qp]*normals[qp]*
                                vector_value*normals[qp];
                    }
                }
            }
      }

      // We have now finished the quadrature point loop,
      // and have therefore applied all the boundary conditions.

      // If this assembly program were to be used on an adaptive mesh,
      // we would have to apply any hanging node constraint equations.
      dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

      // The element matrix and right-hand-side are now built
      // for this element.  Add them to the global matrix and
      // right-hand-side vector.  The  SparseMatrix::add_matrix()
      // and  NumericVector::add_vector() members do this for us.
      matrix.add_matrix (Ke, dof_indices);
      system.rhs->add_vector (Fe, dof_indices);
    }

  // All done!
}
