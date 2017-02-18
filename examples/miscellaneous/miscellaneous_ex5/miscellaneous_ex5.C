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



// <h1>Miscellaneous Example 5 - Interior Penalty Discontinuous Galerkin</h1>
// \author Lorenzo Botti
// \date 2010
//
// This example is based on Adaptivity Example 3, but uses an
// Interior Penalty Discontinuous Galerkin formulation.

#include <iostream>

// LibMesh include files.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_modification.h"
#include "libmesh/elem.h"
#include "libmesh/transient_system.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/fe_interface.h"
#include "libmesh/getpot.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/error_vector.h"
#include "libmesh/kelly_error_estimator.h"
#include "libmesh/discontinuity_measure.h"
#include "libmesh/string_to_enum.h"

#include "libmesh/exact_solution.h"
//#define QORDER TWENTYSIXTH

// Bring in everything from the libMesh namespace
using namespace libMesh;

Number exact_solution (const Point & p,
                       const Parameters & parameters,
                       const std::string &,
                       const std::string &)
{
  const Real x = p(0);
  const Real y = p(1);
  const Real z = p(2);

  if (parameters.get<bool>("singularity"))
    {
      Real theta = atan2(y, x);

      if (theta < 0)
        theta += 2. * libMesh::pi;

      return pow(x*x + y*y, 1./3.)*sin(2./3.*theta) + z;
    }
  else
    {
      return cos(x) * exp(y) * (1. - z);
    }
}

// We now define the gradient of the exact solution, again being careful
// to obtain an angle from atan2 in the correct
// quadrant.
Gradient exact_derivative(const Point & p,
                          const Parameters & parameters,  // es parameters
                          const std::string &,            // sys_name, not needed
                          const std::string &)            // unk_name, not needed
{
  // Gradient value to be returned.
  Gradient gradu;

  // x and y coordinates in space
  const Real x = p(0);
  const Real y = p(1);
  const Real z = p(2);

  if (parameters.get<bool>("singularity"))
    {
      libmesh_assert_not_equal_to (x, 0.);

      // For convenience...
      const Real tt = 2./3.;
      const Real ot = 1./3.;

      // The value of the radius, squared
      const Real r2 = x*x + y*y;

      // The boundary value, given by the exact solution,
      // u_exact = r^(2/3)*sin(2*theta/3).
      Real theta = atan2(y, x);

      // Make sure 0 <= theta <= 2*pi
      if (theta < 0)
        theta += 2. * libMesh::pi;

      // du/dx
      gradu(0) = tt*x*pow(r2,-tt)*sin(tt*theta) - pow(r2,ot)*cos(tt*theta)*tt/(1.+y*y/x/x)*y/x/x;
      gradu(1) = tt*y*pow(r2,-tt)*sin(tt*theta) + pow(r2,ot)*cos(tt*theta)*tt/(1.+y*y/x/x)*1./x;
      gradu(2) = 1.;
    }
  else
    {
      gradu(0) = -sin(x) * exp(y) * (1. - z);
      gradu(1) = cos(x) * exp(y) * (1. - z);
      gradu(2) = -cos(x) * exp(y);
    }
  return gradu;
}

// We now define the matrix assembly function for the
// Laplace system.  We need to first compute element volume
// matrices, and then take into account the boundary
// conditions and the flux integrals, which will be handled
// via an interior penalty method.
void assemble_ellipticdg(EquationSystems & es,
                         const std::string & libmesh_dbg_var(system_name))
{
  libMesh::out << " assembling elliptic dg system... ";
  libMesh::out.flush();

  // It is a good idea to make sure we are assembling
  // the proper system.
  libmesh_assert_equal_to (system_name, "EllipticDG");

  // Get a constant reference to the mesh object.
  const MeshBase & mesh = es.get_mesh();
  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to the LinearImplicitSystem we are solving
  LinearImplicitSystem & ellipticdg_system = es.get_system<LinearImplicitSystem> ("EllipticDG");
  // Get some parameters that we need during assembly
  const Real penalty = es.parameters.get<Real> ("penalty");
  std::string refinement_type = es.parameters.get<std::string> ("refinement");

  // A reference to the DofMap object for this system.  The DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.  We will talk more about the DofMap
  const DofMap & dof_map = ellipticdg_system.get_dof_map();

  // Get a constant reference to the Finite Element type
  // for the first (and only) variable in the system.
  FEType fe_type = ellipticdg_system.variable_type(0);

  // Build a Finite Element object of the specified type.  Since the
  // FEBase::build() member dynamically creates memory we will
  // store the object as a UniquePtr<FEBase>.  This can be thought
  // of as a pointer that will clean up after itself.
  UniquePtr<FEBase> fe  (FEBase::build(dim, fe_type));
  UniquePtr<FEBase> fe_elem_face(FEBase::build(dim, fe_type));
  UniquePtr<FEBase> fe_neighbor_face(FEBase::build(dim, fe_type));

  // Quadrature rules for numerical integration.
#ifdef QORDER
  QGauss qrule (dim, QORDER);
#else
  QGauss qrule (dim, fe_type.default_quadrature_order());
#endif
  fe->attach_quadrature_rule (&qrule);

#ifdef QORDER
  QGauss qface(dim-1, QORDER);
#else
  QGauss qface(dim-1, fe_type.default_quadrature_order());
#endif

  // Tell the finite element object to use our quadrature rule.
  fe_elem_face->attach_quadrature_rule(&qface);
  fe_neighbor_face->attach_quadrature_rule(&qface);

  // Here we define some references to cell-specific data that
  // will be used to assemble the linear system.
  // Data for interior volume integrals
  const std::vector<Real> & JxW = fe->get_JxW();
  const std::vector<std::vector<RealGradient> > & dphi = fe->get_dphi();

  // Data for surface integrals on the element boundary
  const std::vector<std::vector<Real> > &  phi_face = fe_elem_face->get_phi();
  const std::vector<std::vector<RealGradient> > & dphi_face = fe_elem_face->get_dphi();
  const std::vector<Real> & JxW_face = fe_elem_face->get_JxW();
  const std::vector<Point> & qface_normals = fe_elem_face->get_normals();
  const std::vector<Point> & qface_points = fe_elem_face->get_xyz();

  // Data for surface integrals on the neighbor boundary
  const std::vector<std::vector<Real> > &  phi_neighbor_face = fe_neighbor_face->get_phi();
  const std::vector<std::vector<RealGradient> > & dphi_neighbor_face = fe_neighbor_face->get_dphi();

  // Define data structures to contain the element interior matrix
  // and right-hand-side vector contribution.  Following
  // basic finite element terminology we will denote these
  // "Ke" and "Fe".
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  // Data structures to contain the element and neighbor boundary matrix
  // contribution. This matrices will do the coupling beetwen the dofs of
  // the element and those of his neighbors.
  // Ken: matrix coupling elem and neighbor dofs
  DenseMatrix<Number> Kne;
  DenseMatrix<Number> Ken;
  DenseMatrix<Number> Kee;
  DenseMatrix<Number> Knn;

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;

  // Now we will loop over all the elements in the mesh.  We will
  // compute first the element interior matrix and right-hand-side contribution
  // and then the element and neighbors boundary matrix contributions.
  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
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
      const unsigned int n_dofs = dof_indices.size();

      // Compute the element-specific data for the current
      // element.  This involves computing the location of the
      // quadrature points (q_point) and the shape functions
      // (phi, dphi) for the current element.
      fe->reinit (elem);

      // Zero the element matrix and right-hand side before
      // summing them.  We use the resize member here because
      // the number of degrees of freedom might have changed from
      // the last element.
      Ke.resize (n_dofs, n_dofs);
      Fe.resize (n_dofs);

      // Now we will build the element interior matrix.  This involves
      // a double loop to integrate the test funcions (i) against
      // the trial functions (j).
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        for (unsigned int i=0; i<n_dofs; i++)
          for (unsigned int j=0; j<n_dofs; j++)
            Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);

      // Now we adress boundary conditions.
      // We consider Dirichlet bc imposed via the interior penalty method
      // The following loops over the sides of the element.
      // If the element has no neighbor on a side then that
      // side MUST live on a boundary of the domain.
      for (unsigned int side=0; side<elem->n_sides(); side++)
        {
          if (elem->neighbor_ptr(side) == libmesh_nullptr)
            {
              // Pointer to the element face
              fe_elem_face->reinit(elem, side);

              UniquePtr<const Elem> elem_side (elem->build_side_ptr(side));
              // h elemet dimension to compute the interior penalty penalty parameter
              const unsigned int elem_b_order = static_cast<unsigned int> (fe_elem_face->get_order());
              const double h_elem = elem->volume()/elem_side->volume() * 1./pow(elem_b_order, 2.);

              for (unsigned int qp=0; qp<qface.n_points(); qp++)
                {
                  Number bc_value = exact_solution(qface_points[qp], es.parameters, "null", "void");
                  for (unsigned int i=0; i<n_dofs; i++)
                    {
                      // Matrix contribution
                      for (unsigned int j=0; j<n_dofs; j++)
                        {
                          // stability
                          Ke(i,j) += JxW_face[qp] * penalty/h_elem * phi_face[i][qp] * phi_face[j][qp];

                          // consistency
                          Ke(i,j) -=
                            JxW_face[qp] *
                            (phi_face[i][qp] * (dphi_face[j][qp]*qface_normals[qp]) +
                             phi_face[j][qp] * (dphi_face[i][qp]*qface_normals[qp]));
                        }

                      // RHS contributions

                      // stability
                      Fe(i) += JxW_face[qp] * bc_value * penalty/h_elem * phi_face[i][qp];

                      // consistency
                      Fe(i) -= JxW_face[qp] * dphi_face[i][qp] * (bc_value*qface_normals[qp]);
                    }
                }
            }

          // If the element is not on a boundary of the domain
          // we loop over his neighbors to compute the element
          // and neighbor boundary matrix contributions
          else
            {
              // Store a pointer to the neighbor we are currently
              // working on.
              const Elem * neighbor = elem->neighbor_ptr(side);

              // Get the global id of the element and the neighbor
              const unsigned int elem_id = elem->id();
              const unsigned int neighbor_id = neighbor->id();

              // If the neighbor has the same h level and is active
              // perform integration only if our global id is bigger than our neighbor id.
              // We don't want to compute twice the same contributions.
              // If the neighbor has a different h level perform integration
              // only if the neighbor is at a lower level.
              if ((neighbor->active() &&
                   (neighbor->level() == elem->level()) &&
                   (elem_id < neighbor_id)) ||
                  (neighbor->level() < elem->level()))
                {
                  // Pointer to the element side
                  UniquePtr<const Elem> elem_side (elem->build_side_ptr(side));

                  // h dimension to compute the interior penalty penalty parameter
                  const unsigned int elem_b_order = static_cast<unsigned int>(fe_elem_face->get_order());
                  const unsigned int neighbor_b_order = static_cast<unsigned int>(fe_neighbor_face->get_order());
                  const double side_order = (elem_b_order + neighbor_b_order)/2.;
                  const double h_elem = (elem->volume()/elem_side->volume()) * 1./pow(side_order,2.);

                  // The quadrature point locations on the neighbor side
                  std::vector<Point> qface_neighbor_point;

                  // The quadrature point locations on the element side
                  std::vector<Point > qface_point;

                  // Reinitialize shape functions on the element side
                  fe_elem_face->reinit(elem, side);

                  // Get the physical locations of the element quadrature points
                  qface_point = fe_elem_face->get_xyz();

                  // Find their locations on the neighbor
                  unsigned int side_neighbor = neighbor->which_neighbor_am_i(elem);
                  if (refinement_type == "p")
                    fe_neighbor_face->side_map (neighbor,
                                                elem_side.get(),
                                                side_neighbor,
                                                qface.get_points(),
                                                qface_neighbor_point);
                  else
                    FEInterface::inverse_map (elem->dim(),
                                              fe->get_fe_type(),
                                              neighbor,
                                              qface_point,
                                              qface_neighbor_point);

                  // Calculate the neighbor element shape functions at those locations
                  fe_neighbor_face->reinit(neighbor, &qface_neighbor_point);

                  // Get the degree of freedom indices for the
                  // neighbor.  These define where in the global
                  // matrix this neighbor will contribute to.
                  std::vector<dof_id_type> neighbor_dof_indices;
                  dof_map.dof_indices (neighbor, neighbor_dof_indices);
                  const unsigned int n_neighbor_dofs = neighbor_dof_indices.size();

                  // Zero the element and neighbor side matrix before
                  // summing them.  We use the resize member here because
                  // the number of degrees of freedom might have changed from
                  // the last element or neighbor.
                  // Note that Kne and Ken are not square matrices if neighbor
                  // and element have a different p level
                  Kne.resize (n_neighbor_dofs, n_dofs);
                  Ken.resize (n_dofs, n_neighbor_dofs);
                  Kee.resize (n_dofs, n_dofs);
                  Knn.resize (n_neighbor_dofs, n_neighbor_dofs);

                  // Now we will build the element and neighbor
                  // boundary matrices.  This involves
                  // a double loop to integrate the test funcions
                  // (i) against the trial functions (j).
                  for (unsigned int qp=0; qp<qface.n_points(); qp++)
                    {
                      // Kee Matrix. Integrate the element test function i
                      // against the element test function j
                      for (unsigned int i=0; i<n_dofs; i++)
                        {
                          for (unsigned int j=0; j<n_dofs; j++)
                            {
                              // consistency
                              Kee(i,j) -=
                                0.5 * JxW_face[qp] *
                                (phi_face[j][qp]*(qface_normals[qp]*dphi_face[i][qp]) +
                                 phi_face[i][qp]*(qface_normals[qp]*dphi_face[j][qp]));

                              // stability
                              Kee(i,j) += JxW_face[qp] * penalty/h_elem * phi_face[j][qp]*phi_face[i][qp];
                            }
                        }

                      // Knn Matrix. Integrate the neighbor test function i
                      // against the neighbor test function j
                      for (unsigned int i=0; i<n_neighbor_dofs; i++)
                        {
                          for (unsigned int j=0; j<n_neighbor_dofs; j++)
                            {
                              // consistency
                              Knn(i,j) +=
                                0.5 * JxW_face[qp] *
                                (phi_neighbor_face[j][qp]*(qface_normals[qp]*dphi_neighbor_face[i][qp]) +
                                 phi_neighbor_face[i][qp]*(qface_normals[qp]*dphi_neighbor_face[j][qp]));

                              // stability
                              Knn(i,j) +=
                                JxW_face[qp] * penalty/h_elem * phi_neighbor_face[j][qp]*phi_neighbor_face[i][qp];
                            }
                        }

                      // Kne Matrix. Integrate the neighbor test function i
                      // against the element test function j
                      for (unsigned int i=0; i<n_neighbor_dofs; i++)
                        {
                          for (unsigned int j=0; j<n_dofs; j++)
                            {
                              // consistency
                              Kne(i,j) +=
                                0.5 * JxW_face[qp] *
                                (phi_neighbor_face[i][qp]*(qface_normals[qp]*dphi_face[j][qp]) -
                                 phi_face[j][qp]*(qface_normals[qp]*dphi_neighbor_face[i][qp]));

                              // stability
                              Kne(i,j) -= JxW_face[qp] * penalty/h_elem * phi_face[j][qp]*phi_neighbor_face[i][qp];
                            }
                        }

                      // Ken Matrix. Integrate the element test function i
                      // against the neighbor test function j
                      for (unsigned int i=0; i<n_dofs; i++)
                        {
                          for (unsigned int j=0; j<n_neighbor_dofs; j++)
                            {
                              // consistency
                              Ken(i,j) +=
                                0.5 * JxW_face[qp] *
                                (phi_neighbor_face[j][qp]*(qface_normals[qp]*dphi_face[i][qp]) -
                                 phi_face[i][qp]*(qface_normals[qp]*dphi_neighbor_face[j][qp]));

                              // stability
                              Ken(i,j) -= JxW_face[qp] * penalty/h_elem * phi_face[i][qp]*phi_neighbor_face[j][qp];
                            }
                        }
                    }

                  // The element and neighbor boundary matrix are now built
                  // for this side.  Add them to the global matrix
                  // The SparseMatrix::add_matrix() members do this for us.
                  ellipticdg_system.matrix->add_matrix(Kne, neighbor_dof_indices, dof_indices);
                  ellipticdg_system.matrix->add_matrix(Ken, dof_indices, neighbor_dof_indices);
                  ellipticdg_system.matrix->add_matrix(Kee, dof_indices);
                  ellipticdg_system.matrix->add_matrix(Knn, neighbor_dof_indices);
                }
            }
        }
      // The element interior matrix and right-hand-side are now built
      // for this element.  Add them to the global matrix and
      // right-hand-side vector.  The SparseMatrix::add_matrix()
      // and NumericVector::add_vector() members do this for us.
      ellipticdg_system.matrix->add_matrix(Ke, dof_indices);
      ellipticdg_system.rhs->add_vector(Fe, dof_indices);
    }

  libMesh::out << "done" << std::endl;
}



int main (int argc, char** argv)
{
  LibMeshInit init(argc, argv);

  // Skip adaptive examples on a non-adaptive libMesh build
#ifndef LIBMESH_ENABLE_AMR
  libmesh_example_requires(false, "--enable-amr");
#else

  //Parse the input file
  GetPot input_file("miscellaneous_ex5.in");

  //Read in parameters from the input file
  const unsigned int adaptive_refinement_steps = input_file("max_adaptive_r_steps", 3);
  const unsigned int uniform_refinement_steps  = input_file("uniform_h_r_steps", 3);
  const Real refine_fraction                   = input_file("refine_fraction", 0.5);
  const Real coarsen_fraction                  = input_file("coarsen_fraction", 0.);
  const unsigned int max_h_level               = input_file("max_h_level", 10);
  const std::string refinement_type            = input_file("refinement_type","p");
  Order p_order                                = static_cast<Order>(input_file("p_order", 1));
  const std::string element_type               = input_file("element_type", "tensor");
  const Real penalty                           = input_file("ip_penalty", 10.);
  const bool singularity                       = input_file("singularity", true);
  const unsigned int dim                       = input_file("dimension", 3);

  // Skip higher-dimensional examples on a lower-dimensional libMesh build
  libmesh_example_requires(dim <= LIBMESH_DIM, "2D/3D support");


  // Create a mesh, with dimension to be overridden later, distributed
  // across the default MPI communicator.
  Mesh mesh(init.comm());

  if (dim == 1)
    MeshTools::Generation::build_line(mesh, 1, -1., 0.);
  else if (dim == 2)
    mesh.read("lshaped.xda");
  else
    mesh.read ("lshaped3D.xda");

  // Use triangles if the config file says so
  if (element_type == "simplex")
    MeshTools::Modification::all_tri(mesh);

  // Mesh Refinement object
  MeshRefinement mesh_refinement(mesh);
  mesh_refinement.refine_fraction() = refine_fraction;
  mesh_refinement.coarsen_fraction() = coarsen_fraction;
  mesh_refinement.max_h_level() = max_h_level;

  // Do uniform refinement
  for (unsigned int rstep=0; rstep<uniform_refinement_steps; rstep++)
    mesh_refinement.uniformly_refine(1);

  // Crate an equation system object
  EquationSystems equation_system (mesh);

  // Set parameters for the equation system and the solver
  equation_system.parameters.set<Real>("linear solver tolerance") = TOLERANCE * TOLERANCE;
  equation_system.parameters.set<unsigned int>("linear solver maximum iterations") = 1000;
  equation_system.parameters.set<Real>("penalty") = penalty;
  equation_system.parameters.set<bool>("singularity") = singularity;
  equation_system.parameters.set<std::string>("refinement") = refinement_type;

  // Create a system named ellipticdg
  LinearImplicitSystem & ellipticdg_system = equation_system.add_system<LinearImplicitSystem> ("EllipticDG");

  // Add a variable "u" to "ellipticdg" using the p_order specified in the config file
  if (on_command_line("element_type"))
    {
      std::string fe_str =
        command_line_value(std::string("element_type"),
                           std::string("MONOMIAL"));

      if (fe_str != "MONOMIAL" || fe_str != "XYZ")
        libmesh_error_msg("Error: This example must be run with MONOMIAL or XYZ element types.");

      ellipticdg_system.add_variable ("u", p_order, Utility::string_to_enum<FEFamily>(fe_str));
    }
  else
    ellipticdg_system.add_variable ("u", p_order, MONOMIAL);

  // Give the system a pointer to the matrix assembly function
  ellipticdg_system.attach_assemble_function (assemble_ellipticdg);

  // Initialize the data structures for the equation system
  equation_system.init();

  // Construct ExactSolution object and attach solution functions
  ExactSolution exact_sol(equation_system);
  exact_sol.attach_exact_value(exact_solution);
  exact_sol.attach_exact_deriv(exact_derivative);

  // A refinement loop.
  for (unsigned int rstep=0; rstep<adaptive_refinement_steps; ++rstep)
    {
      libMesh::out << "  Beginning Solve " << rstep << std::endl;
      libMesh::out << "Number of elements: " << mesh.n_elem() << std::endl;

      // Solve the system
      ellipticdg_system.solve();

      libMesh::out << "System has: "
                   << equation_system.n_active_dofs()
                   << " degrees of freedom."
                   << std::endl;

      libMesh::out << "Linear solver converged at step: "
                   << ellipticdg_system.n_linear_iterations()
                   << ", final residual: "
                   << ellipticdg_system.final_linear_residual()
                   << std::endl;

      // Compute the error
      exact_sol.compute_error("EllipticDG", "u");

      // Print out the error values
      libMesh::out << "L2-Error is: "
                   << exact_sol.l2_error("EllipticDG", "u")
                   << std::endl;

      // Possibly refine the mesh
      if (rstep+1 < adaptive_refinement_steps)
        {
          // The ErrorVector is a particular StatisticsVector
          // for computing error information on a finite element mesh.
          ErrorVector error;

          // The discontinuity error estimator
          // evaluate the jump of the solution
          // on elements faces
          DiscontinuityMeasure error_estimator;
          error_estimator.estimate_error(ellipticdg_system,error);

          // Take the error in error and decide which elements will be coarsened or refined
          mesh_refinement.flag_elements_by_error_fraction(error);
          if (refinement_type == "p")
            mesh_refinement.switch_h_to_p_refinement();
          if (refinement_type == "hp")
            mesh_refinement.add_p_to_h_refinement();

          // Refine and coarsen the flagged elements
          mesh_refinement.refine_and_coarsen_elements();
          equation_system.reinit();
        }
    }

  // Write out the solution
  // After solving the system write the solution
  // to a ExodusII-formatted plot file.
#ifdef LIBMESH_HAVE_EXODUS_API
  ExodusII_IO (mesh).write_discontinuous_exodusII("lshaped_dg.e", equation_system);
#endif

#endif // #ifndef LIBMESH_ENABLE_AMR

  // All done.
  return 0;
}
