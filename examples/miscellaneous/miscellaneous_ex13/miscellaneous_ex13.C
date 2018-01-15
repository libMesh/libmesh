// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// <h1>Miscellaneous Example 12 - Quad8 Shell Elements</h1>
// \author Sylvain Vallaghe
// \date 2017
//
// This example implements shell elements using 8-noded serendipity quadrilaterals and reduced integration.
// Full integration Quad8 shell elements are prone to locking issues.
// The reduced integration version "eliminates locking in most situations
// although it introduces two spureous mechanisms. Fortunately,
// these mechanisms disappear after assembly of the global stiffness matrix and the element
// can be considered safe for practical purposes".
// (source: Onate, Structural Analysis with the Finite Element Method).
// The "pinched cylinder" problem is solved and the solution
// is compared to analytical values at selected points.

// C++ include files that we need
#include <iostream>

// LibMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature.h"
#include "libmesh/node.h"
#include "libmesh/elem.h"
#include "libmesh/dof_map.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/vector_value.h"
#include "libmesh/tensor_value.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dense_subvector.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/zero_function.h"
#include "libmesh/linear_solver.h"
#include "libmesh/libmesh_nullptr.h"
#include "libmesh/getpot.h"

// Eigen includes
#ifdef LIBMESH_HAVE_EIGEN
#include "libmesh/ignore_warnings.h"
# include <Eigen/Dense>
#include "libmesh/restore_warnings.h"
#endif

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Function prototype.  This is the function that will assemble
// the stiffness matrix and the right-hand-side vector ready
// for solution.
void assemble_shell (EquationSystems & es,
                     const std::string & system_name);

// Begin the main program.
int main (int argc, char ** argv)
{
  // Initialize libMesh.
  LibMeshInit init (argc, argv);

  // Skip this 3D example if libMesh was compiled as 1D/2D-only.
  libmesh_example_requires (3 == LIBMESH_DIM, "3D support");

  // Our input mesh here is in ExodusII format
#ifndef LIBMESH_HAVE_EXODUS_API
  libmesh_example_requires (false, "ExodusII support");
#endif

#ifndef LIBMESH_ENABLE_SECOND_DERIVATIVES
  libmesh_example_requires (false, "second derivatives enabled");
#endif

  // This example does a bunch of linear algebra during assembly, and
  // therefore requires Eigen.
#ifndef LIBMESH_HAVE_EIGEN
  libmesh_example_requires(false, "--enable-eigen");
#endif

  // This example converts between ExodusII and XDR files, therefore
  // it requires XDR support in libmesh.
#ifndef LIBMESH_HAVE_XDR
  libmesh_example_requires (false, "XDR support");
#endif

  // Read the "distributed_load" flag from the command line
  GetPot command_line (argc, argv);
  int distributed_load = 0;
  if (command_line.search(1, "-distributed_load"))
    distributed_load = command_line.next(distributed_load);

  {
    Mesh mesh (init.comm(), 3);

    // To confirm that both ExodusII and Xdr formats work for shell
    // meshes, we read in cylinder.exo, then write out cylinder.xdr,
    // then read in cylinder.exo again below and use that for the rest
    // of the example.
    mesh.read("cylinder.exo");
    mesh.write("cylinder.xdr");
  }
  Mesh mesh (init.comm(), 3);
  mesh.read("cylinder.xdr");

  // Print information about the mesh to the screen.
  mesh.print_info();

  // Create an equation systems object.
  EquationSystems equation_systems (mesh);

  // Declare the system and its variables.
  // Create a linear implicit system named "Shell".
  LinearImplicitSystem & system = equation_systems.add_system<LinearImplicitSystem> ("Shell");

  // Add the three displacement variables "u", "v", "w",
  // and the three rotational variables "theta_x", "theta_y", "theta_z".
  // All variables are second order.
  system.add_variable ("u",SECOND,LAGRANGE);
  system.add_variable ("v",SECOND,LAGRANGE);
  system.add_variable ("w",SECOND,LAGRANGE);
  system.add_variable ("theta_x",SECOND,LAGRANGE);
  system.add_variable ("theta_y",SECOND,LAGRANGE);
  system.add_variable ("theta_z",SECOND,LAGRANGE);

  // Give the system a pointer to the matrix and rhs assembly
  // function.
  system.attach_assemble_function (assemble_shell);

  // Use the parameters of the equation systems object to
  // tell the shell system about the material properties, the
  // shell thickness, and the external load.
  const Real h  = 0.03;
  const Real E  = 3e10;
  const Real nu = 0.3;
  const Real q  = 1;
  equation_systems.parameters.set<Real> ("thickness")       = h;
  equation_systems.parameters.set<Real> ("young's modulus") = E;
  equation_systems.parameters.set<Real> ("poisson ratio")   = nu;
  equation_systems.parameters.set<Real> ("point load")    = q;
  equation_systems.parameters.set<bool>("distributed load")  = (distributed_load != 0);

  // Dirichlet conditions for the pinched cylinder problem.
  // Only one 8th of the cylinder is considered using symmetry considerations.
  // The cylinder longitudinal axis is the y-axis.
  // The four corners of the surface are named A(3,0,0), B(3,3,0), C(0,3,3), D(0,0,3).
  // The point load (pinch) is applied at C in the -z direction.
  // Edge AD is the actual edge of the cylinder and is rigid in the xz-plane.
  // Other edges have symmetric boundary conditions.

  // AB w, theta_x, theta_y
  {
    std::set<boundary_id_type> boundary_ids;
    boundary_ids.insert(7);
    unsigned int variables[] = {2, 3, 4};
    ZeroFunction<> zf;

    // Most DirichletBoundary users will want to supply a "locally
    // indexed" functor
    DirichletBoundary dirichlet_bc 
      (boundary_ids,
       std::vector<unsigned int>(variables, variables+3), zf,
       LOCAL_VARIABLE_ORDER);
    system.get_dof_map().add_dirichlet_boundary(dirichlet_bc);
  }
  // BC v, theta_x, theta_z
  {
    std::set<boundary_id_type> boundary_ids;
    boundary_ids.insert(8);
    unsigned int variables[] = {1, 3, 5};
    ZeroFunction<> zf;

    DirichletBoundary dirichlet_bc
      (boundary_ids,
       std::vector<unsigned int>(variables, variables+3), zf,
       LOCAL_VARIABLE_ORDER);
    system.get_dof_map().add_dirichlet_boundary(dirichlet_bc);
  }
  // CD u, theta_y, theta_z
  {
    std::set<boundary_id_type> boundary_ids;
    boundary_ids.insert(9);
    unsigned int variables[] = {0, 4, 5};
    ZeroFunction<> zf;

    DirichletBoundary dirichlet_bc
      (boundary_ids,
       std::vector<unsigned int>(variables, variables+3), zf,
       LOCAL_VARIABLE_ORDER);
    system.get_dof_map().add_dirichlet_boundary(dirichlet_bc);
  }
  // AD u, w, theta_y
  {
    std::set<boundary_id_type> boundary_ids;
    boundary_ids.insert(10);
    unsigned int variables[] = {0, 2, 4};
    ZeroFunction<> zf;

    DirichletBoundary dirichlet_bc
      (boundary_ids,
       std::vector<unsigned int>(variables, variables+3), zf,
       LOCAL_VARIABLE_ORDER);
    system.get_dof_map().add_dirichlet_boundary(dirichlet_bc);
  }

  // Initialize the data structures for the equation system.
  equation_systems.init();

  // Print information about the system to the screen.
  equation_systems.print_info();

  // This example can be run with EigenSparseLinearSolvers, but it
  // only works with either the CG or SPARSELU types, and SparseLU
  // turns out to be faster.
  if (libMesh::default_solver_package() == EIGEN_SOLVERS)
    system.get_linear_solver()->set_solver_type(SPARSELU);

  // Solve the linear system.
  system.solve();

  // After solving the system, write the solution to an
  // ExodusII output file ready for import in, e.g.,
  // Paraview.
  ExodusII_IO(mesh).write_equation_systems ("out.e", equation_systems);

  // Compare with analytical solution for point load
  if (distributed_load==0)
    {
      // Find the node nearest point C.
      Node * node_C = libmesh_nullptr;
      Point point_C(0, 3, 3);
      {
        Real nearest_dist_sq = std::numeric_limits<Real>::max();

        // Find the closest local node.  On a DistributedMesh we may
        // not even know about the existence of closer non-local
        // nodes.
        for (const auto & node : mesh.local_node_ptr_range())
          {
            const Real dist_sq = (*node - point_C).norm_sq();
            if (dist_sq < nearest_dist_sq)
              {
                nearest_dist_sq = dist_sq;
                node_C = node;
              }
          }

        // Check with other processors to see if any found a closer node
        unsigned int minrank = 0;
        system.comm().minloc(nearest_dist_sq, minrank);

        // Broadcast the ID of the closest node, so every processor can
        // see for certain whether they have it or not.
        dof_id_type nearest_node_id;
        if (system.processor_id() == minrank)
          nearest_node_id = node_C->id();
        system.comm().broadcast(nearest_node_id, minrank);
        node_C = mesh.query_node_ptr(nearest_node_id);
      }

      // Evaluate the z-displacement "w" at the node nearest C.
      Number w = 0;

      // If we know about the closest node, and if we also own the DoFs
      // on that node, then we can evaluate the solution at that node.
      if (node_C)
        {
          const unsigned int w_var = system.variable_number ("w");
          dof_id_type w_dof = node_C->dof_number (system.number(), w_var, 0);
          if (w_dof >= system.get_dof_map().first_dof() &&
              w_dof < system.get_dof_map().end_dof())
            w = system.current_solution(w_dof);
        }
      system.comm().sum(w);


      Number w_C_bar = -E*h*w/q;
      const Real w_C_bar_analytic = 164.24;

      // Print the finite element solution and the analytic
      // prediction to the screen.
      libMesh::out << "z-displacement of the point C: " << w_C_bar << std::endl;
      libMesh::out << "Analytic solution: " << w_C_bar_analytic << std::endl;

      // Evaluate the y-displacement "v" at point D.  This time we'll
      // evaluate at the exact point, not just the closest node.
      Point point_D(0, 0, 3);
      const unsigned int v_var = system.variable_number ("v");
      Number v = system.point_value(v_var, point_D);

      Number v_D_bar = E*h*v/q;
      const Real v_D_bar_analytic = 4.114;

      // Print the finite element solution and the analytic
      // prediction to the screen.
      libMesh::out << "y-displacement of the point D: " << v_D_bar << std::endl;
      libMesh::out << "Analytic solution: " << v_D_bar_analytic << std::endl;
    }

  // All done.
  return 0;
}



// We now define the matrix and rhs vector assembly function
// for the shell system.
void assemble_shell (EquationSystems & es,
                     const std::string & system_name)
{
  // This example requires Eigen to actually work, but we should still
  // let it compile and throw a runtime error if you don't.
#ifdef LIBMESH_HAVE_EIGEN
  // It is a good idea to make sure we are assembling
  // the proper system.
  libmesh_assert_equal_to (system_name, "Shell");

  // Get a constant reference to the mesh object.
  const MeshBase & mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to the shell system object.
  LinearImplicitSystem & system = es.get_system<LinearImplicitSystem> (system_name);

  // Get the shell parameters that we need during assembly.
  const Real h  = es.parameters.get<Real> ("thickness");
  const Real E  = es.parameters.get<Real> ("young's modulus");
  const Real nu = es.parameters.get<Real> ("poisson ratio");
  const Real q  = es.parameters.get<Real> ("point load");
  const bool distributed_load  = es.parameters.get<bool> ("distributed load");

  // The membrane elastic matrix.
  Eigen::Matrix3d Hm;
  Hm <<
    1., nu, 0.,
    nu, 1., 0.,
    0., 0., 0.5 * (1-nu);
  Hm *= h * E/(1-nu*nu);

  // The bending elastic matrix.
  Eigen::Matrix3d Hf;
  Hf <<
    1., nu, 0.,
    nu, 1., 0.,
    0., 0., 0.5 * (1-nu);
  Hf *= h*h*h/12 * E/(1-nu*nu);

  // The shear elastic matrices.
  Eigen::Matrix2d Hc0 = Eigen::Matrix2d::Identity();
  Hc0 *= h * 5./6*E/(2*(1+nu));

  Eigen::Matrix2d Hc1 = Eigen::Matrix2d::Identity();
  Hc1 *= h*h*h/12 * 5./6*E/(2*(1+nu));

  // Get the Finite Element type, this will be
  // the same for all variables.
  FEType fe_type = system.variable_type (0);

  std::unique_ptr<FEBase> fe (FEBase::build(dim, fe_type));
  QGauss qrule (dim, SECOND); // Reduced integration, 2x2 gauss points (instead of 3x3 for full integration)
  fe->attach_quadrature_rule (&qrule);

  // The element Jacobian * quadrature weight at each integration point.
  const std::vector<Real> & JxW = fe->get_JxW();

  // The element shape function and its derivatives evaluated at the
  // quadrature points.
  const std::vector<RealGradient> & dxyzdxi = fe->get_dxyzdxi();
  const std::vector<RealGradient> & dxyzdeta = fe->get_dxyzdeta();
  const std::vector<RealGradient> & d2xyzdxi2 = fe->get_d2xyzdxi2();
  const std::vector<RealGradient> & d2xyzdeta2 = fe->get_d2xyzdeta2();
  const std::vector<RealGradient> & d2xyzdxideta = fe->get_d2xyzdxideta();
  const std::vector<std::vector<Real>> & dphidxi = fe->get_dphidxi();
  const std::vector<std::vector<Real>> & dphideta = fe->get_dphideta();
  const std::vector<std::vector<Real>> & phi = fe->get_phi();

  // A reference to the DofMap object for this system.  The DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.
  const DofMap & dof_map = system.get_dof_map();

  // Define data structures to contain the element stiffness matrix.
  DenseMatrix<Number> Ke;
  DenseSubMatrix<Number> Ke_var[6][6] =
    {
      {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke),
       DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)},
      {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke),
       DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)},
      {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke),
       DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)},
      {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke),
       DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)},
      {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke),
       DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)},
      {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke),
       DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)}
    };

  // Define data structures to contain the element rhs vector.
  DenseVector<Number> Fe;
  DenseSubVector<Number> Fe_w(Fe);

  std::vector<dof_id_type> dof_indices;
  std::vector<std::vector<dof_id_type>> dof_indices_var(6);

  // Now we will loop over all the elements in the mesh.  We will
  // compute the element matrix and right-hand-side contribution.
  for (const auto & elem : mesh.active_local_element_ptr_range())
    {
      dof_map.dof_indices (elem, dof_indices);
      for (unsigned int var=0; var<6; var++)
        dof_map.dof_indices (elem, dof_indices_var[var], var);

      const unsigned int n_dofs   = dof_indices.size();
      const unsigned int n_var_dofs = dof_indices_var[0].size();

      // First compute element data at the nodes
      std::vector<Point> nodes;
      for (unsigned int i=0; i<elem->n_nodes(); ++i)
        nodes.push_back(elem->reference_elem()->node_ref(i));
      fe->reinit (elem, &nodes);

      //Store local orthonormal basis at the nodes
      std::vector<Eigen::Matrix3d> Qnode;
      for (unsigned int i=0; i<elem->n_nodes(); ++i)
        {
          Eigen::Vector3d a1;
          a1 << dxyzdxi[i](0), dxyzdxi[i](1), dxyzdxi[i](2);
          Eigen::Vector3d a2;
          a2 << dxyzdeta[i](0), dxyzdeta[i](1), dxyzdeta[i](2);
          Eigen::Vector3d n;
          n = a1.cross(a2);
          n /= n.norm();

          Real nx = n(0);
          Real ny = n(1);
          Real C  = n(2);
          if (std::abs(1.+C)<1e-6)
            {
              Eigen::Matrix3d Q;
              Q <<
                1, 0, 0,
                0, -1, 0,
                0, 0, -1;
              Qnode.push_back(Q);
            }
          else
            {
              Eigen::Matrix3d Q;
              Q <<
                C+1./(1+C)*ny*ny, -1./(1+C)*nx*ny, nx,
                -1./(1+C)*nx*ny, C+1./(1+C)*nx*nx, ny,
                -nx,             -ny,              C;
              Qnode.push_back(Q);
            }
        }

      Ke.resize (n_dofs, n_dofs);
      for (unsigned int var_i=0; var_i<6; var_i++)
        for (unsigned int var_j=0; var_j<6; var_j++)
          Ke_var[var_i][var_j].reposition (var_i*n_var_dofs, var_j*n_var_dofs, n_var_dofs, n_var_dofs);

      Fe.resize(n_dofs);
      Fe_w.reposition(2*n_var_dofs,n_var_dofs);

      // Reinit element data at the regular Gauss quadrature points
      fe->reinit (elem);

      // Now we will build the element matrix and right-hand-side.
      for (unsigned int qp=0; qp<qrule.n_points(); ++qp)
        {

          //Covariant basis at the quadrature point
          Eigen::Vector3d a1;
          a1 << dxyzdxi[qp](0), dxyzdxi[qp](1), dxyzdxi[qp](2);
          Eigen::Vector3d a2;
          a2 << dxyzdeta[qp](0), dxyzdeta[qp](1), dxyzdeta[qp](2);
          Eigen::Vector3d n;
          n = a1.cross(a2);
          n /= n.norm();
          Eigen::Matrix3d F0;
          F0 <<
            a1(0), a2(0), n(0),
            a1(1), a2(1), n(1),
            a1(2), a2(2), n(2);

          //Contravariant basis
          Eigen::Matrix3d F0it;
          F0it = F0.inverse().transpose();

          //Local orthonormal basis at the quadrature point
          Real nx = n(0);
          Real ny = n(1);
          Real C  = n(2);
          Eigen::Matrix3d Q;
          if (std::abs(1.+C) < 1e-6)
            {
              Q <<
                1, 0, 0,
                0, -1, 0,
                0, 0, -1;
            }
          else
            {
              Q <<
                C+1./(1+C)*ny*ny, -1./(1+C)*nx*ny, nx,
                -1./(1+C)*nx*ny, C+1./(1+C)*nx*nx, ny,
                -nx,             -ny,              C;
            }

          Eigen::Matrix2d C0;
          C0 = F0it.block<3,2>(0,0).transpose()*Q.block<3,2>(0,0);

          // Normal derivatives in reference coordinates
          Eigen::Vector3d d2Xdxi2(d2xyzdxi2[qp](0), d2xyzdxi2[qp](1), d2xyzdxi2[qp](2));
          Eigen::Vector3d d2Xdeta2(d2xyzdeta2[qp](0), d2xyzdeta2[qp](1), d2xyzdeta2[qp](2));
          Eigen::Vector3d d2Xdxideta(d2xyzdxideta[qp](0), d2xyzdxideta[qp](1), d2xyzdxideta[qp](2));

          Eigen::Matrix2d b;
          b <<
            n.dot(d2Xdxi2), n.dot(d2Xdxideta),
            n.dot(d2Xdxideta), n.dot(d2Xdeta2);

          Eigen::Vector3d dndxi = -b(0,0)*F0it.col(0) - b(0,1)*F0it.col(1);
          Eigen::Vector3d dndeta = -b(1,0)*F0it.col(0) - b(1,1)*F0it.col(1);

          Eigen::Matrix2d bhat;
          bhat <<
            F0it.col(1).dot(dndeta), -F0it.col(0).dot(dndeta),
            -F0it.col(1).dot(dndxi), F0it.col(0).dot(dndxi);

          Eigen::Matrix2d bc;
          bc = bhat*C0;

          // Mean curvature
          Real H = 0.5*(dndxi.dot(F0it.col(0))+dndeta.dot(F0it.col(1)));

          // Loop over all pairs of nodes I,J.
          for (unsigned int i=0; i<n_var_dofs; ++i)
            {
              // Matrix B0, zeroth order (through thickness) membrane-bending strain
              Real C1i = dphidxi[i][qp]*C0(0,0) + dphideta[i][qp]*C0(1,0);
              Real C2i = dphidxi[i][qp]*C0(0,1) + dphideta[i][qp]*C0(1,1);

              Eigen::MatrixXd B0I(3, 5);
              B0I = Eigen::MatrixXd::Zero(3, 5);
              B0I.block<1,3>(0,0) = C1i*Q.col(0).transpose();
              B0I.block<1,3>(1,0) = C2i*Q.col(1).transpose();
              B0I.block<1,3>(2,0) = C2i*Q.col(0).transpose()+C1i*Q.col(1).transpose();

              // Matrix B1, first order membrane-bending strain
              Real bc1i = dphidxi[i][qp]*bc(0,0) + dphideta[i][qp]*bc(1,0);
              Real bc2i = dphidxi[i][qp]*bc(0,1) + dphideta[i][qp]*bc(1,1);

              Eigen::Vector2d V1i(-Q.col(0).dot(Qnode[i].col(1)),
                                  Q.col(0).dot(Qnode[i].col(0)));

              Eigen::Vector2d V2i(-Q.col(1).dot(Qnode[i].col(1)),
                                  Q.col(1).dot(Qnode[i].col(0)));

              Eigen::MatrixXd B1I(3,5);
              B1I = Eigen::MatrixXd::Zero(3,5);
              B1I.block<1,3>(0,0) = bc1i*Q.col(0).transpose();
              B1I.block<1,3>(1,0) = bc2i*Q.col(1).transpose();
              B1I.block<1,3>(2,0) = bc2i*Q.col(0).transpose()+bc1i*Q.col(1).transpose();

              B1I.block<1,2>(0,3) = C1i*V1i.transpose();
              B1I.block<1,2>(1,3) = C2i*V2i.transpose();
              B1I.block<1,2>(2,3) = C2i*V1i.transpose()+C1i*V2i.transpose();

              // Matrix B2, second order membrane-bending strain
              Eigen::MatrixXd B2I(3,5);
              B2I = Eigen::MatrixXd::Zero(3,5);

              B2I.block<1,2>(0,3) = bc1i*V1i.transpose();
              B2I.block<1,2>(1,3) = bc2i*V2i.transpose();
              B2I.block<1,2>(2,3) = bc2i*V1i.transpose()+bc1i*V2i.transpose();

              // Matrix Bc0, zeroth order shear strain
              Eigen::MatrixXd Bc0I(2,5);
              Bc0I = Eigen::MatrixXd::Zero(2,5);
              Bc0I.block<1,3>(0,0) = C1i*Q.col(2).transpose();
              Bc0I.block<1,3>(1,0) = C2i*Q.col(2).transpose();
              Bc0I.block<1,2>(0,3) = phi[i][qp]*V1i.transpose();
              Bc0I.block<1,2>(1,3) = phi[i][qp]*V2i.transpose();

              // Matrix Bc1, first order shear strain
              Eigen::MatrixXd Bc1I(2,5);
              Bc1I = Eigen::MatrixXd::Zero(2,5);
              Bc1I.block<1,3>(0,0) = bc1i*Q.col(2).transpose();
              Bc1I.block<1,3>(1,0) = bc2i*Q.col(2).transpose();

              // Drilling dof (in-plane rotation)
              Eigen::Vector2d BdxiI(dphidxi[i][qp],dphideta[i][qp]);
              Eigen::Vector2d BdI = C0.transpose()*BdxiI;

              for (unsigned int j=0; j<n_var_dofs; ++j)
                {

                  // Matrix B0, zeroth order membrane-bending strain
                  Real C1j = dphidxi[j][qp]*C0(0,0) + dphideta[j][qp]*C0(1,0);
                  Real C2j = dphidxi[j][qp]*C0(0,1) + dphideta[j][qp]*C0(1,1);

                  Eigen::MatrixXd B0J(3,5);
                  B0J = Eigen::MatrixXd::Zero(3,5);
                  B0J.block<1,3>(0,0) = C1j*Q.col(0).transpose();
                  B0J.block<1,3>(1,0) = C2j*Q.col(1).transpose();
                  B0J.block<1,3>(2,0) = C2j*Q.col(0).transpose()+C1j*Q.col(1).transpose();

                  // Matrix B1, first order membrane-bending strain
                  Real bc1j = dphidxi[j][qp]*bc(0,0) + dphideta[j][qp]*bc(1,0);
                  Real bc2j = dphidxi[j][qp]*bc(0,1) + dphideta[j][qp]*bc(1,1);

                  Eigen::Vector2d V1j(-Q.col(0).dot(Qnode[j].col(1)),
                                      Q.col(0).dot(Qnode[j].col(0)));

                  Eigen::Vector2d V2j(-Q.col(1).dot(Qnode[j].col(1)),
                                      Q.col(1).dot(Qnode[j].col(0)));

                  Eigen::MatrixXd B1J(3,5);
                  B1J = Eigen::MatrixXd::Zero(3,5);
                  B1J.block<1,3>(0,0) = bc1j*Q.col(0).transpose();
                  B1J.block<1,3>(1,0) = bc2j*Q.col(1).transpose();
                  B1J.block<1,3>(2,0) = bc2j*Q.col(0).transpose()+bc1j*Q.col(1).transpose();

                  B1J.block<1,2>(0,3) = C1j*V1j.transpose();
                  B1J.block<1,2>(1,3) = C2j*V2j.transpose();
                  B1J.block<1,2>(2,3) = C2j*V1j.transpose()+C1j*V2j.transpose();

                  // Matrix B2, second order membrane-bending strain
                  Eigen::MatrixXd B2J(3,5);
                  B2J = Eigen::MatrixXd::Zero(3,5);

                  B2J.block<1,2>(0,3) = bc1j*V1j.transpose();
                  B2J.block<1,2>(1,3) = bc2j*V2j.transpose();
                  B2J.block<1,2>(2,3) = bc2j*V1j.transpose()+bc1j*V2j.transpose();

                  // Matrix Bc0, zeroth order shear strain
                  Eigen::MatrixXd Bc0J(2, 5);
                  Bc0J = Eigen::MatrixXd::Zero(2,5);
                  Bc0J.block<1,3>(0,0) = C1j*Q.col(2).transpose();
                  Bc0J.block<1,3>(1,0) = C2j*Q.col(2).transpose();
                  Bc0J.block<1,2>(0,3) = phi[j][qp]*V1j.transpose();
                  Bc0J.block<1,2>(1,3) = phi[j][qp]*V2j.transpose();

                  // Matrix Bc1, first order shear strain
                  Eigen::MatrixXd Bc1J(2, 5);
                  Bc1J = Eigen::MatrixXd::Zero(2,5);
                  Bc1J.block<1,3>(0,0) = bc1j*Q.col(2).transpose();
                  Bc1J.block<1,3>(1,0) = bc2j*Q.col(2).transpose();

                  // Drilling dof
                  Eigen::Vector2d BdxiJ(dphidxi[j][qp], dphideta[j][qp]);
                  Eigen::Vector2d BdJ = C0.transpose()*BdxiJ;

                  // The total stiffness matrix coupling the nodes
                  // I and J is a sum of membrane, bending and shear contributions.
                  Eigen::MatrixXd local_KIJ(5, 5);
                  local_KIJ = JxW[qp] * (
                                         B0I.transpose() * Hm * B0J
                                         +  B2I.transpose() * Hf * B0J
                                         +  B0I.transpose() * Hf * B2J
                                         +  B1I.transpose() * Hf * B1J
                                         +  2*H * B0I.transpose() * Hf * B1J
                                         +  2*H * B1I.transpose() * Hf * B0J
                                         +  Bc0I.transpose() * Hc0 * Bc0J
                                         +  Bc1I.transpose() * Hc1 * Bc1J
                                         +  2*H * Bc0I.transpose() * Hc1 * Bc1J
                                         +  2*H * Bc1I.transpose() * Hc1 * Bc0J
                                         );

                  // Going from 5 to 6 dofs to add drilling dof
                  Eigen::MatrixXd full_local_KIJ(6, 6);
                  full_local_KIJ = Eigen::MatrixXd::Zero(6, 6);
                  full_local_KIJ.block<5,5>(0,0)=local_KIJ;

                  // Drilling dof stiffness contribution
                  full_local_KIJ(5,5) = Hf(0,0)*JxW[qp]*BdI.transpose()*BdJ;

                  // Transform the stiffness matrix to global coordinates
                  Eigen::MatrixXd global_KIJ(6,6);
                  Eigen::MatrixXd TI(6,6);
                  TI = Eigen::MatrixXd::Identity(6,6);
                  TI.block<3,3>(3,3) = Qnode[i].transpose();
                  Eigen::MatrixXd TJ(6,6);
                  TJ = Eigen::MatrixXd::Identity(6,6);
                  TJ.block<3,3>(3,3) = Qnode[j].transpose();
                  global_KIJ = TI.transpose()*full_local_KIJ*TJ;

                  // Insert the components of the coupling stiffness
                  // matrix KIJ into the corresponding directional
                  // submatrices.
                  for (unsigned int k=0;k<6;k++)
                    for (unsigned int l=0;l<6;l++)
                      Ke_var[k][l](i,j) += global_KIJ(k,l);
                }
            }

        } // end of the quadrature point qp-loop

      if (distributed_load)
        {
          // Loop on shell faces
          for (unsigned int shellface=0; shellface<2; shellface++)
            {
              std::vector<boundary_id_type> bids;
              mesh.get_boundary_info().shellface_boundary_ids(elem, shellface, bids);

              for (std::size_t k=0; k<bids.size(); k++)
                if (bids[k]==11) // sideset id for surface load
                  for (unsigned int qp=0; qp<qrule.n_points(); ++qp)
                    for (unsigned int i=0; i<n_var_dofs; ++i)
                      Fe_w(i) -= JxW[qp] * phi[i][qp];
            }
        }

      // The element matrix is now built for this element.
      // Add it to the global matrix.

      dof_map.constrain_element_matrix_and_vector (Ke,Fe,dof_indices);

      system.matrix->add_matrix (Ke, dof_indices);
      system.rhs->add_vector    (Fe, dof_indices);

    }

  if (!distributed_load)
    {
      //Adding point load to the RHS

      //Pinch position
      Point C(0, 3, 3);

      //Finish assembling rhs so we can set one value
      system.rhs->close();

      MeshBase::const_node_iterator nodeit = mesh.nodes_begin();
      const MeshBase::const_node_iterator node_end = mesh.nodes_end();

      for ( ; nodeit!=node_end; ++nodeit)
        {
          Node & node = **nodeit;
          if ((node-C).norm() < 1e-3)
            system.rhs->set(node.dof_number(0, 2, 0), -q/4);
        }
    }

#else
  // Avoid compiler warnings
  libmesh_ignore(es);
  libmesh_ignore(system_name);
#endif // LIBMESH_HAVE_EIGEN
}
