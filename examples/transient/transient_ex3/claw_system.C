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

// local includes
#include "claw_system.h"

// LibMesh includes
#include "libmesh/dense_matrix.h"
#include "libmesh/dof_map.h"
#include "libmesh/equation_systems.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/fe.h" // FEBase::buid()
#include "libmesh/fe_interface.h" // FEInterface::inverse_map()
#include "libmesh/getpot.h" // GetPot input file parsing
#include "libmesh/libmesh_logging.h" // LOG_SCOPE
#include "libmesh/linear_solver.h" // LinearSolver::reuse_preconditioner()
#include "libmesh/mesh_base.h" // MeshBase::active_local_element_ptr_range()
#include "libmesh/numeric_vector.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/sparse_matrix.h"

namespace libMesh
{

ClawSystem::ClawSystem (EquationSystems & es,
                        const std::string & name,
                        const unsigned int number)
  : Parent(es, name, number),
    _mass_matrix(SparseMatrix<Number>::build(es.comm())),
    _jump_matrix(SparseMatrix<Number>::build(es.comm())),
    _temporal_discretization_type(ForwardEuler),
    _delta_t(0.),
    _n_time_steps(0),
    _LxF_constant(1.),
    _write_interval(1),
    _time(0.)
{
  // set assemble_before_solve flag to false
  // so that we control matrix assembly.
  assemble_before_solve = false;

  // Allocate SparseMatrices. I could not figure out how to do
  // this in the initialization list, I don't think it's possible.
  for (auto i : make_range(2))
  {
    _advection_matrices.push_back(SparseMatrix<Number>::build(es.comm()));
    _avg_matrices.push_back(SparseMatrix<Number>::build(es.comm()));
    _boundary_condition_matrices.push_back(SparseMatrix<Number>::build(es.comm()));
  }

}


ClawSystem::~ClawSystem() = default;

std::string ClawSystem::system_type () const
{
  return "ClawSystem";
}

ClawSystem::TemporalDiscretizationType
ClawSystem::string_to_enum(std::string string_type)
{
  if (string_type == "ForwardEuler")
    return ForwardEuler;

  if (string_type == "RK4")
    return RK4;

  // If we got here, there must be an error
  libmesh_error_msg("Error: Invalid temporal discretization");
}

std::string
ClawSystem::enum_to_string(ClawSystem::TemporalDiscretizationType enum_type)
{
  switch (enum_type)
  {
    case ForwardEuler:
      return "ForwardEuler";

    case RK4:
      return "RK4";

    default:
      libmesh_error_msg("Error: Invalid temporal discretization");
  }
}

void ClawSystem::solve_conservation_law()
{
  // Set the ImplicitSystem "matrix" equal to the mass matrix. This is
  // the only part of the explicit method that requires a "solve()",
  // however we save the LU factorization from this and apply it
  // during subsequent solve() steps.
  this->matrix->zero();
  this->matrix->close();
  this->matrix->add(1., *_mass_matrix);

  this->set_time(0.);

  // For the first solve, make sure we generate a new preconditioner
  linear_solver->reuse_preconditioner(false);

  auto old_solution = NumericVector<Number>::build(this->comm());
  old_solution->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);

  // Temporary and "stage" vectors - only initialized if doing RK timestepping
  std::vector<std::unique_ptr<NumericVector<Number>>> k(4);
  std::unique_ptr<NumericVector<Number>> temp;
  if (_temporal_discretization_type == RK4)
  {
    temp = NumericVector<Number>::build(this->comm());
    temp->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);

    for (auto & vec : k)
    {
      vec = NumericVector<Number>::build(this->comm());
      vec->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);
    }
  }

  // Begin time stepping loop
  for (unsigned int step=1; step<=_n_time_steps; step++)
  {
    *old_solution = *solution;

    switch (_temporal_discretization_type)
    {
      case ForwardEuler:
      {
        this->assemble_claw_rhs(*old_solution);

        // For explicit timestepping schemes, the solve() step
        // "inverts" the mass matrix and then back substitutes the rhs
        // to update the solution. In practice, the mass matrix
        // factorization is computed only once and then the factors
        // are re-used at each step.
        this->solve();

        // Tell the system to reuse the preconditioner (for the mass matrix)
        if (step == 1)
          linear_solver->reuse_preconditioner(true);

        solution->scale(_delta_t);
        solution->add(*old_solution);

        break;
      }

      case RK4:
      {
        // Compute stage 1 (explicit Euler step)
        this->assemble_claw_rhs(*old_solution);
        this->solve();
        *k[0] = *solution;

        // Tell the system to reuse the preconditioner (for the mass matrix)
        if (step==1)
          linear_solver->reuse_preconditioner(true);

        // Compute stage 2
        *temp = *old_solution;
        temp->add(0.5*_delta_t, *k[0]);
        this->assemble_claw_rhs(*temp);
        this->solve();
        *k[1] = *solution;

        // Compute stage 3
        *temp = *old_solution;
        temp->add(0.5*_delta_t, *k[1]);
        this->assemble_claw_rhs(*temp);
        this->solve();
        *k[2] = *solution;

        // Compute stage 4
        *temp = *old_solution;
        temp->add(_delta_t, *k[2]);
        this->assemble_claw_rhs(*temp);
        this->solve();
        *k[3] = *solution;

        // The updated solution is now a linear combination of the 4 "stages"
        // u_{n+1} = u_n + dt*(k1/6 + k2/3 + k3/3 + k4/6)
        *solution = *old_solution;
        solution->add(_delta_t/6., *k[0]);
        solution->add(_delta_t/3., *k[1]);
        solution->add(_delta_t/3., *k[2]);
        solution->add(_delta_t/6., *k[3]);

        break;
      }

      default:
        libmesh_error_msg("Error: Invalid temporal discretization in ClawSystem::solve_conservation_law");
    }

    this->update();

    this->set_time(_time + _delta_t);

#ifdef LIBMESH_HAVE_EXODUS_API
    if (step % get_write_interval() == 0)
      {
        libMesh::out << std::endl << "plotting time step " << step << ", time = " << _time << std::endl;

        std::ostringstream file_name;
        file_name << "claw_solution."
                  << std::setw(4)
                  << std::right
                  << std::setfill('0')
                  << step
                  << ".e";

        // We are writing to a separate file so we write each step as step 1.
        // It would be better to append each new solution to an existing file...
        ExodusII_IO(get_mesh()).write_timestep(
          file_name.str(), this->get_equation_systems(), /*step=*/1, _time);
      }
#endif
  } // end time stepping loop

  linear_solver->reuse_preconditioner(false);
}

void ClawSystem::process_parameters_file (const std::string& parameters_filename)
{
  // First read in data from parameters_filename
  GetPot infile(parameters_filename);

  const std::string temporal_discretization_type_in = infile("temporal_discretization_type", "ForwardEuler");
  const Real LxF_constant_in           = infile("LxF_constant", 1.);
  const Real delta_t_in                = infile("delta_t", 0.05);
  const unsigned int n_time_steps_in   = infile("n_time_steps", 0);
  const unsigned int write_interval_in = infile("write_interval", 1);

  this->set_temporal_discretization_type(ClawSystem::string_to_enum(temporal_discretization_type_in));
  this->set_LxF_constant(LxF_constant_in);
  this->set_delta_t(delta_t_in);
  this->set_n_time_steps(n_time_steps_in);
  this->set_write_interval(write_interval_in);
}

void ClawSystem::print_info()
{
  // Print out info that describes the current setup
  libMesh::out << std::endl << "==== ClawSystem ====" << std::endl;
  libMesh::out << "system name: " << this->name() << std::endl;
  libMesh::out << "LxF constant: " << get_LxF_constant() << std::endl;
  libMesh::out << "delta_t: " << get_delta_t() << std::endl;
  libMesh::out << "n_time_steps: " << get_n_time_steps() << std::endl;
  libMesh::out << "temporal_discretization_type: " << ClawSystem::enum_to_string(this->get_temporal_discretization_type()) << std::endl;
  libMesh::out << "write_interval: " << get_write_interval() << std::endl;
  libMesh::out << std::endl;
}

void ClawSystem::assemble_all_matrices()
{
  this->assemble_mass_matrix();
  this->assemble_advection_matrices();
  this->assemble_avg_coupling_matrices();
  this->assemble_jump_coupling_matrix();
  this->assemble_boundary_condition_matrices();
}

void ClawSystem::assemble_mass_matrix()
{
  LOG_SCOPE("assemble_mass_matrix()", "ClawSystem");

  // clear the mass matrix before filling
  _mass_matrix->zero();

  // Convenient reference to the Mesh
  const MeshBase & mesh = this->get_mesh();
  const unsigned int mesh_dim = mesh.mesh_dimension();

  // Construct an FE object to be used for assembly
  FEType fe_type = this->variable_type(0);
  std::unique_ptr<FEBase> fe = FEBase::build(mesh_dim, fe_type);

  // Quadrature rules for numerical integration.
  QGauss qrule (mesh_dim, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule (&qrule);

  // Pre-request JxW and shape function values at interior quadrature points.
  const std::vector<Real> & JxW = fe->get_JxW();
  const std::vector<std::vector<Real>> & phi = fe->get_phi();

  // Element mass matrix
  DenseMatrix<Number> Ke;

  // DofMap reference and vector used for getting dof_indices on each Elem
  const DofMap & dof_map = get_dof_map();
  std::vector<dof_id_type> dof_indices;

  for (const auto & elem : mesh.active_local_element_ptr_range())
  {
    // Recompute shape function values, etc. on current Elem
    fe->reinit(elem);

    // Get DOFs for current Elem
    dof_map.dof_indices (elem, dof_indices);

    // The number of degrees of freedom on this element for variable 0.
    const unsigned int n_q1_dofs = dof_indices.size();

    // Resize the Elem mass matrix for the current Elem
    Ke.resize(n_q1_dofs, n_q1_dofs);

    // The number of quadrature points on the current Elem
    unsigned int n_qpoints = qrule.n_points();

    // Assemble element mass matrix. Here we loop over vars but assume that
    // every var has the same number of DOFs as variable 0.
    for (auto qp : make_range(n_qpoints))
      for (auto i : make_range(n_q1_dofs))
        for (auto j : make_range(n_q1_dofs))
          Ke(i,j) += JxW[qp] * phi[j][qp] * phi[i][qp];

    // Apply constraints, e.g. periodic constraints and "stamp" the
    // local matrix into the global matrix
    this->get_dof_map().constrain_element_matrix(Ke, dof_indices);
    _mass_matrix->add_matrix (Ke, dof_indices);
  }

  // Close the mass matrix now that we are done with assembly
  _mass_matrix->close();
}

void ClawSystem::assemble_advection_matrices()
{
  LOG_SCOPE("assemble_advection_matrices()", "ClawSystem");

  // Zero all advection matrices before filling
  for (auto & mat : _advection_matrices)
    mat->zero();

  // Convenient reference to the Mesh
  const MeshBase & mesh = this->get_mesh();
  const unsigned int mesh_dim = mesh.mesh_dimension();

  // Check that the mesh dimension matches the number of advection
  // matrices we are going to assemble.
  libmesh_error_msg_if(
    mesh_dim != _advection_matrices.size(),
    "Error: Expected number of advection matrices to match mesh dimension.");

  // Construct an FE object to be used for assembly
  FEType fe_type = this->variable_type(0);
  std::unique_ptr<FEBase> fe = FEBase::build(mesh_dim, fe_type);

  // Quadrature rules for numerical integration.
  QGauss qrule (mesh_dim, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule (&qrule);

  // Pre-request JxW and shape function values/derivatives at interior quadrature points
  const std::vector<Real> & JxW = fe->get_JxW();
  const std::vector<std::vector<Real>> & phi = fe->get_phi();
  const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();

  // We are going to build a separate element advection matrix for each direction.
  std::vector<DenseMatrix<Number>> K_vec(mesh_dim);

  // DofMap reference and vector used for getting dof_indices on each Elem
  const DofMap & dof_map = get_dof_map();
  std::vector<dof_id_type> dof_indices;

  for (const auto & elem : mesh.active_local_element_ptr_range())
  {
    // Recompute shape function values, etc. on current Elem
    fe->reinit(elem);

    // Get DOFs for current Elem
    dof_map.dof_indices (elem, dof_indices);

    // The number of degrees of freedom on this element
    const unsigned int n_q1_dofs = dof_indices.size();

    // Resize all the Elem advection matrices for the current Elem
    for (auto & Ke : K_vec)
      Ke.resize(n_q1_dofs, n_q1_dofs);

    // The number of quadrature points on the current Elem
    unsigned int n_qpoints = qrule.n_points();

    // Assemble the convection operators. There are two separate
    // operators in 2D (one for x-direction and one for y-direction).
    // Note: in the FV (CONSTANT, MONOMIAL) case, these interior
    // contributions are all zero because the spatial derivatives of
    // the single constant "basis function" are all zero.
    for (auto qp : make_range(n_qpoints))
      for (auto i : make_range(n_q1_dofs))
        for (auto j : make_range(n_q1_dofs))
          for (auto dim : make_range(mesh_dim))
            K_vec[dim](i,j) += JxW[qp] * phi[j][qp] * dphi[i][qp](dim);

    // Apply constraints to Elem matrices and "stamp" into the global matrices
    for (auto dim : make_range(mesh_dim))
    {
      dof_map.constrain_element_matrix(K_vec[dim], dof_indices);
      _advection_matrices[dim]->add_matrix(K_vec[dim], dof_indices);
    }
  }

  // Close all advection matrices now that we are done assembling
  for (auto & mat : _advection_matrices)
    mat->close();
}

void ClawSystem::assemble_avg_coupling_matrices()
{
  LOG_SCOPE("assemble_avg_coupling_matrices()", "ClawSystem");

  // Zero all average coupling matrices before filling
  for (auto & mat : _avg_matrices)
    mat->zero();

  // Convenient reference to the Mesh
  const MeshBase & mesh = this->get_mesh();
  const unsigned int mesh_dim = mesh.mesh_dimension();

  // Check that the mesh dimension matches the number of average coupling
  // matrices we are going to assemble.
  libmesh_error_msg_if(
    mesh_dim != _avg_matrices.size(),
    "Error: Expected number of average coupling matrices to match mesh dimension.");

  // Construct an side FE objects for both "current" and "neighbor" Elems
  FEType fe_type = this->variable_type(0);
  std::unique_ptr<FEBase> fe_elem_face = FEBase::build(mesh_dim, fe_type);
  std::unique_ptr<FEBase> fe_neighbor_face =FEBase::build(mesh_dim, fe_type);

  // Quadrature rules for numerical integration.
  QGauss qface(mesh_dim-1, fe_type.default_quadrature_order());
  fe_elem_face->attach_quadrature_rule(&qface);
  fe_neighbor_face->attach_quadrature_rule(&qface);

  // Here we define some references to cell-specific data that
  // will be used to assemble the linear system.
  // Data for surface integrals on the element boundary
  const std::vector<std::vector<Real>> & phi_face = fe_elem_face->get_phi();
  const std::vector<Real> & JxW_face = fe_elem_face->get_JxW();
  const std::vector<Point> & qface_normals = fe_elem_face->get_normals();
  const std::vector<Point> & qface_points = fe_elem_face->get_xyz();

  // Data for surface integrals on the neighbor boundary
  const std::vector<std::vector<Real>> & phi_neighbor_face = fe_neighbor_face->get_phi();

  // Data structures to contain the element/neighbor contributions on
  // interior faces. We store these separately for each space dimension.
  std::vector<DenseMatrix<Number>> Kne_vec(mesh_dim);
  std::vector<DenseMatrix<Number>> Ken_vec(mesh_dim);
  std::vector<DenseMatrix<Number>> Kee_vec(mesh_dim);
  std::vector<DenseMatrix<Number>> Knn_vec(mesh_dim);

  // DofMap reference and vectors used for getting dof_indices on each Elem
  const DofMap & dof_map = get_dof_map();
  std::vector<dof_id_type> dof_indices;
  std::vector<dof_id_type> neighbor_dof_indices;

  for (const auto * elem : mesh.active_local_element_ptr_range())
  {
    dof_map.dof_indices (elem, dof_indices);
    const unsigned int n_dofs = dof_indices.size();

    for (auto side : elem->side_index_range())
    {
      // Skip over boundary sides
      if (elem->neighbor_ptr(side) != nullptr)
      {
        const Elem * neighbor = elem->neighbor_ptr(side);

        // Get the global id of the element and the neighbor
        const unsigned int elem_id = elem->id();
        const unsigned int neighbor_id = neighbor->id();

        // Compare element IDs since we don't want to compute the same contribution twice
        if (elem_id < neighbor_id)
        {
          // Pointer to the element side
          std::unique_ptr<const Elem> elem_side = elem->build_side_ptr(side);

          // Reinitialize shape functions on the element side
          fe_elem_face->reinit(elem, side);

          // Find their locations on the neighbor
          std::vector<Point> qface_neighbor_points;
          FEInterface::inverse_map(
            elem->dim(),
            fe_type,
            neighbor,
            qface_points,
            qface_neighbor_points);

          // Calculate the neighbor element shape functions at those locations
          fe_neighbor_face->reinit(neighbor, &qface_neighbor_points);

          // Get the degree of freedom indices for the
          // neighbor.  These define where in the global
          // matrix this neighbor will contribute to.
          dof_map.dof_indices (neighbor, neighbor_dof_indices);
          const unsigned int n_neighbor_dofs = neighbor_dof_indices.size();

          // Zero the element matrices summing into them.  We use the
          // resize member here because the number of degrees of
          // freedom might have changed from the last element or
          // neighbor.  Note that Kne and Ken are not square matrices
          // if neighbor and element have a different p level
          for (auto dim : make_range(mesh_dim))
          {
            Kne_vec[dim].resize(n_neighbor_dofs, n_dofs);
            Ken_vec[dim].resize(n_dofs, n_neighbor_dofs);
            Kee_vec[dim].resize(n_dofs, n_dofs);
            Knn_vec[dim].resize(n_neighbor_dofs, n_neighbor_dofs);
          }

          // Now we will build the "average" contribution due to the
          // Lax-Friedrichs numerical flux term:
          //
          // "{{f(q)}}[v] \cdot n_e"
          //
          // For the advection problem, this is term is given by:
          // 0.5 * ((\vec{u} q)_e + (\vec{u} q)_n) \cdot \hat{n}_e * [v]
          //
          // The assumption below is that the advective velocity is
          // constant on the entire mesh, so it is pulled outside of
          // the integration loop, hence the reason it does not appear
          // explicitly below.
          for (auto qp : make_range(qface.n_points()))
          {
            // Kee Matrix.
            // This is the "q_e" contribution with v = phi^e_i.
            for (auto i : make_range(n_dofs))
              for (auto j : make_range(n_dofs))
                for (auto dim : make_range(mesh_dim))
                  Kee_vec[dim](i,j) += 0.5 * JxW_face[qp] * phi_face[i][qp] * phi_face[j][qp] * qface_normals[qp](dim);

            // Knn Matrix.
            // This is the "q_n" contribution with v = -phi^n_i.
            // Hence, there is a negative sign.
            for (auto i : make_range(n_neighbor_dofs))
              for (auto j : make_range(n_neighbor_dofs))
                for (auto dim : make_range(mesh_dim))
                  Knn_vec[dim](i,j) += -0.5 * JxW_face[qp] * phi_neighbor_face[i][qp] * phi_neighbor_face[j][qp] * qface_normals[qp](dim);

            // Kne Matrix.
            // This is the "q_e" contribution with v = -phi^n_i.
            // i goes over neighbor DOFs
            // j goes over self DOFs
            // Hence, there is a negative sign.
            for (auto i : make_range(n_neighbor_dofs))
              for (auto j : make_range(n_dofs))
                for (auto dim : make_range(mesh_dim))
                  Kne_vec[dim](i,j) += -0.5 * JxW_face[qp] * phi_neighbor_face[i][qp] * phi_face[j][qp] * qface_normals[qp](dim);

            // Ken Matrix.
            // This is the "q_n" contribution with v = phi^e_i.
            // i goes over self DOFs
            // j goes over neighbor DOFs
            for (auto i : make_range(n_dofs))
              for (auto j : make_range(n_neighbor_dofs))
                for (auto dim : make_range(mesh_dim))
                  Ken_vec[dim](i,j) += 0.5 * JxW_face[qp] * phi_face[i][qp] * phi_neighbor_face[j][qp] * qface_normals[qp](dim);
          }

          // The element and neighbor boundary matrix are now built
          // for this side.  Add them to the appropriate global matrices now.
          for (auto dim : make_range(mesh_dim))
          {
            auto & mat = _avg_matrices[dim];
            mat->add_matrix(Kne_vec[dim], neighbor_dof_indices, dof_indices);
            mat->add_matrix(Ken_vec[dim], dof_indices, neighbor_dof_indices);
            mat->add_matrix(Kee_vec[dim], dof_indices);
            mat->add_matrix(Knn_vec[dim], neighbor_dof_indices);
          }
        }
      }
    }
  }

  // Close all average coupling matrices now that we are done assembling
  for (auto & mat : _avg_matrices)
    mat->close();
}

void ClawSystem::assemble_jump_coupling_matrix()
{
  LOG_SCOPE("assemble_jump_coupling_matrix()", "ClawSystem");

  // clear the matrix
  _jump_matrix->zero();

  const MeshBase & mesh = this->get_mesh();
  const unsigned int mesh_dim = mesh.mesh_dimension();

  // Construct face FE objects for both "current" and "neighbor" Elems
  FEType fe_type = variable_type(0);
  std::unique_ptr<FEBase> fe_elem_face = FEBase::build(mesh_dim, fe_type);
  std::unique_ptr<FEBase> fe_neighbor_face = FEBase::build(mesh_dim, fe_type);

  // Quadrature rules for numerical integration.
  QGauss qface(mesh_dim-1, fe_type.default_quadrature_order());
  fe_elem_face->attach_quadrature_rule(&qface);
  fe_neighbor_face->attach_quadrature_rule(&qface);

  // Here we define some references to cell-specific data that
  // will be used to assemble the linear system.
  // Data for surface integrals on the element boundary
  const std::vector<std::vector<Real>> & phi_face = fe_elem_face->get_phi();
  const std::vector<Real> & JxW_face = fe_elem_face->get_JxW();
  const std::vector<Point> & qface_normals = fe_elem_face->get_normals();
  const std::vector<Point> & qface_points = fe_elem_face->get_xyz();

  // Data for surface integrals on the neighbor boundary
  const std::vector<std::vector<Real>> & phi_neighbor_face = fe_neighbor_face->get_phi();

  // Data structures to contain the element and neighbor boundary matrix
  // contribution. This matrices will do the coupling beetwen the dofs of
  // the element and those of his neighbors.
  // Ken: matrix coupling elem and neighbor dofs
  DenseMatrix<Number> Kne;
  DenseMatrix<Number> Ken;
  DenseMatrix<Number> Kee;
  DenseMatrix<Number> Knn;

  // DofMap reference and vectors used for storing dof_indices
  const DofMap & dof_map = get_dof_map();
  std::vector<dof_id_type> dof_indices;
  std::vector<dof_id_type> neighbor_dof_indices;

  for (const auto & elem : mesh.active_local_element_ptr_range())
  {
    dof_map.dof_indices (elem, dof_indices);
    const unsigned int n_dofs = dof_indices.size();

    for (auto side : elem->side_index_range())
    {
      // Skip over boundary sides
      if (elem->neighbor_ptr(side) != nullptr)
      {
        const Elem * neighbor = elem->neighbor_ptr(side);

        // Get the global id of the element and the neighbor
        const unsigned int elem_id = elem->id();
        const unsigned int neighbor_id = neighbor->id();

        // Compare element IDs since we don't want to compute the same contribution twice
        if (elem_id < neighbor_id)
        {
          // Pointer to the element side
          std::unique_ptr<const Elem> elem_side = elem->build_side_ptr(side);

          // Reinitialize shape functions on the element side
          fe_elem_face->reinit(elem, side);

          // Find side QP locations in neighbor side reference space
          std::vector<Point> qface_neighbor_points;
          FEInterface::inverse_map(
            elem->dim(),
            fe_type,
            neighbor,
            qface_points,
            qface_neighbor_points);

          // Calculate the neighbor element shape functions at those locations
          fe_neighbor_face->reinit(neighbor, &qface_neighbor_points);

          // Get the degree of freedom indices for the
          // neighbor.  These define where in the global
          // matrix this neighbor will contribute to.
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
          for (auto qp : make_range(qface.n_points()))
          {
            // Kee Matrix. This is the "q_e" contribution from the [q][v] term with v = phi^e_i
            // (q_e - q_n) * phi^e_i
            for (auto i : make_range(n_dofs))
              for (auto j : make_range(n_dofs))
                 Kee(i,j) += 0.5 * JxW_face[qp] * phi_face[i][qp] * phi_face[j][qp];

            // Knn Matrix. This is the "q_n" contribution from the [q][v] term with v = -phi^n_i
            // (q_e - q_n) * -phi^n_i
            for (auto i : make_range(n_neighbor_dofs))
              for (auto j : make_range(n_neighbor_dofs))
                 Knn(i,j) += 0.5 * JxW_face[qp] * phi_neighbor_face[i][qp] * phi_neighbor_face[j][qp];

            // Kne Matrix. This is the "q_e" contribution from the [q][v] term with v = -phi^n_i,
            // so there is a negative sign.
            // (q_e - q_n) * -phi^n_i
            for (auto i : make_range(n_neighbor_dofs))
              for (auto j : make_range(n_dofs))
                 Kne(i,j) -= 0.5 * JxW_face[qp] * phi_neighbor_face[i][qp] * phi_face[j][qp];

            // Ken Matrix. This is the "q_n" contribution from the [q][v] term with v = phi^e_i,
            // so there is a negative sign.
            // (q_e - q_n) * phi^e_i
            for (auto i : make_range(n_dofs))
              for (auto j : make_range(n_neighbor_dofs))
                 Ken(i,j) -= 0.5 * JxW_face[qp] * phi_face[i][qp] * phi_neighbor_face[j][qp];
          }

          // The element and neighbor boundary matrix are now built
          // for this side.  Add them to the global matrix
          // The \p SparseMatrix::add_matrix() members do this for us.
          _jump_matrix->add_matrix(Kne, neighbor_dof_indices, dof_indices);
          _jump_matrix->add_matrix(Ken, dof_indices, neighbor_dof_indices);
          _jump_matrix->add_matrix(Kee, dof_indices);
          _jump_matrix->add_matrix(Knn, neighbor_dof_indices);
        }
      }
    }
  }

  // Close the matrix now that we are done assembling it
  _jump_matrix->close();
}

void ClawSystem::assemble_boundary_condition_matrices()
{
  LOG_SCOPE("assemble_boundary_condition_matrices()", "ClawSystem");

  // Zero all boundary condition matrices before filling
  for (auto & mat : _boundary_condition_matrices)
    mat->zero();

  // Construct FE objects for the cell and side domains
  const MeshBase & mesh = this->get_mesh();
  const unsigned int mesh_dim = mesh.mesh_dimension();

  // Check that the mesh dimension matches the number of boundary
  // condition matrices we are going to assemble.
  libmesh_error_msg_if(
    mesh_dim != _boundary_condition_matrices.size(),
    "Error: Expected number of boundary condition matrices to match mesh dimension.");

  // Construct side FE object to be used for assembly
  FEType fe_type = variable_type(0);
  std::unique_ptr<FEBase> fe_elem_face = FEBase::build(mesh_dim, fe_type);
  QGauss qface(mesh_dim-1, fe_type.default_quadrature_order());

  // Tell the face finite element object to use our quadrature rule.
  fe_elem_face->attach_quadrature_rule(&qface);

  // Here we define some references to cell-specific data that
  // will be used to assemble the linear system.
  // Data for surface integrals on the element boundary
  const std::vector<std::vector<Real>> & phi_face = fe_elem_face->get_phi();
  const std::vector<Real> & JxW_face = fe_elem_face->get_JxW();
  const std::vector<Point> & qface_normals = fe_elem_face->get_normals();

  // We are going to build a separate element boundary condition matrix for each direction.
  std::vector<DenseMatrix<Number>> K_vec(mesh_dim);

  // DofMap reference and vector used for getting dof_indices on each Elem
  const DofMap & dof_map = get_dof_map();
  std::vector<dof_id_type> dof_indices;

  for (const auto & elem : mesh.active_local_element_ptr_range())
  {
    dof_map.dof_indices (elem, dof_indices);
    const unsigned int n_dofs = dof_indices.size();

    for (auto side : elem->side_index_range())
    {
      // skip interior sides
      if (elem->neighbor_ptr(side) == nullptr)
      {
        // Reinitialize shape functions on the element side
        fe_elem_face->reinit(elem, side);

        // Resize all the Elem boundary condition matrices for the current Elem
        for (auto & Ke : K_vec)
          Ke.resize(n_dofs, n_dofs);

        // Now we will build the element and neighbor
        // boundary matrices.  This involves
        // a double loop to integrate the test funcions
        // (i) against the trial functions (j).
        for (auto qp : make_range(qface.n_points()))
          for (auto i : make_range(n_dofs))
            for (auto j : make_range(n_dofs))
              for (auto dim : make_range(mesh_dim))
                K_vec[dim](i,j) += JxW_face[qp] * phi_face[i][qp] * phi_face[j][qp] * qface_normals[qp](dim);

        // "Stamp" Elem boundary condition matrices into the global boundary condition matrices
        for (auto dim : make_range(mesh_dim))
          _boundary_condition_matrices[dim]->add_matrix(K_vec[dim], dof_indices);
      }
    }
  }

  // Close all boundary condition matrices now that we are done assembling
  for (auto & mat : _boundary_condition_matrices)
    mat->close();
}

SparseMatrix<Number> &
ClawSystem::get_mass_matrix()
{
  return *_mass_matrix;
}

SparseMatrix<Number> &
ClawSystem::get_advection_matrix(unsigned int dim)
{
  libmesh_error_msg_if(
    dim >= _advection_matrices.size(),
    "Error: Invalid dimension " << dim << " requested.");

  return *_advection_matrices[dim];
}

SparseMatrix<Number> &
ClawSystem::get_avg_matrix(unsigned int dim)
{
  libmesh_error_msg_if(
    dim >= _avg_matrices.size(),
    "Error: Invalid dimension " << dim << " requested.");

  return *_avg_matrices[dim];
}

SparseMatrix<Number> &
ClawSystem::get_jump_matrix()
{
  return *_jump_matrix;
}

SparseMatrix<Number> &
ClawSystem::get_boundary_condition_matrix(unsigned int dim)
{
  libmesh_error_msg_if(
    dim >= _boundary_condition_matrices.size(),
    "Error: Invalid dimension " << dim << " requested.");

  return *_boundary_condition_matrices[dim];
}

void
ClawSystem::set_delta_t(Real delta_t_in)
{
  _delta_t = delta_t_in;
}

Real
ClawSystem::get_delta_t()
{
  return _delta_t;
}

void
ClawSystem::set_n_time_steps(unsigned int n_time_steps_in)
{
  _n_time_steps = n_time_steps_in;
}

unsigned int
ClawSystem::get_n_time_steps()
{
  return _n_time_steps;
}

void
ClawSystem::set_LxF_constant(Real LxF_constant_in)
{
  _LxF_constant = LxF_constant_in;
}

Real
ClawSystem::get_LxF_constant()
{
  return _LxF_constant;
}

void
ClawSystem::set_temporal_discretization_type(TemporalDiscretizationType td_in)
{
  _temporal_discretization_type = td_in;
}

ClawSystem::TemporalDiscretizationType
ClawSystem::get_temporal_discretization_type()
{
  return _temporal_discretization_type;
}

void
ClawSystem::set_write_interval(unsigned int write_interval_in)
{
  _write_interval = write_interval_in;
}

unsigned int
ClawSystem::get_write_interval()
{
  return _write_interval;
}

void
ClawSystem::set_time(Real time_in)
{
  _time = time_in;
}

Real
ClawSystem::get_time()
{
  return _time;
}

void ClawSystem::write_out_discretization_matrices()
{
  _mass_matrix->print_matlab("mass_matrix.m");
  _jump_matrix->print_matlab("jump_matrix.m");

  for (auto i : index_range(_avg_matrices))
    _avg_matrices[i]->print_matlab("avg_matrix_" + std::to_string(i+1) + ".m");

  for (auto i : index_range(_advection_matrices))
    _advection_matrices[i]->print_matlab("advection_matrix_" + std::to_string(i+1) + ".m");

  for (auto i : index_range(_boundary_condition_matrices))
    _boundary_condition_matrices[i]->print_matlab("boundary_condition_matrix_" + std::to_string(i+1) + ".m");
}

void ClawSystem::init_data ()
{
  Parent::init_data();

  DofMap & dof_map = this->get_dof_map();

  // Helper lambda function that attaches the input matrix to the
  // DofMap, initializes it, and finally zeros it.
  auto prepare_matrix = [&](SparseMatrix<Number> & mat)
  {
    dof_map.attach_matrix(mat);
    mat.init();
    mat.zero();
  };

  prepare_matrix(*_mass_matrix);

  for (auto & mat : _advection_matrices)
    prepare_matrix(*mat);

  for (auto & mat : _avg_matrices)
    prepare_matrix(*mat);

  prepare_matrix(*_jump_matrix);

  for (auto & mat : _boundary_condition_matrices)
    prepare_matrix(*mat);
}

} // namespace libMesh
