// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


// <h1> Systems Example 7 - Large deformation elasticity (St. Venant-Kirchoff material) </h1>
// \author Lorenzo Zanon
// \author David Knezevic
// \date 2014
//
// In this example, we consider an elastic cantilever beam modeled as a St. Venant-Kirchoff
// material (which is an extension of the linear elastic material model to the nonlinear
// regime). The implementation presented here uses NonlinearImplicitSystem.
//
// We formulate the PDE on the reference geometry (\Omega) as opposed to the deformed
// geometry (\Omega^deformed). As a result (e.g. see Ciarlet's 3D elasticity book,
// Theorem 2.6-2) the PDE is given as follows:
//
//     \int_\Omega F_im Sigma_mj v_i,j = \int_\Omega f_i v_i + \int_\Gamma g_i v_i ds
//
// where:
//  * F is the deformation gradient, F = I + du/dx (x here refers to reference coordinates).
//  * Sigma is the second Piola-Kirchoff stress, which for the St. Venant Kirchoff model is
//    given by Sigma_ij = C_ijkl E_kl, where E_kl is the strain,
//    E_kl = 0.5 * (u_k,l + u_l,k + u_m,k u_m,l).
//  * f is a body load.
//  * g is a surface traction on the surface \Gamma.
//
// In this example we only consider a body load (e.g. gravity), hence we set g = 0.

// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <cmath>

// Various include files needed for the mesh & solver functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/elem.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/getpot.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/zero_function.h"
#include "libmesh/enum_solver_package.h"

// The nonlinear solver and system we will be using
#include "libmesh/nonlinear_solver.h"
#include "libmesh/nonlinear_implicit_system.h"

#define BOUNDARY_ID_MIN_Z 0
#define BOUNDARY_ID_MIN_Y 1
#define BOUNDARY_ID_MAX_X 2
#define BOUNDARY_ID_MAX_Y 3
#define BOUNDARY_ID_MIN_X 4
#define BOUNDARY_ID_MAX_Z 5

using namespace libMesh;



class LargeDeformationElasticity : public NonlinearImplicitSystem::ComputeResidual,
                                   public NonlinearImplicitSystem::ComputeJacobian
{
private:
  EquationSystems & es;

public:

  LargeDeformationElasticity (EquationSystems & es_in) :
    es(es_in)
  {}

  /**
   * Kronecker delta function.
   */
  Real kronecker_delta(unsigned int i,
                       unsigned int j)
  {
    return i == j ? 1. : 0.;
  }

  /**
   * Evaluate the fourth order tensor (C_ijkl) that relates stress to strain.
   */
  Real elasticity_tensor(Real young_modulus,
                         Real poisson_ratio,
                         unsigned int i,
                         unsigned int j,
                         unsigned int k,
                         unsigned int l)
  {
    // Define the Lame constants
    const Real lambda_1 = (young_modulus*poisson_ratio)/((1.+poisson_ratio)*(1.-2.*poisson_ratio));
    const Real lambda_2 = young_modulus/(2.*(1.+poisson_ratio));

    return lambda_1 * kronecker_delta(i,j) * kronecker_delta(k,l) +
      lambda_2 * (kronecker_delta(i,k) * kronecker_delta(j,l) + kronecker_delta(i,l) * kronecker_delta(j,k));
  }


  /**
   * Evaluate the Jacobian of the nonlinear system.
   */
  virtual void jacobian (const NumericVector<Number> & soln,
                         SparseMatrix<Number> & jacobian,
                         NonlinearImplicitSystem & /*sys*/)
  {
    const Real young_modulus = es.parameters.get<Real>("young_modulus");
    const Real poisson_ratio = es.parameters.get<Real>("poisson_ratio");

    const MeshBase & mesh = es.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    NonlinearImplicitSystem & system =
      es.get_system<NonlinearImplicitSystem>("NonlinearElasticity");

    const unsigned int u_var = system.variable_number ("u");

    const DofMap & dof_map = system.get_dof_map();

    FEType fe_type = dof_map.variable_type(u_var);
    std::unique_ptr<FEBase> fe (FEBase::build(dim, fe_type));
    QGauss qrule (dim, fe_type.default_quadrature_order());
    fe->attach_quadrature_rule (&qrule);

    std::unique_ptr<FEBase> fe_face (FEBase::build(dim, fe_type));
    QGauss qface (dim-1, fe_type.default_quadrature_order());
    fe_face->attach_quadrature_rule (&qface);

    const std::vector<Real> & JxW = fe->get_JxW();
    const std::vector<std::vector<Real>> & phi = fe->get_phi();
    const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();

    DenseMatrix<Number> Ke;
    DenseSubMatrix<Number> Ke_var[3][3] =
      {
        {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)},
        {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)},
        {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)}
      };

    std::vector<dof_id_type> dof_indices;
    std::vector<std::vector<dof_id_type>> dof_indices_var(3);

    jacobian.zero();

    for (const auto & elem : mesh.active_local_element_ptr_range())
      {
        dof_map.dof_indices (elem, dof_indices);
        for (unsigned int var=0; var<3; var++)
          dof_map.dof_indices (elem, dof_indices_var[var], var);

        const unsigned int n_dofs = dof_indices.size();
        const unsigned int n_var_dofs = dof_indices_var[0].size();

        fe->reinit (elem);

        Ke.resize (n_dofs, n_dofs);
        for (unsigned int var_i=0; var_i<3; var_i++)
          for (unsigned int var_j=0; var_j<3; var_j++)
            Ke_var[var_i][var_j].reposition (var_i*n_var_dofs, var_j*n_var_dofs, n_var_dofs, n_var_dofs);

        for (unsigned int qp=0; qp<qrule.n_points(); qp++)
          {
            DenseVector<Number> u_vec(3);
            DenseMatrix<Number> grad_u(3, 3);
            for (unsigned int var_i=0; var_i<3; var_i++)
              {
                for (unsigned int j=0; j<n_var_dofs; j++)
                  u_vec(var_i) += phi[j][qp]*soln(dof_indices_var[var_i][j]);

                // Row is variable u1, u2, or u3, column is x, y, or z
                for (unsigned int var_j=0; var_j<3; var_j++)
                  for (unsigned int j=0; j<n_var_dofs; j++)
                    grad_u(var_i,var_j) += dphi[j][qp](var_j)*soln(dof_indices_var[var_i][j]);
              }

            DenseMatrix<Number> strain_tensor(3, 3);
            for (unsigned int i=0; i<3; i++)
              for (unsigned int j=0; j<3; j++)
                {
                  strain_tensor(i,j) += 0.5 * (grad_u(i,j) + grad_u(j,i));

                  for (unsigned int k=0; k<3; k++)
                    strain_tensor(i,j) += 0.5 * grad_u(k,i)*grad_u(k,j);
                }

            // Define the deformation gradient
            DenseMatrix<Number> F(3, 3);
            F = grad_u;
            for (unsigned int var=0; var<3; var++)
              F(var, var) += 1.;

            DenseMatrix<Number> stress_tensor(3, 3);

            for (unsigned int i=0; i<3; i++)
              for (unsigned int j=0; j<3; j++)
                for (unsigned int k=0; k<3; k++)
                  for (unsigned int l=0; l<3; l++)
                    stress_tensor(i,j) +=
                      elasticity_tensor(young_modulus, poisson_ratio, i, j, k, l) * strain_tensor(k, l);

            for (unsigned int dof_i=0; dof_i<n_var_dofs; dof_i++)
              for (unsigned int dof_j=0; dof_j<n_var_dofs; dof_j++)
                {
                  for (unsigned int i=0; i<3; i++)
                    for (unsigned int j=0; j<3; j++)
                      for (unsigned int m=0; m<3; m++)
                        Ke_var[i][i](dof_i,dof_j) += JxW[qp] *
                          (-dphi[dof_j][qp](m) * stress_tensor(m,j) * dphi[dof_i][qp](j));

                  for (unsigned int i=0; i<3; i++)
                    for (unsigned int j=0; j<3; j++)
                      for (unsigned int k=0; k<3; k++)
                        for (unsigned int l=0; l<3; l++)
                          {
                            Number FxC_ijkl = 0.;
                            for (unsigned int m=0; m<3; m++)
                              FxC_ijkl += F(i,m) * elasticity_tensor(young_modulus, poisson_ratio, m, j, k, l);

                            Ke_var[i][k](dof_i,dof_j) += JxW[qp] *
                              (-0.5 * FxC_ijkl * dphi[dof_j][qp](l) * dphi[dof_i][qp](j));

                            Ke_var[i][l](dof_i,dof_j) += JxW[qp] *
                              (-0.5 * FxC_ijkl * dphi[dof_j][qp](k) * dphi[dof_i][qp](j));

                            for (unsigned int n=0; n<3; n++)
                              Ke_var[i][n](dof_i,dof_j) += JxW[qp] *
                                (-0.5 * FxC_ijkl * (dphi[dof_j][qp](k) * grad_u(n,l) + dphi[dof_j][qp](l) * grad_u(n,k)) * dphi[dof_i][qp](j));
                          }
                }
          }

        dof_map.constrain_element_matrix (Ke, dof_indices);
        jacobian.add_matrix (Ke, dof_indices);
      }
  }

  /**
   * Evaluate the residual of the nonlinear system.
   */
  virtual void residual (const NumericVector<Number> & soln,
                         NumericVector<Number> & residual,
                         NonlinearImplicitSystem & /*sys*/)
  {
    const Real young_modulus = es.parameters.get<Real>("young_modulus");
    const Real poisson_ratio = es.parameters.get<Real>("poisson_ratio");
    const Real forcing_magnitude = es.parameters.get<Real>("forcing_magnitude");

    const MeshBase & mesh = es.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    NonlinearImplicitSystem & system =
      es.get_system<NonlinearImplicitSystem>("NonlinearElasticity");

    const unsigned int u_var = system.variable_number ("u");

    const DofMap & dof_map = system.get_dof_map();

    FEType fe_type = dof_map.variable_type(u_var);
    std::unique_ptr<FEBase> fe (FEBase::build(dim, fe_type));
    QGauss qrule (dim, fe_type.default_quadrature_order());
    fe->attach_quadrature_rule (&qrule);

    std::unique_ptr<FEBase> fe_face (FEBase::build(dim, fe_type));
    QGauss qface (dim-1, fe_type.default_quadrature_order());
    fe_face->attach_quadrature_rule (&qface);

    const std::vector<Real> & JxW = fe->get_JxW();
    const std::vector<std::vector<Real>> & phi = fe->get_phi();
    const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();

    DenseVector<Number> Re;

    DenseSubVector<Number> Re_var[3] =
      {DenseSubVector<Number>(Re),
       DenseSubVector<Number>(Re),
       DenseSubVector<Number>(Re)};

    std::vector<dof_id_type> dof_indices;
    std::vector<std::vector<dof_id_type>> dof_indices_var(3);

    residual.zero();

    for (const auto & elem : mesh.active_local_element_ptr_range())
      {
        dof_map.dof_indices (elem, dof_indices);
        for (unsigned int var=0; var<3; var++)
          dof_map.dof_indices (elem, dof_indices_var[var], var);

        const unsigned int n_dofs = dof_indices.size();
        const unsigned int n_var_dofs = dof_indices_var[0].size();

        fe->reinit (elem);

        Re.resize (n_dofs);
        for (unsigned int var=0; var<3; var++)
          Re_var[var].reposition (var*n_var_dofs, n_var_dofs);

        for (unsigned int qp=0; qp<qrule.n_points(); qp++)
          {
            DenseVector<Number> u_vec(3);
            DenseMatrix<Number> grad_u(3, 3);
            for (unsigned int var_i=0; var_i<3; var_i++)
              {
                for (unsigned int j=0; j<n_var_dofs; j++)
                  u_vec(var_i) += phi[j][qp]*soln(dof_indices_var[var_i][j]);

                // Row is variable u, v, or w column is x, y, or z
                for (unsigned int var_j=0; var_j<3; var_j++)
                  for (unsigned int j=0; j<n_var_dofs; j++)
                    grad_u(var_i,var_j) += dphi[j][qp](var_j)*soln(dof_indices_var[var_i][j]);
              }

            DenseMatrix<Number> strain_tensor(3, 3);
            for (unsigned int i=0; i<3; i++)
              for (unsigned int j=0; j<3; j++)
                {
                  strain_tensor(i,j) += 0.5 * (grad_u(i,j) + grad_u(j,i));

                  for (unsigned int k=0; k<3; k++)
                    strain_tensor(i,j) += 0.5 * grad_u(k,i)*grad_u(k,j);
                }

            // Define the deformation gradient
            DenseMatrix<Number> F(3, 3);
            F = grad_u;
            for (unsigned int var=0; var<3; var++)
              F(var, var) += 1.;

            DenseMatrix<Number> stress_tensor(3, 3);

            for (unsigned int i=0; i<3; i++)
              for (unsigned int j=0; j<3; j++)
                for (unsigned int k=0; k<3; k++)
                  for (unsigned int l=0; l<3; l++)
                    stress_tensor(i,j) +=
                      elasticity_tensor(young_modulus, poisson_ratio, i, j, k, l) * strain_tensor(k,l);

            DenseVector<Number> f_vec(3);
            f_vec(0) = 0.;
            f_vec(1) = 0.;
            f_vec(2) = -forcing_magnitude;

            for (unsigned int dof_i=0; dof_i<n_var_dofs; dof_i++)
              for (unsigned int i=0; i<3; i++)
                {
                  for (unsigned int j=0; j<3; j++)
                    {
                      Number FxStress_ij = 0.;
                      for (unsigned int m=0; m<3; m++)
                        FxStress_ij += F(i,m) * stress_tensor(m,j);

                      Re_var[i](dof_i) += JxW[qp] * (-FxStress_ij * dphi[dof_i][qp](j));
                    }

                  Re_var[i](dof_i) += JxW[qp] * (f_vec(i) * phi[dof_i][qp]);
                }
          }

        dof_map.constrain_element_vector (Re, dof_indices);
        residual.add_vector (Re, dof_indices);
      }
  }

  /**
   * Compute the Cauchy stress for the current solution.
   */
  void compute_stresses()
  {
    const Real young_modulus = es.parameters.get<Real>("young_modulus");
    const Real poisson_ratio = es.parameters.get<Real>("poisson_ratio");

    const MeshBase & mesh = es.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    NonlinearImplicitSystem & system =
      es.get_system<NonlinearImplicitSystem>("NonlinearElasticity");

    unsigned int displacement_vars[3];
    displacement_vars[0] = system.variable_number ("u");
    displacement_vars[1] = system.variable_number ("v");
    displacement_vars[2] = system.variable_number ("w");
    const unsigned int u_var = system.variable_number ("u");

    const DofMap & dof_map = system.get_dof_map();
    FEType fe_type = dof_map.variable_type(u_var);
    std::unique_ptr<FEBase> fe (FEBase::build(dim, fe_type));
    QGauss qrule (dim, fe_type.default_quadrature_order());
    fe->attach_quadrature_rule (&qrule);

    const std::vector<Real> & JxW = fe->get_JxW();
    const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();

    // Also, get a reference to the ExplicitSystem
    ExplicitSystem & stress_system = es.get_system<ExplicitSystem>("StressSystem");
    const DofMap & stress_dof_map = stress_system.get_dof_map();
    unsigned int sigma_vars[6];
    sigma_vars[0] = stress_system.variable_number ("sigma_00");
    sigma_vars[1] = stress_system.variable_number ("sigma_01");
    sigma_vars[2] = stress_system.variable_number ("sigma_02");
    sigma_vars[3] = stress_system.variable_number ("sigma_11");
    sigma_vars[4] = stress_system.variable_number ("sigma_12");
    sigma_vars[5] = stress_system.variable_number ("sigma_22");

    // Storage for the stress dof indices on each element
    std::vector<std::vector<dof_id_type>> dof_indices_var(system.n_vars());
    std::vector<dof_id_type> stress_dof_indices_var;

    // To store the stress tensor on each element
    DenseMatrix<Number> elem_avg_stress_tensor(3, 3);

    for (const auto & elem : mesh.active_local_element_ptr_range())
      {
        for (unsigned int var=0; var<3; var++)
          dof_map.dof_indices (elem, dof_indices_var[var], displacement_vars[var]);

        const unsigned int n_var_dofs = dof_indices_var[0].size();

        fe->reinit (elem);

        // clear the stress tensor
        elem_avg_stress_tensor.resize(3, 3);

        for (unsigned int qp=0; qp<qrule.n_points(); qp++)
          {
            DenseMatrix<Number> grad_u(3, 3);
            // Row is variable u1, u2, or u3, column is x, y, or z
            for (unsigned int var_i=0; var_i<3; var_i++)
              for (unsigned int var_j=0; var_j<3; var_j++)
                for (unsigned int j=0; j<n_var_dofs; j++)
                  grad_u(var_i,var_j) += dphi[j][qp](var_j) * system.current_solution(dof_indices_var[var_i][j]);

            DenseMatrix<Number> strain_tensor(3, 3);
            for (unsigned int i=0; i<3; i++)
              for (unsigned int j=0; j<3; j++)
                {
                  strain_tensor(i,j) += 0.5 * (grad_u(i,j) + grad_u(j,i));

                  for (unsigned int k=0; k<3; k++)
                    strain_tensor(i,j) += 0.5 * grad_u(k,i)*grad_u(k,j);
                }

            // Define the deformation gradient
            DenseMatrix<Number> F(3, 3);
            F = grad_u;
            for (unsigned int var=0; var<3; var++)
              F(var, var) += 1.;

            DenseMatrix<Number> stress_tensor(3, 3);
            for (unsigned int i=0; i<3; i++)
              for (unsigned int j=0; j<3; j++)
                for (unsigned int k=0; k<3; k++)
                  for (unsigned int l=0; l<3; l++)
                    stress_tensor(i,j) +=
                      elasticity_tensor(young_modulus, poisson_ratio, i, j, k, l) * strain_tensor(k, l);

            // stress_tensor now holds the second Piola-Kirchoff stress (PK2) at point qp.
            // However, in this example we want to compute the Cauchy stress which is given by
            // 1/det(F) * F * PK2 * F^T, hence we now apply this transformation.
            stress_tensor.scale(1./F.det());
            stress_tensor.left_multiply(F);
            stress_tensor.right_multiply_transpose(F);

            // We want to plot the average Cauchy stress on each element, hence
            // we integrate stress_tensor
            elem_avg_stress_tensor.add(JxW[qp], stress_tensor);
          }

        // Get the average stress per element by dividing by volume
        elem_avg_stress_tensor.scale(1./elem->volume());

        // load elem_sigma data into stress_system
        unsigned int stress_var_index = 0;
        for (unsigned int i=0; i<3; i++)
          for (unsigned int j=i; j<3; j++)
            {
              stress_dof_map.dof_indices (elem, stress_dof_indices_var, sigma_vars[stress_var_index]);

              // We are using CONSTANT MONOMIAL basis functions, hence we only need to get
              // one dof index per variable
              dof_id_type dof_index = stress_dof_indices_var[0];

              if ((stress_system.solution->first_local_index() <= dof_index) &&
                  (dof_index < stress_system.solution->last_local_index()))
                stress_system.solution->set(dof_index, elem_avg_stress_tensor(i,j));

              stress_var_index++;
            }
      }

    // Should call close and update when we set vector entries directly
    stress_system.solution->close();
    stress_system.update();
  }

};


int main (int argc, char ** argv)
{
  LibMeshInit init (argc, argv);

  // This example requires the PETSc nonlinear solvers
  libmesh_example_requires(libMesh::default_solver_package() == PETSC_SOLVERS, "--enable-petsc");

  // We use a 3D domain.
  libmesh_example_requires(LIBMESH_DIM > 2, "--disable-1D-only --disable-2D-only");

  GetPot infile("systems_of_equations_ex7.in");
  const Real x_length = infile("x_length", 0.);
  const Real y_length = infile("y_length", 0.);
  const Real z_length = infile("z_length", 0.);
  const Real n_elem_x = infile("n_elem_x", 0);
  const Real n_elem_y = infile("n_elem_y", 0);
  const Real n_elem_z = infile("n_elem_z", 0);
  const std::string approx_order = infile("approx_order", "FIRST");
  const std::string fe_family = infile("fe_family", "LAGRANGE");

  const Real young_modulus = infile("Young_modulus", 1.0);
  const Real poisson_ratio = infile("poisson_ratio", 0.3);
  const Real forcing_magnitude = infile("forcing_magnitude", 0.001);

  const Real nonlinear_abs_tol = infile("nonlinear_abs_tol", 1.e-8);
  const Real nonlinear_rel_tol = infile("nonlinear_rel_tol", 1.e-8);
  const unsigned int nonlinear_max_its = infile("nonlinear_max_its", 50);

  const unsigned int n_solves = infile("n_solves", 10);
  const Real force_scaling = infile("force_scaling", 5.0);

  Mesh mesh(init.comm());

  MeshTools::Generation::build_cube(mesh,
                                    n_elem_x,
                                    n_elem_y,
                                    n_elem_z,
                                    0., x_length,
                                    0., y_length,
                                    0., z_length,
                                    HEX27);

  mesh.print_info();

  EquationSystems equation_systems (mesh);
  LargeDeformationElasticity lde(equation_systems);

  NonlinearImplicitSystem & system =
    equation_systems.add_system<NonlinearImplicitSystem> ("NonlinearElasticity");

  unsigned int u_var =
    system.add_variable("u",
                        Utility::string_to_enum<Order>   (approx_order),
                        Utility::string_to_enum<FEFamily>(fe_family));

  unsigned int v_var =
    system.add_variable("v",
                        Utility::string_to_enum<Order>   (approx_order),
                        Utility::string_to_enum<FEFamily>(fe_family));

  unsigned int w_var =
    system.add_variable("w",
                        Utility::string_to_enum<Order>   (approx_order),
                        Utility::string_to_enum<FEFamily>(fe_family));

  // Also, initialize an ExplicitSystem to store stresses
  ExplicitSystem & stress_system =
    equation_systems.add_system<ExplicitSystem> ("StressSystem");
  stress_system.add_variable("sigma_00", CONSTANT, MONOMIAL);
  stress_system.add_variable("sigma_01", CONSTANT, MONOMIAL);
  stress_system.add_variable("sigma_02", CONSTANT, MONOMIAL);
  stress_system.add_variable("sigma_11", CONSTANT, MONOMIAL);
  stress_system.add_variable("sigma_12", CONSTANT, MONOMIAL);
  stress_system.add_variable("sigma_22", CONSTANT, MONOMIAL);

  equation_systems.parameters.set<Real>         ("nonlinear solver absolute residual tolerance") = nonlinear_abs_tol;
  equation_systems.parameters.set<Real>         ("nonlinear solver relative residual tolerance") = nonlinear_rel_tol;
  equation_systems.parameters.set<unsigned int> ("nonlinear solver maximum iterations")          = nonlinear_max_its;

  system.nonlinear_solver->residual_object = &lde;
  system.nonlinear_solver->jacobian_object = &lde;

  equation_systems.parameters.set<Real>("young_modulus") = young_modulus;
  equation_systems.parameters.set<Real>("poisson_ratio") = poisson_ratio;
  equation_systems.parameters.set<Real>("forcing_magnitude") = forcing_magnitude;

  // Attach Dirichlet boundary conditions
  std::set<boundary_id_type> clamped_boundaries;
  clamped_boundaries.insert(BOUNDARY_ID_MIN_X);

  std::vector<unsigned int> uvw;
  uvw.push_back(u_var);
  uvw.push_back(v_var);
  uvw.push_back(w_var);

  ZeroFunction<Number> zero;

  // Most DirichletBoundary users will want to supply a "locally
  // indexed" functor
  system.get_dof_map().add_dirichlet_boundary
    (DirichletBoundary (clamped_boundaries, uvw, zero,
                        LOCAL_VARIABLE_ORDER));

  equation_systems.init();
  equation_systems.print_info();

  // Provide a loop here so that we can do a sequence of solves
  // where solve n gives a good starting guess for solve n+1.
  // This "continuation" approach is helpful for solving for
  // large values of "forcing_magnitude".
  // Set n_solves and force_scaling in nonlinear_elasticity.in.
  for (unsigned int count=0; count<n_solves; count++)
    {
      Real previous_forcing_magnitude = equation_systems.parameters.get<Real>("forcing_magnitude");
      equation_systems.parameters.set<Real>("forcing_magnitude") = previous_forcing_magnitude*force_scaling;

      libMesh::out << "Performing solve "
                   << count
                   << ", forcing_magnitude: "
                   << equation_systems.parameters.get<Real>("forcing_magnitude")
                   << std::endl;

      system.solve();

      libMesh::out << "System solved at nonlinear iteration "
                   << system.n_nonlinear_iterations()
                   << " , final nonlinear residual norm: "
                   << system.final_nonlinear_residual()
                   << std::endl
                   << std::endl;

      libMesh::out << "Computing stresses..." << std::endl;

      lde.compute_stresses();

#ifdef LIBMESH_HAVE_EXODUS_API
      std::stringstream filename;
      filename << "solution_" << count << ".exo";
      ExodusII_IO (mesh).write_equation_systems(filename.str(), equation_systems);
#endif
    }

  return 0;
}
