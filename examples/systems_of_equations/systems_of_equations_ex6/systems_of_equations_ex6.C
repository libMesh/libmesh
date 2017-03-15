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



// <h1> Systems Example 6 - 3D Linear Elastic Cantilever </h1>
// \author David Knezevic
// \date 2012
//
// This is a 3D version of systems_of_equations_ex4. The weak form PDE for
// equilibrium elasticity is:
//
//     \int_\Omega Sigma_ij v_i,j = \int_\Omega f_i v_i + \int_\Gamma g_i v_i ds,
//
// for all admissible test functions v, where:
//  * Sigma is the stress tensor, which for linear elasticity is
//    given by Sigma_ij = C_ijkl u_k,l.
//  * f is a body load.
//  * g is a surface traction on the surface \Gamma.


// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>

// libMesh includes
#include "libmesh/libmesh_config.h"
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/gnuplot_io.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dense_subvector.h"
#include "libmesh/perf_log.h"
#include "libmesh/elem.h"
#include "libmesh/boundary_info.h"
#include "libmesh/zero_function.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/getpot.h"
#include "libmesh/solver_configuration.h"
#include "libmesh/petsc_linear_solver.h"
#include "libmesh/petsc_macro.h"

#define x_scaling 1.3

// boundary IDs
#define BOUNDARY_ID_MIN_Z 0
#define BOUNDARY_ID_MIN_Y 1
#define BOUNDARY_ID_MAX_X 2
#define BOUNDARY_ID_MAX_Y 3
#define BOUNDARY_ID_MIN_X 4
#define BOUNDARY_ID_MAX_Z 5
#define NODE_BOUNDARY_ID 10
#define EDGE_BOUNDARY_ID 20

// Bring in everything from the libMesh namespace
using namespace libMesh;

#ifdef LIBMESH_HAVE_PETSC
// This class allows us to set the solver and preconditioner
// to be appropriate for linear elasticity.
class PetscSolverConfiguration : public SolverConfiguration
{
public:

  PetscSolverConfiguration(PetscLinearSolver<Number> & petsc_linear_solver) :
    _petsc_linear_solver(petsc_linear_solver)
  {
  }

  virtual void configure_solver()
  {
    PetscErrorCode ierr = 0;
    ierr = KSPSetType (_petsc_linear_solver.ksp(), const_cast<KSPType>(KSPCG));
    libmesh_assert(ierr == 0);

    ierr = PCSetType (_petsc_linear_solver.pc(), const_cast<PCType>(PCBJACOBI));
    libmesh_assert(ierr == 0);
  }

  // The linear solver object that we are configuring
  PetscLinearSolver<Number> & _petsc_linear_solver;

};
#endif

class LinearElasticity : public System::Assembly
{
private:
  EquationSystems & es;

public:

  LinearElasticity (EquationSystems & es_in) :
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
  Real elasticity_tensor(unsigned int i,
                         unsigned int j,
                         unsigned int k,
                         unsigned int l)
  {
    // Hard code material parameters for the sake of simplicity
    const Real poisson_ratio = 0.3;
    const Real young_modulus = 1.;

    // Define the Lame constants
    const Real lambda_1 = (young_modulus*poisson_ratio)/((1.+poisson_ratio)*(1.-2.*poisson_ratio));
    const Real lambda_2 = young_modulus/(2.*(1.+poisson_ratio));

    return lambda_1 * kronecker_delta(i, j) * kronecker_delta(k, l) +
      lambda_2 * (kronecker_delta(i, k) * kronecker_delta(j, l) + kronecker_delta(i, l) * kronecker_delta(j, k));
  }

  /**
   * Assemble the system matrix and right-hand side vector.
   */
  void assemble()
  {
    const MeshBase & mesh = es.get_mesh();

    const unsigned int dim = mesh.mesh_dimension();

    LinearImplicitSystem & system = es.get_system<LinearImplicitSystem>("Elasticity");

    const unsigned int u_var = system.variable_number ("u");

    const DofMap & dof_map = system.get_dof_map();
    FEType fe_type = dof_map.variable_type(u_var);
    UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
    QGauss qrule (dim, fe_type.default_quadrature_order());
    fe->attach_quadrature_rule (&qrule);

    UniquePtr<FEBase> fe_face (FEBase::build(dim, fe_type));
    QGauss qface(dim-1, fe_type.default_quadrature_order());
    fe_face->attach_quadrature_rule (&qface);

    const std::vector<Real> & JxW = fe->get_JxW();
    const std::vector<std::vector<Real> > & phi = fe->get_phi();
    const std::vector<std::vector<RealGradient> > & dphi = fe->get_dphi();

    DenseMatrix<Number> Ke;
    DenseSubMatrix<Number> Ke_var[3][3] =
      {
        {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)},
        {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)},
        {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)}
      };

    DenseVector<Number> Fe;

    DenseSubVector<Number> Fe_var[3] =
      {DenseSubVector<Number>(Fe),
       DenseSubVector<Number>(Fe),
       DenseSubVector<Number>(Fe)};

    std::vector<dof_id_type> dof_indices;
    std::vector< std::vector<dof_id_type> > dof_indices_var(3);

    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

    for ( ; el != end_el; ++el)
      {
        const Elem * elem = *el;

        dof_map.dof_indices (elem, dof_indices);
        for (unsigned int var=0; var<3; var++)
          dof_map.dof_indices (elem, dof_indices_var[var], var);

        const unsigned int n_dofs   = dof_indices.size();
        const unsigned int n_var_dofs = dof_indices_var[0].size();

        fe->reinit (elem);

        Ke.resize (n_dofs, n_dofs);
        for (unsigned int var_i=0; var_i<3; var_i++)
          for (unsigned int var_j=0; var_j<3; var_j++)
            Ke_var[var_i][var_j].reposition (var_i*n_var_dofs, var_j*n_var_dofs, n_var_dofs, n_var_dofs);

        Fe.resize (n_dofs);
        for (unsigned int var=0; var<3; var++)
          Fe_var[var].reposition (var*n_var_dofs, n_var_dofs);

        for (unsigned int qp=0; qp<qrule.n_points(); qp++)
          {
            // assemble \int_Omega C_ijkl u_k,l v_i,j \dx
            for (unsigned int dof_i=0; dof_i<n_var_dofs; dof_i++)
              for (unsigned int dof_j=0; dof_j<n_var_dofs; dof_j++)
                for (unsigned int i=0; i<3; i++)
                  for (unsigned int j=0; j<3; j++)
                    for (unsigned int k=0; k<3; k++)
                      for (unsigned int l=0; l<3; l++)
                        Ke_var[i][k](dof_i,dof_j) +=
                          JxW[qp] * elasticity_tensor(i,j,k,l) * dphi[dof_j][qp](l) * dphi[dof_i][qp](j);

            // assemble \int_Omega f_i v_i \dx
            DenseVector<Number> f_vec(3);
            f_vec(0) =  0.;
            f_vec(1) =  0.;
            f_vec(2) = -1.;
            for (unsigned int dof_i=0; dof_i<n_var_dofs; dof_i++)
              for (unsigned int i=0; i<3; i++)
                Fe_var[i](dof_i) += JxW[qp] * (f_vec(i) * phi[dof_i][qp]);
          }

        // assemble \int_\Gamma g_i v_i \ds
        DenseVector<Number> g_vec(3);
        g_vec(0) = 0.;
        g_vec(1) = 0.;
        g_vec(2) = -1.;
        {
          for (unsigned int side=0; side<elem->n_sides(); side++)
            if (elem->neighbor_ptr(side) == libmesh_nullptr)
              {
                const std::vector<std::vector<Real> > & phi_face = fe_face->get_phi();
                const std::vector<Real> & JxW_face = fe_face->get_JxW();

                fe_face->reinit(elem, side);

                // Apply a traction
                for (unsigned int qp=0; qp<qface.n_points(); qp++)
                  if (mesh.get_boundary_info().has_boundary_id(elem, side, BOUNDARY_ID_MAX_X))
                    for (unsigned int dof_i=0; dof_i<n_var_dofs; dof_i++)
                      for (unsigned int i=0; i<3; i++)
                        Fe_var[i](dof_i) += JxW_face[qp] * (g_vec(i) * phi_face[dof_i][qp]);
              }
        }

        dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

        system.matrix->add_matrix (Ke, dof_indices);
        system.rhs->add_vector    (Fe, dof_indices);
      }
  }

  // Post-process the solution to compute stresses
  void compute_stresses()
  {
    const MeshBase & mesh = es.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    LinearImplicitSystem & system = es.get_system<LinearImplicitSystem>("Elasticity");

    unsigned int displacement_vars[3];
    displacement_vars[0] = system.variable_number ("u");
    displacement_vars[1] = system.variable_number ("v");
    displacement_vars[2] = system.variable_number ("w");
    const unsigned int u_var = system.variable_number ("u");

    const DofMap & dof_map = system.get_dof_map();
    FEType fe_type = dof_map.variable_type(u_var);
    UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
    QGauss qrule (dim, fe_type.default_quadrature_order());
    fe->attach_quadrature_rule (&qrule);

    const std::vector<Real> & JxW = fe->get_JxW();
    const std::vector<std::vector<Real> > & phi = fe->get_phi();
    const std::vector<std::vector<RealGradient> > & dphi = fe->get_dphi();

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
    unsigned int vonMises_var = stress_system.variable_number ("vonMises");

    // Storage for the stress dof indices on each element
    std::vector< std::vector<dof_id_type> > dof_indices_var(system.n_vars());
    std::vector<dof_id_type> stress_dof_indices_var;
    std::vector<dof_id_type> vonmises_dof_indices_var;

    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

    for ( ; el != end_el; ++el)
      {
        const Elem * elem = *el;

        for (unsigned int var=0; var<3; var++)
          dof_map.dof_indices (elem, dof_indices_var[var], displacement_vars[var]);

        const unsigned int n_var_dofs = dof_indices_var[0].size();

        fe->reinit (elem);

        std::vector< DenseMatrix<Number> > stress_tensor_qp(qrule.n_points());
        for (unsigned int qp=0; qp<qrule.n_points(); qp++)
          {
            stress_tensor_qp[qp].resize(3,3);

            // Row is variable u1, u2, or u3, column is x, y, or z
            DenseMatrix<Number> grad_u(3,3);
            for (unsigned int var_i=0; var_i<3; var_i++)
              for (unsigned int var_j=0; var_j<3; var_j++)
                for (unsigned int j=0; j<n_var_dofs; j++)
                  grad_u(var_i,var_j) += dphi[j][qp](var_j) * system.current_solution(dof_indices_var[var_i][j]);

            for (unsigned int var_i=0; var_i<3; var_i++)
              for (unsigned int var_j=0; var_j<3; var_j++)
                for (unsigned int k=0; k<3; k++)
                  for (unsigned int l=0; l<3; l++)
                    stress_tensor_qp[qp](var_i,var_j) += elasticity_tensor(var_i,var_j,k,l) * grad_u(k,l);
          }

        stress_dof_map.dof_indices (elem, vonmises_dof_indices_var, vonMises_var);
        std::vector< DenseMatrix<Number> > elem_sigma_vec(vonmises_dof_indices_var.size());
        for (std::size_t index=0; index<elem_sigma_vec.size(); index++)
          elem_sigma_vec[index].resize(3,3);

        // Below we project each component of the stress tensor onto a L2_LAGRANGE discretization.
        // Note that this gives a discontinuous stress plot on element boundaries, which is
        // appropriate. We then also get the von Mises stress from the projected stress tensor.
        unsigned int stress_var_index = 0;
        for (unsigned int var_i=0; var_i<3; var_i++)
          for (unsigned int var_j=var_i; var_j<3; var_j++)
            {
              stress_dof_map.dof_indices (elem, stress_dof_indices_var, sigma_vars[stress_var_index]);

              const unsigned int n_proj_dofs = stress_dof_indices_var.size();

              DenseMatrix<Real> Me(n_proj_dofs, n_proj_dofs);
              for (unsigned int qp=0; qp<qrule.n_points(); qp++)
                {
                  for(unsigned int i=0; i<n_proj_dofs; i++)
                    for(unsigned int j=0; j<n_proj_dofs; j++)
                      {
                        Me(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp]);
                      }
                }

              DenseVector<Number> Fe(n_proj_dofs);
              for (unsigned int qp=0; qp<qrule.n_points(); qp++)
                for(unsigned int i=0; i<n_proj_dofs; i++)
                  {
                    Fe(i) += JxW[qp] * stress_tensor_qp[qp](var_i,var_j) * phi[i][qp];
                  }

              DenseVector<Number> projected_data;
              Me.cholesky_solve(Fe, projected_data);

              for(unsigned int index=0; index<n_proj_dofs; index++)
                {
                  dof_id_type dof_index = stress_dof_indices_var[index];
                  if ((stress_system.solution->first_local_index() <= dof_index) &&
                      (dof_index < stress_system.solution->last_local_index()))
                    stress_system.solution->set(dof_index, projected_data(index));

                  elem_sigma_vec[index](var_i,var_j) = projected_data(index);
                }

              stress_var_index++;
            }

        for (std::size_t index=0; index<elem_sigma_vec.size(); index++)
          {
            elem_sigma_vec[index](1,0) = elem_sigma_vec[index](0,1);
            elem_sigma_vec[index](2,0) = elem_sigma_vec[index](0,2);
            elem_sigma_vec[index](2,1) = elem_sigma_vec[index](1,2);

            // Get the von Mises stress from the projected stress tensor
            Number vonMises_value = std::sqrt(0.5*(pow(elem_sigma_vec[index](0,0) - elem_sigma_vec[index](1,1), 2.) +
                                                   pow(elem_sigma_vec[index](1,1) - elem_sigma_vec[index](2,2), 2.) +
                                                   pow(elem_sigma_vec[index](2,2) - elem_sigma_vec[index](0,0), 2.) +
                                                   6.*(pow(elem_sigma_vec[index](0,1), 2.) +
                                                       pow(elem_sigma_vec[index](1,2), 2.) +
                                                       pow(elem_sigma_vec[index](2,0), 2.))));

            dof_id_type dof_index = vonmises_dof_indices_var[index];

            if ((stress_system.solution->first_local_index() <= dof_index) &&
                (dof_index < stress_system.solution->last_local_index()))
              stress_system.solution->set(dof_index, vonMises_value);
          }
      }

    // Should call close and update when we set vector entries directly
    stress_system.solution->close();
    stress_system.update();
  }
};


// Begin the main program.
int main (int argc, char ** argv)
{
  // Initialize libMesh and any dependent libraries
  LibMeshInit init (argc, argv);

  // Initialize the cantilever mesh
  const unsigned int dim = 3;

  // Make sure libMesh was compiled for 3D
  libmesh_example_requires(dim == LIBMESH_DIM, "3D support");

  // Create a 3D mesh distributed across the default MPI communicator.
  Mesh mesh(init.comm(), dim);
  MeshTools::Generation::build_cube (mesh,
                                     32,
                                     8,
                                     4,
                                     0., 1.*x_scaling,
                                     0., 0.3,
                                     0., 0.1,
                                     HEX8);

  // Print information about the mesh to the screen.
  mesh.print_info();

  // Let's add some node and edge boundary conditions
  // Each processor should know about each boundary condition it can
  // see, so we loop over all elements, not just local elements.
  MeshBase::const_element_iterator       el     = mesh.elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.elements_end();
  for ( ; el != end_el; ++el)
    {
      const Elem * elem = *el;

      unsigned int
        side_max_x = 0, side_min_y = 0,
        side_max_y = 0, side_max_z = 0;

      bool
        found_side_max_x = false, found_side_max_y = false,
        found_side_min_y = false, found_side_max_z = false;

      for (unsigned int side=0; side<elem->n_sides(); side++)
        {
          if (mesh.get_boundary_info().has_boundary_id(elem, side, BOUNDARY_ID_MAX_X))
            {
              side_max_x = side;
              found_side_max_x = true;
            }

          if (mesh.get_boundary_info().has_boundary_id(elem, side, BOUNDARY_ID_MIN_Y))
            {
              side_min_y = side;
              found_side_min_y = true;
            }

          if (mesh.get_boundary_info().has_boundary_id(elem, side, BOUNDARY_ID_MAX_Y))
            {
              side_max_y = side;
              found_side_max_y = true;
            }

          if (mesh.get_boundary_info().has_boundary_id(elem, side, BOUNDARY_ID_MAX_Z))
            {
              side_max_z = side;
              found_side_max_z = true;
            }
        }

      // If elem has sides on boundaries
      // BOUNDARY_ID_MAX_X, BOUNDARY_ID_MAX_Y, BOUNDARY_ID_MAX_Z
      // then let's set a node boundary condition
      if (found_side_max_x && found_side_max_y && found_side_max_z)
        for (unsigned int n=0; n<elem->n_nodes(); n++)
          if (elem->is_node_on_side(n, side_max_x) &&
              elem->is_node_on_side(n, side_max_y) &&
              elem->is_node_on_side(n, side_max_z))
            mesh.get_boundary_info().add_node(elem->node_ptr(n), NODE_BOUNDARY_ID);


      // If elem has sides on boundaries
      // BOUNDARY_ID_MAX_X and BOUNDARY_ID_MIN_Y
      // then let's set an edge boundary condition
      if (found_side_max_x && found_side_min_y)
        for (unsigned int e=0; e<elem->n_edges(); e++)
          if (elem->is_edge_on_side(e, side_max_x) &&
              elem->is_edge_on_side(e, side_min_y))
            mesh.get_boundary_info().add_edge(elem, e, EDGE_BOUNDARY_ID);
    }

  // Create an equation systems object.
  EquationSystems equation_systems (mesh);

  // Declare the system and its variables.
  // Create a system named "Elasticity"
  LinearImplicitSystem & system =
    equation_systems.add_system<LinearImplicitSystem> ("Elasticity");

#ifdef LIBMESH_HAVE_PETSC
  // Attach a SolverConfiguration object to system.linear_solver
  PetscLinearSolver<Number> * petsc_linear_solver =
    libmesh_cast_ptr<PetscLinearSolver<Number>*>(system.get_linear_solver());
  libmesh_assert(petsc_linear_solver);
  PetscSolverConfiguration petsc_solver_config(*petsc_linear_solver);
  petsc_linear_solver->set_solver_configuration(petsc_solver_config);
#endif

  // Add three displacement variables, u and v, to the system
  unsigned int u_var = system.add_variable("u", FIRST, LAGRANGE);
  unsigned int v_var = system.add_variable("v", FIRST, LAGRANGE);
  unsigned int w_var = system.add_variable("w", FIRST, LAGRANGE);

  LinearElasticity le(equation_systems);
  system.attach_assemble_object(le);

  std::set<boundary_id_type> boundary_ids;
  boundary_ids.insert(BOUNDARY_ID_MIN_X);
  boundary_ids.insert(NODE_BOUNDARY_ID);
  boundary_ids.insert(EDGE_BOUNDARY_ID);

  // Create a vector storing the variable numbers which the BC applies to
  std::vector<unsigned int> variables;
  variables.push_back(u_var);
  variables.push_back(v_var);
  variables.push_back(w_var);

  // Create a ZeroFunction to initialize dirichlet_bc
  ZeroFunction<> zf;

  // Most DirichletBoundary users will want to supply a "locally
  // indexed" functor
  DirichletBoundary dirichlet_bc(boundary_ids, variables, zf,
                                 LOCAL_VARIABLE_ORDER);

  // We must add the Dirichlet boundary condition _before_
  // we call equation_systems.init()
  system.get_dof_map().add_dirichlet_boundary(dirichlet_bc);

  // Also, initialize an ExplicitSystem to store stresses
  ExplicitSystem & stress_system =
    equation_systems.add_system<ExplicitSystem> ("StressSystem");

  stress_system.add_variable("sigma_00", FIRST, L2_LAGRANGE);
  stress_system.add_variable("sigma_01", FIRST, L2_LAGRANGE);
  stress_system.add_variable("sigma_02", FIRST, L2_LAGRANGE);
  stress_system.add_variable("sigma_11", FIRST, L2_LAGRANGE);
  stress_system.add_variable("sigma_12", FIRST, L2_LAGRANGE);
  stress_system.add_variable("sigma_22", FIRST, L2_LAGRANGE);
  stress_system.add_variable("vonMises", FIRST, L2_LAGRANGE);

  // Initialize the data structures for the equation system.
  equation_systems.init();

  // Print information about the system to the screen.
  equation_systems.print_info();

  // Solve the system
  system.solve();

  // Post-process the solution to compute the stresses
  le.compute_stresses();

  // Plot the solution
#ifdef LIBMESH_HAVE_EXODUS_API

  // Use single precision in this case (reduces the size of the exodus file)
  ExodusII_IO exo_io(mesh, /*single_precision=*/true);
  exo_io.write_discontinuous_exodusII("displacement_and_stress.exo", equation_systems);

#endif // #ifdef LIBMESH_HAVE_EXODUS_API

  // All done.
  return 0;
}
