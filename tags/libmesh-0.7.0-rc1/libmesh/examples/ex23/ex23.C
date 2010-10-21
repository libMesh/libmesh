// rbOOmit: An implementation of the Certified Reduced Basis method.
// Copyright (C) 2009, 2010 David J. Knezevic
//
//     This file is part of rbOOmit.

// rbOOmit is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
  
// rbOOmit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <cmath>
#include <set>

// Basic include file needed for the mesh functionality.
#include "libmesh.h"
#include "mesh.h"
#include "mesh_generation.h"
#include "gmv_io.h"
#include "equation_systems.h"
#include "fe.h"
#include "quadrature_gauss.h"
#include "dof_map.h"
#include "sparse_matrix.h"
#include "numeric_vector.h"
#include "dense_matrix.h"
#include "dense_vector.h"
#include "fe_interface.h"
#include "getpot.h"
#include "o_string_stream.h"
#include "elem.h"
#include "simple_rb.h"

// In this example problem we use the Certified Reduced Basis method
// to solve a steady convection-diffusion problem on the unit square.
// The reduced basis method relies on an expansion of the PDE in the
// form \sum_q=1^Q_a theta_a^q(\mu) * a^q(u,v) = \sum_q=1^Q_f theta_f^q(\mu) f^q(v)
// where theta_a, theta_f are parameter dependent functions and 
// a^q, f^q are parameter independent operators (\mu denotes a parameter).

// We first attach the parameter dependent functions and paramater
// independent operators to the RBSystem. Then in Offline mode, a
// reduced basis space is generated and written out to the directory
// "offline_data". In Online mode, the reduced basis data in "offline_data"
// is read in and used to solve the reduced problem for the parameters
// specified in ex23.in.

// We also attach four outputs to the system which are averages over certain
// subregions of the domain. In Online mode, we print out the values of these
// outputs as well as rigorous error bounds with respect to the output
// associated with the "truth" finite element discretization.


// Populate the sets non-Dirichlet and Dirichlet dofs,
// used to impose boundary conditions
void initialize_dirichlet_dofs(FEMContext &c,
                               System& system,
                               std::set<unsigned int>& dirichlet_dofs_set);

// Expansion of the PDE operator
Number theta_a_0(std::vector<Real>& ) { return 0.05; }
void A0(FEMContext&, System&);
Number theta_a_1(std::vector<Real>& mu) { return mu[0]; }
void A1(FEMContext&, System&);
Number theta_a_2(std::vector<Real>& mu) { return mu[1]; }
void A2(FEMContext&, System&);

// Expansion of the RHS
Number theta_f_0(std::vector<Real>& ) { return 1.; }
void F0(FEMContext&, System&);

// Expansion of the outputs
Number theta_l(std::vector<Real>& ) { return 1.; }
void L0_assembly(FEMContext&, System&);
void L1_assembly(FEMContext&, System&);
void L2_assembly(FEMContext&, System&);
void L3_assembly(FEMContext&, System&);

// The main program.
int main (int argc, char** argv)
{
  // Initialize libMesh.
  LibMeshInit init (argc, argv);

#if !defined(LIBMESH_HAVE_XDR)
  // We need XDR support to write out reduced bases
  libmesh_example_assert(false, "--enable-xdr");
#elif !defined(LIBMESH_HAVE_PETSC)
  // FIXME: This example currently segfaults with Trilinos?
  libmesh_example_assert(false, "--enable-petsc");
#endif

  // Parse the input file (ex23.in) using GetPot
  std::string parameters_filename = "ex23.in";
  GetPot infile(parameters_filename);

  unsigned int n_elem = infile("n_elem", 1);       // Determines the number of elements in the "truth" mesh
  const unsigned int dim = 2;                      // The number of spatial dimensions
  
  bool online_mode = infile("online_mode", false); // Are we in Online mode?

  // Create a one-dimensional mesh.
  Mesh mesh (dim);
  MeshTools::Generation::build_square (mesh,
                                       n_elem, n_elem,
                                       0., 1.,
                                       0., 1.,
                                       QUAD4);

  // Create an equation systems object.
  EquationSystems equation_systems (mesh);

  // We override RBSystem with SimpleRB in order to specialize a few functions.
  SimpleRB & system =
    equation_systems.add_system<SimpleRB> ("RBConvectionDiffusion");

  // Point the systems to the input file defining the problem
  system.parameters_filename = parameters_filename;

  system.attach_dirichlet_dof_initialization(initialize_dirichlet_dofs);
  
  // Attach the expansion of the PDE operator. The third argument
  // here refers to assembly over the boundary of the mesh, but this
  // problem only requires internal assembly and hence it is set to NULL
  system.attach_A_q(theta_a_0, A0, NULL);
  system.attach_A_q(theta_a_1, A1, NULL);
  system.attach_A_q(theta_a_2, A2, NULL);
  
  // Attach the expansion of the RHS
  system.attach_F_q(theta_f_0, F0, NULL);
  
  // Attach output assembly
  std::vector<theta_q_fptr> theta_l_vector(1); theta_l_vector[0] = theta_l;
  std::vector<affine_assembly_fptr> L_bndry_vector(1);  L_bndry_vector[0] = NULL;
  std::vector<affine_assembly_fptr> L0_intrr_vector(1); L0_intrr_vector[0] = L0_assembly;
  std::vector<affine_assembly_fptr> L1_intrr_vector(1); L1_intrr_vector[0] = L1_assembly;
  std::vector<affine_assembly_fptr> L2_intrr_vector(1); L2_intrr_vector[0] = L2_assembly;
  std::vector<affine_assembly_fptr> L3_intrr_vector(1); L3_intrr_vector[0] = L3_assembly;
  system.attach_output(theta_l_vector, L0_intrr_vector, L_bndry_vector);
  system.attach_output(theta_l_vector, L1_intrr_vector, L_bndry_vector);
  system.attach_output(theta_l_vector, L2_intrr_vector, L_bndry_vector);
  system.attach_output(theta_l_vector, L3_intrr_vector, L_bndry_vector);
  
  // We reuse the operator A0 as the inner product matrix
  system.attach_inner_prod_assembly(A0);

  // Initialize the data structures for the equation system.
  equation_systems.init ();

  if(system.initialize_calN_dependent_data)
  {
    // Print out some information about the "truth" discretization
    mesh.print_info();
    equation_systems.print_info();
  }
  
  // Initialize the RB data structures.
  // If we're in Offline Mode (online_mode == false) then
  // also pre-assemble the Dirichlet dofs list, RB operators etc
  system.initialize_RB_system(online_mode);

  if(!online_mode)
  {
    // Compute the reduced basis space by computing "snapshots", i.e.
    // "truth" solves, at well-chosen parameter values and employing
    // these snapshots as basis functions.
    system.train_reduced_basis();
    
    // Write out the reduced basis data
    system.write_offline_data_to_files();
  }
  else
  {
    // Read in the reduced basis data
    system.read_offline_data_from_files();
    
    // Get the parameters at which we do a reduced basis solve
    unsigned int online_N = infile("online_N",1);
    std::vector<Real> online_mu_vector(system.get_n_params());
    for(unsigned int i=0; i<system.get_n_params(); i++)
    {
      online_mu_vector[i] = infile("online_mu", online_mu_vector[i], i);
    }

    // Set the parameters to online_mu_vector
    system.set_current_parameters(online_mu_vector);
    system.print_current_parameters();

    // Now do the Online solve using the precomputed reduced basis
    system.RB_solve(online_N);

    // Print out outputs as well as the corresponding output error bounds.
    std::cout << "output 1, value = " << system.RB_outputs[0]
              << ", bound = " << system.RB_output_error_bounds[0]
              << std::endl;
    std::cout << "output 2, value = " << system.RB_outputs[1]
              << ", bound = " << system.RB_output_error_bounds[1]
              << std::endl;
    std::cout << "output 3, value = " << system.RB_outputs[2]
              << ", bound = " << system.RB_output_error_bounds[2]
              << std::endl;
    std::cout << "output 4, value = " << system.RB_outputs[3]
              << ", bound = " << system.RB_output_error_bounds[3]
              << std::endl;

    // If we stored the basis functions in the offline_data directory, plot the RB solution
    if(system.store_basis_functions)
    {
      system.load_RB_solution();
      GMVIO(mesh).write_equation_systems ("RB_sol.gmv",equation_systems);
      
      system.load_basis_function(0);
      GMVIO(mesh).write_equation_systems ("bf0.gmv",equation_systems);
    }
  }

  return 0;
}

// The Laplacian
void A0(FEMContext &c, System& system)
{
  SimpleRB& rb_system = libmesh_cast_ref<SimpleRB&>(system);
  const unsigned int u_var = rb_system.u_var;

  const std::vector<Real> &JxW =
    c.element_fe_var[u_var]->get_JxW();

  // The velocity shape function gradients at interior
  // quadrature points.
  const std::vector<std::vector<RealGradient> >& dphi =
    c.element_fe_var[u_var]->get_dphi();

  // The number of local degrees of freedom in each variable
  const unsigned int n_u_dofs = c.dof_indices_var[u_var].size();

  // Now we will build the affine operator
  unsigned int n_qpoints = c.element_qrule->n_points();

  for (unsigned int qp=0; qp != n_qpoints; qp++)
    for (unsigned int i=0; i != n_u_dofs; i++)
      for (unsigned int j=0; j != n_u_dofs; j++)
        c.elem_jacobian(i,j) += JxW[qp] * dphi[j][qp]*dphi[i][qp];
}

// Convection in the x-direction
void A1(FEMContext &c, System& system)
{
  SimpleRB& rb_system = libmesh_cast_ref<SimpleRB&>(system);
  const unsigned int u_var = rb_system.u_var;

  const std::vector<Real> &JxW =
    c.element_fe_var[u_var]->get_JxW();

  const std::vector<std::vector<Real> >& phi =
    c.element_fe_var[u_var]->get_phi();

  const std::vector<std::vector<RealGradient> >& dphi =
    c.element_fe_var[u_var]->get_dphi();

  // The number of local degrees of freedom in each variable
  const unsigned int n_u_dofs = c.dof_indices_var[u_var].size();

  // Now we will build the affine operator
  unsigned int n_qpoints = c.element_qrule->n_points();

  for (unsigned int qp=0; qp != n_qpoints; qp++)
    for (unsigned int i=0; i != n_u_dofs; i++)
      for (unsigned int j=0; j != n_u_dofs; j++)
        c.elem_jacobian(i,j) += JxW[qp] *dphi[j][qp](0)*phi[i][qp];
}

// Convection in the y-direction
void A2(FEMContext &c, System& system)
{
  SimpleRB& rb_system = libmesh_cast_ref<SimpleRB&>(system);
  const unsigned int u_var = rb_system.u_var;

  const std::vector<Real> &JxW =
    c.element_fe_var[u_var]->get_JxW();

  const std::vector<std::vector<Real> >& phi =
    c.element_fe_var[u_var]->get_phi();

  const std::vector<std::vector<RealGradient> >& dphi =
    c.element_fe_var[u_var]->get_dphi();

  // The number of local degrees of freedom in each variable
  const unsigned int n_u_dofs = c.dof_indices_var[u_var].size();

  // Now we will build the affine operator
  unsigned int n_qpoints = c.element_qrule->n_points();

  for (unsigned int qp=0; qp != n_qpoints; qp++)
    for (unsigned int i=0; i != n_u_dofs; i++)
      for (unsigned int j=0; j != n_u_dofs; j++)
        c.elem_jacobian(i,j) += JxW[qp] *dphi[j][qp](1)*phi[i][qp];
}

// Source term, 1 throughout the domain
void F0(FEMContext &c, System& system)
{
  SimpleRB& rb_system = libmesh_cast_ref<SimpleRB&>(system);
  const unsigned int u_var = rb_system.u_var;

  const std::vector<Real> &JxW =
    c.element_fe_var[u_var]->get_JxW();

  const std::vector<std::vector<Real> >& phi =
    c.element_fe_var[u_var]->get_phi();

  // The number of local degrees of freedom in each variable
  const unsigned int n_u_dofs = c.dof_indices_var[u_var].size();

  // Now we will build the affine operator
  unsigned int n_qpoints = c.element_qrule->n_points();

  for (unsigned int qp=0; qp != n_qpoints; qp++)
    for (unsigned int i=0; i != n_u_dofs; i++)
      c.elem_residual(i) += JxW[qp] * ( 1.*phi[i][qp] );
}

// Output: Average value over the region [0.7,0.8]x[0.7,0.8]
void L0_assembly(FEMContext &c, System& system)
{
  SimpleRB& rb_system = libmesh_cast_ref<SimpleRB&>(system);
  const unsigned int u_var = rb_system.u_var;

  const std::vector<Real> &JxW =
    c.element_fe_var[u_var]->get_JxW();

  const std::vector<std::vector<Real> >& phi =
    c.element_fe_var[u_var]->get_phi();

  // The number of local degrees of freedom in each variable
  const unsigned int n_u_dofs = c.dof_indices_var[u_var].size();

  // Now we will build the affine operator
  unsigned int n_qpoints = c.element_qrule->n_points();

  Point centroid = c.elem->centroid();
  if( (0.7 <= centroid(0)) && (centroid(0) <= 0.8) &&
      (0.7 <= centroid(1)) && (centroid(1) <= 0.8) )
    for (unsigned int qp=0; qp != n_qpoints; qp++)
      for (unsigned int i=0; i != n_u_dofs; i++)
        c.elem_residual(i) += JxW[qp] * ( 1.*phi[i][qp] ) / 0.01;
}

// Output: Average value over the region [0.2,0.3]x[0.7,0.8]
void L1_assembly(FEMContext &c, System& system)
{
  SimpleRB& rb_system = libmesh_cast_ref<SimpleRB&>(system);
  const unsigned int u_var = rb_system.u_var;

  const std::vector<Real> &JxW =
    c.element_fe_var[u_var]->get_JxW();

  const std::vector<std::vector<Real> >& phi =
    c.element_fe_var[u_var]->get_phi();

  // The number of local degrees of freedom in each variable
  const unsigned int n_u_dofs = c.dof_indices_var[u_var].size();

  // Now we will build the affine operator
  unsigned int n_qpoints = c.element_qrule->n_points();

  Point centroid = c.elem->centroid();
  if( (0.2 <= centroid(0)) && (centroid(0) <= 0.3) &&
      (0.7 <= centroid(1)) && (centroid(1) <= 0.8) )
    for (unsigned int qp=0; qp != n_qpoints; qp++)
      for (unsigned int i=0; i != n_u_dofs; i++)
        c.elem_residual(i) += JxW[qp] * ( 1.*phi[i][qp] ) / 0.01;
}

// Output: Average value over the region [0.2,0.3]x[0.2,0.3]
void L2_assembly(FEMContext &c, System& system)
{
  SimpleRB& rb_system = libmesh_cast_ref<SimpleRB&>(system);
  const unsigned int u_var = rb_system.u_var;

  const std::vector<Real> &JxW =
    c.element_fe_var[u_var]->get_JxW();

  const std::vector<std::vector<Real> >& phi =
    c.element_fe_var[u_var]->get_phi();

  // The number of local degrees of freedom in each variable
  const unsigned int n_u_dofs = c.dof_indices_var[u_var].size();

  // Now we will build the affine operator
  unsigned int n_qpoints = c.element_qrule->n_points();

  Point centroid = c.elem->centroid();
  if( (0.2 <= centroid(0)) && (centroid(0) <= 0.3) &&
      (0.2 <= centroid(1)) && (centroid(1) <= 0.3) )
    for (unsigned int qp=0; qp != n_qpoints; qp++)
      for (unsigned int i=0; i != n_u_dofs; i++)
        c.elem_residual(i) += JxW[qp] * ( 1.*phi[i][qp] ) / 0.01;
}

// Output: Average value over the region [0.7,0.8]x[0.2,0.3]
void L3_assembly(FEMContext &c, System& system)
{
  SimpleRB& rb_system = libmesh_cast_ref<SimpleRB&>(system);
  const unsigned int u_var = rb_system.u_var;

  const std::vector<Real> &JxW =
    c.element_fe_var[u_var]->get_JxW();

  const std::vector<std::vector<Real> >& phi =
    c.element_fe_var[u_var]->get_phi();

  // The number of local degrees of freedom in each variable
  const unsigned int n_u_dofs = c.dof_indices_var[u_var].size();

  // Now we will build the affine operator
  unsigned int n_qpoints = c.element_qrule->n_points();

  Point centroid = c.elem->centroid();
  if( (0.7 <= centroid(0)) && (centroid(0) <= 0.8) &&
      (0.2 <= centroid(1)) && (centroid(1) <= 0.3) )
    for (unsigned int qp=0; qp != n_qpoints; qp++)
      for (unsigned int i=0; i != n_u_dofs; i++)
        c.elem_residual(i) += JxW[qp] * ( 1.*phi[i][qp] ) / 0.01;
}

// Build a list (element-by-element) of the Dirichlet dofs. In this case
// all boundary dofs are Dirichlet.
void initialize_dirichlet_dofs(FEMContext &c, System& system,
                               std::set<unsigned int>& dirichlet_dofs_set)
{
  SimpleRB& rb_system = libmesh_cast_ref<SimpleRB&>(system);
  const unsigned int u_var = rb_system.u_var;

  std::vector<unsigned int> side_dofs;
  FEInterface::dofs_on_side(c.elem, c.dim, c.element_fe_var[u_var]->get_fe_type(),
                            c.side, side_dofs);

  for(unsigned int ii=0; ii<side_dofs.size(); ii++)
    dirichlet_dofs_set.insert(c.dof_indices[side_dofs[ii]]);
}
