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


// <h1>Miscellaneous Example 15 - Laplace equation on unbound domain: Gravitational field computation</h1>
// \author Hubert Weissmann
// \date 2021
//
// This example uses second order infinite elements to solve the Laplace equation to compute the gravitational
// field of three prism-shaped objects with same density. Probably the most important question that you always wanted
// to study!
//
// Here infinite elements are used with frequency set to 0 so they can be used to describe the asymptotic behaviour
// of the solution accurately.
//
#include <iostream>
#include <fstream>
// libMesh include files.
#include "libmesh/getpot.h" // for input-argument parsing
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/elem.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dof_map.h"
#include "libmesh/point_locator_tree.h"
// for infinite elements:
#include "libmesh/inf_fe.h"
#include "libmesh/inf_elem_builder.h"
// for refinement:
#include "libmesh/error_vector.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/kelly_error_estimator.h"
// for distributed mesh editing
#include "libmesh/parallel_ghost_sync.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

void assemble_func(EquationSystems & es, const std::string & system_name)
{
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
  const MeshBase& mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  LinearImplicitSystem & f_system = es.get_system<LinearImplicitSystem> (system_name);
  const DofMap& dof_map = f_system.get_dof_map();
  const FEType fe_type = dof_map.variable_type(0);
  std::unique_ptr<FEBase> fe (FEBase::build(dim, fe_type));
  std::unique_ptr<FEBase> inf_fe (FEBase::build_InfFE(dim, fe_type));

  const std::set<dof_id_type > charged_objects = es.parameters.get<std::set<dof_id_type> >("charged_elem_id");
  const auto qrule=fe_type.default_quadrature_rule(/*dim = */ 3, /*extra order = */ 3);

  fe->attach_quadrature_rule (&*qrule);
  inf_fe->attach_quadrature_rule (&*qrule);

  DenseVector<Number> f;
  DenseMatrix<Number> M;

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;

  // Now we will loop over all the elements in the mesh that
  // live on the local processor.
  for (const auto & elem : mesh.active_local_element_ptr_range())
    {
      // Get the degree of freedom indices for the
      // current element.  These define where in the global
      // matrix and right-hand-side this element will
      // contribute to.
      dof_map.dof_indices (elem, dof_indices);

      // unifyging finite and infinite elements
      FEBase * cfe = libmesh_nullptr;

      if (elem->infinite())
        cfe = inf_fe.get();
      else
        cfe = fe.get();

      // The element Jacobian * quadrature weight at each integration point.
      const std::vector<Real>&                        JxW = cfe->get_JxWxdecay_sq();
      // The element shape functions evaluated at the quadrature points.
      const std::vector<std::vector<Real> >&         phi  = cfe->get_phi_over_decayxR();
      const std::vector<std::vector<RealGradient> >& dphi = cfe->get_dphi_over_decayxR();
      const std::vector<Real>&                      sob_w = cfe->get_Sobolev_weightxR_sq();
      const std::vector<RealGradient>&             dsob_w = cfe->get_Sobolev_dweightxR_sq();

      // Compute the element-specific data for the current
      // element. This involves computing the location of the
      // quadrature points (q_point) and the shape functions
      // (phi, dphi) for the current element.
      cfe->reinit (elem);
      M.resize (dof_indices.size(), dof_indices.size());
      f.resize (dof_indices.size());

      Real rho;
      // check if charged_objects contains elem->id(), we add the corresponding charge to the rhs vector
      if (charged_objects.count(elem->id()))
        rho=1./elem->volume();
      else
        rho = 0;

      const unsigned int max_qp = cfe->n_quadrature_points();
      for (unsigned int qp=0; qp<max_qp; ++qp)
        {
          const unsigned int n_sf =
            FEInterface::n_dofs(cfe->get_fe_type(), elem);
          for (unsigned int i=0; i<n_sf; ++i)
            {
              if (rho > 0)
                f(i) -=4*pi*JxW[qp]*rho*phi[i][qp];

              // store the test functions derivative (which contains the extra sobolev weight)
              // separately; so we don't need to recompute it for each j.
              const RealGradient dphi_i=dphi[i][qp]*sob_w[qp]+phi[i][qp]*dsob_w[qp];
              for (unsigned int j=0; j<n_sf; ++j)
                {
                  M(i,j) += JxW[qp]*dphi_i*dphi[j][qp];
                }
            }
        }
      dof_map.constrain_element_matrix_and_vector(M, f, dof_indices);
      f_system.matrix->add_matrix(M, dof_indices);
      f_system.rhs->add_vector(f, dof_indices);
    }
  f_system.rhs->close();
  f_system.matrix->close();
  /**
   * All done!
   */
#else //ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
  // Avoid compiler warnings
  libmesh_ignore(es, system_name);
#endif
  return;
}

int main (int argc, char** argv)
{
  // Initialize libMesh and the dependent libraries.
  LibMeshInit init (argc, argv);

  // Tell the user what we are doing.
  std::cout << "Running " << argv[0];
  for (int i=1; i<argc; i++)
    std::cout << " " << argv[i];
  std::cout << std::endl << std::endl;

#ifndef LIBMESH_ENABLE_INFINITE_ELEMENTS
  libmesh_example_requires(false, "--enable-ifem");
#else

  // Skip this example if libMesh was compiled with <3 dimensions.
  // INFINITE ELEMENTS ARE IMPLEMENTED ONLY FOR 3 DIMENSIONS AT THE MOMENT.
  libmesh_example_requires(3 <= LIBMESH_DIM, "3D support");

  int dim = 3;
  // creation of an empty mesh-object
  Mesh mesh(init.comm(), dim);

  // We're going to keep track of elements by id, so we don't want
  // even a DistributedMesh to renumber them.
  mesh.allow_renumbering(false);

  // fill the meshes with a spherical grid of type HEX27 with radius r
  MeshTools::Generation::build_cube (mesh, /*nx= */5, /*ny=*/5, /*nz=*/3,
                                     /*xmin=*/ -5.2, /*xmax=*/5.2,
                                     /*ymin=*/ -5.2, /*ymax=*/5.2,
                                     /*zmin=*/ -3.2, /*zmax=*/3.2,
                                     PRISM18 /* Currently, these elements cannot be viewed by paraview.*/
                                     );

  mesh.print_info();

  // The infinite elements are attached to all elements that build the outer surface of the FEM-region.
  InfElemBuilder builder(mesh);
  builder.build_inf_elem(true);
  // find the neighbours; for correct linking the two areas
  mesh.find_neighbors();

  // Reassign subdomain_id() of all infinite elements.
  // and their neighbours. This makes finding the surface
  // between these elemests much easier.
  for (auto & elem : mesh.element_ptr_range())
    if (elem->infinite())
      {
        elem->subdomain_id() = 1;
        // the base elements are always the 0-th neighbor.
        elem->neighbor_ptr(0)->subdomain_id()=2;
      }

  // If we're on a distributed mesh then we might have missed some
  // ghosted base elements with remote infinite elements.
  if (!mesh.is_serial())
    {
      SyncSubdomainIds sync_obj(mesh);
      Parallel::sync_dofobject_data_by_id
        (mesh.comm(), mesh.elements_begin(), mesh.elements_end(),
         sync_obj);
    }

  // Now we set the sources of the field: prism-shaped objects that are
  // determined here by containing certain points:
  auto p_pt_lctr = mesh.sub_point_locator();
  auto & pt_lctr = *p_pt_lctr;
  pt_lctr.enable_out_of_mesh_mode();

  std::set<dof_id_type> charged_elem_ids;
  {
    Point pt_0(-3.,-3.0,-1.5);
    Point pt_1(2.,-2.6,-1.5);
    Point pt_2(2., 3.1, 1.7);
    const Elem * elem_0=pt_lctr(pt_0);
    if (elem_0)
      charged_elem_ids.insert(elem_0->id());
    const Elem * elem_1=pt_lctr(pt_1);
    if (elem_1)
      charged_elem_ids.insert(elem_1->id());
    const Elem * elem_2=pt_lctr(pt_2);
    if (elem_2)
      charged_elem_ids.insert(elem_2->id());

    // On a distributed mesh we might not have every point in a
    // semilocal element on every processor
    if (!mesh.is_serial())
      mesh.comm().set_union(charged_elem_ids);

    // But we should have every point on *some* processor
    libmesh_assert_equal_to(charged_elem_ids.size(), 3);
  }

  // Create an equation systems object
  EquationSystems eq_sys (mesh);

  // This is the only system added here.
  LinearImplicitSystem & eig_sys = eq_sys.add_system<LinearImplicitSystem> ("Poisson");

  eq_sys.parameters.set<std::set<dof_id_type> >("charged_elem_id")=charged_elem_ids;

  //set the complete type of the variable
  FEType fe_type(SECOND, LAGRANGE, FOURTH, JACOBI_20_00, CARTESIAN);

  // Name the variable of interest 'phi' and approximate it as \p fe_type.
  eig_sys.add_variable("phi", fe_type);

  // attach the name of the function that assembles the matrix equation:
  eig_sys.attach_assemble_function (assemble_func);

  // Initialize the data structures for the equation system.
  eq_sys.init();

  // Solve system. This function calls the assemble-functions.
  eig_sys.solve();

#ifdef LIBMESH_ENABLE_AMR
  for (unsigned int i=0; i< 2; ++i)
    {
      ErrorVector error;
      MeshRefinement mesh_refinement(mesh);
      KellyErrorEstimator error_estimator;

      error_estimator.estimate_error(eig_sys, error);
      mesh_refinement.refine_fraction()=0.3;
      mesh_refinement.flag_elements_by_error_fraction(error);
      error_estimator.estimate_error(eig_sys, error);
      mesh_refinement.refine_and_coarsen_elements();
      eq_sys.reinit();

      // in the refined mesh, find the elements that describe the
      // sources of gravitational field: Due to refinement, there are
      // successively more elements to account for the same object.
      std::set<dof_id_type> charged_elem_children;
      for(auto id : charged_elem_ids)
        {
          Elem * charged_elem = mesh.query_elem_ptr(id);
          if (!charged_elem)
            {
              libmesh_assert(!mesh.is_serial());
              continue;
            }

          if (charged_elem->has_children())
            {
              for(auto & child : charged_elem->child_ref_range())
                if (!child.is_remote())
                  charged_elem_children.insert(child.id());
            }
          else
            charged_elem_children.insert(charged_elem->id());
        }

      charged_elem_ids=charged_elem_children;
      if (!mesh.is_serial())
        mesh.comm().set_union(charged_elem_ids);

      eq_sys.parameters.set<std::set<dof_id_type> >("charged_elem_id")=charged_elem_ids;

      // re-assemble and than solve again.
      eig_sys.solve();
    }
#endif // LIBMESH_ENABLE_AMR

  // All done.
#endif // LIBMESH_ENABLE_INFINITE_ELEMENTS
  return 0;
}

