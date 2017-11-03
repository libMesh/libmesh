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

// C++ includes
#include <limits> // for std::numeric_limits::max
#include <math.h>    // for sqrt


// Local Includes
#include "libmesh/hp_coarsentest.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dof_map.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_interface.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/elem.h"
#include "libmesh/error_vector.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/quadrature.h"
#include "libmesh/system.h"
#include "libmesh/tensor_value.h"

#ifdef LIBMESH_ENABLE_AMR

namespace libMesh
{

//-----------------------------------------------------------------
// HPCoarsenTest implementations

void HPCoarsenTest::add_projection(const System & system,
                                   const Elem * elem,
                                   unsigned int var)
{
  // If we have children, we need to add their projections instead
  if (!elem->active())
    {
      libmesh_assert(!elem->subactive());
      for (auto & child : elem->child_ref_range())
        this->add_projection(system, &child, var);
      return;
    }

  // The DofMap for this system
  const DofMap & dof_map = system.get_dof_map();

  // The type of finite element to use for this variable
  const FEType & fe_type = dof_map.variable_type (var);

  const FEContinuity cont = fe->get_continuity();

  fe->reinit(elem);

  dof_map.dof_indices(elem, dof_indices, var);

  const unsigned int n_dofs =
    cast_int<unsigned int>(dof_indices.size());

  FEInterface::inverse_map (system.get_mesh().mesh_dimension(),
                            fe_type, coarse, *xyz_values, coarse_qpoints);

  fe_coarse->reinit(coarse, &coarse_qpoints);

  const unsigned int n_coarse_dofs =
    cast_int<unsigned int>(phi_coarse->size());

  if (Uc.size() == 0)
    {
      Ke.resize(n_coarse_dofs, n_coarse_dofs);
      Ke.zero();
      Fe.resize(n_coarse_dofs);
      Fe.zero();
      Uc.resize(n_coarse_dofs);
      Uc.zero();
    }
  libmesh_assert_equal_to (Uc.size(), phi_coarse->size());

  // Loop over the quadrature points
  for (unsigned int qp=0; qp<qrule->n_points(); qp++)
    {
      // The solution value at the quadrature point
      Number val = libMesh::zero;
      Gradient grad;
      Tensor hess;

      for (unsigned int i=0; i != n_dofs; i++)
        {
          dof_id_type dof_num = dof_indices[i];
          val += (*phi)[i][qp] *
            system.current_solution(dof_num);
          if (cont == C_ZERO || cont == C_ONE)
            grad.add_scaled((*dphi)[i][qp],system.current_solution(dof_num));
          // grad += (*dphi)[i][qp] *
          //  system.current_solution(dof_num);
          if (cont == C_ONE)
            hess.add_scaled((*d2phi)[i][qp], system.current_solution(dof_num));
          // hess += (*d2phi)[i][qp] *
          //  system.current_solution(dof_num);
        }

      // The projection matrix and vector
      for (unsigned int i=0; i != Fe.size(); ++i)
        {
          Fe(i) += (*JxW)[qp] *
            (*phi_coarse)[i][qp]*val;
          if (cont == C_ZERO || cont == C_ONE)
            Fe(i) += (*JxW)[qp] *
              (grad*(*dphi_coarse)[i][qp]);
          if (cont == C_ONE)
            Fe(i) += (*JxW)[qp] *
              hess.contract((*d2phi_coarse)[i][qp]);
          // Fe(i) += (*JxW)[qp] *
          //  (*d2phi_coarse)[i][qp].contract(hess);

          for (unsigned int j=0; j != Fe.size(); ++j)
            {
              Ke(i,j) += (*JxW)[qp] *
                (*phi_coarse)[i][qp]*(*phi_coarse)[j][qp];
              if (cont == C_ZERO || cont == C_ONE)
                Ke(i,j) += (*JxW)[qp] *
                  (*dphi_coarse)[i][qp]*(*dphi_coarse)[j][qp];
              if (cont == C_ONE)
                Ke(i,j) += (*JxW)[qp] *
                  ((*d2phi_coarse)[i][qp].contract((*d2phi_coarse)[j][qp]));
            }
        }
    }
}

void HPCoarsenTest::select_refinement (System & system)
{
  LOG_SCOPE("select_refinement()", "HPCoarsenTest");

  // The current mesh
  MeshBase & mesh = system.get_mesh();

  // The dimensionality of the mesh
  const unsigned int dim = mesh.mesh_dimension();

  // The number of variables in the system
  const unsigned int n_vars = system.n_vars();

  // The DofMap for this system
  const DofMap & dof_map = system.get_dof_map();

  // The system number (for doing bad hackery)
  const unsigned int sys_num = system.number();

  // Check for a valid component_scale
  if (!component_scale.empty())
    {
      if (component_scale.size() != n_vars)
        libmesh_error_msg("ERROR: component_scale is the wrong size:\n" \
                          << " component_scale.size()=" \
                          << component_scale.size()     \
                          << "\n n_vars=" \
                          << n_vars);
    }
  else
    {
      // No specified scaling.  Scale all variables by one.
      component_scale.resize (n_vars, 1.0);
    }

  // Resize the error_per_cell vectors to handle
  // the number of elements, initialize them to 0.
  std::vector<ErrorVectorReal> h_error_per_cell(mesh.max_elem_id(), 0.);
  std::vector<ErrorVectorReal> p_error_per_cell(mesh.max_elem_id(), 0.);

  // Loop over all the variables in the system
  for (unsigned int var=0; var<n_vars; var++)
    {
      // Possibly skip this variable
      if (!component_scale.empty())
        if (component_scale[var] == 0.0) continue;

      // The type of finite element to use for this variable
      const FEType & fe_type = dof_map.variable_type (var);

      // Finite element objects for a fine (and probably a coarse)
      // element will be needed
      fe = FEBase::build (dim, fe_type);
      fe_coarse = FEBase::build (dim, fe_type);

      // Any cached coarse element results have expired
      coarse = libmesh_nullptr;
      unsigned int cached_coarse_p_level = 0;

      const FEContinuity cont = fe->get_continuity();
      libmesh_assert (cont == DISCONTINUOUS || cont == C_ZERO ||
                      cont == C_ONE);

      // Build an appropriate quadrature rule
      qrule = fe_type.default_quadrature_rule(dim);

      // Tell the refined finite element about the quadrature
      // rule.  The coarse finite element need not know about it
      fe->attach_quadrature_rule (qrule.get());

      // We will always do the integration
      // on the fine elements.  Get their Jacobian values, etc..
      JxW = &(fe->get_JxW());
      xyz_values = &(fe->get_xyz());

      // The shape functions
      phi = &(fe->get_phi());
      phi_coarse = &(fe_coarse->get_phi());

      // The shape function derivatives
      if (cont == C_ZERO || cont == C_ONE)
        {
          dphi = &(fe->get_dphi());
          dphi_coarse = &(fe_coarse->get_dphi());
        }

      // The shape function second derivatives
      if (cont == C_ONE)
        {
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
          d2phi = &(fe->get_d2phi());
          d2phi_coarse = &(fe_coarse->get_d2phi());
#else
          libmesh_error_msg("Minimization of H2 error without second derivatives is not possible.");
#endif
        }

      // Iterate over all the active elements in the mesh
      // that live on this processor.
      for (const auto & elem : mesh.active_local_element_ptr_range())
        {
          // We're only checking elements that are already flagged for h
          // refinement
          if (elem->refinement_flag() != Elem::REFINE)
            continue;

          const dof_id_type e_id = elem->id();

          // Find the projection onto the parent element,
          // if necessary
          if (elem->parent() &&
              (coarse != elem->parent() ||
               cached_coarse_p_level != elem->p_level()))
            {
              Uc.resize(0);

              coarse = elem->parent();
              cached_coarse_p_level = elem->p_level();

              unsigned int old_parent_level = coarse->p_level();
              (const_cast<Elem *>(coarse))->hack_p_level(elem->p_level());

              this->add_projection(system, coarse, var);

              (const_cast<Elem *>(coarse))->hack_p_level(old_parent_level);

              // Solve the h-coarsening projection problem
              Ke.cholesky_solve(Fe, Uc);
            }

          fe->reinit(elem);

          // Get the DOF indices for the fine element
          dof_map.dof_indices (elem, dof_indices, var);

          // The number of quadrature points
          const unsigned int n_qp = qrule->n_points();

          // The number of DOFS on the fine element
          const unsigned int n_dofs =
            cast_int<unsigned int>(dof_indices.size());

          // The number of nodes on the fine element
          const unsigned int n_nodes = elem->n_nodes();

          // The average element value (used as an ugly hack
          // when we have nothing p-coarsened to compare to)
          // Real average_val = 0.;
          Number average_val = 0.;

          // Calculate this variable's contribution to the p
          // refinement error

          if (elem->p_level() == 0)
            {
              unsigned int n_vertices = 0;
              for (unsigned int n = 0; n != n_nodes; ++n)
                if (elem->is_vertex(n))
                  {
                    n_vertices++;
                    const Node & node = elem->node_ref(n);
                    average_val += system.current_solution
                      (node.dof_number(sys_num,var,0));
                  }
              average_val /= n_vertices;
            }
          else
            {
              unsigned int old_elem_level = elem->p_level();
              (const_cast<Elem *>(elem))->hack_p_level(old_elem_level - 1);

              fe_coarse->reinit(elem, &(qrule->get_points()));

              const unsigned int n_coarse_dofs =
                cast_int<unsigned int>(phi_coarse->size());

              (const_cast<Elem *>(elem))->hack_p_level(old_elem_level);

              Ke.resize(n_coarse_dofs, n_coarse_dofs);
              Ke.zero();
              Fe.resize(n_coarse_dofs);
              Fe.zero();

              // Loop over the quadrature points
              for (unsigned int qp=0; qp<qrule->n_points(); qp++)
                {
                  // The solution value at the quadrature point
                  Number val = libMesh::zero;
                  Gradient grad;
                  Tensor hess;

                  for (unsigned int i=0; i != n_dofs; i++)
                    {
                      dof_id_type dof_num = dof_indices[i];
                      val += (*phi)[i][qp] *
                        system.current_solution(dof_num);
                      if (cont == C_ZERO || cont == C_ONE)
                        grad.add_scaled((*dphi)[i][qp], system.current_solution(dof_num));
                      // grad += (*dphi)[i][qp] *
                      //  system.current_solution(dof_num);
                      if (cont == C_ONE)
                        hess.add_scaled((*d2phi)[i][qp], system.current_solution(dof_num));
                      // hess += (*d2phi)[i][qp] *
                      //  system.current_solution(dof_num);
                    }

                  // The projection matrix and vector
                  for (unsigned int i=0; i != Fe.size(); ++i)
                    {
                      Fe(i) += (*JxW)[qp] *
                        (*phi_coarse)[i][qp]*val;
                      if (cont == C_ZERO || cont == C_ONE)
                        Fe(i) += (*JxW)[qp] *
                          grad * (*dphi_coarse)[i][qp];
                      if (cont == C_ONE)
                        Fe(i) += (*JxW)[qp] *
                          hess.contract((*d2phi_coarse)[i][qp]);

                      for (unsigned int j=0; j != Fe.size(); ++j)
                        {
                          Ke(i,j) += (*JxW)[qp] *
                            (*phi_coarse)[i][qp]*(*phi_coarse)[j][qp];
                          if (cont == C_ZERO || cont == C_ONE)
                            Ke(i,j) += (*JxW)[qp] *
                              (*dphi_coarse)[i][qp]*(*dphi_coarse)[j][qp];
                          if (cont == C_ONE)
                            Ke(i,j) += (*JxW)[qp] *
                              ((*d2phi_coarse)[i][qp].contract((*d2phi_coarse)[j][qp]));
                        }
                    }
                }

              // Solve the p-coarsening projection problem
              Ke.cholesky_solve(Fe, Up);
            }

          // loop over the integration points on the fine element
          for (unsigned int qp=0; qp<n_qp; qp++)
            {
              Number value_error = 0.;
              Gradient grad_error;
              Tensor hessian_error;
              for (unsigned int i=0; i<n_dofs; i++)
                {
                  const dof_id_type dof_num = dof_indices[i];
                  value_error += (*phi)[i][qp] *
                    system.current_solution(dof_num);
                  if (cont == C_ZERO || cont == C_ONE)
                    grad_error.add_scaled((*dphi)[i][qp], system.current_solution(dof_num));
                  // grad_error += (*dphi)[i][qp] *
                  //  system.current_solution(dof_num);
                  if (cont == C_ONE)
                    hessian_error.add_scaled((*d2phi)[i][qp], system.current_solution(dof_num));
                  // hessian_error += (*d2phi)[i][qp] *
                  //    system.current_solution(dof_num);
                }
              if (elem->p_level() == 0)
                {
                  value_error -= average_val;
                }
              else
                {
                  for (unsigned int i=0; i<Up.size(); i++)
                    {
                      value_error -= (*phi_coarse)[i][qp] * Up(i);
                      if (cont == C_ZERO || cont == C_ONE)
                        grad_error.subtract_scaled((*dphi_coarse)[i][qp], Up(i));
                      // grad_error -= (*dphi_coarse)[i][qp] * Up(i);
                      if (cont == C_ONE)
                        hessian_error.subtract_scaled((*d2phi_coarse)[i][qp], Up(i));
                      // hessian_error -= (*d2phi_coarse)[i][qp] * Up(i);
                    }
                }

              p_error_per_cell[e_id] += static_cast<ErrorVectorReal>
                (component_scale[var] *
                 (*JxW)[qp] * TensorTools::norm_sq(value_error));
              if (cont == C_ZERO || cont == C_ONE)
                p_error_per_cell[e_id] += static_cast<ErrorVectorReal>
                  (component_scale[var] *
                   (*JxW)[qp] * grad_error.norm_sq());
              if (cont == C_ONE)
                p_error_per_cell[e_id] += static_cast<ErrorVectorReal>
                  (component_scale[var] *
                   (*JxW)[qp] * hessian_error.norm_sq());
            }

          // Calculate this variable's contribution to the h
          // refinement error

          if (!elem->parent())
            {
              // For now, we'll always start with an h refinement
              h_error_per_cell[e_id] =
                std::numeric_limits<ErrorVectorReal>::max() / 2;
            }
          else
            {
              FEInterface::inverse_map (dim, fe_type, coarse,
                                        *xyz_values, coarse_qpoints);

              unsigned int old_parent_level = coarse->p_level();
              (const_cast<Elem *>(coarse))->hack_p_level(elem->p_level());

              fe_coarse->reinit(coarse, &coarse_qpoints);

              (const_cast<Elem *>(coarse))->hack_p_level(old_parent_level);

              // The number of DOFS on the coarse element
              unsigned int n_coarse_dofs =
                cast_int<unsigned int>(phi_coarse->size());

              // Loop over the quadrature points
              for (unsigned int qp=0; qp<n_qp; qp++)
                {
                  // The solution difference at the quadrature point
                  Number value_error = libMesh::zero;
                  Gradient grad_error;
                  Tensor hessian_error;

                  for (unsigned int i=0; i != n_dofs; ++i)
                    {
                      const dof_id_type dof_num = dof_indices[i];
                      value_error += (*phi)[i][qp] *
                        system.current_solution(dof_num);
                      if (cont == C_ZERO || cont == C_ONE)
                        grad_error.add_scaled((*dphi)[i][qp], system.current_solution(dof_num));
                      // grad_error += (*dphi)[i][qp] *
                      //  system.current_solution(dof_num);
                      if (cont == C_ONE)
                        hessian_error.add_scaled((*d2phi)[i][qp], system.current_solution(dof_num));
                      // hessian_error += (*d2phi)[i][qp] *
                      //  system.current_solution(dof_num);
                    }

                  for (unsigned int i=0; i != n_coarse_dofs; ++i)
                    {
                      value_error -= (*phi_coarse)[i][qp] * Uc(i);
                      if (cont == C_ZERO || cont == C_ONE)
                        // grad_error -= (*dphi_coarse)[i][qp] * Uc(i);
                        grad_error.subtract_scaled((*dphi_coarse)[i][qp], Uc(i));
                      if (cont == C_ONE)
                        hessian_error.subtract_scaled((*d2phi_coarse)[i][qp], Uc(i));
                      // hessian_error -= (*d2phi_coarse)[i][qp] * Uc(i);
                    }

                  h_error_per_cell[e_id] += static_cast<ErrorVectorReal>
                    (component_scale[var] *
                     (*JxW)[qp] * TensorTools::norm_sq(value_error));
                  if (cont == C_ZERO || cont == C_ONE)
                    h_error_per_cell[e_id] += static_cast<ErrorVectorReal>
                      (component_scale[var] *
                       (*JxW)[qp] * grad_error.norm_sq());
                  if (cont == C_ONE)
                    h_error_per_cell[e_id] += static_cast<ErrorVectorReal>
                      (component_scale[var] *
                       (*JxW)[qp] * hessian_error.norm_sq());
                }

            }
        }
    }

  // Now that we've got our approximations for p_error and h_error, let's see
  // if we want to switch any h refinement flags to p refinement

  // Iterate over all the active elements in the mesh
  // that live on this processor.
  for (auto & elem : mesh.active_local_element_ptr_range())
    {
      // We're only checking elements that are already flagged for h
      // refinement
      if (elem->refinement_flag() != Elem::REFINE)
        continue;

      const dof_id_type e_id = elem->id();

      unsigned int dofs_per_elem = 0, dofs_per_p_elem = 0;

      // Loop over all the variables in the system
      for (unsigned int var=0; var<n_vars; var++)
        {
          // The type of finite element to use for this variable
          const FEType & fe_type = dof_map.variable_type (var);

          // FIXME: we're overestimating the number of DOFs added by h
          // refinement
          FEType elem_fe_type = fe_type;
          elem_fe_type.order =
            static_cast<Order>(fe_type.order + elem->p_level());
          dofs_per_elem +=
            FEInterface::n_dofs(dim, elem_fe_type, elem->type());

          elem_fe_type.order =
            static_cast<Order>(fe_type.order + elem->p_level() + 1);
          dofs_per_p_elem +=
            FEInterface::n_dofs(dim, elem_fe_type, elem->type());
        }

      const unsigned int new_h_dofs = dofs_per_elem *
        (elem->n_children() - 1);

      const unsigned int new_p_dofs = dofs_per_p_elem -
        dofs_per_elem;

      /*
        libMesh::err << "Cell " << e_id << ": h = " << elem->hmax()
        << ", p = " << elem->p_level() + 1 << "," << std::endl
        << "     h_error = " << h_error_per_cell[e_id]
        << ", p_error = " << p_error_per_cell[e_id] << std::endl
        << "     new_h_dofs = " << new_h_dofs
        << ", new_p_dofs = " << new_p_dofs << std::endl;
      */
      const Real p_value =
        std::sqrt(p_error_per_cell[e_id]) * p_weight / new_p_dofs;
      const Real h_value =
        std::sqrt(h_error_per_cell[e_id]) /
        static_cast<Real>(new_h_dofs);
      if (p_value > h_value)
        {
          elem->set_p_refinement_flag(Elem::REFINE);
          elem->set_refinement_flag(Elem::DO_NOTHING);
        }
    }

  // libMesh::MeshRefinement will now assume that users have set
  // refinement flags consistently on all processors, but all we've
  // done so far is set refinement flags on local elements.  Let's
  // make sure that flags on geometrically ghosted elements are all
  // set to whatever their owners decided.
  MeshRefinement(mesh).make_flags_parallel_consistent();
}

} // namespace libMesh

#endif // #ifdef LIBMESH_ENABLE_AMR
