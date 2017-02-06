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
#include <algorithm> // for std::fill
#include <cstdlib> // *must* precede <cmath> for proper std:abs() on PGI, Sun Studio CC
#include <cmath>    // for sqrt


// Local Includes
#include "libmesh/libmesh_common.h"
#include "libmesh/exact_error_estimator.h"
#include "libmesh/dof_map.h"
#include "libmesh/equation_systems.h"
#include "libmesh/error_vector.h"
#include "libmesh/fe_base.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/elem.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_function.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/quadrature.h"
#include "libmesh/system.h"
#include "libmesh/tensor_tools.h"

namespace libMesh
{

//-----------------------------------------------------------------
// ErrorEstimator implementations
void ExactErrorEstimator::attach_exact_value (Number fptr(const Point & p,
                                                          const Parameters & parameters,
                                                          const std::string & sys_name,
                                                          const std::string & unknown_name))
{
  libmesh_assert(fptr);
  _exact_value = fptr;

  // We're not using a fine grid solution
  _equation_systems_fine = libmesh_nullptr;

  // We're not using user-provided functors
  this->clear_functors();
}


void ExactErrorEstimator::attach_exact_values (std::vector<FunctionBase<Number> *> f)
{
  // Clear out any previous _exact_values entries, then add a new
  // entry for each system.
  for (std::size_t i=0; i != _exact_values.size(); ++i)
    delete (_exact_values[i]);

  _exact_values.clear();
  _exact_values.resize(f.size(), libmesh_nullptr);

  // We use clone() to get non-sliced copies of FunctionBase
  // subclasses, but we don't currently put the resulting UniquePtrs
  // into an STL container.
  for (std::size_t i=0; i != f.size(); ++i)
    if (f[i])
      _exact_values[i] = f[i]->clone().release();
}


void ExactErrorEstimator::attach_exact_value (unsigned int sys_num,
                                              FunctionBase<Number> * f)
{
  if (_exact_values.size() <= sys_num)
    _exact_values.resize(sys_num+1, libmesh_nullptr);

  if (f)
    _exact_values[sys_num] = f->clone().release();
}


void ExactErrorEstimator::attach_exact_deriv (Gradient gptr(const Point & p,
                                                            const Parameters & parameters,
                                                            const std::string & sys_name,
                                                            const std::string & unknown_name))
{
  libmesh_assert(gptr);
  _exact_deriv = gptr;

  // We're not using a fine grid solution
  _equation_systems_fine = libmesh_nullptr;

  // We're not using user-provided functors
  this->clear_functors();
}


void ExactErrorEstimator::attach_exact_derivs (std::vector<FunctionBase<Gradient> *> g)
{
  // Clear out any previous _exact_derivs entries, then add a new
  // entry for each system.
  for (std::size_t i=0; i != _exact_derivs.size(); ++i)
    delete (_exact_derivs[i]);

  _exact_derivs.clear();
  _exact_derivs.resize(g.size(), libmesh_nullptr);

  // We use clone() to get non-sliced copies of FunctionBase
  // subclasses, but we don't currently put the resulting UniquePtrs
  // into an STL container.
  for (std::size_t i=0; i != g.size(); ++i)
    if (g[i])
      _exact_derivs[i] = g[i]->clone().release();
}


void ExactErrorEstimator::attach_exact_deriv (unsigned int sys_num,
                                              FunctionBase<Gradient> * g)
{
  if (_exact_derivs.size() <= sys_num)
    _exact_derivs.resize(sys_num+1, libmesh_nullptr);

  if (g)
    _exact_derivs[sys_num] = g->clone().release();
}




void ExactErrorEstimator::attach_exact_hessian (Tensor hptr(const Point & p,
                                                            const Parameters & parameters,
                                                            const std::string & sys_name,
                                                            const std::string & unknown_name))
{
  libmesh_assert(hptr);
  _exact_hessian = hptr;

  // We're not using a fine grid solution
  _equation_systems_fine = libmesh_nullptr;

  // We're not using user-provided functors
  this->clear_functors();
}


void ExactErrorEstimator::attach_exact_hessians (std::vector<FunctionBase<Tensor> *> h)
{
  // Clear out any previous _exact_hessians entries, then add a new
  // entry for each system.
  for (std::size_t i=0; i != _exact_hessians.size(); ++i)
    delete (_exact_hessians[i]);

  _exact_hessians.clear();
  _exact_hessians.resize(h.size(), libmesh_nullptr);

  // We use clone() to get non-sliced copies of FunctionBase
  // subclasses, but we don't currently put the resulting UniquePtrs
  // into an STL container.
  for (std::size_t i=0; i != h.size(); ++i)
    if (h[i])
      _exact_hessians[i] = h[i]->clone().release();
}


void ExactErrorEstimator::attach_exact_hessian (unsigned int sys_num,
                                                FunctionBase<Tensor> * h)
{
  if (_exact_hessians.size() <= sys_num)
    _exact_hessians.resize(sys_num+1, libmesh_nullptr);

  if (h)
    _exact_hessians[sys_num] = h->clone().release();
}


void ExactErrorEstimator::attach_reference_solution (EquationSystems * es_fine)
{
  libmesh_assert(es_fine);
  _equation_systems_fine = es_fine;

  // If we're using a fine grid solution, we're not using exact value
  // function pointers or functors.
  _exact_value = libmesh_nullptr;
  _exact_deriv = libmesh_nullptr;
  _exact_hessian = libmesh_nullptr;

  this->clear_functors();
}

void ExactErrorEstimator::estimate_error (const System & system,
                                          ErrorVector & error_per_cell,
                                          const NumericVector<Number> * solution_vector,
                                          bool estimate_parent_error)
{
  // Ignore the fact that this variable is unused when !LIBMESH_ENABLE_AMR
  libmesh_ignore(estimate_parent_error);

  // The current mesh
  const MeshBase & mesh = system.get_mesh();

  // The dimensionality of the mesh
  const unsigned int dim = mesh.mesh_dimension();

  // The number of variables in the system
  const unsigned int n_vars = system.n_vars();

  // The DofMap for this system
  const DofMap & dof_map = system.get_dof_map();

  // Resize the error_per_cell vector to be
  // the number of elements, initialize it to 0.
  error_per_cell.resize (mesh.max_elem_id());
  std::fill (error_per_cell.begin(), error_per_cell.end(), 0.);

  // Prepare current_local_solution to localize a non-standard
  // solution vector if necessary
  if (solution_vector && solution_vector != system.solution.get())
    {
      NumericVector<Number> * newsol =
        const_cast<NumericVector<Number> *>(solution_vector);
      System & sys = const_cast<System &>(system);
      newsol->swap(*sys.solution);
      sys.update();
    }

  // Loop over all the variables in the system
  for (unsigned int var=0; var<n_vars; var++)
    {
      // Possibly skip this variable
      if (error_norm.weight(var) == 0.0) continue;

      // The (string) name of this variable
      const std::string & var_name = system.variable_name(var);

      // The type of finite element to use for this variable
      const FEType & fe_type = dof_map.variable_type (var);

      UniquePtr<FEBase> fe (FEBase::build (dim, fe_type));

      // Build an appropriate Gaussian quadrature rule
      UniquePtr<QBase> qrule =
        fe_type.default_quadrature_rule (dim,
                                         _extra_order);

      fe->attach_quadrature_rule (qrule.get());

      // Prepare a global solution and a MeshFunction of the fine system if we need one
      UniquePtr<MeshFunction> fine_values;
      UniquePtr<NumericVector<Number> > fine_soln = NumericVector<Number>::build(system.comm());
      if (_equation_systems_fine)
        {
          const System & fine_system = _equation_systems_fine->get_system(system.name());

          std::vector<Number> global_soln;
          // FIXME - we're assuming that the fine system solution gets
          // used even when a different vector is used for the coarse
          // system
          fine_system.update_global_solution(global_soln);
          fine_soln->init
            (cast_int<numeric_index_type>(global_soln.size()), true,
             SERIAL);
          (*fine_soln) = global_soln;

          fine_values = UniquePtr<MeshFunction>
            (new MeshFunction(*_equation_systems_fine,
                              *fine_soln,
                              fine_system.get_dof_map(),
                              fine_system.variable_number(var_name)));
          fine_values->init();
        } else {
        // Initialize functors if we're using them
        for (std::size_t i=0; i != _exact_values.size(); ++i)
          if (_exact_values[i])
            _exact_values[i]->init();

        for (std::size_t i=0; i != _exact_derivs.size(); ++i)
          if (_exact_derivs[i])
            _exact_derivs[i]->init();

        for (std::size_t i=0; i != _exact_hessians.size(); ++i)
          if (_exact_hessians[i])
            _exact_hessians[i]->init();
      }

      // Request the data we'll need to compute with
      fe->get_JxW();
      fe->get_phi();
      fe->get_dphi();
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
      fe->get_d2phi();
#endif
      fe->get_xyz();

#ifdef LIBMESH_ENABLE_AMR
      // If we compute on parent elements, we'll want to do so only
      // once on each, so we need to keep track of which we've done.
      std::vector<bool> computed_var_on_parent;

      if (estimate_parent_error)
        computed_var_on_parent.resize(error_per_cell.size(), false);
#endif

      // TODO: this ought to be threaded (and using subordinate
      // MeshFunction objects in each thread rather than a single
      // master)

      // Iterate over all the active elements in the mesh
      // that live on this processor.
      MeshBase::const_element_iterator
        elem_it  = mesh.active_local_elements_begin();
      const MeshBase::const_element_iterator
        elem_end = mesh.active_local_elements_end();

      for (;elem_it != elem_end; ++elem_it)
        {
          // e is necessarily an active element on the local processor
          const Elem * elem = *elem_it;
          const dof_id_type e_id = elem->id();

#ifdef LIBMESH_ENABLE_AMR
          // See if the parent of element e has been examined yet;
          // if not, we may want to compute the estimator on it
          const Elem * parent = elem->parent();

          // We only can compute and only need to compute on
          // parents with all active children
          bool compute_on_parent = true;
          if (!parent || !estimate_parent_error)
            compute_on_parent = false;
          else
            for (unsigned int c=0; c != parent->n_children(); ++c)
              if (!parent->child_ptr(c)->active())
                compute_on_parent = false;

          if (compute_on_parent &&
              !computed_var_on_parent[parent->id()])
            {
              computed_var_on_parent[parent->id()] = true;

              // Compute a projection onto the parent
              DenseVector<Number> Uparent;
              FEBase::coarsened_dof_values(*(system.current_local_solution),
                                           dof_map, parent, Uparent,
                                           var, false);

              error_per_cell[parent->id()] +=
                static_cast<ErrorVectorReal>
                (find_squared_element_error(system, var_name,
                                            parent, Uparent,
                                            fe.get(),
                                            fine_values.get()));
            }
#endif

          // Get the local to global degree of freedom maps
          std::vector<dof_id_type> dof_indices;
          dof_map.dof_indices (elem, dof_indices, var);
          const unsigned int n_dofs =
            cast_int<unsigned int>(dof_indices.size());
          DenseVector<Number> Uelem(n_dofs);
          for (unsigned int i=0; i != n_dofs; ++i)
            Uelem(i) = system.current_solution(dof_indices[i]);

          error_per_cell[e_id] +=
            static_cast<ErrorVectorReal>
            (find_squared_element_error(system, var_name, elem,
                                        Uelem, fe.get(),
                                        fine_values.get()));

        } // End loop over active local elements
    } // End loop over variables



  // Each processor has now computed the error contribuions
  // for its local elements.  We need to sum the vector
  // and then take the square-root of each component.  Note
  // that we only need to sum if we are running on multiple
  // processors, and we only need to take the square-root
  // if the value is nonzero.  There will in general be many
  // zeros for the inactive elements.

  // First sum the vector of estimated error values
  this->reduce_error(error_per_cell, system.comm());

  // Compute the square-root of each component.
  {
    LOG_SCOPE("std::sqrt()", "ExactErrorEstimator");
    for (std::size_t i=0; i<error_per_cell.size(); i++)
      if (error_per_cell[i] != 0.)
        {
          libmesh_assert_greater (error_per_cell[i], 0.);
          error_per_cell[i] = std::sqrt(error_per_cell[i]);
        }
  }

  // If we used a non-standard solution before, now is the time to fix
  // the current_local_solution
  if (solution_vector && solution_vector != system.solution.get())
    {
      NumericVector<Number> * newsol =
        const_cast<NumericVector<Number> *>(solution_vector);
      System & sys = const_cast<System &>(system);
      newsol->swap(*sys.solution);
      sys.update();
    }
}



Real ExactErrorEstimator::find_squared_element_error(const System & system,
                                                     const std::string & var_name,
                                                     const Elem * elem,
                                                     const DenseVector<Number> & Uelem,
                                                     FEBase * fe,
                                                     MeshFunction * fine_values) const
{
  // The (string) name of this system
  const std::string & sys_name = system.name();
  const unsigned int sys_num = system.number();

  const unsigned int var = system.variable_number(var_name);
  const unsigned int var_component =
    system.variable_scalar_number(var, 0);

  const Parameters & parameters = system.get_equation_systems().parameters;

  // reinitialize the element-specific data
  // for the current element
  fe->reinit (elem);

  // Get the data we need to compute with
  const std::vector<Real> &                       JxW          = fe->get_JxW();
  const std::vector<std::vector<Real> > &         phi_values   = fe->get_phi();
  const std::vector<std::vector<RealGradient> > & dphi_values  = fe->get_dphi();
  const std::vector<Point> &                      q_point      = fe->get_xyz();
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  const std::vector<std::vector<RealTensor> > &   d2phi_values = fe->get_d2phi();
#endif

  // The number of shape functions
  const unsigned int n_sf =
    cast_int<unsigned int>(Uelem.size());

  // The number of quadrature points
  const unsigned int n_qp =
    cast_int<unsigned int>(JxW.size());

  Real error_val = 0;

  // Begin the loop over the Quadrature points.
  //
  for (unsigned int qp=0; qp<n_qp; qp++)
    {
      // Real u_h = 0.;
      // RealGradient grad_u_h;

      Number u_h = 0.;

      Gradient grad_u_h;
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
      Tensor grad2_u_h;
#endif

      // Compute solution values at the current
      // quadrature point.  This reqiures a sum
      // over all the shape functions evaluated
      // at the quadrature point.
      for (unsigned int i=0; i<n_sf; i++)
        {
          // Values from current solution.
          u_h      += phi_values[i][qp]*Uelem(i);
          grad_u_h += dphi_values[i][qp]*Uelem(i);
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
          grad2_u_h += d2phi_values[i][qp]*Uelem(i);
#endif
        }

      // Compute the value of the error at this quadrature point
      if (error_norm.type(var) == L2 ||
          error_norm.type(var) == H1 ||
          error_norm.type(var) == H2)
        {
          Number val_error = u_h;
          if (_exact_value)
            val_error -= _exact_value(q_point[qp],parameters,sys_name,var_name);
          else if (_exact_values.size() > sys_num && _exact_values[sys_num])
            val_error -= _exact_values[sys_num]->
              component(var_component, q_point[qp], system.time);
          else if (_equation_systems_fine)
            val_error -= (*fine_values)(q_point[qp]);

          // Add the squares of the error to each contribution
          error_val += JxW[qp]*TensorTools::norm_sq(val_error);
        }

      // Compute the value of the error in the gradient at this
      // quadrature point
      if (error_norm.type(var) == H1 ||
          error_norm.type(var) == H1_SEMINORM ||
          error_norm.type(var) == H2)
        {
          Gradient grad_error = grad_u_h;
          if(_exact_deriv)
            grad_error -= _exact_deriv(q_point[qp],parameters,sys_name,var_name);
          else if (_exact_derivs.size() > sys_num && _exact_derivs[sys_num])
            grad_error -= _exact_derivs[sys_num]->
              component(var_component, q_point[qp], system.time);
          else if(_equation_systems_fine)
            grad_error -= fine_values->gradient(q_point[qp]);

          error_val += JxW[qp]*grad_error.norm_sq();
        }


#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
      // Compute the value of the error in the hessian at this
      // quadrature point
      if ((error_norm.type(var) == H2_SEMINORM ||
           error_norm.type(var) == H2))
        {
          Tensor grad2_error = grad2_u_h;
          if(_exact_hessian)
            grad2_error -= _exact_hessian(q_point[qp],parameters,sys_name,var_name);
          else if (_exact_hessians.size() > sys_num && _exact_hessians[sys_num])
            grad2_error -= _exact_hessians[sys_num]->
              component(var_component, q_point[qp], system.time);
          else if (_equation_systems_fine)
            grad2_error -= fine_values->hessian(q_point[qp]);

          error_val += JxW[qp]*grad2_error.norm_sq();
        }
#endif

    } // end qp loop

  libmesh_assert_greater_equal (error_val, 0.);

  return error_val;
}



void ExactErrorEstimator::clear_functors()
{
  // delete will clean up any cloned functors and no-op on any NULL
  // pointers

  for (std::size_t i=0; i != _exact_values.size(); ++i)
    delete (_exact_values[i]);

  for (std::size_t i=0; i != _exact_derivs.size(); ++i)
    delete (_exact_derivs[i]);

  for (std::size_t i=0; i != _exact_hessians.size(); ++i)
    delete (_exact_hessians[i]);
}



} // namespace libMesh
