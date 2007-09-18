// $Id: exact_error_estimator.C,v 1.10 2007-09-18 22:13:40 roystgnr Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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
#include <cmath>    // for sqrt


// Local Includes
#include "libmesh_common.h"
#include "exact_error_estimator.h"
#include "dof_map.h"
#include "equation_systems.h"
#include "error_vector.h"
#include "fe.h"
#include "fe_interface.h"
#include "libmesh_logging.h"
#include "elem.h"
#include "mesh.h"
#include "mesh_function.h"
#include "numeric_vector.h"
#include "quadrature.h"
#include "system.h"

//-----------------------------------------------------------------
// ErrorEstimator implementations
void ExactErrorEstimator::attach_exact_value (Number fptr(const Point& p,
                                                          const Parameters& parameters,
                                                          const std::string& sys_name,
                                                          const std::string& unknown_name))
{
  assert (fptr != NULL);
  _exact_value = fptr;

  // If we're not using a fine grid solution
  _equation_systems_fine = NULL;

}


void ExactErrorEstimator::attach_exact_deriv (Gradient fptr(const Point& p,
                                                            const Parameters& parameters,
                                                            const std::string& sys_name,
                                                            const std::string& unknown_name))
{
  assert (fptr != NULL);
  _exact_deriv = fptr;

  // If we're not using a fine grid solution
  _equation_systems_fine = NULL;

}


void ExactErrorEstimator::attach_exact_hessian (Tensor fptr(const Point& p,
                                                            const Parameters& parameters,
                                                            const std::string& sys_name,
                                                            const std::string& unknown_name))
{
  assert (fptr != NULL);
  _exact_hessian = fptr;

  // If we're not using a fine grid solution
  _equation_systems_fine = NULL;
}

void ExactErrorEstimator::attach_reference_solution (EquationSystems* es_fine)
{
  assert (es_fine != NULL);
  _equation_systems_fine = es_fine;

  // If we're using a fine grid solution, we're not using exact values
  _exact_value = NULL;
  _exact_deriv = NULL;
  _exact_hessian = NULL;
}

void ExactErrorEstimator::estimate_error (const System& system,
					  ErrorVector& error_per_cell,
					  bool estimate_parent_error)
{
  // The current mesh
  const Mesh& mesh = system.get_mesh();

  // The dimensionality of the mesh
  const unsigned int dim = mesh.mesh_dimension();
  
  // The number of variables in the system
  const unsigned int n_vars = system.n_vars();
  
  // The DofMap for this system
  const DofMap& dof_map = system.get_dof_map();

  // Resize the error_per_cell vector to be
  // the number of elements, initialize it to 0.
  error_per_cell.resize (mesh.n_elem());
  std::fill (error_per_cell.begin(), error_per_cell.end(), 0.);

  // Check for the use of component_mask
  this->convert_component_mask_to_scale();
  
  // Check for a valid component_scale
  if (!component_scale.empty())
    {
      if (component_scale.size() != n_vars)
	{
	  std::cerr << "ERROR: component_scale is the wrong size:"
		    << std::endl
		    << " component_scale.scale()=" << component_scale.size()
		    << std::endl
		    << " n_vars=" << n_vars
		    << std::endl;
	  error();
	}
    }
  else
    {
      // No specified scaling.  Scale all variables by one.
      component_scale.resize (n_vars);
      std::fill (component_scale.begin(), component_scale.end(), 1.0);
    }

  // Loop over all the variables in the system
  for (unsigned int var=0; var<n_vars; var++)
    {
      // Possibly skip this variable
      if (!component_scale.empty())
	if (component_scale[var] == 0.0) continue;

      // The (string) name of this variable
      const std::string& var_name = system.variable_name(var);
      
      // The type of finite element to use for this variable
      const FEType& fe_type = dof_map.variable_type (var);
      
      AutoPtr<FEBase> fe (FEBase::build (dim, fe_type));

      // Build an appropriate Gaussian quadrature rule
      AutoPtr<QBase> qrule =
        fe_type.default_quadrature_rule (mesh.mesh_dimension(),
                                         _extra_order);

      fe->attach_quadrature_rule (qrule.get());

      // Prepare a global solution and a MeshFunction of the fine system if we need one
      AutoPtr<MeshFunction> fine_values;
      AutoPtr<NumericVector<Number> > fine_soln = NumericVector<Number>::build();
      if (_equation_systems_fine)
      {
	const System& fine_system = _equation_systems_fine->get_system(system.name());

	std::vector<Number> global_soln;
	fine_system.update_global_solution(global_soln);
	fine_soln->init(fine_system.solution->size(),fine_system.solution->size());
	(*fine_soln) = global_soln;

	fine_values = AutoPtr<MeshFunction>
	  (new MeshFunction(*_equation_systems_fine,
			    *fine_soln,
			    fine_system.get_dof_map(),
			    fine_system.variable_number(var_name)));
	fine_values->init();
      }
      
      // Request the data we'll need to compute with
      fe->get_JxW();
      fe->get_phi();
      fe->get_dphi();
#ifdef ENABLE_SECOND_DERIVATIVES
      fe->get_d2phi();
#endif
      fe->get_xyz();

      // The global DOF indices for elements e & f
      std::vector<unsigned int> dof_indices;

      // Iterate over all the active elements in the mesh
      // that live on this processor.
      MeshBase::const_element_iterator
        elem_it  = mesh.active_local_elements_begin();
      const MeshBase::const_element_iterator
        elem_end = mesh.active_local_elements_end(); 

	{
	  // e is necessarily an active element on the local processor
	  const Elem* elem = *elem_it;
	  const unsigned int e_id = elem->id();

#ifdef ENABLE_AMR
          // See if the parent of element e has been examined yet;
          // if not, we may want to compute the estimator on it
          const Elem* parent = elem->parent();

          // We only can compute and only need to compute on
          // parents with all active children
          bool compute_on_parent = true;
          if (!parent || !estimate_parent_error)
            compute_on_parent = false;
          else
            for (unsigned int c=0; c != parent->n_children(); ++c)
              if (!parent->child(c)->active())
                compute_on_parent = false;

          if (compute_on_parent &&
              !error_per_cell[parent->id()])
            {
              // Compute a projection onto the parent
              DenseVector<Number> Uparent;
              FEBase::coarsened_dof_values(*(system.solution),
                                           dof_map, parent, Uparent,
                                           var, false);

              error_per_cell[parent->id()] =
		find_squared_element_error(system, var_name, parent, Uparent,
                                           fe.get(), fine_values.get());
            }
#endif

          // Get the local to global degree of freedom maps
          dof_map.dof_indices (elem, dof_indices, var);
          DenseVector<Number> Uelem(dof_indices.size());
          for (unsigned int i=0; i != dof_indices.size(); ++i)
            Uelem(i) = system.current_solution(dof_indices[i]);

          error_per_cell[e_id] =
            find_squared_element_error(system, var_name, elem, Uelem,
                                       fe.get(), fine_values.get());

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
  this->reduce_error(error_per_cell);

  // Compute the square-root of each component.
  START_LOG("std::sqrt()", "ExactErrorEstimator");
  for (unsigned int i=0; i<error_per_cell.size(); i++)
    {
      
      if (error_per_cell[i] != 0.)
	{
	  assert (error_per_cell[i] > 0.);
	  error_per_cell[i] = std::sqrt(error_per_cell[i]);
	}
      
      
    }
  STOP_LOG("std::sqrt()", "ExactErrorEstimator");

  //  STOP_LOG("flux_jumps()", "KellyErrorEstimator");
}



Real ExactErrorEstimator::find_squared_element_error(const System& system,
                                                     const std::string& var_name,
                                                     const Elem *elem,
                                                     const DenseVector<Number> &Uelem,
                                                     FEBase *fe,
                                                     MeshFunction *fine_values) const
{
  // The (string) name of this system
  const std::string& sys_name = system.name();
  
  const Parameters& parameters = system.get_equation_systems().parameters;

  // reinitialize the element-specific data
  // for the current element
  fe->reinit (elem);

  // Get the data we need to compute with
  const std::vector<Real> &                      JxW          = fe->get_JxW();
  const std::vector<std::vector<Real> >&         phi_values   = fe->get_phi();
  const std::vector<std::vector<RealGradient> >& dphi_values  = fe->get_dphi();
  const std::vector<Point>&                      q_point      = fe->get_xyz();
#ifdef ENABLE_SECOND_DERIVATIVES
  const std::vector<std::vector<RealTensor> >&   d2phi_values = fe->get_d2phi();
#endif

  // The number of shape functions
  const unsigned int n_sf = Uelem.size();

  // The number of quadrature points
  const unsigned int n_qp = JxW.size();

  double L2normsq = 0., H1seminormsq = 0., H2seminormsq = 0.;
  
  // Begin the loop over the Quadrature points.
  //
  for (unsigned int qp=0; qp<n_qp; qp++)
    {
      // Real u_h = 0.;
      // RealGradient grad_u_h;

      Number u_h = 0.;

#ifndef USE_COMPLEX_NUMBERS
      RealGradient grad_u_h;
#else
      // Gradient     grad_u_h;
      RealGradient grad_u_h_re;
      RealGradient grad_u_h_im;
#endif
#ifdef ENABLE_SECOND_DERIVATIVES
  #ifndef USE_COMPLEX_NUMBERS
      RealTensor grad2_u_h;
  #else
      RealTensor grad2_u_h_re;
      RealTensor grad2_u_h_im;
  #endif
#endif

      // Compute solution values at the current
      // quadrature point.  This reqiures a sum
      // over all the shape functions evaluated
      // at the quadrature point.
      for (unsigned int i=0; i<n_sf; i++)
        {
          // Values from current solution.
          u_h      += phi_values[i][qp]*Uelem(i);
#ifndef USE_COMPLEX_NUMBERS
          grad_u_h += dphi_values[i][qp]*Uelem(i);
#else
          grad_u_h_re += dphi_values[i][qp]*Uelem(i).real();
          grad_u_h_im += dphi_values[i][qp]*Uelem(i).imag();
#endif
#ifdef ENABLE_SECOND_DERIVATIVES
  #ifndef USE_COMPLEX_NUMBERS
          grad2_u_h += d2phi_values[i][qp]*Uelem(i);
  #else
          grad2_u_h_re += d2phi_values[i][qp]*Uelem(i).real();
          grad2_u_h_im += d2phi_values[i][qp]*Uelem(i).imag();
  #endif
#endif
        }

#ifdef USE_COMPLEX_NUMBERS
      Gradient grad_u_h (grad_u_h_re, grad_u_h_im);
  #ifdef ENABLE_SECOND_DERIVATIVES
      Tensor grad2_u_h (grad2_u_h_re, grad2_u_h_im);
  #endif
#endif

      // Compute the value of the error at this quadrature point
      Number val_error = 0;
      if(_exact_value)
	val_error = u_h - _exact_value(q_point[qp],parameters,sys_name,var_name);
      else if(_equation_systems_fine)
	val_error = u_h - (*fine_values)(q_point[qp]);

      // Add the squares of the error to each contribution
#ifndef USE_COMPLEX_NUMBERS
      L2normsq += JxW[qp]*(val_error*val_error);
#else
      L2normsq += JxW[qp]*std::norm(val_error);
#endif

      // Compute the value of the error in the gradient at this
      // quadrature point
      if ((_exact_deriv || _equation_systems_fine) && _sobolev_order > 0)
        {
          Gradient grad_error;
	  if(_exact_deriv)
	    grad_error = grad_u_h - _exact_deriv(q_point[qp],parameters,sys_name,var_name);
	  else if(_equation_systems_fine)
	    grad_error = grad_u_h - fine_values->gradient(q_point[qp]);

#ifndef USE_COMPLEX_NUMBERS
          H1seminormsq += JxW[qp]*(grad_error*grad_error);
#else
          H1seminormsq += JxW[qp]*std::abs(grad_error*grad_error);
#endif
        }


#ifdef ENABLE_SECOND_DERIVATIVES
      // Compute the value of the error in the hessian at this
      // quadrature point
      if ((_exact_hessian || _equation_systems_fine) && _sobolev_order > 1)
        {
	  Tensor grad2_error;
	  if(_exact_hessian)
	    grad2_error = grad2_u_h - _exact_hessian(q_point[qp],parameters,sys_name,var_name);
	  else if (_equation_systems_fine)
	    grad2_error = grad2_u_h - fine_values->hessian(q_point[qp]);

          H2seminormsq += JxW[qp]*std::abs(grad2_error.contract(grad2_error));
        }
#endif

    } // end qp loop

  assert (L2normsq     >= 0.);
  assert (H1seminormsq >= 0.);
	  
  Real error_val = L2normsq;
  if (_sobolev_order > 0)
    error_val += H1seminormsq;
  if (_sobolev_order > 1)
    error_val += H2seminormsq;
}
