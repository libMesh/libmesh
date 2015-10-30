// $Id: exact_error_estimator.C,v 1.3 2006-06-01 22:35:21 roystgnr Exp $

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
#include "fe.h"
#include "fe_interface.h"
#include "libmesh_logging.h"
#include "elem.h"
#include "mesh.h"
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
}


void ExactErrorEstimator::attach_exact_deriv (Gradient fptr(const Point& p,
                                                            const Parameters& parameters,
                                                            const std::string& sys_name,
                                                            const std::string& unknown_name))
{
  assert (fptr != NULL);
  _exact_deriv = fptr;
}


void ExactErrorEstimator::attach_exact_hessian (Tensor fptr(const Point& p,
                                                            const Parameters& parameters,
                                                            const std::string& sys_name,
                                                            const std::string& unknown_name))
{
  assert (fptr != NULL);
  _exact_hessian = fptr;
}


void ExactErrorEstimator::estimate_error (const System& system,
					  std::vector<float>& error_per_cell)
{
  // The current mesh
  const Mesh& mesh = system.get_mesh();

  // The dimensionality of the mesh
  const unsigned int dim = mesh.mesh_dimension();
  
  // The number of variables in the system
  const unsigned int n_vars = system.n_vars();
  
  // The DofMap for this system
  const DofMap& dof_map = system.get_dof_map();

  const Parameters& parameters = system.get_equation_systems().parameters;

  // Resize the error_per_cell vector to be
  // the number of elements, initialize it to 0.
  error_per_cell.resize (mesh.n_elem());
  std::fill (error_per_cell.begin(), error_per_cell.end(), 0.);

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

      // The (string) name of this system
      const std::string& sys_name = system.name();
      
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
      
      const std::vector<Real> &
        JxW          = fe->get_JxW();
      const std::vector<std::vector<Real> >&
        phi_values   = fe->get_phi();
      const std::vector<std::vector<RealGradient> >&
        dphi_values  = fe->get_dphi();
#ifdef ENABLE_SECOND_DERIVATIVES
      const std::vector<std::vector<RealTensor> >&
        d2phi_values = fe->get_d2phi();
#endif
      const std::vector<Point>& q_point = fe->get_xyz();

      // The global DOF indices for elements e & f
      std::vector<unsigned int> dof_indices;

      // Iterate over all the active elements in the mesh
      // that live on this processor.
      MeshBase::const_element_iterator
        elem_it  = mesh.active_local_elements_begin();
      const MeshBase::const_element_iterator
        elem_end = mesh.active_local_elements_end(); 

      for (; elem_it != elem_end; ++elem_it)
	{
	  // e is necessarily an active element on the local processor
	  const Elem* elem = *elem_it;
	  const unsigned int e_id = elem->id();
          double L2normsq = 0., H1seminormsq = 0., H2seminormsq = 0.;

          // reinitialize the element-specific data
          // for the current element
          fe->reinit (elem);

          // Get the local to global degree of freedom maps
          dof_map.dof_indices (elem, dof_indices, var);

          // The number of quadrature points
          const unsigned int n_qp = qrule->n_points();

          // The number of shape functions
          const unsigned int n_sf = dof_indices.size();

          //
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
                  u_h      += phi_values[i][qp]*system.current_solution (dof_indices[i]);
#ifndef USE_COMPLEX_NUMBERS
                  grad_u_h += dphi_values[i][qp]*system.current_solution (dof_indices[i]);
#else
                  grad_u_h_re += dphi_values[i][qp]*system.current_solution (dof_indices[i]).real();
                  grad_u_h_im += dphi_values[i][qp]*system.current_solution (dof_indices[i]).imag();
#endif
#ifdef ENABLE_SECOND_DERIVATIVES
  #ifndef USE_COMPLEX_NUMBERS
                  grad2_u_h += d2phi_values[i][qp]*system.current_solution (dof_indices[i]);
  #else
                  grad2_u_h_re += d2phi_values[i][qp]*system.current_solution (dof_indices[i]).real();
                  grad2_u_h_im += d2phi_values[i][qp]*system.current_solution (dof_indices[i]).imag();
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
              const Number val_error = (u_h - _exact_value(q_point[qp],
                                                           parameters,
                                                           sys_name,
                                                           var_name));

              // Add the squares of the error to each contribution
#ifndef USE_COMPLEX_NUMBERS
              L2normsq += JxW[qp]*(val_error*val_error);
#else
	      L2normsq += JxW[qp]*std::norm(val_error);
#endif

	      // Compute the value of the error in the gradient at this
              // quadrature point
              if (_exact_deriv != NULL && _sobolev_order > 0)
                {
                  Gradient grad_error = (grad_u_h - _exact_deriv(q_point[qp],
                                                                 parameters,
                                                                 sys_name,
                                                                 var_name));

#ifndef USE_COMPLEX_NUMBERS
                  H1seminormsq += JxW[qp]*(grad_error*grad_error);
#else
                  H1seminormsq += JxW[qp]*std::abs(grad_error*grad_error);
#endif
                }


#ifdef ENABLE_SECOND_DERIVATIVES
	      // Compute the value of the error in the hessian at this
              // quadrature point
              if (_exact_hessian != NULL && _sobolev_order > 1)
                {
                  Tensor grad2_error = (grad2_u_h - _exact_hessian(q_point[qp],
                                                                   parameters,
                                                                   sys_name,
                                                                   var_name));

                  H2seminormsq += JxW[qp]*(grad2_error.contract(grad2_error));
                }
#endif

            } // end qp loop

          error_per_cell[e_id] = L2normsq;
          if (_sobolev_order > 0)
            error_per_cell[e_id] += H1seminormsq;
          if (_sobolev_order > 1)
            error_per_cell[e_id] += H2seminormsq;

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
    if (error_per_cell[i] != 0.)
      error_per_cell[i] = std::sqrt(error_per_cell[i]);
  STOP_LOG("std::sqrt()", "ExactErrorEstimator");

  //  STOP_LOG("flux_jumps()", "KellyErrorEstimator");
}
