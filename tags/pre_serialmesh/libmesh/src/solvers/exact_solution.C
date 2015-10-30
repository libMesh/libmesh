// $Id: exact_solution.C,v 1.31 2007-06-15 22:34:34 roystgnr Exp $

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


// Local includes
#include "dof_map.h"
#include "elem.h"
#include "exact_solution.h"
#include "equation_systems.h"
#include "fe.h"
#include "fe_interface.h"
#include "mesh.h"
#include "mesh_function.h"
#include "numeric_vector.h"
#include "parallel.h"
#include "quadrature.h"
#include "tensor_value.h"
#include "vector_value.h"

ExactSolution::ExactSolution(EquationSystems& es) :
  _exact_value (NULL),
  _exact_deriv (NULL),
  _exact_hessian (NULL),
  _equation_systems(es),
  _equation_systems_fine(NULL),
  _extra_order(0)
{
  // Initialize the _errors data structure which holds all
  // the eventual values of the error.
  for (unsigned int sys=0; sys<_equation_systems.n_systems(); ++sys)
    {
      // Reference to the system
      const System& system = _equation_systems.get_system(sys);
      
      // The name of the system
      const std::string& sys_name = system.name();

      // The SystemErrorMap to be inserted
      ExactSolution::SystemErrorMap sem;

      for (unsigned int var=0; var<system.n_vars(); ++var)
	{
	  // The name of this variable
	  const std::string& var_name = system.variable_name(var);
	  sem[var_name] = std::vector<Number>(3, 0.);
	}
      
      _errors[sys_name] = sem;
    }
}


void ExactSolution::attach_reference_solution (EquationSystems* es_fine)
{
  assert (es_fine != NULL);
  _equation_systems_fine = es_fine;

  // If we're using a fine grid solution, we're not using exact values
  _exact_value = NULL;
  _exact_deriv = NULL;
  _exact_hessian = NULL;
}


void ExactSolution::attach_exact_value (Number fptr(const Point& p,
						    const Parameters& parameters,
						    const std::string& sys_name,
						    const std::string& unknown_name))
{
  assert (fptr != NULL);
  _exact_value = fptr;

  // If we're using exact values, we're not using a fine grid solution
  _equation_systems_fine = NULL;
}


void ExactSolution::attach_exact_deriv (Gradient fptr(const Point& p,
						      const Parameters& parameters,
						      const std::string& sys_name,
						      const std::string& unknown_name))
{
  assert (fptr != NULL);
  _exact_deriv = fptr;

  // If we're using exact values, we're not using a fine grid solution
  _equation_systems_fine = NULL;
}


void ExactSolution::attach_exact_hessian (Tensor fptr(const Point& p,
						      const Parameters& parameters,
						      const std::string& sys_name,
						      const std::string& unknown_name))
{
  assert (fptr != NULL);
  _exact_hessian = fptr;

  // If we're using exact values, we're not using a fine grid solution
  _equation_systems_fine = NULL;
}




std::vector<Number>& ExactSolution::_check_inputs(const std::string& sys_name,
						  const std::string& unknown_name)
{
  // If no exact solution function or fine grid solution has been
  // attached, we now just compute the solution norm (i.e. the
  // difference from an "exact solution" of zero
/*
  if (_exact_value == NULL)
    {
      std::cerr << "Cannot compute error, you must provide a "
		<< "function which computes the exact solution."
		<< std::endl;
      error();
    }
*/
  
  // Make sure the requested sys_name exists.
  std::map<std::string, SystemErrorMap>::iterator sys_iter =
    _errors.find(sys_name);
  
  if (sys_iter == _errors.end())
    {
      std::cerr << "Sorry, couldn't find the requested system '"
                << sys_name << "'."
		<< std::endl;
      error();
    }
  
  // Make sure the requested unknown_name exists.
  SystemErrorMap::iterator var_iter = (*sys_iter).second.find(unknown_name);

  if (var_iter == (*sys_iter).second.end())
    {
      std::cerr << "Sorry, couldn't find the requested variable '"
                << unknown_name << "'."
		<< std::endl;
      error();
    }

  // Return a reference to the proper error entry
  return (*var_iter).second;
}






void ExactSolution::compute_error(const std::string& sys_name,
				  const std::string& unknown_name)
{
  // Check the inputs for validity, and get a reference
  // to the proper location to store the error
  std::vector<Number>& error_vals = this->_check_inputs(sys_name,
							unknown_name);
  this->_compute_error(sys_name,
		       unknown_name,
		       error_vals);
}





Number ExactSolution::l2_error(const std::string& sys_name,
			       const std::string& unknown_name)
{
  
  // Check the inputs for validity, and get a reference
  // to the proper location to store the error
  std::vector<Number>& error_vals = this->_check_inputs(sys_name,
							unknown_name);
  
  // Return the square root of the first component of the
  // computed error.
  return std::sqrt(error_vals[0]);
}







Number ExactSolution::h1_error(const std::string& sys_name,
			       const std::string& unknown_name)
{
  // If the user has supplied no exact derivative function, we
  // just integrate the H1 norm of the solution; i.e. its
  // difference from an "exact solution" of zero.
/*
  if (_exact_deriv == NULL)
    {
      std::cerr << "Cannot compute H1 error, you must provide a "
		<< "function which computes the gradient of the exact solution."
		<< std::endl;
      error();
    }
*/
  
  // Check the inputs for validity, and get a reference
  // to the proper location to store the error
  std::vector<Number>& error_vals = this->_check_inputs(sys_name,
							unknown_name);
  
  // Return the square root of the sum of the computed errors.
  return std::sqrt(error_vals[0] + error_vals[1]);
}







Number ExactSolution::h2_error(const std::string& sys_name,
			       const std::string& unknown_name)
{
  // If the user has supplied no exact derivative functions, we
  // just integrate the H1 norm of the solution; i.e. its
  // difference from an "exact solution" of zero.
/*
  if (_exact_deriv == NULL || _exact_hessian == NULL)
    {
      std::cerr << "Cannot compute H2 error, you must provide functions "
		<< "which computes the gradient and hessian of the "
                << "exact solution."
		<< std::endl;
      error();
    }
*/
  
  // Check the inputs for validity, and get a reference
  // to the proper location to store the error
  std::vector<Number>& error_vals = this->_check_inputs(sys_name,
						      unknown_name);
  
  // Return the square root of the sum of the computed errors.
  return std::sqrt(error_vals[0] + error_vals[1] + error_vals[2]);
}








void ExactSolution::_compute_error(const std::string& sys_name,
				   const std::string& unknown_name,
				   std::vector<Number>& error_vals)
{
  // Make sure we aren't "overconfigured"
  assert (!(_exact_value && _equation_systems_fine));

  // Get a reference to the system whose error is being computed.
  // If we have a fine grid, however, we'll integrate on that instead
  // for more accuracy.
  const System& computed_system = _equation_systems_fine ?
    _equation_systems_fine->get_system(sys_name) :
    _equation_systems.get_system (sys_name);

  // Prepare a global solution and a MeshFunction of the coarse system if we need one
  AutoPtr<MeshFunction> coarse_values;
  AutoPtr<NumericVector<Number> > comparison_soln = NumericVector<Number>::build();
  if (_equation_systems_fine)
    {
      const System& comparison_system
	= _equation_systems.get_system(sys_name);

      std::vector<Number> global_soln;
      comparison_system.update_global_solution(global_soln);
      comparison_soln->init(comparison_system.solution->size(),
                            comparison_system.solution->size());
      (*comparison_soln) = global_soln;

      coarse_values = AutoPtr<MeshFunction>
	(new MeshFunction(_equation_systems,
			  *comparison_soln,
			  comparison_system.get_dof_map(),
			  comparison_system.variable_number(unknown_name)));
      coarse_values->init();
    }

  // Get a reference to the dofmap and mesh for that system
  const DofMap& computed_dof_map = computed_system.get_dof_map();

  const MeshBase& _mesh = computed_system.get_mesh();

  // Zero the error before summation
  error_vals = std::vector<Number>(3, 0.);

  // get the EquationSystems parameters for use with _exact_value
  const Parameters& parameters = this->_equation_systems.parameters;

  // Get the current time, in case the exact solution depends on it.
  // Steady systems of equations do not have a time parameter, so this
  // routine needs to take that into account.
  // FIXME!!!
  // const Real time = 0.;//_equation_systems.parameter("time");  

  // Construct Quadrature rule based on default quadrature order
  const unsigned int var = computed_system.variable_number(unknown_name);
  const FEType& fe_type  = computed_dof_map.variable_type(var);

  AutoPtr<QBase> qrule =
    fe_type.default_quadrature_rule (_mesh.mesh_dimension(),
                                     _extra_order);

  // Construct finite element object
  
  AutoPtr<FEBase> fe(FEBase::build(_mesh.mesh_dimension(), fe_type));

  // Attach quadrature rule to FE object
  fe->attach_quadrature_rule (qrule.get());
  
  // The Jacobian*weight at the quadrature points.
  const std::vector<Real>& JxW                               = fe->get_JxW();
  
  // The value of the shape functions at the quadrature points
  // i.e. phi(i) = phi_values[i][qp] 
  const std::vector<std::vector<Real> >&  phi_values         = fe->get_phi();
  
  // The value of the shape function gradients at the quadrature points
  const std::vector<std::vector<RealGradient> >& dphi_values = fe->get_dphi();

#ifdef ENABLE_SECOND_DERIVATIVES
  // The value of the shape function second derivatives at the quadrature points
  const std::vector<std::vector<RealTensor> >& d2phi_values = fe->get_d2phi();
#endif
  
  // The XYZ locations (in physical space) of the quadrature points
  const std::vector<Point>& q_point                          = fe->get_xyz();
	    
  // The global degree of freedom indices associated
  // with the local degrees of freedom.
  std::vector<unsigned int> dof_indices;


  //
  // Begin the loop over the elements
  //
  MeshBase::const_element_iterator       el     = _mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = _mesh.active_local_elements_end(); 

  for ( ; el != end_el; ++el)
    {
      // Store a pointer to the element we are currently
      // working on.  This allows for nicer syntax later.
      const Elem* elem = *el;

      // reinitialize the element-specific data
      // for the current element
      fe->reinit (elem);

      // Get the local to global degree of freedom maps
      computed_dof_map.dof_indices    (elem, dof_indices, var);
      
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
	      u_h      += phi_values[i][qp]*computed_system.current_solution  (dof_indices[i]);
#ifndef USE_COMPLEX_NUMBERS
	      grad_u_h += dphi_values[i][qp]*computed_system.current_solution (dof_indices[i]);
#else
	      grad_u_h_re += dphi_values[i][qp]*computed_system.current_solution (dof_indices[i]).real();
	      grad_u_h_im += dphi_values[i][qp]*computed_system.current_solution (dof_indices[i]).imag();
#endif
#ifdef ENABLE_SECOND_DERIVATIVES
  #ifndef USE_COMPLEX_NUMBERS
	      grad2_u_h += d2phi_values[i][qp]*computed_system.current_solution (dof_indices[i]);
  #else
	      grad2_u_h_re += d2phi_values[i][qp]*computed_system.current_solution (dof_indices[i]).real();
	      grad2_u_h_im += d2phi_values[i][qp]*computed_system.current_solution (dof_indices[i]).imag();
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
	  Number exact_val = 0.0;
          if (_exact_value)
	    exact_val = _exact_value(q_point[qp], parameters,
				     sys_name, unknown_name);
	  else if (_equation_systems_fine)
	    exact_val = (*coarse_values)(q_point[qp]);

	  const Number val_error = u_h - exact_val;

	  // Add the squares of the error to each contribution
	  error_vals[0] += JxW[qp]*(val_error*val_error);


	  // Compute the value of the error in the gradient at this
	  // quadrature point
          Gradient exact_grad;
	  if (_exact_deriv)
	    exact_grad = _exact_deriv(q_point[qp], parameters,
				      sys_name, unknown_name);
	  else if (_equation_systems_fine)
	    exact_grad = coarse_values->gradient(q_point[qp]);

	  const Gradient grad_error = grad_u_h - exact_grad;

	  error_vals[1] += JxW[qp]*(grad_error*grad_error);


#ifdef ENABLE_SECOND_DERIVATIVES
	  // Compute the value of the error in the hessian at this
	  // quadrature point
          Tensor exact_hess;
	  if (_exact_hessian)
	    exact_hess = _exact_hessian(q_point[qp], parameters,
					sys_name, unknown_name);
	  else if (_equation_systems_fine)
	    exact_hess = coarse_values->hessian(q_point[qp]);

	  const Tensor grad2_error = grad2_u_h - exact_hess;

	  error_vals[2] += JxW[qp]*(grad2_error.contract(grad2_error));
#endif
	  
	} // end qp loop
    } // end element loop

  // Add up the error values on all processors
  Parallel::sum(error_vals);
}
