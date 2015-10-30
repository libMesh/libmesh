// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "fe_base.h"
#include "function_base.h"
#include "mesh_base.h"
#include "mesh_function.h"
#include "numeric_vector.h"
#include "parallel.h"
#include "quadrature.h"
#include "tensor_value.h"
#include "vector_value.h"
#include "wrapped_function.h"

namespace libMesh
{

ExactSolution::ExactSolution(const EquationSystems& es) :
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
	  sem[var_name] = std::vector<Real>(5, 0.);
	}

      _errors[sys_name] = sem;
    }
}


ExactSolution::~ExactSolution()
{
  // delete will clean up any cloned functors and no-op on any NULL
  // pointers

  for (unsigned int i=0; i != _exact_values.size(); ++i)
    delete (_exact_values[i]);

  for (unsigned int i=0; i != _exact_derivs.size(); ++i)
    delete (_exact_derivs[i]);

  for (unsigned int i=0; i != _exact_hessians.size(); ++i)
    delete (_exact_hessians[i]);
}


void ExactSolution::attach_reference_solution (const EquationSystems* es_fine)
{
  libmesh_assert (es_fine != NULL);
  _equation_systems_fine = es_fine;

  // If we're using a fine grid solution, we're not using exact values
  _exact_values.clear();
  _exact_derivs.clear();
  _exact_hessians.clear();
}


void ExactSolution::attach_exact_value (Number fptr(const Point& p,
						    const Parameters& parameters,
						    const std::string& sys_name,
						    const std::string& unknown_name))
{
  libmesh_assert (fptr != NULL);

  // Clear out any previous _exact_values entries, then add a new
  // entry for each system.
  _exact_values.clear();

  for (unsigned int sys=0; sys<_equation_systems.n_systems(); ++sys)
    {
      const System& system = _equation_systems.get_system(sys);
      _exact_values.push_back
	(new WrappedFunction<Number>
          (system, fptr, &_equation_systems.parameters));
    }

  // If we're using exact values, we're not using a fine grid solution
  _equation_systems_fine = NULL;
}


void ExactSolution::attach_exact_values (std::vector<FunctionBase<Number> *> f)
{
  // Clear out any previous _exact_values entries, then add a new
  // entry for each system.
  for (unsigned int i=0; i != _exact_values.size(); ++i)
    delete (_exact_values[i]);

  _exact_values.clear();
  _exact_values.resize(f.size(), NULL);

  // We use clone() to get non-sliced copies of FunctionBase
  // subclasses, but we can't put the resulting AutoPtrs into an STL
  // container.
  for (unsigned int i=0; i != f.size(); ++i)
    if (f[i])
      _exact_values[i] = f[i]->clone().release();
}


void ExactSolution::attach_exact_value (unsigned int sys_num,
                                        FunctionBase<Number> * f)
{
  if (_exact_values.size() <= sys_num)
  _exact_values.resize(sys_num+1, NULL);

  if (f)
    _exact_values[sys_num] = f->clone().release();
}


void ExactSolution::attach_exact_deriv (Gradient gptr(const Point& p,
						      const Parameters& parameters,
						      const std::string& sys_name,
						      const std::string& unknown_name))
{
  libmesh_assert (gptr != NULL);

  // Clear out any previous _exact_derivs entries, then add a new
  // entry for each system.
  _exact_derivs.clear();

  for (unsigned int sys=0; sys<_equation_systems.n_systems(); ++sys)
    {
      const System& system = _equation_systems.get_system(sys);
      _exact_derivs.push_back
	(new WrappedFunction<Gradient>
          (system, gptr, &_equation_systems.parameters));
    }

  // If we're using exact values, we're not using a fine grid solution
  _equation_systems_fine = NULL;
}


void ExactSolution::attach_exact_derivs (std::vector<FunctionBase<Gradient> *> g)
{
  // Clear out any previous _exact_derivs entries, then add a new
  // entry for each system.
  for (unsigned int i=0; i != _exact_derivs.size(); ++i)
    delete (_exact_derivs[i]);

  _exact_derivs.clear();
  _exact_derivs.resize(g.size(), NULL);

  // We use clone() to get non-sliced copies of FunctionBase
  // subclasses, but we can't put the resulting AutoPtrs into an STL
  // container.
  for (unsigned int i=0; i != g.size(); ++i)
    if (g[i])
      _exact_derivs[i] = g[i]->clone().release();
}


void ExactSolution::attach_exact_deriv (unsigned int sys_num,
                                        FunctionBase<Gradient>* g)
{
  if (_exact_derivs.size() <= sys_num)
  _exact_derivs.resize(sys_num+1, NULL);

  if (g)
    _exact_derivs[sys_num] = g->clone().release();
}


void ExactSolution::attach_exact_hessian (Tensor hptr(const Point& p,
						      const Parameters& parameters,
						      const std::string& sys_name,
						      const std::string& unknown_name))
{
  libmesh_assert (hptr != NULL);

  // Clear out any previous _exact_hessians entries, then add a new
  // entry for each system.
  _exact_hessians.clear();

  for (unsigned int sys=0; sys<_equation_systems.n_systems(); ++sys)
    {
      const System& system = _equation_systems.get_system(sys);
      _exact_hessians.push_back
	(new WrappedFunction<Tensor>
          (system, hptr, &_equation_systems.parameters));
    }

  // If we're using exact values, we're not using a fine grid solution
  _equation_systems_fine = NULL;
}


void ExactSolution::attach_exact_hessians (std::vector<FunctionBase<Tensor> *> h)
{
  // Clear out any previous _exact_hessians entries, then add a new
  // entry for each system.
  for (unsigned int i=0; i != _exact_hessians.size(); ++i)
    delete (_exact_hessians[i]);

  _exact_hessians.clear();
  _exact_hessians.resize(h.size(), NULL);

  // We use clone() to get non-sliced copies of FunctionBase
  // subclasses, but we can't put the resulting AutoPtrs into an STL
  // container.
  for (unsigned int i=0; i != h.size(); ++i)
    if (h[i])
      _exact_hessians[i] = h[i]->clone().release();
}


void ExactSolution::attach_exact_hessian (unsigned int sys_num,
                                          FunctionBase<Tensor>* h)
{
  if (_exact_hessians.size() <= sys_num)
  _exact_hessians.resize(sys_num+1, NULL);

  if (h)
    _exact_hessians[sys_num] = h->clone().release();
}


std::vector<Real>& ExactSolution::_check_inputs(const std::string& sys_name,
                                                const std::string& unknown_name)
{
  // If no exact solution function or fine grid solution has been
  // attached, we now just compute the solution norm (i.e. the
  // difference from an "exact solution" of zero

  // Make sure the requested sys_name exists.
  std::map<std::string, SystemErrorMap>::iterator sys_iter =
    _errors.find(sys_name);

  if (sys_iter == _errors.end())
    {
      libMesh::err << "Sorry, couldn't find the requested system '"
                    << sys_name << "'."
		    << std::endl;
      libmesh_error();
    }

  // Make sure the requested unknown_name exists.
  SystemErrorMap::iterator var_iter = (*sys_iter).second.find(unknown_name);

  if (var_iter == (*sys_iter).second.end())
    {
      libMesh::err << "Sorry, couldn't find the requested variable '"
                    << unknown_name << "'."
		    << std::endl;
      libmesh_error();
    }

  // Return a reference to the proper error entry
  return (*var_iter).second;
}






void ExactSolution::compute_error(const std::string& sys_name,
				  const std::string& unknown_name)
{
  // Check the inputs for validity, and get a reference
  // to the proper location to store the error
  std::vector<Real>& error_vals = this->_check_inputs(sys_name,
                                                      unknown_name);
  this->_compute_error(sys_name,
		       unknown_name,
		       error_vals);
}





Real ExactSolution::error_norm(const std::string& sys_name,
                               const std::string& unknown_name,
                               const FEMNormType& norm)
{
  // Check the inputs for validity, and get a reference
  // to the proper location to store the error
  std::vector<Real>& error_vals = this->_check_inputs(sys_name,
                                                      unknown_name);

  switch (norm)
    {
    case L2:
      return std::sqrt(error_vals[0]);
    case H1:
      return std::sqrt(error_vals[0] + error_vals[1]);
    case H2:
      return std::sqrt(error_vals[0] + error_vals[1] + error_vals[2]);
    case H1_SEMINORM:
      return std::sqrt(error_vals[1]);
    case H2_SEMINORM:
      return std::sqrt(error_vals[2]);
    case L1:
      return error_vals[3];
    case L_INF:
      return error_vals[4];
    // Currently only Sobolev norms/seminorms are supported
    default:
      libmesh_error();
    }
}







Real ExactSolution::l2_error(const std::string& sys_name,
                             const std::string& unknown_name)
{

  // Check the inputs for validity, and get a reference
  // to the proper location to store the error
  std::vector<Real>& error_vals = this->_check_inputs(sys_name,
						      unknown_name);

  // Return the square root of the first component of the
  // computed error.
  return std::sqrt(error_vals[0]);
}







Real ExactSolution::l1_error(const std::string& sys_name,
                             const std::string& unknown_name)
{

  // Check the inputs for validity, and get a reference
  // to the proper location to store the error
  std::vector<Real>& error_vals = this->_check_inputs(sys_name,
                                                      unknown_name);

  // Return the square root of the first component of the
  // computed error.
  return error_vals[3];
}







Real ExactSolution::l_inf_error(const std::string& sys_name,
                                const std::string& unknown_name)
{

  // Check the inputs for validity, and get a reference
  // to the proper location to store the error
  std::vector<Real>& error_vals = this->_check_inputs(sys_name,
                                                      unknown_name);

  // Return the square root of the first component of the
  // computed error.
  return error_vals[4];
}







Real ExactSolution::h1_error(const std::string& sys_name,
                             const std::string& unknown_name)
{
  // If the user has supplied no exact derivative function, we
  // just integrate the H1 norm of the solution; i.e. its
  // difference from an "exact solution" of zero.

  // Check the inputs for validity, and get a reference
  // to the proper location to store the error
  std::vector<Real>& error_vals = this->_check_inputs(sys_name,
                                                      unknown_name);

  // Return the square root of the sum of the computed errors.
  return std::sqrt(error_vals[0] + error_vals[1]);
}







Real ExactSolution::h2_error(const std::string& sys_name,
                             const std::string& unknown_name)
{
  // If the user has supplied no exact derivative functions, we
  // just integrate the H2 norm of the solution; i.e. its
  // difference from an "exact solution" of zero.

  // Check the inputs for validity, and get a reference
  // to the proper location to store the error
  std::vector<Real>& error_vals = this->_check_inputs(sys_name,
						      unknown_name);

  // Return the square root of the sum of the computed errors.
  return std::sqrt(error_vals[0] + error_vals[1] + error_vals[2]);
}








void ExactSolution::_compute_error(const std::string& sys_name,
				   const std::string& unknown_name,
				   std::vector<Real>& error_vals)
{
  // This function must be run on all processors at once
  parallel_only();

  // Make sure we aren't "overconfigured"
  libmesh_assert (!(_exact_values.size() && _equation_systems_fine));

  // Get a reference to the system whose error is being computed.
  // If we have a fine grid, however, we'll integrate on that instead
  // for more accuracy.
  const System& computed_system = _equation_systems_fine ?
    _equation_systems_fine->get_system(sys_name) :
    _equation_systems.get_system (sys_name);

  const Real time = _equation_systems.get_system(sys_name).time;

  const unsigned int sys_num = computed_system.number();
  const unsigned int var = computed_system.variable_number(unknown_name);
  const unsigned int var_component =
    computed_system.variable_scalar_number(var, 0);

  // Prepare a global solution and a MeshFunction of the coarse system if we need one
  AutoPtr<MeshFunction> coarse_values;
  AutoPtr<NumericVector<Number> > comparison_soln = NumericVector<Number>::build();
  if (_equation_systems_fine)
    {
      const System& comparison_system
	= _equation_systems.get_system(sys_name);

      std::vector<Number> global_soln;
      comparison_system.update_global_solution(global_soln);
      comparison_soln->init(comparison_system.solution->size(), true, SERIAL);
      (*comparison_soln) = global_soln;

      coarse_values = AutoPtr<MeshFunction>
	(new MeshFunction(_equation_systems,
			  *comparison_soln,
			  comparison_system.get_dof_map(),
			  comparison_system.variable_number(unknown_name)));
      coarse_values->init();
    }

  // Initialize any functors we're going to use
  for (unsigned int i=0; i != _exact_values.size(); ++i)
    if (_exact_values[i])
      _exact_values[i]->init();

  for (unsigned int i=0; i != _exact_derivs.size(); ++i)
    if (_exact_derivs[i])
      _exact_derivs[i]->init();

  for (unsigned int i=0; i != _exact_hessians.size(); ++i)
    if (_exact_hessians[i])
      _exact_hessians[i]->init();

  // Get a reference to the dofmap and mesh for that system
  const DofMap& computed_dof_map = computed_system.get_dof_map();

  const MeshBase& _mesh = computed_system.get_mesh();

  // Zero the error before summation
  error_vals = std::vector<Real>(5, 0.);

  // Construct Quadrature rule based on default quadrature order
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

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
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
  // TODO: this ought to be threaded (and using subordinate
  // MeshFunction objects in each thread rather than a single
  // master)
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
	      u_h      += phi_values[i][qp]*computed_system.current_solution  (dof_indices[i]);
	      grad_u_h += dphi_values[i][qp]*computed_system.current_solution (dof_indices[i]);
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
	      grad2_u_h += d2phi_values[i][qp]*computed_system.current_solution (dof_indices[i]);
#endif
	    }

	  // Compute the value of the error at this quadrature point
	  Number exact_val = 0.0;
          if (_exact_values.size() > sys_num && _exact_values[sys_num])
	    exact_val =
              _exact_values[sys_num]->
		component(var_component, q_point[qp], time);
	  else if (_equation_systems_fine)
	    exact_val = (*coarse_values)(q_point[qp]);

	  const Number val_error = u_h - exact_val;

	  // Add the squares of the error to each contribution
	  error_vals[0] += JxW[qp]*libmesh_norm(val_error);
	  Real norm = sqrt(libmesh_norm(val_error));
	  error_vals[3] += JxW[qp]*norm;
	  if(error_vals[4]<norm) { error_vals[4] = norm; }

	  // Compute the value of the error in the gradient at this
	  // quadrature point
          Gradient exact_grad;
	  if (_exact_derivs.size() > sys_num && _exact_derivs[sys_num])
	    exact_grad =
               _exact_derivs[sys_num]->
                 component(var_component, q_point[qp], time);
	  else if (_equation_systems_fine)
	    exact_grad = coarse_values->gradient(q_point[qp]);

	  const Gradient grad_error = grad_u_h - exact_grad;

	  error_vals[1] += JxW[qp]*grad_error.size_sq();


#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
	  // Compute the value of the error in the hessian at this
	  // quadrature point
          Tensor exact_hess;
	  if (_exact_hessians.size() > sys_num && _exact_hessians[sys_num])
	    exact_hess =
              _exact_hessians[sys_num]->
                component(var_component, q_point[qp], time);
	  else if (_equation_systems_fine)
	    exact_hess = coarse_values->hessian(q_point[qp]);

	  const Tensor grad2_error = grad2_u_h - exact_hess;

	  error_vals[2] += JxW[qp]*grad2_error.size_sq();
#endif

	} // end qp loop
    } // end element loop

  // Add up the error values on all processors, except for the L-infty
  // norm, for which the maximum is computed.
  Real l_infty_norm = error_vals[4];
  Parallel::max(l_infty_norm);
  Parallel::sum(error_vals);
  error_vals[4] = l_infty_norm;
}

} // namespace libMesh
