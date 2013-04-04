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
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/exact_solution.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe_base.h"
#include "libmesh/function_base.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_function.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/parallel.h"
#include "libmesh/quadrature.h"
#include "libmesh/wrapped_function.h"
#include "libmesh/fe_interface.h"
#include "libmesh/raw_accessor.h"
#include "libmesh/tensor_tools.h"

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
  libmesh_assert(es_fine);
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
  libmesh_assert(fptr);

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
  libmesh_assert(gptr);

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
  libmesh_assert(hptr);

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

  libmesh_assert( _equation_systems.has_system(sys_name) );
  const System& sys = _equation_systems.get_system<System>( sys_name );

  libmesh_assert( sys.has_variable( unknown_name ) );
  switch( FEInterface::field_type(sys.variable_type( unknown_name )) )
    {
    case TYPE_SCALAR:
      {
	this->_compute_error<Real>(sys_name,
				   unknown_name,
				   error_vals);
	break;
      }
    case TYPE_VECTOR:
      {
	this->_compute_error<RealGradient>(sys_name,
					   unknown_name,
					   error_vals);
	break;
      }
    default:
      libmesh_error();
    }

  return;
}





Real ExactSolution::error_norm(const std::string& sys_name,
                               const std::string& unknown_name,
                               const FEMNormType& norm)
{
  // Check the inputs for validity, and get a reference
  // to the proper location to store the error
  std::vector<Real>& error_vals = this->_check_inputs(sys_name,
                                                      unknown_name);

  libmesh_assert(_equation_systems.has_system(sys_name));
  libmesh_assert(_equation_systems.get_system(sys_name).has_variable( unknown_name ));
  const FEType& fe_type = _equation_systems.get_system(sys_name).variable_type(unknown_name);

  switch (norm)
    {
    case L2:
      return std::sqrt(error_vals[0]);
    case H1:
      return std::sqrt(error_vals[0] + error_vals[1]);
    case H2:
      return std::sqrt(error_vals[0] + error_vals[1] + error_vals[2]);
    case HCURL:
      {
	if(FEInterface::field_type(fe_type) == TYPE_SCALAR)
	  {
	    libMesh::err << "Cannot compute HCurl error norm of scalar-valued variables!" << std::endl;
	    libmesh_error();
	  }
	else
	  return std::sqrt(error_vals[0] + error_vals[5]);
      }
    case HDIV:
      {
	if(FEInterface::field_type(fe_type) == TYPE_SCALAR)
	  {
	    libMesh::err << "Cannot compute HDiv error norm of scalar-valued variables!" << std::endl;
	    libmesh_error();
	  }
	else
	  return std::sqrt(error_vals[0] + error_vals[6]);
      }
    case H1_SEMINORM:
      return std::sqrt(error_vals[1]);
    case H2_SEMINORM:
      return std::sqrt(error_vals[2]);
    case HCURL_SEMINORM:
      {
	if(FEInterface::field_type(fe_type) == TYPE_SCALAR)
	  {
	    libMesh::err << "Cannot compute HCurl error seminorm of scalar-valued variables!" << std::endl;
	    libmesh_error();
	  }
	else
	  return std::sqrt(error_vals[5]);
      }
    case HDIV_SEMINORM:
      {
	if(FEInterface::field_type(fe_type) == TYPE_SCALAR)
	  {
	    libMesh::err << "Cannot compute HDiv error seminorm of scalar-valued variables!" << std::endl;
	    libmesh_error();
	  }
	else
	  return std::sqrt(error_vals[6]);
      }
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


Real ExactSolution::hcurl_error(const std::string& sys_name,
				const std::string& unknown_name)
{
  return this->error_norm(sys_name,unknown_name,HCURL);
}


Real ExactSolution::hdiv_error(const std::string& sys_name,
			       const std::string& unknown_name)
{
  return this->error_norm(sys_name,unknown_name,HDIV);
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







template< typename OutputShape>
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
  AutoPtr<NumericVector<Number> > comparison_soln = NumericVector<Number>::build(libMesh::default_solver_package(), _equation_systems.communicator());
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

  const unsigned int dim = _mesh.mesh_dimension();

  // Zero the error before summation
  // 0 - sum of square of function error (L2)
  // 1 - sum of square of gradient error (H1 semi)
  // 2 - sum of square of Hessian error (H2 semi)
  // 3 - sum of sqrt(square of function error) (L1)
  // 4 - max of sqrt(square of function error) (Linfty)
  // 5 - sum of square of curl error (HCurl semi)
  // 6 - sum of square of div error (HDiv semi)
  error_vals = std::vector<Real>(7, 0.);

  // Construct Quadrature rule based on default quadrature order
  const FEType& fe_type  = computed_dof_map.variable_type(var);

  unsigned int n_vec_dim = FEInterface::n_vec_dim( _mesh, fe_type );

  // FIXME: MeshFunction needs to be updated to support vector-valued
  //        elements before we can use a reference solution.
  if( (n_vec_dim > 1) && _equation_systems_fine )
    {
      libMesh::err << "Error calculation using reference solution not yet\n"
		   << "supported for vector-valued elements."
		   << std::endl;
      libmesh_not_implemented();
    }

  AutoPtr<QBase> qrule =
    fe_type.default_quadrature_rule (dim,
                                     _extra_order);

  // Construct finite element object

  AutoPtr<FEGenericBase<OutputShape> > fe(FEGenericBase<OutputShape>::build(dim, fe_type));

  // Attach quadrature rule to FE object
  fe->attach_quadrature_rule (qrule.get());

  // The Jacobian*weight at the quadrature points.
  const std::vector<Real>& JxW                               = fe->get_JxW();

  // The value of the shape functions at the quadrature points
  // i.e. phi(i) = phi_values[i][qp]
  const std::vector<std::vector<OutputShape> >&  phi_values         = fe->get_phi();

  // The value of the shape function gradients at the quadrature points
  const std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputGradient> >&
		    dphi_values = fe->get_dphi();

  // The value of the shape function curls at the quadrature points
  // Only computed for vector-valued elements
  const std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputShape> >* curl_values = NULL;

  // The value of the shape function divergences at the quadrature points
  // Only computed for vector-valued elements
  const std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputDivergence> >* div_values = NULL;

  if( FEInterface::field_type(fe_type) == TYPE_VECTOR )
    {
      curl_values = &fe->get_curl_phi();
      div_values = &fe->get_div_phi();
    }

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  // The value of the shape function second derivatives at the quadrature points
  const std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputTensor> >&
		    d2phi_values = fe->get_d2phi();
#endif

  // The XYZ locations (in physical space) of the quadrature points
  const std::vector<Point>& q_point                          = fe->get_xyz();

  // The global degree of freedom indices associated
  // with the local degrees of freedom.
  std::vector<dof_id_type> dof_indices;


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
      const unsigned int n_sf =
	libmesh_cast_int<unsigned int>(dof_indices.size());

      //
      // Begin the loop over the Quadrature points.
      //
      for (unsigned int qp=0; qp<n_qp; qp++)
	{
	  // Real u_h = 0.;
	  // RealGradient grad_u_h;

	  typename FEGenericBase<OutputShape>::OutputNumber u_h = 0.;

	  typename FEGenericBase<OutputShape>::OutputNumberGradient grad_u_h;
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
	  typename FEGenericBase<OutputShape>::OutputNumberTensor grad2_u_h;
#endif
	  typename FEGenericBase<OutputShape>::OutputNumber curl_u_h = 0.0;
	  typename FEGenericBase<OutputShape>::OutputNumberDivergence div_u_h = 0.0;

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
	      if( FEInterface::field_type(fe_type) == TYPE_VECTOR )
		{
		  curl_u_h += (*curl_values)[i][qp]*computed_system.current_solution (dof_indices[i]);
		  div_u_h += (*div_values)[i][qp]*computed_system.current_solution (dof_indices[i]);
		}
	    }

	  // Compute the value of the error at this quadrature point
	  typename FEGenericBase<OutputShape>::OutputNumber exact_val = 0;
	  RawAccessor<typename FEGenericBase<OutputShape>::OutputNumber> exact_val_accessor( exact_val, dim );
          if (_exact_values.size() > sys_num && _exact_values[sys_num])
	    {
	      for( unsigned int c = 0; c < n_vec_dim; c++)
		exact_val_accessor(c) =
		  _exact_values[sys_num]->
		  component(var_component+c, q_point[qp], time);
	    }
	  else if (_equation_systems_fine)
	    {
	      // FIXME: Needs to be updated for vector-valued elements
	      exact_val = (*coarse_values)(q_point[qp]);
	    }
	  const typename FEGenericBase<OutputShape>::OutputNumber val_error = u_h - exact_val;

	  // Add the squares of the error to each contribution
	  Real error_sq = TensorTools::norm_sq(val_error);
	  error_vals[0] += JxW[qp]*error_sq;

	  Real norm = sqrt(error_sq);
	  error_vals[3] += JxW[qp]*norm;

	  if(error_vals[4]<norm) { error_vals[4] = norm; }

	  // Compute the value of the error in the gradient at this
	  // quadrature point
	  typename FEGenericBase<OutputShape>::OutputNumberGradient exact_grad;
	  RawAccessor<typename FEGenericBase<OutputShape>::OutputNumberGradient> exact_grad_accessor( exact_grad, _mesh.spatial_dimension() );
	  if (_exact_derivs.size() > sys_num && _exact_derivs[sys_num])
	    {
	      for( unsigned int c = 0; c < n_vec_dim; c++)
		for( unsigned int d = 0; d < _mesh.spatial_dimension(); d++ )
		  exact_grad_accessor(d + c*_mesh.spatial_dimension() ) =
		    _exact_derivs[sys_num]->
		    component(var_component+c, q_point[qp], time)(d);
	    }
	  else if (_equation_systems_fine)
	    {
	      // FIXME: Needs to be updated for vector-valued elements
	      exact_grad = coarse_values->gradient(q_point[qp]);
	    }

	  const typename FEGenericBase<OutputShape>::OutputNumberGradient grad_error = grad_u_h - exact_grad;

	  error_vals[1] += JxW[qp]*grad_error.size_sq();


	  if( FEInterface::field_type(fe_type) == TYPE_VECTOR )
	    {
	      // Compute the value of the error in the curl at this
	      // quadrature point
	      typename FEGenericBase<OutputShape>::OutputNumber exact_curl = 0.0;
	      if (_exact_derivs.size() > sys_num && _exact_derivs[sys_num])
		{
		  exact_curl = TensorTools::curl_from_grad( exact_grad );
		}
	      else if (_equation_systems_fine)
		{
		  // FIXME: Need to implement curl for MeshFunction and support reference
		  //        solution for vector-valued elements
		}

	      const typename FEGenericBase<OutputShape>::OutputNumber curl_error = curl_u_h - exact_curl;

	      error_vals[5] += JxW[qp]*TensorTools::norm_sq(curl_error);

	      // Compute the value of the error in the divergence at this
	      // quadrature point
	      typename FEGenericBase<OutputShape>::OutputNumberDivergence exact_div = 0.0;
	      if (_exact_derivs.size() > sys_num && _exact_derivs[sys_num])
		{
		  exact_div = TensorTools::div_from_grad( exact_grad );
		}
	      else if (_equation_systems_fine)
		{
		  // FIXME: Need to implement div for MeshFunction and support reference
		  //        solution for vector-valued elements
		}

	      const typename FEGenericBase<OutputShape>::OutputNumberDivergence div_error = div_u_h - exact_div;

	      error_vals[6] += JxW[qp]*TensorTools::norm_sq(div_error);
	    }

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
	  // Compute the value of the error in the hessian at this
	  // quadrature point
	  typename FEGenericBase<OutputShape>::OutputNumberTensor exact_hess;
	  RawAccessor<typename FEGenericBase<OutputShape>::OutputNumberTensor> exact_hess_accessor( exact_hess, dim );
	  if (_exact_hessians.size() > sys_num && _exact_hessians[sys_num])
	    {
	      //FIXME: This needs to be implemented to support rank 3 tensors
	      //       which can't happen until type_n_tensor is fully implemented
	      //       and a RawAccessor<TypeNTensor> is fully implemented
	      if( FEInterface::field_type(fe_type) == TYPE_VECTOR )
		libmesh_not_implemented();

	      for( unsigned int c = 0; c < n_vec_dim; c++)
		for( unsigned int d = 0; d < dim; d++ )
		  for( unsigned int e =0; e < dim; e++ )
		    exact_hess_accessor(d + e*dim + c*dim*dim) =
		      _exact_hessians[sys_num]->
		      component(var_component+c, q_point[qp], time)(d,e);
	    }
	  else if (_equation_systems_fine)
	    {
	      // FIXME: Needs to be updated for vector-valued elements
	      exact_hess = coarse_values->hessian(q_point[qp]);
	    }

	  const typename FEGenericBase<OutputShape>::OutputNumberTensor grad2_error = grad2_u_h - exact_hess;

	  // FIXME: PB: Is this what we want for rank 3 tensors?
	  error_vals[2] += JxW[qp]*grad2_error.size_sq();
#endif

	} // end qp loop
    } // end element loop

  // Add up the error values on all processors, except for the L-infty
  // norm, for which the maximum is computed.
  Real l_infty_norm = error_vals[4];
  CommWorld.max(l_infty_norm);
  CommWorld.sum(error_vals);
  error_vals[4] = l_infty_norm;
}

  // Explicit instantiations of templated member functions
  template void ExactSolution::_compute_error<Real>(const std::string&, const std::string&, std::vector<Real>&);
  template void ExactSolution::_compute_error<RealGradient>(const std::string&, const std::string&, std::vector<Real>&);

} // namespace libMesh
