// $Id: mesh_function.C,v 1.9 2006-03-29 19:16:40 roystgnr Exp $

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


// Local Includes
#include "mesh_function.h"
#include "equation_systems.h"
#include "numeric_vector.h"
#include "dof_map.h"
#include "point_locator_base.h"
#include "fe_interface.h"
#include "fe_compute_data.h"
#include "mesh.h"
#include "point.h"


//------------------------------------------------------------------
// MeshFunction methods
MeshFunction::MeshFunction (const EquationSystems& eqn_systems,
			    const NumericVector<Number>& vec,
			    const DofMap& dof_map,
			    const std::vector<unsigned int>& vars,
			    const FunctionBase* master) :
  FunctionBase   (master),
  _eqn_systems   (eqn_systems),
  _vector        (vec),
  _dof_map       (dof_map),
  _system_vars   (vars),
  _point_locator (NULL)
{
}



MeshFunction::MeshFunction (const EquationSystems& eqn_systems,
			    const NumericVector<Number>& vec,
			    const DofMap& dof_map,
			    const unsigned int var,
			    const FunctionBase* master) :
  FunctionBase   (master),
  _eqn_systems   (eqn_systems),
  _vector        (vec),
  _dof_map       (dof_map),
  _system_vars   (var),
  _point_locator (NULL)
{
//   std::vector<unsigned int> buf (1);
//   buf[0] = var;
//   _system_vars (buf);
}







MeshFunction::~MeshFunction ()
{
  // only delete the point locator when we are the master
  if ((this->_point_locator != NULL) && (this->_master == NULL))
    delete this->_point_locator;
}




void MeshFunction::init ()
{
  // are indices of the desired variable(s) provided?
  assert (this->_system_vars.size() > 0);

  // Don't do twice...
  if (this->_initialized)
    {
      assert(this->_point_locator != NULL);
      return;
    }

  /*
   * set up the PointLocator: either someone else
   * is the master (go and get the address of his
   * point locator) or this object is the master 
   * (build the point locator  on our own).
   */
  if (this->_master != NULL)
    {
      // we aren't the master
      const MeshFunction* master =
	dynamic_cast<const MeshFunction*>(this->_master);
      
      if (master->_point_locator == NULL)
        {
	  std::cerr << "ERROR: When the master-servant concept is used,"
		    << std::endl
		    << " the master has to be initialized first!"
		    << std::endl;
	  error();
	}
      else
        {
	  this->_point_locator = master->_point_locator;
	}
    }
  else
    {
      // we are the master: build the point locator

      // constant reference to the other mesh
      const Mesh& mesh = this->_eqn_systems.get_mesh();

      // build the point locator.  Only \p TREE version available
      AutoPtr<PointLocatorBase> ap (PointLocatorBase::build (TREE, mesh));
      this->_point_locator = ap.release();

      // initialize the point locator, so that it is ready for use
      this->_point_locator->init();
    }


  // ready for use
  this->_initialized = true;
}




Number MeshFunction::operator() (const Point& p, 
				 const Real time)
{
  assert (this->initialized());
  // At the moment the function we call ignores the time
  assert (time == 0.);
  
  DenseVector<Number> buf (1);
  this->operator() (p, time, buf);
  return buf(0);
}




void MeshFunction::operator() (const Point& p,
			       const Real,
			       DenseVector<Number>& output)
{
  assert (this->initialized());

  // locate the point in the other mesh
  const Elem* element = this->_point_locator->operator()(p);

  // resize the output vector to the number of output values
  // that the user told us
  output.resize (this->_system_vars.size());


  {
    const unsigned int dim = this->_eqn_systems.get_mesh().mesh_dimension();


    /*
     * Get local coordinates to feed these into compute_data().  
     * Note that the fe_type can safely be used from the 0-variable,
     * since the inverse mapping is the same for all FEFamilies
     */
    const Point mapped_point (FEInterface::inverse_map (dim, 
							this->_dof_map.variable_type(0),
							element, 
							p));


    // loop over all vars
    for (unsigned int index=0; index < this->_system_vars.size(); index++)
      {
	/*
	 * the data for this variable
	 */
	const unsigned int var = _system_vars[index];
	const FEType& fe_type = this->_dof_map.variable_type(var);

	/**
	 * Build an FEComputeData that contains both input and output data
	 * for the specific compute_data method.
	 */
	{
	  FEComputeData data (this->_eqn_systems, mapped_point);
    
	  FEInterface::compute_data (dim, fe_type, element, data);

	  // where the solution values for the var-th variable are stored
	  std::vector<unsigned int> dof_indices;
	  this->_dof_map.dof_indices (element, dof_indices, var);
	
	  // interpolate the solution
	  {
	    Number value = 0.;

	    for (unsigned int i=0; i<dof_indices.size(); i++)
	      value += this->_vector(dof_indices[i]) * data.shape[i];

	    output(index) = value;
	  }

	}

	// next variable
      }
  }

  // all done
  return;
}




