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


// Local Includes
#include "libmesh/mesh_function.h"
#include "libmesh/dense_vector.h"
#include "libmesh/equation_systems.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dof_map.h"
#include "libmesh/point_locator_tree.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fe_compute_data.h"
#include "libmesh/mesh_base.h"
#include "libmesh/point.h"

namespace libMesh
{


//------------------------------------------------------------------
// MeshFunction methods
MeshFunction::MeshFunction (const EquationSystems& eqn_systems,
			    const NumericVector<Number>& vec,
			    const DofMap& dof_map,
			    const std::vector<unsigned int>& vars,
			    const FunctionBase<Number>* master) :
  FunctionBase<Number> (master),
  ParallelObject       (eqn_systems),
  _eqn_systems         (eqn_systems),
  _vector              (vec),
  _dof_map             (dof_map),
  _system_vars         (vars),
  _point_locator       (NULL),
  _out_of_mesh_mode    (false),
  _out_of_mesh_value   ()
{
}



MeshFunction::MeshFunction (const EquationSystems& eqn_systems,
			    const NumericVector<Number>& vec,
			    const DofMap& dof_map,
			    const unsigned int var,
			    const FunctionBase<Number>* master) :
  FunctionBase<Number> (master),
  ParallelObject       (eqn_systems),
  _eqn_systems         (eqn_systems),
  _vector              (vec),
  _dof_map             (dof_map),
  _system_vars         (1,var),
  _point_locator       (NULL),
  _out_of_mesh_mode    (false),
  _out_of_mesh_value   ()
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




void MeshFunction::init (const Trees::BuildType /*point_locator_build_type*/)
{
  // are indices of the desired variable(s) provided?
  libmesh_assert_greater (this->_system_vars.size(), 0);

  // Don't do twice...
  if (this->_initialized)
    {
      libmesh_assert(this->_point_locator);
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
	libmesh_cast_ptr<const MeshFunction*>(this->_master);

      if (master->_point_locator == NULL)
        {
	  libMesh::err << "ERROR: When the master-servant concept is used,"
		        << std::endl
		        << " the master has to be initialized first!"
		        << std::endl;
	  libmesh_error();
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
      const MeshBase& mesh = this->_eqn_systems.get_mesh();

      // build the point locator.  Only \p TREE version available
      //AutoPtr<PointLocatorBase> ap (PointLocatorBase::build (TREE, mesh));
      //this->_point_locator = ap.release();
      // this->_point_locator = new PointLocatorTree (mesh, point_locator_build_type);
       this->_point_locator = mesh.sub_point_locator().release();

      // Point locator no longer needs to be initialized.
      //      this->_point_locator->init();
    }


  // ready for use
  this->_initialized = true;
}


void
MeshFunction::clear ()
{
  // only delete the point locator when we are the master
  if ((this->_point_locator != NULL) && (this->_master == NULL))
    {
      delete this->_point_locator;
      this->_point_locator = NULL;
    }
  this->_initialized = false;
}



AutoPtr<FunctionBase<Number> > MeshFunction::clone () const
{
  return AutoPtr<FunctionBase<Number> >
    (new MeshFunction
      (_eqn_systems, _vector, _dof_map, _system_vars, this));
}



Number MeshFunction::operator() (const Point& p,
				 const Real time)
{
  libmesh_assert (this->initialized());

  DenseVector<Number> buf (1);
  this->operator() (p, time, buf);
  return buf(0);
}



Gradient MeshFunction::gradient (const Point& p,
				 const Real time)
{
  libmesh_assert (this->initialized());

  std::vector<Gradient> buf (1);
  this->gradient(p, time, buf);
  return buf[0];
}



#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
Tensor MeshFunction::hessian (const Point& p,
			      const Real time)
{
  libmesh_assert (this->initialized());

  std::vector<Tensor> buf (1);
  this->hessian(p, time, buf);
  return buf[0];
}
#endif



void MeshFunction::operator() (const Point& p,
			       const Real,
			       DenseVector<Number>& output)
{
  libmesh_assert (this->initialized());

  /* Ensure that in the case of a master mesh function, the
     out-of-mesh mode is enabled either for both or for none.  This is
     important because the out-of-mesh mode is also communicated to
     the point locator.  Since this is time consuming, enable it only
     in debug mode.  */
#ifdef DEBUG
  if (this->_master != NULL)
    {
      const MeshFunction* master =
	libmesh_cast_ptr<const MeshFunction*>(this->_master);
      if(_out_of_mesh_mode!=master->_out_of_mesh_mode)
	{
	  libMesh::err << "ERROR: If you use out-of-mesh-mode in connection with master mesh functions, you must enable out-of-mesh mode for both the master and the slave mesh function." << std::endl;
	  libmesh_error();
	}
    }
#endif

  // locate the point in the other mesh
  const Elem* element = this->_point_locator->operator()(p);

  // If we have an element, but it's not a local element, then we
  // either need to have a serialized vector or we need to find a
  // local element sharing the same point.
  if (element &&
     (element->processor_id() != this->processor_id()) &&
     _vector.type() != SERIAL)
    {
      // look for a local element containing the point
      std::set<const Elem*> point_neighbors;
      element->find_point_neighbors(p, point_neighbors);
      element = NULL;
      std::set<const Elem*>::const_iterator       it  = point_neighbors.begin();
      const std::set<const Elem*>::const_iterator end = point_neighbors.end();
      for (; it != end; ++it)
        {
          const Elem* elem = *it;
          if (elem->processor_id() == this->processor_id())
            {
              element = elem;
              break;
            }
        }
    }

  if (!element)
    {
      output = _out_of_mesh_value;
    }
  else
    {
      // resize the output vector to the number of output values
      // that the user told us
      output.resize (libmesh_cast_int<unsigned int>
		     (this->_system_vars.size()));


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
	      std::vector<dof_id_type> dof_indices;
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
    }

  // all done
  return;
}



void MeshFunction::gradient (const Point& p,
			     const Real,
			     std::vector<Gradient>& output)
{
  libmesh_assert (this->initialized());

  /* Ensure that in the case of a master mesh function, the
     out-of-mesh mode is enabled either for both or for none.  This is
     important because the out-of-mesh mode is also communicated to
     the point locator.  Since this is time consuming, enable it only
     in debug mode.  */
#ifdef DEBUG
  if (this->_master != NULL)
    {
      const MeshFunction* master =
	libmesh_cast_ptr<const MeshFunction*>(this->_master);
      if(_out_of_mesh_mode!=master->_out_of_mesh_mode)
	{
	  libMesh::err << "ERROR: If you use out-of-mesh-mode in connection with master mesh functions, you must enable out-of-mesh mode for both the master and the slave mesh function." << std::endl;
	  libmesh_error();
	}
    }
#endif

  // locate the point in the other mesh
  const Elem* element = this->_point_locator->operator()(p);

  // If we have an element, but it's not a local element, then we
  // either need to have a serialized vector or we need to find a
  // local element sharing the same point.
  if (element &&
     (element->processor_id() != this->processor_id()) &&
     _vector.type() != SERIAL)
    {
      // look for a local element containing the point
      std::set<const Elem*> point_neighbors;
      element->find_point_neighbors(p, point_neighbors);
      element = NULL;
      std::set<const Elem*>::const_iterator       it  = point_neighbors.begin();
      const std::set<const Elem*>::const_iterator end = point_neighbors.end();
      for (; it != end; ++it)
        {
          const Elem* elem = *it;
          if (elem->processor_id() == this->processor_id())
            {
              element = elem;
              break;
            }
        }
    }

  if (!element)
    {
      output.resize(0);
    }
  else
    {
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

        std::vector<Point> point_list (1, mapped_point);

	// loop over all vars
	for (unsigned int index=0; index < this->_system_vars.size(); index++)
	  {
	    /*
	     * the data for this variable
	     */
	    const unsigned int var = _system_vars[index];
	    const FEType& fe_type = this->_dof_map.variable_type(var);

            AutoPtr<FEBase> point_fe (FEBase::build(dim, fe_type));
            const std::vector<std::vector<RealGradient> >& dphi = point_fe->get_dphi();
            point_fe->reinit(element, &point_list);

	    // where the solution values for the var-th variable are stored
	    std::vector<dof_id_type> dof_indices;
	    this->_dof_map.dof_indices (element, dof_indices, var);

	    // interpolate the solution
	    Gradient grad(0.);

	    for (unsigned int i=0; i<dof_indices.size(); i++)
	      grad.add_scaled(dphi[i][0], this->_vector(dof_indices[i]));

	    output[index] = grad;
	  }
      }
    }

  // all done
  return;
}



#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
void MeshFunction::hessian (const Point& p,
			    const Real,
			    std::vector<Tensor>& output)
{
  libmesh_assert (this->initialized());

  /* Ensure that in the case of a master mesh function, the
     out-of-mesh mode is enabled either for both or for none.  This is
     important because the out-of-mesh mode is also communicated to
     the point locator.  Since this is time consuming, enable it only
     in debug mode.  */
#ifdef DEBUG
  if (this->_master != NULL)
    {
      const MeshFunction* master =
	libmesh_cast_ptr<const MeshFunction*>(this->_master);
      if(_out_of_mesh_mode!=master->_out_of_mesh_mode)
	{
	  libMesh::err << "ERROR: If you use out-of-mesh-mode in connection with master mesh functions, you must enable out-of-mesh mode for both the master and the slave mesh function." << std::endl;
	  libmesh_error();
	}
    }
#endif

  // locate the point in the other mesh
  const Elem* element = this->_point_locator->operator()(p);

  // If we have an element, but it's not a local element, then we
  // either need to have a serialized vector or we need to find a
  // local element sharing the same point.
  if (element &&
     (element->processor_id() != this->processor_id()) &&
     _vector.type() != SERIAL)
    {
      // look for a local element containing the point
      std::set<const Elem*> point_neighbors;
      element->find_point_neighbors(p, point_neighbors);
      element = NULL;
      std::set<const Elem*>::const_iterator       it  = point_neighbors.begin();
      const std::set<const Elem*>::const_iterator end = point_neighbors.end();
      for (; it != end; ++it)
        {
          const Elem* elem = *it;
          if (elem->processor_id() == this->processor_id())
            {
              element = elem;
              break;
            }
        }
    }

  if (!element)
    {
      output.resize(0);
    }
  else
    {
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

        std::vector<Point> point_list (1, mapped_point);

	// loop over all vars
	for (unsigned int index=0; index < this->_system_vars.size(); index++)
	  {
	    /*
	     * the data for this variable
	     */
	    const unsigned int var = _system_vars[index];
	    const FEType& fe_type = this->_dof_map.variable_type(var);

            AutoPtr<FEBase> point_fe (FEBase::build(dim, fe_type));
            const std::vector<std::vector<RealTensor> >& d2phi =
			    point_fe->get_d2phi();
            point_fe->reinit(element, &point_list);

	    // where the solution values for the var-th variable are stored
	    std::vector<dof_id_type> dof_indices;
	    this->_dof_map.dof_indices (element, dof_indices, var);

	    // interpolate the solution
	    Tensor hess;

	    for (unsigned int i=0; i<dof_indices.size(); i++)
	      hess.add_scaled(d2phi[i][0], this->_vector(dof_indices[i]));

	    output[index] = hess;
	  }
      }
    }

  // all done
  return;
}
#endif



const PointLocatorBase& MeshFunction::get_point_locator (void) const
{
  libmesh_assert (this->initialized());
  return *_point_locator;
}

void MeshFunction::enable_out_of_mesh_mode(const DenseVector<Number>& value)
{
  libmesh_assert (this->initialized());
  _point_locator->enable_out_of_mesh_mode();
  _out_of_mesh_mode = true;
  _out_of_mesh_value = value;
}

void MeshFunction::enable_out_of_mesh_mode(const Number& value)
{
  DenseVector<Number> v(1);
  v(0) = value;
  this->enable_out_of_mesh_mode(v);
}

void MeshFunction::disable_out_of_mesh_mode(void)
{
  libmesh_assert (this->initialized());
  _point_locator->disable_out_of_mesh_mode();
  _out_of_mesh_mode = false;
}

} // namespace libMesh
