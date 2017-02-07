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
#include "libmesh/elem.h"

namespace libMesh
{


//------------------------------------------------------------------
// MeshFunction methods
MeshFunction::MeshFunction (const EquationSystems & eqn_systems,
                            const NumericVector<Number> & vec,
                            const DofMap & dof_map,
                            const std::vector<unsigned int> & vars,
                            const FunctionBase<Number> * master) :
  FunctionBase<Number> (master),
  ParallelObject       (eqn_systems),
  _eqn_systems         (eqn_systems),
  _vector              (vec),
  _dof_map             (dof_map),
  _system_vars         (vars),
  _point_locator       (libmesh_nullptr),
  _out_of_mesh_mode    (false),
  _out_of_mesh_value   ()
{
}



MeshFunction::MeshFunction (const EquationSystems & eqn_systems,
                            const NumericVector<Number> & vec,
                            const DofMap & dof_map,
                            const unsigned int var,
                            const FunctionBase<Number> * master) :
  FunctionBase<Number> (master),
  ParallelObject       (eqn_systems),
  _eqn_systems         (eqn_systems),
  _vector              (vec),
  _dof_map             (dof_map),
  _system_vars         (1,var),
  _point_locator       (libmesh_nullptr),
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
  if (this->_master == libmesh_nullptr)
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
  if (this->_master != libmesh_nullptr)
    {
      // we aren't the master
      const MeshFunction * master =
        cast_ptr<const MeshFunction *>(this->_master);

      if (master->_point_locator == libmesh_nullptr)
        libmesh_error_msg("ERROR: When the master-servant concept is used, the master has to be initialized first!");

      else
        {
          this->_point_locator = master->_point_locator;
        }
    }
  else
    {
      // we are the master: build the point locator

      // constant reference to the other mesh
      const MeshBase & mesh = this->_eqn_systems.get_mesh();

      // build the point locator.  Only \p TREE version available
      //UniquePtr<PointLocatorBase> ap (PointLocatorBase::build (TREE, mesh));
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
  if ((this->_point_locator != libmesh_nullptr) && (this->_master == libmesh_nullptr))
    {
      delete this->_point_locator;
      this->_point_locator = libmesh_nullptr;
    }
  this->_initialized = false;
}



UniquePtr<FunctionBase<Number> > MeshFunction::clone () const
{
  FunctionBase<Number> * mf_clone =
    new MeshFunction(_eqn_systems, _vector, _dof_map, _system_vars, this);

  if(this->initialized())
    mf_clone->init();

  return UniquePtr< FunctionBase<Number> >(mf_clone);
}



Number MeshFunction::operator() (const Point & p,
                                 const Real time)
{
  libmesh_assert (this->initialized());

  DenseVector<Number> buf (1);
  this->operator() (p, time, buf);
  return buf(0);
}

std::map<const Elem *, Number> MeshFunction::discontinuous_value (const Point & p,
                                                                  const Real time)
{
  libmesh_assert (this->initialized());

  std::map<const Elem *, DenseVector<Number> > buffer;
  this->discontinuous_value (p, time, buffer);
  std::map<const Elem *, Number> return_value;
  for (std::map<const Elem *, DenseVector<Number> >::const_iterator it = buffer.begin(); it != buffer.end(); ++it)
    return_value[it->first] = it->second(0);
  // NOTE: If no suitable element is found, then the map return_value is empty. This
  // puts burden on the user of this function but I don't really see a better way.
  return return_value;
}

Gradient MeshFunction::gradient (const Point & p,
                                 const Real time)
{
  libmesh_assert (this->initialized());

  std::vector<Gradient> buf (1);
  this->gradient(p, time, buf);
  return buf[0];
}

std::map<const Elem *, Gradient> MeshFunction::discontinuous_gradient (const Point & p,
                                                                       const Real time)
{
  libmesh_assert (this->initialized());

  std::map<const Elem *, std::vector<Gradient> > buffer;
  this->discontinuous_gradient (p, time, buffer);
  std::map<const Elem *, Gradient> return_value;
  for (std::map<const Elem *, std::vector<Gradient> >::const_iterator it = buffer.begin(); it != buffer.end(); ++it)
    return_value[it->first] = it->second[0];
  // NOTE: If no suitable element is found, then the map return_value is empty. This
  // puts burden on the user of this function but I don't really see a better way.
  return return_value;
}

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
Tensor MeshFunction::hessian (const Point & p,
                              const Real time)
{
  libmesh_assert (this->initialized());

  std::vector<Tensor> buf (1);
  this->hessian(p, time, buf);
  return buf[0];
}
#endif

void MeshFunction::operator() (const Point & p,
                               const Real time,
                               DenseVector<Number> & output)
{
  this->operator() (p, time, output, libmesh_nullptr);
}

void MeshFunction::operator() (const Point & p,
                               const Real,
                               DenseVector<Number> & output,
                               const std::set<subdomain_id_type> * subdomain_ids)
{
  libmesh_assert (this->initialized());

  const Elem * element = this->find_element(p,subdomain_ids);

  if (!element)
    {
      output = _out_of_mesh_value;
    }
  else
    {
      // resize the output vector to the number of output values
      // that the user told us
      output.resize (cast_int<unsigned int>
                     (this->_system_vars.size()));


      {
        const unsigned int dim = element->dim();


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
        for (std::size_t index=0; index < this->_system_vars.size(); index++)
          {
            /*
             * the data for this variable
             */
            const unsigned int var = _system_vars[index];
            const FEType & fe_type = this->_dof_map.variable_type(var);

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

                for (std::size_t i=0; i<dof_indices.size(); i++)
                  value += this->_vector(dof_indices[i]) * data.shape[i];

                output(index) = value;
              }

            }

            // next variable
          }
      }
    }
}


void MeshFunction::discontinuous_value (const Point & p,
                                        const Real time,
                                        std::map<const Elem *, DenseVector<Number> > & output)
{
  this->discontinuous_value (p, time, output, libmesh_nullptr);
}



void MeshFunction::discontinuous_value (const Point & p,
                                        const Real,
                                        std::map<const Elem *, DenseVector<Number> > & output,
                                        const std::set<subdomain_id_type> * subdomain_ids)
{
  libmesh_assert (this->initialized());

  // clear the output map
  output.clear();

  // get the candidate elements
  std::set<const Elem *> candidate_element = this->find_elements(p,subdomain_ids);

  // loop through all candidates, if the set is empty this function will return an
  // empty map
  for (std::set<const Elem *>::const_iterator it = candidate_element.begin(); it != candidate_element.end(); ++it)
    {
      const Elem * element = (*it);

      const unsigned int dim = element->dim();

      // define a temporary vector to store all values
      DenseVector<Number> temp_output (cast_int<unsigned int>(this->_system_vars.size()));

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
      for (std::size_t index=0; index < this->_system_vars.size(); index++)
        {
          /*
           * the data for this variable
           */
          const unsigned int var = _system_vars[index];
          const FEType & fe_type = this->_dof_map.variable_type(var);

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

              for (std::size_t i=0; i<dof_indices.size(); i++)
                value += this->_vector(dof_indices[i]) * data.shape[i];

              temp_output(index) = value;
            }

          }

          // next variable
        }

      // Insert temp_output into output
      output[element] = temp_output;
    }
}



void MeshFunction::gradient (const Point & p,
                             const Real,
                             std::vector<Gradient> & output,
                             const std::set<subdomain_id_type> * subdomain_ids)
{
  libmesh_assert (this->initialized());

  const Elem * element = this->find_element(p,subdomain_ids);

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
        const unsigned int dim = element->dim();


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
        for (std::size_t index=0; index < this->_system_vars.size(); index++)
          {
            /*
             * the data for this variable
             */
            const unsigned int var = _system_vars[index];
            const FEType & fe_type = this->_dof_map.variable_type(var);

            UniquePtr<FEBase> point_fe (FEBase::build(dim, fe_type));
            const std::vector<std::vector<RealGradient> > & dphi = point_fe->get_dphi();
            point_fe->reinit(element, &point_list);

            // where the solution values for the var-th variable are stored
            std::vector<dof_id_type> dof_indices;
            this->_dof_map.dof_indices (element, dof_indices, var);

            // interpolate the solution
            Gradient grad(0.);

            for (std::size_t i=0; i<dof_indices.size(); i++)
              grad.add_scaled(dphi[i][0], this->_vector(dof_indices[i]));

            output[index] = grad;
          }
      }
    }
}


void MeshFunction::discontinuous_gradient (const Point & p,
                                           const Real time,
                                           std::map<const Elem *, std::vector<Gradient> > & output)
{
  this->discontinuous_gradient (p, time, output, libmesh_nullptr);
}



void MeshFunction::discontinuous_gradient (const Point & p,
                                           const Real,
                                           std::map<const Elem *, std::vector<Gradient> > & output,
                                           const std::set<subdomain_id_type> * subdomain_ids)
{
  libmesh_assert (this->initialized());

  // clear the output map
  output.clear();

  // get the candidate elements
  std::set<const Elem *> candidate_element = this->find_elements(p,subdomain_ids);

  // loop through all candidates, if the set is empty this function will return an
  // empty map
  for (std::set<const Elem *>::const_iterator it = candidate_element.begin(); it != candidate_element.end(); ++it)
    {
      const Elem * element = (*it);

      const unsigned int dim = element->dim();

      // define a temporary vector to store all values
      std::vector<Gradient> temp_output (cast_int<unsigned int>(this->_system_vars.size()));

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
      std::vector<Point> point_list (1, mapped_point);
      for (std::size_t index = 0 ; index < this->_system_vars.size(); ++index)
        {
          /*
           * the data for this variable
           */
          const unsigned int var = _system_vars[index];
          const FEType & fe_type = this->_dof_map.variable_type(var);

          UniquePtr<FEBase> point_fe (FEBase::build(dim, fe_type));
          const std::vector<std::vector<RealGradient> > & dphi = point_fe->get_dphi();
          point_fe->reinit(element, &point_list);

          // where the solution values for the var-th variable are stored
          std::vector<dof_id_type> dof_indices;
          this->_dof_map.dof_indices (element, dof_indices, var);

          Gradient grad(0.);

          for (std::size_t i = 0; i < dof_indices.size(); ++i)
            grad.add_scaled(dphi[i][0], this->_vector(dof_indices[i]));

          temp_output[index] = grad;

          // next variable
        }

      // Insert temp_output into output
      output[element] = temp_output;
    }
}



#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
void MeshFunction::hessian (const Point & p,
                            const Real,
                            std::vector<Tensor> & output,
                            const std::set<subdomain_id_type> * subdomain_ids)
{
  libmesh_assert (this->initialized());

  const Elem * element = this->find_element(p,subdomain_ids);

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
        const unsigned int dim = element->dim();


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
        for (std::size_t index=0; index < this->_system_vars.size(); index++)
          {
            /*
             * the data for this variable
             */
            const unsigned int var = _system_vars[index];
            const FEType & fe_type = this->_dof_map.variable_type(var);

            UniquePtr<FEBase> point_fe (FEBase::build(dim, fe_type));
            const std::vector<std::vector<RealTensor> > & d2phi =
              point_fe->get_d2phi();
            point_fe->reinit(element, &point_list);

            // where the solution values for the var-th variable are stored
            std::vector<dof_id_type> dof_indices;
            this->_dof_map.dof_indices (element, dof_indices, var);

            // interpolate the solution
            Tensor hess;

            for (std::size_t i=0; i<dof_indices.size(); i++)
              hess.add_scaled(d2phi[i][0], this->_vector(dof_indices[i]));

            output[index] = hess;
          }
      }
    }
}
#endif

const Elem * MeshFunction::find_element(const Point & p,
                                        const std::set<subdomain_id_type> * subdomain_ids) const
{
  /* Ensure that in the case of a master mesh function, the
     out-of-mesh mode is enabled either for both or for none.  This is
     important because the out-of-mesh mode is also communicated to
     the point locator.  Since this is time consuming, enable it only
     in debug mode.  */
#ifdef DEBUG
  if (this->_master != libmesh_nullptr)
    {
      const MeshFunction * master =
        cast_ptr<const MeshFunction *>(this->_master);
      if(_out_of_mesh_mode!=master->_out_of_mesh_mode)
        libmesh_error_msg("ERROR: If you use out-of-mesh-mode in connection with master mesh " \
                          << "functions, you must enable out-of-mesh mode for both the master and the slave mesh function.");
    }
#endif

  // locate the point in the other mesh
  const Elem * element = this->_point_locator->operator()(p,subdomain_ids);

  // If we have an element, but it's not a local element, then we
  // either need to have a serialized vector or we need to find a
  // local element sharing the same point.
  if (element &&
      (element->processor_id() != this->processor_id()) &&
      _vector.type() != SERIAL)
    {
      // look for a local element containing the point
      std::set<const Elem *> point_neighbors;
      element->find_point_neighbors(p, point_neighbors);
      element = libmesh_nullptr;
      std::set<const Elem *>::const_iterator       it  = point_neighbors.begin();
      const std::set<const Elem *>::const_iterator end = point_neighbors.end();
      for (; it != end; ++it)
        {
          const Elem * elem = *it;
          if (elem->processor_id() == this->processor_id())
            {
              element = elem;
              break;
            }
        }
    }

  return element;
}

std::set<const Elem *> MeshFunction::find_elements(const Point & p,
                                                   const std::set<subdomain_id_type> * subdomain_ids) const
{
  /* Ensure that in the case of a master mesh function, the
     out-of-mesh mode is enabled either for both or for none.  This is
     important because the out-of-mesh mode is also communicated to
     the point locator.  Since this is time consuming, enable it only
     in debug mode.  */
#ifdef DEBUG
  if (this->_master != libmesh_nullptr)
    {
      const MeshFunction * master =
        cast_ptr<const MeshFunction *>(this->_master);
      if(_out_of_mesh_mode!=master->_out_of_mesh_mode)
        libmesh_error_msg("ERROR: If you use out-of-mesh-mode in connection with master mesh " \
                          << "functions, you must enable out-of-mesh mode for both the master and the slave mesh function.");
    }
#endif

  // locate the point in the other mesh
  std::set<const Elem *> candidate_elements;
  std::set<const Elem *> final_candidate_elements;
  this->_point_locator->operator()(p,candidate_elements,subdomain_ids);
  for (std::set<const Elem *>::const_iterator elem_it = candidate_elements.begin(); elem_it != candidate_elements.end(); ++elem_it)
    {
      const Elem * element = (*elem_it);
      // If we have an element, but it's not a local element, then we
      // either need to have a serialized vector or we need to find a
      // local element sharing the same point.
      if (element &&
          (element->processor_id() != this->processor_id()) &&
          _vector.type() != SERIAL)
        {
          // look for a local element containing the point
          std::set<const Elem *> point_neighbors;
          element->find_point_neighbors(p, point_neighbors);
          std::set<const Elem *>::const_iterator       it  = point_neighbors.begin();
          const std::set<const Elem *>::const_iterator end = point_neighbors.end();
          for (; it != end; ++it)
            {
              const Elem * elem = *it;
              if (elem->processor_id() == this->processor_id())
                {
                  final_candidate_elements.insert(elem);
                  break;
                }
            }
        }
      else
        final_candidate_elements.insert(element);
    }

  return final_candidate_elements;
}

const PointLocatorBase & MeshFunction::get_point_locator (void) const
{
  libmesh_assert (this->initialized());
  return *_point_locator;
}

void MeshFunction::enable_out_of_mesh_mode(const DenseVector<Number> & value)
{
  libmesh_assert (this->initialized());
  _point_locator->enable_out_of_mesh_mode();
  _out_of_mesh_mode = true;
  _out_of_mesh_value = value;
}

void MeshFunction::enable_out_of_mesh_mode(const Number & value)
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

void MeshFunction::set_point_locator_tolerance(Real tol)
{
  _point_locator->set_close_to_point_tol(tol);
}

void MeshFunction::unset_point_locator_tolerance()
{
  _point_locator->unset_close_to_point_tol();
}

} // namespace libMesh
