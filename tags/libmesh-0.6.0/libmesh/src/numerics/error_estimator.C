// $Id: error_estimator.C,v 1.25 2007-04-11 23:06:42 roystgnr Exp $

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


// Local include files
#include "libmesh_common.h"
#include "error_estimator.h"
#include "error_vector.h"
#include "equation_systems.h"




// ErrorEstimator functions

#ifdef HAVE_MPI
void ErrorEstimator::reduce_error (std::vector<ErrorVectorReal>& error_per_cell) const
{
  // Each processor has now computed the error contribuions
  // for its local elements.  We need to sum the vector to
  // recover the error for each element.  Note that we only
  // need to sum if we are running on multiple processors
  if (libMesh::n_processors() > 1)
    {
      // Allreduce requires 2 buffers.  Copy the
      // error_per_cell vector into the epc vector
      std::vector<ErrorVectorReal> epc (error_per_cell);
      
      MPI_Allreduce (&epc[0], &error_per_cell[0],
		     error_per_cell.size(),
		     MPI_ERRORVECTORREAL, MPI_SUM, libMesh::COMM_WORLD);
    }  
}
#else
void ErrorEstimator::reduce_error (std::vector<ErrorVectorReal>&) const
{
}
#endif



void ErrorEstimator::convert_component_mask_to_scale()
{
  if (!component_mask.empty())
    {
      // component_mask has been replaced by component_scale,
      // and will be removed in future libMesh versions
      deprecated();

      component_scale.resize(component_mask.size(),0.0);
      for (unsigned int i=0; i != component_mask.size(); ++i)
        {
          if (component_mask[i])
            component_scale[i] = 1.0;
        }
    }
}



void ErrorEstimator::estimate_errors(const EquationSystems& equation_systems,
                                     ErrorVector& error_per_cell,
                                     std::map<const System*, std::vector<float> >& component_scales,
                                     bool estimate_parent_error)
{
  // This is a brand-new function; if you're using it you should
  // already have stopped using component_mask
  assert(component_mask.empty());

  std::vector<float> old_component_scale = this->component_scale;

  // Sum the error values from each system
  for (unsigned int s = 0; s != equation_systems.n_systems(); ++s)
    {
      ErrorVector system_error_per_cell;
      const System &sys = equation_systems.get_system(s);
      if (component_scales.find(&sys) == component_scales.end())
        this->component_scale = old_component_scale;
      else
        this->component_scale = component_scales[&sys];

      this->estimate_error(sys, system_error_per_cell, estimate_parent_error);

      if (s)
        {
          assert(error_per_cell.size() == system_error_per_cell.size());
          for (unsigned int i=0; i != error_per_cell.size(); ++i)
            error_per_cell[i] += system_error_per_cell[i];
        }
      else
        error_per_cell = system_error_per_cell;
    }

  // Restore our old state before returning
  this->component_scale = old_component_scale;
}



/**
 * FIXME: This is a default implementation - derived classes should
 * reimplement it for efficiency.
 */
void ErrorEstimator::estimate_errors(const EquationSystems& equation_systems,
                                     ErrorMap& errors_per_cell,
                                     bool estimate_parent_error)
{
  // This is a brand-new function; if you're using it you should
  // already have stopped using component_mask
  assert(component_mask.empty());

  std::vector<float> old_component_scale = this->component_scale;

  // Find the requested error values from each system
  for (unsigned int s = 0; s != equation_systems.n_systems(); ++s)
    {
      const System &sys = equation_systems.get_system(s);

      unsigned int n_vars = sys.n_vars();

      // Calculate error in only one variable
      this->component_scale.clear();
      this->component_scale.resize(n_vars, 0.0);

      for (unsigned int v = 0; v != n_vars; ++v)
        {
          // Only fill in ErrorVectors the user asks for
          if (errors_per_cell.find(std::make_pair(&sys, v)) ==
              errors_per_cell.end())
            continue;

          this->component_scale[v] = 1.0;

          this->estimate_error
            (sys, *errors_per_cell[std::make_pair(&sys, v)],
             estimate_parent_error);

          this->component_scale[v] = 0.0;
        }
    }

  // Restore our old state before returning
  this->component_scale = old_component_scale;
}
