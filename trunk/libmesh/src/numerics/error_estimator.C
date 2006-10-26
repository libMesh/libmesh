// $Id: error_estimator.C,v 1.22 2006-10-26 04:22:48 roystgnr Exp $

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
