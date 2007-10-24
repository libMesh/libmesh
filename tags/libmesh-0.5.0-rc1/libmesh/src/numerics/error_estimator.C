// $Id: error_estimator.C,v 1.19 2005-05-11 23:11:59 benkirk Exp $

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




// ErrorEstimator functions
void ErrorEstimator::reduce_error (std::vector<float>& error_per_cell) const
{
  // Each processor has now computed the error contribuions
  // for its local elements.  We need to sum the vector to
  // recover the error for each element.  Note that we only
  // need to sum if we are running on multiple processors
#ifdef HAVE_MPI
  if (libMesh::n_processors() > 1)
    {
      // Allreduce requires 2 buffers.  Copy the
      // error_per_cell vector into the epc vector
      std::vector<float> epc (error_per_cell);
      
      MPI_Allreduce (&epc[0], &error_per_cell[0],
		     error_per_cell.size(),
		     MPI_FLOAT, MPI_SUM, libMesh::COMM_WORLD);
    }  
#endif
}
