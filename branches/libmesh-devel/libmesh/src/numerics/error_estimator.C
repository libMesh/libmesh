// $Id: error_estimator.C,v 1.1.2.1 2003-05-14 22:29:35 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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
#include <algorithm>

// Local Includes
#include "error_estimator.h"
#include "equation_systems.h"
#include "system_base.h"
#include "dof_map.h"
#include "mesh.h"
#include "fe.h"




//-----------------------------------------------------------------
// ErrorEstimator implementations
void ErrorEstimator::flux_jump (const EquationSystems& es,
				const std::string& name,
				std::vector<Real>& error_per_cell)
{
  /*

  - e & f are global element ids
  
  Case (1.) Elements are at the same level, e<f
            Compute the flux jump on the face and
	    add it as a contribution to error_per_cell[e]
	    and error_per_cell[f]
  
                   ----------------------
		  |           |          |
		  |           |    f     |
		  |           |          |
		  |    e      |---> n    | 
		  |           |          |
		  |           |          |
                   ----------------------


   Case (2.) The neighbor is at a higher level.
             Compute the flux jump on e's face and
	     add it as a contribution to error_per_cell[e]
	     and error_per_cell[f]

                   ----------------------
		  |     |     |          |
		  |     |  e  |---> n    |
		  |     |     |          |
		  |-----------|    f     | 
		  |     |     |          |
		  |     |     |          |
		  |     |     |          |
                   ----------------------
*/
   
  // The current mesh
  const Mesh& mesh = es.get_mesh();

  // The dimensionality of the mesh
  const unsigned int dim = mesh.mesh_dimension();
  
  // The System object to estimate the error for
  const SystemBase& system = es(name);

  // The DofMap for this system
  const DofMap& dof_map = system.get_dof_map();

  // Resize the error_per_cell vector to be
  // the number of elements, initialize it to 0.
  error_per_cell.resize (mesh.n_elem());
  std::fill (error_per_cell.begin(), error_per_cell.end(), 0.);

  

  
  // Loop over all the variables in the system
  for (unsigned int var=0; var<system.n_vars(); var++)
    {
      const FEType& fe_type = dof_map.variable_type (var);

      // Finite element objects for the same face from
      // different sides
      AutoPtr<FEBase> fe_face_e (FEBase::build (dim, fe_type));
      AutoPtr<FEBase> fe_face_f (FEBase::build (dim, fe_type));

      // By convention we will always do the integration
      // on the face of element e.  Get its Jacobian values, etc..
      const std::vector<Real>&  JxW_face    = fe_face_e->get_JxW();
      const std::vector<Point>& qface_point = fe_face_e->get_xyz();
      
      
    }
}
