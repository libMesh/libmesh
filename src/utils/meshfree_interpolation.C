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
#include "libmesh/meshfree_interpolation.h"



namespace libMesh
{
  
  void MeshfreeInterpolation::clear ()
  {
    _names.clear();
    _src_pts.clear();
    _src_vals.clear();
  }



  void MeshfreeInterpolation::add_field_data (const std::vector<std::string> &field_names,
					      const std::vector<Point>  &pts,
					      const std::vector<Number> &vals)
  {
    libmesh_here();
    libmesh_assert_equal_to (field_names.size()*pts.size(), vals.size());

    // If we already have field variables, we assume we are appending.
    // that means the names and ordering better be identical!
    if (!_names.empty())
      {
	if (_names.size() != field_names.size())
	  {
	    libMesh::err << "ERROR:  when adding field data to an existing list the\n"
			 << "varaible list must be the same!\n";
	    libmesh_error();	      	      
	  }
	for (unsigned int v=0; v<_names.size(); v++)
	  if (_names[v] != field_names[v])
	    {
	      libMesh::err << "ERROR:  when adding field data to an existing list the\n"
			   << "varaible list must be the same!\n";
	      libmesh_error();	      	      
	    }
	    
      }
    
    // otherwise copy the names
    else
      _names = field_names;

    // append the data
  }

  
}

