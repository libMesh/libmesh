// $Id: mesh_data_unv_support.C,v 1.2 2003-05-15 23:34:35 benkirk Exp $

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
#include <fstream>

// Local includes
#include "mesh_data.h"
#include "mesh_base.h"



//------------------------------------------------------
// MeshData UNV support functions
void MeshData::read_unv(const std::string& /* name */)
{
  /*
   * When reading data, make sure the id maps are ok
   */
  assert (_node_id_map_closed);
  assert (_elem_id_map_closed);



  std::cerr << "ERROR: not yet implemented." << std::endl;
  error();




  /*
   * finished reading.  Ready for use
   */
  this->_node_data_closed = true;
  this->_elem_data_closed = true;
}




void MeshData::write_unv (const std::string& /* name */)
{
  /*
   * make sure the id maps are ready
   * and that we have data to write
   */
  assert (_node_id_map_closed);
  assert (_elem_id_map_closed);

  assert (_node_data_closed);
  assert (_elem_data_closed);
  



  std::cerr << "ERROR: not yet implemented." << std::endl;
  error();
}

