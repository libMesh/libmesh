// $Id: boundary_data.C,v 1.1 2003-05-14 11:54:37 ddreyer Exp $

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

// Local includes
#include "boundary_data.h"
#include "mesh_base.h"



//------------------------------------------------------
// BoundaryData functions
BoundaryData::BoundaryData(const MeshBase& m) :
  _mesh               (m),
  _node_id_map_closed (false),
  _node_initialized   (false),
  _elem_id_map_closed (false),
  _elem_initialized   (false)
{
}




BoundaryData::~BoundaryData()
{
  clear();
}




void BoundaryData::clear()
{
  _node_data.clear();
  _elem_data.clear();
  _node_initialized = false;
  _elem_initialized = false;
}





void  BoundaryData::slim ()
{
  assert (_node_id_map_closed);
  assert (_elem_id_map_closed);

  _node_id.clear();
  _elem_id.clear();

  _id_node.clear();
  _id_elem.clear();
}





void BoundaryData::close_node_map ()
{
  assert (!_node_id.empty());
  _node_id_map_closed = true;
}




void BoundaryData::close_elem_map ()
{
  assert (!_elem_id.empty());
  _elem_id_map_closed = true;
}




void BoundaryData::read (const std::string& name)
{
  START_LOG("read()", "BoundaryData");

  // at least one of the maps should be closed.
  assert (_node_id_map_closed || _elem_id_map_closed);


  // Read the file based on extension
  {
    if (name.rfind(".xda") < name.size())
      read_xdr (name);

    else if (name.rfind(".xdr")  < name.size())
      read_xdr_binary (name);

    else if (name.rfind(".unv") < name.size())
      read_unv (name);

    else
      {
	std::cerr << " ERROR: Unrecognized file extension: " << name
		  << "\n   I understand the following:\n\n"
		  << "     *.xda  -- Internal ASCII format\n"
		  << "     *.xdr  -- Internal binary format,\n"
		  << "               compatible with XdrMGF\n"
		  << "     *.unv  -- I-deas format\n"
		  << std::endl;
	error();

      }    
  }
  STOP_LOG("read()", "BoundaryData");
}






void BoundaryData::write (const std::string& name)
{
  START_LOG("write()", "BoundaryData");

  assert ((this->_node_initialized) || (this->_elem_initialized));
  
  // Read the file based on extension
  {
    if (name.rfind(".xda") < name.size())
      write_xdr (name);

    else if (name.rfind(".xdr")  < name.size())
      write_xdr_binary (name);

    else if (name.rfind(".unv") < name.size())
      write_unv (name);

    else
      {
	std::cerr << " ERROR: Unrecognized file extension: " << name
		  << "\n   I understand the following:\n\n"
		  << "     *.xda  -- Internal ASCII format\n"
		  << "     *.xdr  -- Internal binary format,\n"
		  << "               compatible with XdrMGF\n"
		  << "     *.unv  -- I-deas format\n"
		  << std::endl;
	error();

      }    
  }
  STOP_LOG("write()", "BoundaryData");
}

