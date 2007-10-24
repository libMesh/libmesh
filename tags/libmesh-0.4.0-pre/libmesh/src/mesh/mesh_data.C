// $Id: mesh_data.C,v 1.1 2003-05-15 19:43:34 ddreyer Exp $

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
#include "mesh_data.h"
#include "mesh_base.h"



//------------------------------------------------------
// MeshData functions
MeshData::MeshData(const MeshBase& m) :
  _mesh               (m),
  _data_descriptor    (""),
  _n_val_per_node     (0),
  _n_val_per_elem     (0),
  _node_id_map_closed (false),
  _node_data_closed   (false),
  _elem_id_map_closed (false),
  _elem_data_closed   (false),
  _active             (false)
{
}




MeshData::~MeshData()
{
  clear();
}





void MeshData::activate (const std::string& descriptor)
{
  _active = true;
  _data_descriptor = descriptor;
}






void MeshData::clear()
{
  _data_descriptor    = "";
  _node_data.clear();
  _elem_data.clear();
  _node_data_closed   = false;
  _elem_data_closed   = false;
}





void  MeshData::slim (const bool node_id_map,
		      const bool elem_id_map)
{
  assert (_active);

  if (node_id_map)
    {
      // dumb check
      assert (_node_id_map_closed);

      _node_id_map_closed = false;
      _node_id.clear();
      _id_node.clear();
    }

  if (elem_id_map)
    {
      // dumb check
      assert (_elem_id_map_closed);

      _elem_id_map_closed = false;
      _elem_id.clear();
      _id_elem.clear();
    }
}





void MeshData::close_foreign_id_maps ()
{
  if (_active)
    {
      assert (!_elem_id.empty());
      assert (!_id_elem.empty());
      assert (!_node_id.empty());
      assert (!_id_node.empty());

      _elem_id_map_closed = true;
      _node_id_map_closed = true;
    }
}





void MeshData::read (const std::string& name)
{
  START_LOG("read()", "MeshData");

  assert (_active);

  // the id maps have to be closed before reading
  assert (_elem_id_map_closed && _node_id_map_closed);

  // Read the file based on extension
  {
    if (name.rfind(".xta") < name.size())
      read_xdr (name, true);

    else if (name.rfind(".xtr")  < name.size())
      read_xdr (name, false);

    else if (name.rfind(".unv") < name.size())
      read_unv (name);

    else
      {
	std::cerr << " ERROR: Unrecognized file extension: " << name
		  << "\n   I understand the following:\n\n"
		  << "     *.xta  -- Internal ASCII data format\n"
		  << "     *.xtr  -- Internal binary data format\n"
		  << "     *.unv  -- I-deas format\n"
		  << std::endl;
	error();

      }    
  }
  STOP_LOG("read()", "MeshData");
}






void MeshData::write (const std::string& name)
{
  START_LOG("write()", "MeshData");

  assert (_active);

  // the id maps have to be closed before reading
  assert (_elem_id_map_closed && _node_id_map_closed);
  
  // Read the file based on extension
  {
    if (name.rfind(".xta") < name.size())
      write_xdr (name, true);

    else if (name.rfind(".xtr")  < name.size())
      write_xdr (name, false);

    else if (name.rfind(".unv") < name.size())
      write_unv (name);

    else
      {
	std::cerr << " ERROR: Unrecognized file extension: " << name
		  << "\n   I understand the following:\n\n"
		  << "     *.xta  -- Internal ASCII data format\n"
		  << "     *.xtr  -- Internal binary data format\n"
		  << "     *.unv  -- I-deas format\n"
		  << std::endl;
	error();

      }    
  }
  STOP_LOG("write()", "MeshData");
}







const Node* MeshData::foreign_id_to_node (const unsigned int fid) const
{
  assert (_active);
  assert (_node_id_map_closed);

  std::map<const unsigned int,
           const Node*>::const_iterator pos = _id_node.find(fid);

  if (pos == _id_node.end())
    {
      std::cerr << "ERROR: Have no Node* associated with the foreign id = "
		<< fid
		<< std::endl;
      error();
      return NULL;
    }
  else
    return (*pos).second;
}





const unsigned int MeshData::node_to_foreign_id (const Node* n) const
{
  assert (_active);
  assert (_node_id_map_closed);

  assert (n != NULL);

  /*
   * look it up in the map
   */
  std::map<const Node*,
	   const unsigned int>::const_iterator pos = _node_id.find(n);

  if (pos == _node_id.end())
    {
      std::cerr << "ERROR: No foreign id stored for the node "
		<< "with the libMesh id = "
		<< n->id()
		<< std::endl;
      error();
      return 0;
    }
  else
    return (*pos).second;
}








const Elem* MeshData::foreign_id_to_elem (const unsigned int fid) const
{
  assert (_active);
  assert (_elem_id_map_closed);

  std::map<const unsigned int,
           const Elem*>::const_iterator pos = _id_elem.find(fid);

  if (pos == _id_elem.end())
    {
      std::cerr << "ERROR: Have no Elem* associated with the foreign id = "
		<< fid
		<< std::endl;
      error();
      return NULL;
    }
  else
    return (*pos).second;
}





const unsigned int MeshData::elem_to_foreign_id (const Elem* e) const
{
  assert (_active);
  assert (_elem_id_map_closed);

  assert (e != NULL);

  /*
   * look it up in the map
   */
  std::map<const Elem*,
	   const unsigned int>::const_iterator pos = _elem_id.find(e);

  if (pos == _elem_id.end())
    {
      std::cerr << "ERROR: No foreign id stored for the element "
		<< "with the libMesh id = "
		<< e->id()
		<< std::endl;
      error();
      return 0;
    }
  else
    return (*pos).second;
}

