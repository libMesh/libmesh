// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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
#include <sstream>

// Local includes
#include "mesh_data.h"
#include "mesh_base.h"
#include "libmesh_logging.h"
#include "elem.h"


//------------------------------------------------------
// MeshData functions
MeshData::MeshData(const MeshBase& m) :
  _mesh               (m),
  _data_descriptor    (""),
  _node_id_map_closed (false),
  _node_data_closed   (false),
  _elem_id_map_closed (false),
  _elem_data_closed   (false),
  _active             (false),
  _compatibility_mode (false),
  _unv_header         (NULL)
{
}




MeshData::~MeshData()
{
  clear();
}





void MeshData::activate (const std::string& descriptor)
{
#ifdef DEBUG
  if (_compatibility_mode)
      std::cerr << "WARNING: MeshData was in compatibility mode, now being activated."
		<< std::endl;
#endif

  _compatibility_mode = false;
  _active = true;
  _data_descriptor = descriptor;
}





void MeshData::enable_compatibility_mode (const std::string& descriptor)
{
  if (!_active)
    {
      _compatibility_mode = true;
      _active = false;
      // do as if the id maps are already closed
      _node_id_map_closed = true;
      _elem_id_map_closed = true;
      _data_descriptor = descriptor;
      // we can safely clear the id maps
      _node_id.clear();
      _id_node.clear();
      _elem_id.clear();
      _id_elem.clear();
    }
#ifdef DEBUG
  else
      std::cerr << "WARNING: MeshData was in compatibility mode, now being activated."
		<< std::endl;
#endif
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
  if (this->active())
    {

      if (node_id_map)
        {
	  // dumb check
	  libmesh_assert (_node_id_map_closed);

	  _node_id_map_closed = false;
	  _node_id.clear();
	  _id_node.clear();
	}

      if (elem_id_map)
        {
	  // dumb check
	  libmesh_assert (_elem_id_map_closed);

	  _elem_id_map_closed = false;
	  _elem_id.clear();
	  _id_elem.clear();
	}
    }

#ifdef DEBUG
  else if (this->compatibility_mode())
    {
      std::cerr << "WARNING: No need for MeshData::slim() in compatibility mode." << std::endl;
    }
#endif
}




void MeshData::translate (const MeshBase& out_mesh,
			  std::vector<Number>& values,
			  std::vector<std::string>& names) const
{
  libmesh_assert (_active || _compatibility_mode);

  START_LOG("translate()", "MeshData");

  const unsigned int n_comp = this->n_val_per_node();

  // transfer our nodal data to a vector
  // that may be written concurrently
  // with the \p out_mesh.
  {
    // reserve memory for the nodal data
    values.reserve(n_comp*out_mesh.n_nodes());
    
    // iterate over the mesh's nodes
    MeshBase::const_node_iterator       nodes_it  = out_mesh.nodes_begin();
    const MeshBase::const_node_iterator nodes_end = out_mesh.nodes_end();

    // Do not use the \p get_data() method, but the operator()
    // method, since this returns by default a zero value,
    // when there is no nodal data.
    for (; nodes_it != nodes_end; ++nodes_it)
      {
	const Node* node = *nodes_it;
	
	for (unsigned int c= 0; c<n_comp; c++)
	  values.push_back(this->operator()(node, c));
      }
  }


  
  // Now we have the data, nicely stored in \p values.
  // It remains to give names to the data, then write to
  // file.
  {
    names.reserve(n_comp);
    
    // this naming scheme only works up to n_comp=100
    // (at least for gmv-accepted variable names)
    libmesh_assert(n_comp < 100);

    for (unsigned int n=0; n<n_comp; n++)
      {
	std::ostringstream name_buf;
	name_buf << "bc_" << n;
	names.push_back(name_buf.str());
      }
  }

  STOP_LOG("translate()", "MeshData");
}




void MeshData::close_foreign_id_maps ()
{
  if (_active)
    {
      libmesh_assert (!_elem_id.empty());
      libmesh_assert (!_id_elem.empty());
      libmesh_assert (!_node_id.empty());
      libmesh_assert (!_id_node.empty());

      _elem_id_map_closed = true;
      _node_id_map_closed = true;
    }
}





void MeshData::read (const std::string& name)
{
  START_LOG("read()", "MeshData");

  libmesh_assert (_active || _compatibility_mode);

  // the id maps have to be closed before reading
  // (note that in compatibility mode these are also true)
  libmesh_assert (_elem_id_map_closed && _node_id_map_closed);

#ifdef DEBUG
  if (this->compatibility_mode())
      std::cerr << "WARNING: MeshData in compatibility mode, node and element ids" << std::endl
		<< "         stored in file may be totally different from libMesh ids!" << std::endl;
#endif

  // Read the file based on extension.  We let all processors read the
  // data because it would be inaccurate to let only one processor
  // have it and we're too lazy to code up a proper parallel read or
  // read+broadcast right now.

  //if (libMesh::processor_id() == 0)
  //  {
      if (name.rfind(".xta") < name.size())
	this->read_xdr (name, READ);
      
      else if (name.rfind(".xtr")  < name.size())
	this->read_xdr (name, DECODE);
      
      else if (name.rfind(".unv") < name.size())
	this->read_unv (name);
      
      else if ((name.rfind(".node") < name.size()) ||
	       (name.rfind(".ele") < name.size()))
	this->read_tetgen (name);
      
      else
	{
	  std::cerr << " ERROR: Unrecognized file extension: " << name
		    << "\n   I understand the following:\n\n"
		    << "     *.xta  -- Internal ASCII data format\n"
		    << "     *.xtr  -- Internal binary data format\n"
		    << "     *.unv  -- I-deas format\n"
		    << std::endl;
	  libmesh_error();
	  
	}    
    //}
  STOP_LOG("read()", "MeshData");
}






void MeshData::write (const std::string& name)
{
  START_LOG("write()", "MeshData");

  libmesh_assert (_active || _compatibility_mode);

  // the id maps have to be closed before writing
  // (note that in compatibility mode these are also true)
  libmesh_assert (_elem_id_map_closed && _node_id_map_closed);
  
#ifdef DEBUG
  if (this->compatibility_mode())
      std::cerr << "WARNING: MeshData in compatibility mode.  Node and element ids" << std::endl
		<< "         written to file may differ from libMesh numbering" << std::endl
		<< "         next time this file is read!" << std::endl;
#endif

  // Read the file based on extension
  {
    if (name.rfind(".xta") < name.size())
      write_xdr (name, WRITE);

    else if (name.rfind(".xtr")  < name.size())
      write_xdr (name, ENCODE);

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
	libmesh_error();

      }    
  }
  STOP_LOG("write()", "MeshData");
}




std::string MeshData::get_info() const
{
  std::ostringstream out;

  if (this->active() || this->compatibility_mode())
    {
      out << " MeshData Information:\n";
      if (this->active())
	  out << "  object activated.\n";
      if (this->compatibility_mode())
	  out << "  object in compatibility mode.\n";
      if (this->_data_descriptor != "")
	  out << "  descriptor=" << this->_data_descriptor << '\n';
      if (this->elem_initialized())
	  out << "  Element associated data initialized.\n"
	      << "   n_val_per_elem()=" << this->n_val_per_elem() << '\n'
	      << "   n_elem_data()=" << this->n_elem_data() << '\n';
      if (this->node_initialized())
	  out << "  Node associated data initialized.\n"
	      << "   n_val_per_node()=" << this->n_val_per_node() << '\n'
	      << "   n_node_data()=" << this->n_node_data() << '\n';
    }
  else
      out << " MeshData neither active nor in compatibility mode.\n";

  return out.str();
}




void MeshData::print_info(std::ostream& os) const
{
  os << this->get_info()
     << std::endl;
}


std::ostream& operator << (std::ostream& os, const MeshData& m)
{
  m.print_info(os);
  return os;
}




const Node* MeshData::foreign_id_to_node (const unsigned int fid) const
{
  if (_active)
    {
      // when active, use our _id_node map
      libmesh_assert (_node_id_map_closed);

      std::map<unsigned int,
	       const Node*>::const_iterator pos = _id_node.find(fid);

      if (pos == _id_node.end())
        {
	  std::cerr << "ERROR: Have no Node* associated with the foreign id = "
		    << fid
		    << std::endl;
	  libmesh_error();
	  return NULL;
	}
      else
	  return pos->second;
    }
  else if (_compatibility_mode)
      // when only in compatibility mode, 
      // return the node stored in the MeshBase 
      // under its current id
      return this->_mesh.node_ptr(fid);

  // should never get here
  libmesh_error();
  return NULL;
}





unsigned int MeshData::node_to_foreign_id (const Node* n) const
{
  libmesh_assert (n != NULL);

  if (_active)
    {
      // when active, use our _node_id map
      libmesh_assert (_node_id_map_closed);

      // look it up in the map
      std::map<const Node*,
	       unsigned int>::const_iterator pos = _node_id.find(n);
      
      if (pos == _node_id.end())
        {
	  std::cerr << "ERROR: No foreign id stored for the node "
		    << "with the libMesh id = "
		    << n->id()
		    << std::endl;
	  libmesh_error();
	  return 0;
	}
      else
	  return pos->second;
    }
  else if (_compatibility_mode)
    // when only in compatibility mode, 
    // return libMesh's node id
    return n->id();

  // should never get here
  libmesh_error();
  return 0;
}








const Elem* MeshData::foreign_id_to_elem (const unsigned int fid) const
{
  if (_active)
    {
      // when active, use our _id_elem map
      libmesh_assert (_elem_id_map_closed);
      
      std::map<unsigned int,
	       const Elem*>::const_iterator pos = _id_elem.find(fid);
      
      if (pos == _id_elem.end())
        {
	  std::cerr << "ERROR: Have no Elem* associated with the foreign id = "
		    << fid
		    << std::endl;
	  libmesh_error();
	  return NULL;
	}
      else
	  return pos->second;
    }
  else if (_compatibility_mode)
    // when only in compatibility mode, 
    // return element using the libMesh id
    return this->_mesh.elem(fid);

  // should never get here
  libmesh_error();
  return NULL;
}





unsigned int MeshData::elem_to_foreign_id (const Elem* e) const
{
  libmesh_assert (e != NULL);

  if (_active)
    {
      // when active, use our _id_elem map
      libmesh_assert (_elem_id_map_closed);

             // look it up in the map
             std::map<const Elem*,
	       unsigned int>::const_iterator pos = _elem_id.find(e);

      if (pos == _elem_id.end())
        {
	  std::cerr << "ERROR: No foreign id stored for the element "
		    << "with the libMesh id = "
		    << e->id()
		    << std::endl;
	  libmesh_error();
	  return 0;
	}
      else
	  return pos->second;
    }
  else if (_compatibility_mode)
    // when only in compatibility mode, 
    // return libMesh's element id
    return e->id();

  // should never get here
  libmesh_error();
  return 0;
}







void MeshData::insert_node_data (std::map<const Node*,
				 std::vector<Number> >& nd,
				 const bool close_elem_data)
{
  libmesh_assert (this->_active || this->_compatibility_mode);
  // these are also true in compatibility mode
  libmesh_assert (this->_node_id_map_closed);

  if (this->_node_data_closed)
    {
      std::cerr << "ERROR: Nodal data already closed!  Use clear() first!"
		<< std::endl;
      libmesh_error();
    }

  libmesh_assert (this->_node_data.empty());

#ifdef DEBUG
  std::map<const Node*, 
           std::vector<Number> >::const_iterator nd_pos = nd.begin();
  std::map<const Node*, 
           std::vector<Number> >::const_iterator nd_end = nd.end();

  // Compare entity-by-entity that the
  // sizes of the std::vector's are identical.
  // For this, simply take the length of the 0th
  // entry as reference length, and compare this
  // with the length of the 1st, 2nd...
  libmesh_assert (nd_pos != nd_end);
  const unsigned int reference_length = (*nd_pos).second.size();

  // advance, so that we compare with the 1st
  ++nd_pos;

  for (; nd_pos != nd_end; ++nd_pos)
    if ( (*nd_pos).second.size() != reference_length) 
      {
	std::cerr << "ERROR: Size mismatch."
		  << std::endl;
	libmesh_error();
      }
#endif

  // copy over
  _node_data = nd;

  // we may freely trash the nd
  nd.clear();

  // close node data
  this->_node_data_closed = true;

  // if user wants to, then close elem data, too
  if (close_elem_data)
    {
      libmesh_assert((this->_elem_id_map_closed));
      this->_elem_data_closed = true;
    }
}





void MeshData::insert_elem_data (std::map<const Elem*,
				 std::vector<Number> >& ed,
				 const bool close_node_data)
{
  libmesh_assert (this->_active || this->_compatibility_mode);
  // these are also true in compatibility mode
  libmesh_assert (this->_elem_id_map_closed);

  if (this->_elem_data_closed)
    {
      std::cerr << "ERROR: Element data already closed!  Use clear() first!"
		<< std::endl;
      libmesh_error();
    }

  libmesh_assert (this->_elem_data.empty());

#ifdef DEBUG
  std::map<const Elem*, 
           std::vector<Number> >::const_iterator ed_pos = ed.begin();
  std::map<const Elem*, 
           std::vector<Number> >::const_iterator ed_end = ed.end();

  // Compare entity-by-entity that the
  // sizes of the std::vector's are identical.
  const unsigned int reference_length = (*ed_pos).second.size();
  ++ed_pos;

  for (; ed_pos != ed_end; ++ed_pos)
    if ( (*ed_pos).second.size() != reference_length) 
      {
	std::cerr << "ERROR: Size mismatch."
		  << std::endl;
	libmesh_error();
      }
#endif

  // copy over
  _elem_data = ed;

  // we may freely trash the ed
  ed.clear();

  // close elem data
  this->_elem_data_closed = true;

  // if user wants to, then close node data, too
  if (close_node_data)
    {
      libmesh_assert((this->_node_id_map_closed));
      this->_node_data_closed = true;
    }
}





unsigned int MeshData::n_val_per_node () const
{
  libmesh_assert (this->_active || this->_compatibility_mode);
  libmesh_assert (this->_node_data_closed);

  if (!this->_node_data.empty())
    {
      std::map<const Node*, 
	       std::vector<Number> >::const_iterator pos = _node_data.begin();
      libmesh_assert (pos != _node_data.end());
      return (pos->second.size());
    }
  else
      return 0;
}




unsigned int MeshData::n_node_data () const
{
  libmesh_assert (this->_active || this->_compatibility_mode);
  libmesh_assert (this->_node_data_closed);

  return this->_node_data.size();
}




unsigned int MeshData::n_val_per_elem () const
{
  libmesh_assert (this->_active || this->_compatibility_mode);
  libmesh_assert (this->_elem_data_closed);

  if (!_elem_data.empty())
    {
      std::map<const Elem*, 
	       std::vector<Number> >::const_iterator pos = _elem_data.begin();
      libmesh_assert (pos != _elem_data.end());
      return (pos->second.size());
    }
  else
      return 0;
}




unsigned int MeshData::n_elem_data () const
{
  libmesh_assert (this->_active || this->_compatibility_mode);
  libmesh_assert (this->_elem_data_closed);

  return _elem_data.size();
}




void MeshData::assign (const MeshData& omd)
{
  this->_data_descriptor    = omd._data_descriptor;
  this->_node_id_map_closed = omd._node_id_map_closed;
  this->_node_data_closed   = omd._node_data_closed;

  // we have to be able to modify our elem id maps
  libmesh_assert (!this->_elem_id_map_closed);

  this->_elem_data_closed   = omd._elem_data_closed;
  this->_active             = omd._active;
  this->_compatibility_mode = omd._compatibility_mode;

  // this is ok because we do not manage the UnvHeader
  // in terms of memory, but only hold a pointer to it...
  this->_unv_header         = omd._unv_header;

  // Now copy the foreign id maps -- but only for the 
  // nodes.  The nodes of the boundary mesh are actually
  // nodes of the volume mesh.
  this->_node_id = omd._node_id;
  this->_id_node = omd._id_node;

  // The element vector of the boundary mesh contains elements
  // that are new, and there _cannot_ be any associated
  // foreign id in the maps.  Therefore, fill the maps with
  // the libMesh id's.  But only when the other MeshData
  // has element ids.
  if ((this->_active) && (omd._elem_id.size() != 0))
    {

      MeshBase::const_element_iterator       elem_it  = _mesh.elements_begin();
      const MeshBase::const_element_iterator elem_end = _mesh.elements_end();

      for (; elem_it != elem_end; ++elem_it)
        {
	  const Elem* elem = *elem_it;  
	  this->add_foreign_elem_id(elem, elem->id());
	}
    }

  // now we can safely assign omd's value
  this->_elem_id_map_closed   = omd._elem_id_map_closed;
  

  // and finally the node- and element-associated data
  this->_node_data = omd._node_data;
  this->_elem_data = omd._elem_data;
}
