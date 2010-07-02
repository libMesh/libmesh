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
#include <fstream>

// Local includes
#include "mesh_data.h"
#include "mesh_base.h"
#include "xdr_cxx.h"
#include "elem.h"

namespace libMesh
{


//------------------------------------------------------
// MeshData functions
void MeshData::read_xdr (const std::string& name,
			 const XdrMODE mode)
{
  /**
   * This code implements the output of the MeshData
   * object in XDR format.  This warrants some documentation.
   * The output consists of 8 sections:
   *
   *   1.) The name of the data stored, if provided (string)
   *
   *   2.) A switch whether real or complex data is stored (string)
   *
   *   3.) The number of nodes for which values are stored
   *       (unsigned int)
   *
   *   4.) The number of elements for which values are stored
   *       (unsigned int)
   *
   *   for each node
   *     
   *     5.) The foreign node id (unsigned int)
   *     
   *     6.) The actual values (vector of real/complex)
   * 
   *   end node loop
   *
   *   for each element
   *     
   *     7.) The foreign element id (unsigned int)
   *     
   *     8.) The actual values (vector of real/complex)
   * 
   *   end node loop
   *
   * Note that the actual IO is handled through the Xdr class 
   * (to be renamed later?) which provides a uniform interface to 
   * both the XDR (eXternal Data Representation) interface and standard
   * ASCII output.  Thus this one section of code will write XDR or ASCII
   * files with no changes.
   */ 


  // we should better be active or in compatibility mode
  libmesh_assert (_active || _compatibility_mode);


  // make sure the id maps are ready
  libmesh_assert (_node_id_map_closed);
  libmesh_assert (_elem_id_map_closed);


  /**
   * clear the data, but keep the id maps
   */
  this->clear();


  Xdr io(name, mode);


  /*
   * all processors read the data in the same format,
   * but only the processor that owns the element stores
   * element-associated data.  For nodes, i haven't come
   * up with such asmart idea, yet... :-P
   */
  const unsigned int proc_id = _mesh.processor_id();



  /**
   * 1.)
   *
   * Read the descriptive name
   */
  {
    std::string desc = "";
    io.data (desc);
    this->_data_descriptor = desc;
  }



  /**
   * 2.)
   *
   * Read: either real or complex
   */
  {
    std::string vtype="";
    io.data (vtype);
#ifdef LIBMESH_USE_COMPLEX_NUMBERS
    if (vtype != "COMPLEX")
      {
	libMesh::err << "ERROR: File does not contain complex-valued data!"
		      << std::endl;
	libmesh_error();
      }
#elif LIBMESH_USE_REAL_NUMBERS
    if (vtype != "REAL")
      {
	libMesh::err << "ERROR: File does not contain real-valued data!"
		      << std::endl;
	libmesh_error();
      }
#else
    /*
     * What number type is this?
     */
    libmesh_error();
#endif
  }



  /**
   * 3.)  
   *
   * Read the number of nodes for which data is there
   */
  unsigned int n_node = 0;
  io.data (n_node);


  /**
   * 4.)  
   *
   * Read the number of elements for which data is there
   */
  unsigned int n_elem = 0;
  io.data (n_elem);

  unsigned int previous_values_size = 0;

  for (unsigned int n_cnt=0; n_cnt < n_node; n_cnt++)
    {
      /**
       * 5.)  
       *
       * Read the foreign node id, locate the
       * Node* associated with this foreign id
       */
      unsigned int f_id = 0;
      io.data (f_id);

      const Node* node = foreign_id_to_node(f_id);


      /**
       * 6.)  
       *
       * the actual values for this node, Xdr knows
       * the length
       */
      {
        std::vector<Number> values;
	io.data (values);


#ifdef DEBUG
	/*
	 * make sure the size of the values vectors
	 * are identical for all nodes
	 */
	if (n_cnt == 0)
	    previous_values_size = values.size();
	else
	  {
	    if (previous_values_size != values.size())
	      {
		libMesh::err << "ERROR: Size mismatch for n_cnt = " << n_cnt << std::endl;
		libmesh_error();
	      }
	  }
#endif


	/**
	 * insert this node and the values in the _node_data 
	 */
	_node_data.insert (std::make_pair(node, values));
      }
    }



  previous_values_size = 0;

  for (unsigned int n_cnt=0; n_cnt < n_elem; n_cnt++)
    {
      /**
       * 7.)  
       *
       * Read the foreign elem id, locate the Elem*
       */
      unsigned int f_id = 0;
      io.data (f_id);

      const Elem* elem = foreign_id_to_elem(f_id);


      /**
       * 8.)  
       *
       * the actual values for this elem, Xdr knows
       * how many
       */
      {
        std::vector<Number> values;
	io.data (values);


#ifdef DEBUG
	/*
	 * make sure the size of the values vectors
	 * are identical for all elements
	 */
	if (n_cnt == 0)
	    previous_values_size = values.size();
	else
	  {
	    if (previous_values_size != values.size())
	      {
		libMesh::err << "ERROR: Size mismatch for n_cnt = " << n_cnt << std::endl;
		libmesh_error();
	      }
	  }
#endif


	/**
	 * insert this elem and the values in our _elem_data 
	 * @e only when we own this element!
	 */
	if (elem->processor_id() == proc_id)
	  _elem_data.insert (std::make_pair(elem, values));
      }
    }


  /*
   * finished reading.  Now ready for use, provided
   * there was any data contained in the file.
   */
  libmesh_assert ((this->_node_data.size() != 0) || (this->_elem_data.size() != 0));

  this->_node_data_closed = true;
  this->_elem_data_closed = true;
}






void MeshData::write_xdr (const std::string& name,
			  const XdrMODE mode)
{
  /**
   * This code implements the output of the MeshData
   * object in XDR format.  This warrants some documentation.
   * The output consists of 8 sections:
   *
   *   1.) The name of the data stored, if provided (string)
   *
   *   2.) A switch whether real or complex data is stored (string)
   *
   *   3.) The number of nodes for which values are stored
   *       (unsigned int)
   *
   *   4.) The number of elements for which values are stored
   *       (unsigned int)
   *
   *   for each node
   *     
   *     5.) The foreign node id (unsigned int)
   *     
   *     6.) The actual values (vector of real/complex)
   * 
   *   end node loop
   *
   *   for each element
   *     
   *     7.) The foreign element id (unsigned int)
   *     
   *     8.) The actual values (vector of real/complex)
   * 
   *   end node loop
   *
   * Note that the actual IO is handled through the Xdr class 
   * (to be renamed later?) which provides a uniform interface to 
   * both the XDR (eXternal Data Representation) interface and standard
   * ASCII output.  Thus this one section of code will write XDR or ASCII
   * files with no changes.
   */ 

  /*
   * make sure the id maps are ready
   * and that we have data to write
   */
  libmesh_assert (_node_id_map_closed);
  libmesh_assert (_elem_id_map_closed);

  libmesh_assert (_node_data_closed);
  libmesh_assert (_elem_data_closed);
  

  Xdr io(name, mode);


  // all processors write the data in the same format
  //const unsigned int proc_id = _mesh.processor_id();

  /**
   * 1.)
   *
   * Write the descriptive name
   */
  {
    std::string desc = this->_data_descriptor;
    io.data (desc, "# Data description");
  }



  /**
   * 2.)
   *
   * Write: either real or complex
   */
  {
#ifdef LIBMESH_USE_COMPLEX_NUMBERS
    std::string desc = "COMPLEX";
#elif LIBMESH_USE_REAL_NUMBERS
    std::string desc = "REAL";
#else
better_you_choke_this...
#endif
    io.data (desc, "# type of values");
  }



  /**
   * 3.)  
   *
   * Write the number of nodes for which data is there
   */
  {
    unsigned int n_node = this->_node_data.size();
    io.data (n_node, "# No. of nodes for which data is stored");
  }


  /**
   * 4.)  
   *
   * Write the number of elements for which data is there
   */
  {
    unsigned int n_elem = this->_elem_data.size();
    io.data (n_elem, "# No. of elements for which data is stored");
  }




  std::map<const Node*, 
           std::vector<Number> >::const_iterator nit = _node_data.begin ();

  for (; nit != _node_data.end(); ++nit)
    {
      const Node* node = (*nit).first;

      /**
       * 5.)  
       *
       * Write the foreign node id
       */
      {
	unsigned int f_id = node_to_foreign_id (node);
	io.data (f_id, "# Foreign node id");
      }


      /**
       * 6.)  
       *
       * the actual values for this node
       */
      {
	/* 
	 * since we are iterating over our @e own 
	 * map, this libmesh_assert should never break...
	 */
	libmesh_assert (this->has_data(node));

	const std::vector<Number>& values = this->get_data(node);
	
	/*
	 * copy the data to a local buf, since
	 * the Xdr class needs write access, even
	 * when only reading data
	 */
	std::vector<Number> buf = values;
	io.data (buf, "# Values");
      }
    }







  std::map<const Elem*, 
           std::vector<Number> >::const_iterator eit = _elem_data.begin ();

  for (; eit != _elem_data.end(); ++eit)
    {
      const Elem* elem = (*eit).first;

      /**
       * 7.)  
       *
       * Write the foreign element id
       */
      {
	unsigned int f_id = elem_to_foreign_id (elem);
	io.data (f_id, "# Foreign element id");
      }


      /**
       * 8.)  
       *
       * the actual values for this element
       */
      {
	/* 
	 * since we are iterating over our @e own 
	 * map, this libmesh_assert should never break...
	 */
	libmesh_assert (this->has_data(elem));

	const std::vector<Number>& values = this->get_data(elem);
	
	/*
	 * copy the data to a local buf, since
	 * the Xdr class needs write access, even
	 * when only reading data
	 */
	std::vector<Number> buf = values;
	io.data (buf, "# Values");
      }
    }
}

} // namespace libMesh


