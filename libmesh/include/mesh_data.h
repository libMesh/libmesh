// $Id: mesh_data.h,v 1.6 2003-06-07 14:36:08 ddreyer Exp $

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



#ifndef __mesh_data_h__
#define __mesh_data_h__

// C++ includes
#include <map>
#include <vector>
#include <string>



// Local Includes
#include "libmesh.h"
#include "node.h"
#include "elem.h"
#include "enum_xdr_mode.h"



// Forward Declarations
class MeshBase;
class UnvMeshInterface;
class XdrInterface;



/**
 * Yet another Mesh-something class...  What's this good for: 
 * \p MeshData handles actual data on entities (nodes, elements)
 * of meshes.  \p MeshBase owns a \p MeshData for dealing with files
 * that contain nodal or element-oriented data, numbered in the same 
 * format as the corresponding mesh file.
 *
 * @author Daniel Dreyer, 2003
 */

// ------------------------------------------------------------
// MeshData class definition
class MeshData 
{
public:

  /**
   * Default Constructor.  Takes const reference
   * to the mesh it belongs to.
   */
  MeshData (const MeshBase& m);

  /**
   * Destructor.
   */
  ~MeshData ();

  /**
   * When \p MeshData should be used, it has to be activated
   * first, @e prior to reading in a mesh with the \p Mesh::read()
   * methods.  Optionally takes a string that should help the user
   * in identifying the data later on.
   */
  void activate (const std::string& descriptor="");

  /**
   * Clears the data fields, but leaves the id maps
   * untouched.  Useful for clearing data for a new
   * data file.  Use \p slim() to delete the maps.
   */
  void clear ();

  /**
   * Once the data is properly read from file, the id 
   * maps can safely be cleared.  However, if this object
   * should remain able to @e write nodal or element oriented 
   * data to file, this method should better @e not be used.
   * Use the appropriate \p bool to select the id map that
   * should be cleared.  By default, both id maps are deleted.
   */
  void slim (const bool node_id_map = true,
	     const bool elem_id_map = true);

  /**
   * Read mesh data from file named \p name.  
   * Guess format from the file extension.  Note that
   * prior to this you have to at least either
   * \p close_node_map() or \p close_elem_map().
   */
  void read (const std::string& name);

  /**
   * Write mesh data to file named \p name.  
   * Guess format from the file extension.
   */
  void write (const std::string& name);

  /**
   * @returns the \f$ i^{th} \f$ value (defaults to 0) associated 
   * with node \p node.  Returns \p libMesh::zero when there
   * is no such \p node in the map.
   */
  Number operator() (const Node* node, 
		     const unsigned int i=0) const;

  /**
   * @returns the \f$ i^{th} \f$ value associated with the 
   * node \p node.  Other than before, this function terminates
   * with an error when data for a non-existent \p Node* is
   * requested.
   */
  Number safe (const Node* node, 
	       const unsigned int i=0) const;

  /**
   * Stores @e all data associated with the node \p node 
   * in the vector \p data.  May resize \p data.
   */
  void operator() (const Node* node,
		   std::vector<Number>& data) const;

  /**
   * @returns the \f$ i^{th} \f$ value (defaults to 0) associated 
   * with element \p elem.  Returns \p libMesh::zero when there
   * is no data for \p elem in the map.
   */
  Number operator() (const Elem* elem, 
		     const unsigned int i=0) const;

  /**
   * @returns the \f$ i^{th} \f$ value associated with the 
   * element \p elem.  Other than before, this function terminates
   * with an error when data for a non-existent \p Elem* is
   * requested.
   */
  Number safe (const Elem* elem, 
	       const unsigned int i=0) const;

  /**
   * Stores @e all data associated with the element \p elem
   * in the vector \p data.  May resize \p data.
   */
  void operator() (const Elem* elem,
		   std::vector<Number>& data) const;

  /**
   * @returns \p true when this object is active and working.
   * Use \p activate() to bring this object alive.
   */
  bool active () const;

  /**
   * @returns \p true when this object is properly initialized
   * and ready for use for @e element associated data, \p false 
   * otherwise.
   */
  bool elem_initialized () const;

  /**
   * @returns \p true when this object is properly initialized
   * and ready for use for @e nodal data, \p false otherwise.
   */
  bool node_initialized () const;


  //----------------------------------------------------------
  // Methods for accessing the node and element maps.
  // Heavily used by the \p read() and \p write() methods.
  /**
   * @returns the \p Node* that this foreign id maps to.
   */
  const Node* foreign_id_to_node (const unsigned int fid) const;

  /**
   * @returns the \p Elem* that this foreign id maps to.
   */
  const Elem* foreign_id_to_elem (const unsigned int fid) const;

  /**
   * @returns the foreign id this \p Node* maps to.
   */
  unsigned int node_to_foreign_id (const Node* n) const;

  /**
   * @returns the foreign id this \p Elem* maps to.
   */
  unsigned int elem_to_foreign_id (const Elem* n) const;


protected:


  //----------------------------------------------------------
  // Methods used by mesh importes to communicate node/element
  // labels to this \p MeshData
  /**
   * In general, \p MeshData gathers nodal data
   * from a file, but it needs to relate this data
   * with the \p Node* of the current mesh.  Mesh
   * importers simply use this method to add such
   * a map.
   */
  void add_foreign_node_id (const Node* node, 
			    const unsigned int foreign_node_id);

  /**
   * In general, \p MeshData gathers element-associated
   * data from file, but it needs to relate this data
   * with the \p Elem* of the current mesh.  Mesh
   * importers simply use this method to add such
   * a map.
   */
  void add_foreign_elem_id (const Elem* elem, 
			    const unsigned int foreign_elem_id);

  /**
   * Signal to this object that the mesh importer finished
   * adding node and element foreign-id maps.
   */
  void close_foreign_id_maps ();



  //----------------------------------------------------------
  // read/write Methods
  /**
   * Read nodal/element oriented data in UNV format.
   */
  void read_unv (const std::string& name);

  /**
   * Write nodal/element oriented data in UNV format.
   */
  void write_unv (const std::string& name);

  /**
   * Read nodal/element oriented data using the
   * \p Xdr class that enables both ASCII and
   * binary format through the same interface.  
   * By default uses ASCII format, but may easily
   * be changed setting \p mode to \p DECODE.
   */
  void read_xdr (const std::string& name,
		 const XdrMODE mode = READ);

  /**
   * Write nodal data in format comparable to
   * the XDR format already known from \p Mesh.
   * By default uses ASCII format, but may easily
   * be changed setting \p mode to \p ENCODE.
   */
  void write_xdr (const std::string& name,
		  const XdrMODE mode = WRITE);


  /**
   * The mesh this object belongs to
   */
  const MeshBase& _mesh;

  /**
   * Some name the user gave to the data when this
   * object got activated
   */
  std::string _data_descriptor;


  //--------------------------------------------------
  // node associated data & maps
  /**
   * The map containing pointers to nodes in the mesh
   * and the corresponding data.
   */
  std::map<const Node*,
           std::vector<Number> > _node_data;

  /**
   * Maps node pointers to node numbers in the @e foreign
   * format.  
   */
  std::map<const Node*,
           unsigned int> _node_id;

  /**
   * Maps @e foreign node ids to node pointers of the
   * current mesh.
   */
  std::map<unsigned int,
           const Node*> _id_node;

  /**
   * The number of real values that are stored per node
   */
  unsigned int _n_val_per_node;



  //--------------------------------------------------
  // element associated data & maps
  /**
   * Maps element pointers to the element-associated data
   */
  std::map<const Elem*,
           std::vector<Number> > _elem_data;

  /**
   * Maps element pointers to element labels in the @e foreign
   * format.  
   */
  std::map<const Elem*,
           unsigned int> _elem_id;
  /**
   * Maps @e foreign element labels to element pointers of the
   * current mesh.
   */
  std::map<unsigned int,
           const Elem*> _id_elem;

  /**
   * The number of real values that are stored per element
   */
  unsigned int _n_val_per_elem;



  //--------------------------------------------------------
  /**
   * \p true when the mesh importer finished adding
   * node id maps.
   */
  bool _node_id_map_closed;

  /**
   * \p true when the nodal data are properly initialized,
   * false otherwise.
   */
  bool _node_data_closed;


  //--------------------------------------------------------
  /**
   * \p true when the mesh importer finished adding
   * element id maps.
   */
  bool _elem_id_map_closed;

  /**
   * \p true when the element based data are properly initialized,
   * false otherwise.
   */
  bool _elem_data_closed;


  //--------------------------------------------------------
  /**
   * \p true when this object is set active (to gather data
   * during mesh import).
   */
  bool _active;

  /**
   * Make the mesh importer class \p UnvInterface friend, so
   * that it can communicate foreign node ids to this class.
   */
  friend class UnvMeshInterface;

  /**
   * Make the mesh importer class \p XdrInterface friend, so
   * that it can communicate foreign node ids to this class.
   */
  friend class XdrInterface;

};


// ------------------------------------------------------------
// MeshData inline methods
inline
Number MeshData::operator() (const Node* node, 
			     const unsigned int i) const
{
  assert (_active);
  assert (_node_id_map_closed);
  assert (_node_data_closed);

  std::map<const Node*, 
           std::vector<Number> >::const_iterator pos = _node_data.find(node);

  if (pos == _node_data.end())
      return libMesh::zero;
  else
    {
      assert (i < (*pos).second.size());
      return (*pos).second[i];
    }
}



inline
Number MeshData::safe (const Node* node, 
		       const unsigned int i) const
{
  assert (_active);
  assert (_node_id_map_closed);
  assert (_node_data_closed);

  std::map<const Node*, 
           std::vector<Number> >::const_iterator pos = _node_data.find(node);

  if (pos == _node_data.end())
    {
#ifdef DEBUG
      std::cerr << "ERROR: No data stored for the node with id" 
		<< node->id() << std::endl;
#endif
      error();
    }


  assert (i < (*pos).second.size());
  return (*pos).second[i];
}



inline
void MeshData::operator() (const Node* node,
			   std::vector<Number>& data) const
{
  assert (_active);
  assert (_node_id_map_closed);
  assert (_node_data_closed);

  std::map<const Node*, 
           std::vector<Number> >::const_iterator pos = _node_data.find(node);

  if (pos == _node_data.end())
    {
      data.clear();
    }
  else
    {
      data = (*pos).second;
    }
  return;
}



inline
Number MeshData::operator() (const Elem* elem, 
			     const unsigned int i) const
{
  assert (_active);
  assert (_elem_id_map_closed);
  assert (_elem_data_closed);

  std::map<const Elem*, 
           std::vector<Number> >::const_iterator pos = _elem_data.find(elem);

  if (pos == _elem_data.end())
    return libMesh::zero;
  
  
  assert (i < (*pos).second.size());
  return (*pos).second[i];
}



inline
Number MeshData::safe (const Elem* elem, 
		       const unsigned int i) const
{
  assert (_active);
  assert (_elem_id_map_closed);
  assert (_elem_data_closed);

  std::map<const Elem*, 
           std::vector<Number> >::const_iterator pos = _elem_data.find(elem);

  if (pos == _elem_data.end())
    {
#ifdef DEBUG
      std::cerr << "ERROR: No data stored for the element with id " 
		<< elem->id() << std::endl;
#endif
      error();
    }


  assert (i < (*pos).second.size());
  return (*pos).second[i];
}



inline
void MeshData::operator() (const Elem* elem,
			   std::vector<Number>& data) const
{
  assert (_active);
  assert (_elem_id_map_closed);
  assert (_elem_data_closed);

  std::map<const Elem*, 
           std::vector<Number> >::const_iterator pos = _elem_data.find(elem);

  if (pos == _elem_data.end())
    {
      data.clear();
    }
  else
    {
      data = (*pos).second;
    }
  return;
}



inline
bool MeshData::active() const
{
  return _active;
}



inline
bool MeshData::elem_initialized() const
{
  return (_active && _elem_data_closed);
}



inline
bool MeshData::node_initialized() const
{
  return (_active && _node_data_closed);
}



inline  
void MeshData::add_foreign_node_id (const Node* node, 
				    const unsigned int foreign_node_id)
{
  if (_active)
    {
      assert (!_node_id_map_closed);
      assert (node                             != NULL);
      assert (_node_id.find(node)              == _node_id.end());
      assert (_id_node.find(foreign_node_id)   == _id_node.end());

      /*
       * _always_ insert in _id_node and _node_id.  If we would 
       * use the mesh.node(unsigned int) method or the node.id()
       * to get Node* and unsigned int, respectively, we would not
       * be safe any more when the mesh gets refined or re-numbered
       * within libMesh. And we could get in big trouble that would
       * be hard to find when importing data _after_ having refined...
       */
      _node_id.insert(std::make_pair(node, foreign_node_id));
      _id_node.insert(std::make_pair(foreign_node_id, node));
    }
}



inline  
void MeshData::add_foreign_elem_id (const Elem* elem, 
				    const unsigned int foreign_elem_id)
{
  if (_active)
    {
      assert (!_elem_id_map_closed);
      assert (elem                             != NULL);
      assert (_elem_id.find(elem)              == _elem_id.end());
      assert (_id_elem.find(foreign_elem_id)   == _id_elem.end());

      _elem_id.insert(std::make_pair(elem, foreign_elem_id));
      _id_elem.insert(std::make_pair(foreign_elem_id, elem));
    }
}




#endif
