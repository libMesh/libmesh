// $Id: boundary_data.h,v 1.1 2003-05-14 11:54:36 ddreyer Exp $

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



#ifndef __boundary_data_h__
#define __boundary_data_h__

// C++ includes
#include <map>
#include <vector>
#include <string>


// Local Includes
#include "libmesh.h"
#include "node.h"
#include "elem.h"


// Forward Declarations
class MeshBase;



/**
 * Yet another Boundary-something class...  What's this good for: 
 *  \p BoundaryInfo handles @e information about boundaries: which id
 * is associated with which node etc.  Now \p BoundaryData handles actual
 * data on boundaries of meshes.  In general, a \p BoundaryData 
 * belongs to a \p BoundaryInfo for dealing with files that 
 * contain nodal or element-oriented data, numbered in the same 
 * format as the corresponding mesh file.  It makes also sense
 * to separate \p BoundaryData from \p Mesh, since generally there
 * are multiple data files for different load cases associated with 
 * one mesh file.
 *
 * @author Daniel Dreyer, 2003
 */

// ------------------------------------------------------------
// BoundaryData class definition
class BoundaryData 
{
public:

  /**
   * Default Constructor.  Takes const reference
   * to the mesh it belongs to.
   */
  BoundaryData (const MeshBase& m);

  /**
   * Destructor.
   */
  ~BoundaryData ();

  /**
   * Clears the data fields, but leaves the id maps
   * untouched.  Useful for clearing data for a new
   * data file.  Use \p slim() to delete the maps.
   */
  void clear ();

  /**
   * Once the data is properly read from file, the id 
   * maps can safely be cleared.  However, if this object
   * should remain able to write nodal or element oriented 
   * data to file, this method should @e not be used.
   */
  void slim ();

  /**
   * In general, \p BoundaryData gathers nodal data
   * from a file, but it needs to relate this data
   * with the \p Node* of the current mesh.  Mesh
   * importers simply use this method to add such
   * a map.
   */
  void add_foreign_node_id (const Node* node, 
			    const unsigned int foreign_node_id);

  /**
   * In general, \p BoundaryData gathers element-associated
   * data from file, but it needs to relate this data
   * with the \p Elem* of the current mesh.  Mesh
   * importers simply use this method to add such
   * a map.
   */
  void add_foreign_elem_id (const Elem* elem, 
			    const unsigned int foreign_elem_id);

  /**
   * Signal to this object that the mesh importer finished
   * adding node id maps.
   */
  void close_node_map ();

  /**
   * Signal to this object that the mesh importer finished
   * adding elem id maps.
   */
  void close_elem_map ();

  /**
   * Read boundary data from file named \p name.  
   * Guess format from the file extension.  Note that
   * prior to this you have to at least either
   * \p close_node_map() or \p close_elem_map().
   */
  void read (const std::string& name);

  /**
   * Write boundary data to file named \p name.  
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
   * is no such \p elem in the map.
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
   * @returns \p true when this object is properly initialized
   * and ready for use, \p false otherwise.
   */
  bool initialized () const;


protected:

  /**
   * Read nodal/element oriented data in UNV format.
   */
  void read_unv (const std::string& name);

  /**
   * Write nodal/element oriented data in UNV format.
   */
  void write_unv (const std::string& name);

  /**
   * Read nodal/element oriented data in format comparable
   * to the XDR format already known from \p Mesh.
   * This method actually expects an ASCII-file.
   */
  void read_xdr (const std::string& name);

  /**
   * Read nodal/element oriented data in format comparable
   * to the XDR format already known from \p Mesh.
   * Expects a true XDR-Encoded binary file.
   */
  void read_xdr_binary (const std::string& name);

  /**
   * Write nodal data in format comparable to
   * the XDR format already known from \p Mesh.
   * This method actually expects an ASCII-file.
   */
  void write_xdr (const std::string& name);

  /**
   * Same, but expects a true XDR-Encoded binary file.
   */
  void write_xdr_binary (const std::string& name);

  /**
   * The mesh this object belongs to
   */
  const MeshBase& _mesh;

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
   * \p true when the mesh importer finished adding
   * node id maps.
   */
  bool _node_id_map_closed;

  /**
   * \p true when the nodal data are properly initialized,
   * false otherwise.
   */
  bool _node_initialized;

  /**
   * \p true when the mesh importer finished adding
   * element id maps.
   */
  bool _elem_id_map_closed;

  /**
   * \p true when the element based data are properly initialized,
   * false otherwise.
   */
  bool _elem_initialized;

};


// ------------------------------------------------------------
// BoundaryData inline methods
inline  
void BoundaryData::add_foreign_node_id (const Node* node, 
					const unsigned int foreign_node_id)
{
  assert (node != NULL);
  assert (_node_id.find(node) == _node_id.end());
  assert (!_node_id_map_closed);

  this->_node_id.insert(std::make_pair(node, foreign_node_id));
  this->_id_node.insert(std::make_pair(foreign_node_id, node));

}



inline  
void BoundaryData::add_foreign_elem_id (const Elem* elem, 
					const unsigned int foreign_elem_id)
{
  assert (elem != NULL);
  assert (_elem_id.find(elem) == _elem_id.end());
  assert (!_elem_id_map_closed);

  this->_elem_id.insert(std::make_pair(elem, foreign_elem_id));
  this->_id_elem.insert(std::make_pair(foreign_elem_id, elem));
}



inline
Number BoundaryData::operator() (const Node* node, 
				 const unsigned int i) const
{
  assert (this->_node_initialized);

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
Number BoundaryData::safe (const Node* node, 
			   const unsigned int i) const
{
  assert (this->_node_initialized);

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
  else
    {
      assert (i < (*pos).second.size());
      return (*pos).second[i];
    }
}



inline
void BoundaryData::operator() (const Node* node,
			       std::vector<Number>& data) const
{
  assert (this->_node_initialized);

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
Number BoundaryData::operator() (const Elem* elem, 
				 const unsigned int i) const
{
  assert (this->_elem_initialized);

  std::map<const Elem*, 
           std::vector<Number> >::const_iterator pos = _elem_data.find(elem);

  if (pos == _elem_data.end())
      return libMesh::zero;
  else
    {
      assert (i < (*pos).second.size());
      return (*pos).second[i];
    }
}



inline
Number BoundaryData::safe (const Elem* elem, 
			   const unsigned int i) const
{
  assert (this->_elem_initialized);

  std::map<const Elem*, 
           std::vector<Number> >::const_iterator pos = _elem_data.find(elem);

  if (pos == _elem_data.end())
    {
#ifdef DEBUG
      std::cerr << "ERROR: No data stored for the element with the node ids:" 
		<< std::endl;
      for (unsigned int n=0; n<elem->n_nodes(); n++)
	  std::cerr << elem->get_node(n)->id() << std::endl;
#endif
      error();
    }
  else
    {
      assert (i < (*pos).second.size());
      return (*pos).second[i];
    }
}



inline
void BoundaryData::operator() (const Elem* elem,
			       std::vector<Number>& data) const
{
  assert (this->_elem_initialized);

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
bool BoundaryData::initialized() const
{
  /*
   * when there is an element id map, then
   * both _elem_initialized _and_
   * _node_initialized have to be true for
   * the whole object to be initialized.
   * Otherwise, only the _node_initialized
   * has to be true.
   */
  if (_elem_id.empty() && (!_elem_id_map_closed))
    return (_elem_initialized && _node_initialized);
  else
    return (_node_initialized);
}




#endif
