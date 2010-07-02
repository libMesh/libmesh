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



#ifndef __boundary_info_h__
#define __boundary_info_h__

// C++ includes
#include <vector>
#include <map>
#include <set>

// Local includes
#include "libmesh_common.h"

namespace libMesh
{


// Forward declarations
class Elem;
class Node;
class MeshBase;
class BoundaryMesh;
class MeshData;


/**
 * The \p BoundaryInfo class contains information relevant
 * to boundary conditions: it does not hold actual boundary 
 * condition data (check \p MeshData for that), but can mark
 * element faces and nodes with ids useful for identifying the 
 * type of boundary condtion.  It can also build a mesh that
 * just includes boundary elements/faces.
 *
 * TODO[JWP]: Generalize this to work with MeshBase again.
 */

//------------------------------------------------------
// BoundaryInfo class definition
class BoundaryInfo
{
protected:
  friend class MeshBase;
  
  /**
   * Constructor.  Takes a reference to the mesh.
   * The BoundaryInfo class is only used internally
   * by the Mesh class.  A user should never instantiate
   * this class.  Therefore the constructor is protected.
   */ 
  BoundaryInfo (const MeshBase& m);

public:

  /**
   * Destructor.  Not much to do.
   */ 
  ~BoundaryInfo ();

  /**
   * Clears the underlying data structures.
   * Returns the object to a pristine state
   * with no data stored.
   */
  void clear ();
  
  /**
   * Close the data structures and prepare for use.
   * Synchronizes the \p boundary_mesh
   * data structures with the \p mesh data structures.
   * Allows the \p boundary_mesh to be used like any other mesh.
   * Before this is called the \p boundary_mesh data structure is
   * empty.
   *
   * If you are using a MeshData class with this Mesh, you can
   * pass a pointer to both the boundary_mesh's MeshData object,
   * and the MeshData object used for this mesh.
   */
  void sync (BoundaryMesh& boundary_mesh,
	     MeshData* boundary_mesh_data=NULL,
	     MeshData* this_mesh_data=NULL);
  
  /**
   * Close the data structures and prepare for use.
   * Synchronizes the \p boundary_mesh
   * data structures with the \p mesh data structures.
   * Allows the \p boundary_mesh to be used like any other mesh.
   * Before this is called the \p boundary_mesh data structure is
   * empty.  Only boundary elements with the specified ids are
   * extracted.
   *
   * If you are using a MeshData class with this Mesh, you can
   * pass a pointer to both the boundary_mesh's MeshData object,
   * and the MeshData object used for this mesh.
   */
  void sync (const std::set<short int> &requested_boundary_ids,
	     BoundaryMesh& boundary_mesh);
  
  /**
   * Add \p Node \p node with boundary id \p id to the boundary
   * information data structures.
   */ 
  void add_node (const Node* node,
		 const short int id);

  /**
   * Add node number \p node with boundary id \p id to the boundary
   * information data structures.
   */ 
  void add_node (const unsigned int node,
		 const short int id);

  /**
   * Add side \p side of element number \p elem with boundary id \p id
   * to the boundary information data structure.
   */
  void add_side (const unsigned int elem,
		 const unsigned short int side,
		 const short int id);

  /**
   * Add side \p side of element \p elem with boundary id \p id
   * to the boundary information data structure.
   */
  void add_side (const Elem* elem,
		 const unsigned short int side,
		 const short int id);

  /**
   * Removes the boundary conditions associated with node \p node,
   * if any exist.
   */
  void remove (const Node* node);

  /**
   * Removes the boundary conditions associated with element \p elem,
   * if any exist.
   */
  void remove (const Elem* elem);

  /**
   * Removes all boundary conditions associated with side \p side of
   * element \p elem, if any exist.
   */
  void remove_side (const Elem* elem,
                    const unsigned short int side);

  /**
   * Removes the boundary id \p id from side \p side of element \p
   * elem, if it exists.
   */
  void remove_side (const Elem* elem,
                    const unsigned short int side,
                    const short int id);

  /**
   * Returns the number of user-specified boundary ids.
   */
  unsigned int n_boundary_ids () const { return _boundary_ids.size(); }

  /**
   * Returns the boundary ids associated with \p Node \p node.
   */
  std::vector<short int> boundary_ids (const Node* node) const;

  /**
   * Returns the boundary id associated with the \p side side of
   * element \p elem.  Note that only one id per side is allowed,
   * however multiple sides per element are allowed.  Returns \p invalid_id
   * if the \p side does not have an associated boundary id, hence
   * \p invalid_id can be used as the default boundary id.
   */
  short int boundary_id (const Elem* const elem,
			 const unsigned short int side) const;

  /**
   * Returns a side of element \p elem whose associated boundary id is 
   * \p boundary_id if such a side exists.
   * If multiple sides of \p elem have the same id, only the lowest numbered
   * such side is returned.
   *
   * Returns \p invalid_uint if no side has the requested boundary id.
   */
  unsigned int side_with_boundary_id(const Elem* const elem,
				     const unsigned short int boundary_id) const;

  /**
   * Builds the list of unique node boundary ids.
   */
  void build_node_boundary_ids(std::vector<short int> &b_ids);

  /**
   * Builds the list of unique side boundary ids.
   */
  void build_side_boundary_ids(std::vector<short int> &b_ids);

  /**
   * @returns the number of element-based boundary conditions.
   */
  unsigned int n_boundary_conds () const;
  
  /**
   * Creates a list of nodes and ids for those nodes.
   */
  void build_node_list (std::vector<unsigned int>& nl,
			std::vector<short int>&    il) const;

  /**
   * Adds nodes with boundary ids based on the side's boundary
   * ids they are connected to.
   */
  void build_node_list_from_side_list();

  /**
   * Creates a list of element numbers, sides, and  and ids for those sides.
   */
  void build_side_list (std::vector<unsigned int>&       el,
			std::vector<unsigned short int>& sl,
			std::vector<short int>&          il) const;

  /**
   * @returns the user-specified boundary ids.
   */
  const std::set<short int>& get_boundary_ids () const
  { return _boundary_ids; }

  /**
   * Print the boundary information data structure.
   */
  void print_info (std::ostream& out=libMesh::out) const;
  
  /**
   * Number used for internal use. This is the return value
   * if a boundary condition is not specified.
   */
  static const short int invalid_id;


 private:


  /**
   * The Mesh this boundary info pertains to.
   */
  const MeshBase& _mesh;

  /**
   * Data structure that maps nodes in the mesh
   * to boundary ids.
   */  
  std::multimap<const Node*,
                short int> _boundary_node_id;

  /**
   * Data structure that maps sides of elements
   * to boundary ids.
   */
  std::multimap<const Elem*,
                std::pair<unsigned short int, short int> >
                                             _boundary_side_id;

  /**
   * A collection of user-specified boundary ids.
   */
  std::set<short int> _boundary_ids;



//   /**
//    * Functor class for printing a single node's info
//    * To be used with "for_each".
//    */
//   class PrintNodeInfo 
//   {
//   public:
//     inline
//     void operator() (const std::pair<const Node*, short int>& np) const
//     {
//       std::cout << "  (" << np.first->id()
// 		<< ", "  << np.second
// 		<< ")"  << std::endl;
//     }
//   };


//   /**
//    * Functor class for printing a single side's info.
//    * To be used with "for_each".
//    */
//   class PrintSideInfo
//   {
//   public:
//     PrintSideInfo() {}
//     inline
//     void operator() (const std::pair<const Elem*, std::pair<unsigned short int,short int> >& sp) const
//     {
//       std::cout << "  (" << sp.first->id()
// 		<< ", "  << sp.second.first
// 		<< ", "  << sp.second.second 
// 		<< ")"   << std::endl;
//     }
//   };


  
  /**
   * Functor class for initializing a map.
   * The entries being added to the map
   * increase by exactly one each time.
   * The desctructor also inserts the
   * invalid_id entry.
   */
  class Fill
  {
  public:
    Fill(std::map<short int, unsigned int>& im) : id_map(im), cnt(0) {}

    ~Fill()
    {
      id_map[invalid_id] = cnt;
    }
    
    inline
    void operator() (const short int & pos)
    {
      id_map[pos] = cnt++;
    }
    
  private:
    std::map<short int, unsigned int>& id_map;
    unsigned int cnt;
  };
  
};






// ------------------------------------------------------------
// BoundaryInfo inline methods
inline
void BoundaryInfo::remove (const Node* node)
{
  libmesh_assert (node != NULL);
  
  // Erase everything associated with node
  _boundary_node_id.erase (node);
}



inline
void BoundaryInfo::remove (const Elem* elem)
{
  libmesh_assert (elem != NULL);
  
  // Erase everything associated with elem
  _boundary_side_id.erase (elem);
}

} // namespace libMesh

#endif
