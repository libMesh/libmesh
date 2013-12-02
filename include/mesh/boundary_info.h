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



#ifndef LIBMESH_BOUNDARY_INFO_H
#define LIBMESH_BOUNDARY_INFO_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/id_types.h"
#include "libmesh/parallel_object.h"

// C++ includes
#include <cstddef>
#include <map>
#include <set>
#include <vector>

namespace libMesh
{


// Forward declarations
class Elem;
class Node;
class MeshBase;
class UnstructuredMesh;
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
  class BoundaryInfo : public ParallelObject
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
   * Actual copying operation.
   *
   * Note that it does not copy the mesh over (for obvious reasons).
   */
  BoundaryInfo& operator=(const BoundaryInfo& other_boundary_info);


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
  void sync (UnstructuredMesh& boundary_mesh,
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
  void sync (const std::set<boundary_id_type> &requested_boundary_ids,
	     UnstructuredMesh& boundary_mesh,
	     MeshData* boundary_mesh_data=NULL,
	     MeshData* this_mesh_data=NULL);

  /**
   * Add \p Node \p node with boundary id \p id to the boundary
   * information data structures.
   */
  void add_node (const Node* node,
		 const boundary_id_type id);

  /**
   * Add node number \p node with boundary id \p id to the boundary
   * information data structures.
   */
  void add_node (const dof_id_type node,
		 const boundary_id_type id);

  /**
   * Add \p Node \p node with boundary ids \p ids to the boundary
   * information data structure.
   */
  void add_node (const Node* node,
		 const std::vector<boundary_id_type>& ids);

  /**
   * Clears all the boundary information from all of the nodes in the mesh
   */
  void clear_boundary_node_ids();

  /**
   * Add edge \p edge of element number \p elem with boundary id \p id
   * to the boundary information data structure.
   * Edge-based boundary IDs should only be used in 3D.
   */
  void add_edge (const dof_id_type elem,
		 const unsigned short int edge,
		 const boundary_id_type id);

  /**
   * Add edge \p edge of element \p elem with boundary id \p id
   * to the boundary information data structure.
   * Edge-based boundary IDs should only be used in 3D.
   */
  void add_edge (const Elem* elem,
		 const unsigned short int edge,
		 const boundary_id_type id);

  /**
   * Add edge \p edge of element \p elem with boundary ids \p ids
   * to the boundary information data structure.
   * Edge-based boundary IDs should only be used in 3D.
   */
  void add_edge (const Elem* elem,
		 const unsigned short int edge,
		 const std::vector<boundary_id_type>& ids);

  /**
   * Add side \p side of element number \p elem with boundary id \p id
   * to the boundary information data structure.
   */
  void add_side (const dof_id_type elem,
		 const unsigned short int side,
		 const boundary_id_type id);

  /**
   * Add side \p side of element \p elem with boundary id \p id
   * to the boundary information data structure.
   */
  void add_side (const Elem* elem,
		 const unsigned short int side,
		 const boundary_id_type id);

  /**
   * Add side \p side of element \p elem with boundary ids \p ids
   * to the boundary information data structure.
   */
  void add_side (const Elem* elem,
		 const unsigned short int side,
		 const std::vector<boundary_id_type>& ids);

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
   * Removes all boundary conditions associated with edge \p edge of
   * element \p elem, if any exist.
   */
  void remove_edge (const Elem* elem,
                    const unsigned short int edge);

  /**
   * Removes the boundary id \p id from edge \p edge of element \p
   * elem, if it exists.
   */
  void remove_edge (const Elem* elem,
                    const unsigned short int edge,
                    const boundary_id_type id);

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
                    const boundary_id_type id);

  /**
   * Returns the number of user-specified boundary ids.
   */
  std::size_t n_boundary_ids () const { return _boundary_ids.size(); }

  /**
   * Returns true iff the given node is associated with the given id.
   */
  bool has_boundary_id (const Node* const node,
                        const boundary_id_type id) const;

  /**
   * Returns the boundary ids associated with \p Node \p node.
   */
  std::vector<boundary_id_type> boundary_ids (const Node* node) const;

  /**
   * Returns the number of boundary ids associated with \p Node \p node.
   */
  unsigned int n_boundary_ids (const Node* node) const;

  /**
   * Returns the number of boundary ids associated with the \p edge
   * edge of element \p elem.
   * Edge-based boundary IDs should only be used in 3D.
   */
  unsigned int n_edge_boundary_ids (const Elem* const elem,
                                    const unsigned short int edge) const;

  /**
   * Returns the list of boundary ids associated with the \p edge edge of
   * element \p elem.
   * Edge-based boundary IDs should only be used in 3D.
   */
  std::vector<boundary_id_type> edge_boundary_ids (const Elem* const elem,
                                                   const unsigned short int edge) const;

  /**
   * Returns the list of raw boundary ids associated with the \p edge
   * edge of element \p elem.  These ids are ``raw'' because they
   * exclude ids which are implicit, such as a child's inheritance of
   * its ancestors' boundary id.
   * Edge-based boundary IDs should only be used in 3D.
   */
  std::vector<boundary_id_type> raw_edge_boundary_ids (const Elem* const elem,
                                                       const unsigned short int edge) const;

  /**
   * Returns true iff the given side of the given element is
   * associated with the given id.
   */
  bool has_boundary_id (const Elem* const elem,
                        const unsigned short int side,
                        const boundary_id_type id) const;

  /**
   * Returns the boundary id associated with the \p side side of
   * element \p elem.  Note that only one id per side is allowed,
   * however multiple sides per element are allowed.  Returns \p invalid_id
   * if the \p side does not have an associated boundary id, hence
   * \p invalid_id can be used as the default boundary id.
   */
  boundary_id_type boundary_id (const Elem* const elem,
			        const unsigned short int side) const;

  /**
   * Returns the number of boundary ids associated with the \p side
   * side of element \p elem.
   */
  unsigned int n_boundary_ids (const Elem* const elem,
                               const unsigned short int side) const;

  /**
   * Returns the list of boundary ids associated with the \p side side of
   * element \p elem.
   */
  std::vector<boundary_id_type> boundary_ids (const Elem* const elem,
                                              const unsigned short int side) const;

  /**
   * Returns the list of raw boundary ids associated with the \p side
   * side of element \p elem.  These ids are ``raw'' because they
   * exclude ids which are implicit, such as a child's inheritance of
   * its ancestors' boundary id.
   */
  std::vector<boundary_id_type> raw_boundary_ids (const Elem* const elem,
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
				     const boundary_id_type boundary_id) const;

  /**
   * Builds the list of unique node boundary ids.
   */
  void build_node_boundary_ids(std::vector<boundary_id_type> &b_ids);

  /**
   * Builds the list of unique side boundary ids.
   */
  void build_side_boundary_ids(std::vector<boundary_id_type> &b_ids);

  /**
   * @returns the number of element-side-based boundary conditions.
   */
  std::size_t n_boundary_conds () const;

  /**
   * @returns the number of edge-based boundary conditions.
   * Edge-based boundary IDs should only be used in 3D.
   */
  std::size_t n_edge_conds () const;

  /**
   * @returns the number of node-based boundary conditions.
   */
  std::size_t n_nodeset_conds () const;

  /**
   * Creates a list of nodes and ids for those nodes.
   */
  void build_node_list (std::vector<dof_id_type>&      node_id_list,
			std::vector<boundary_id_type>& bc_id_list) const;

  /**
   * Adds nodes with boundary ids based on the side's boundary
   * ids they are connected to.
   */
  void build_node_list_from_side_list();

  /**
   * Adds sides to a sideset if every node on that side are in the same
   * sideset
   */
  void build_side_list_from_node_list();

  /**
   * Creates a list of element numbers, sides, and  and ids for those sides.
   */
  void build_side_list (std::vector<dof_id_type>&        element_id_list,
			std::vector<unsigned short int>& side_list,
			std::vector<boundary_id_type>&   bc_id_list) const;

  /**
   * @returns the user-specified boundary ids.
   */
  const std::set<boundary_id_type>& get_boundary_ids () const
  { return _boundary_ids; }

  /**
   * Returns a reference to the set of all boundary IDs
   * specified on sides.
   */
  const std::set<boundary_id_type>& get_side_boundary_ids () const
  { return _side_boundary_ids; }

  /**
   * Returns a reference to the set of all boundary IDs
   * specified on edges.
   * Edge-based boundary IDs should only be used in 3D.
   */
  const std::set<boundary_id_type>& get_edge_boundary_ids () const
  { return _edge_boundary_ids; }

  /**
   * Returns a reference to the set of all boundary IDs
   * specified on nodes.
   */
  const std::set<boundary_id_type>& get_node_boundary_ids () const
  { return _node_boundary_ids; }


  /**
   * Print the boundary information data structure.
   */
  void print_info (std::ostream& out=libMesh::out) const;

  /**
   * Print a summary of the boundary information.
   */
  void print_summary (std::ostream& out=libMesh::out) const;

  /**
   * Returns a writable reference for getting/setting an optional
   * name for a sideset name.
   */
  std::string& sideset_name(boundary_id_type id);

  /**
   * Returns a writable reference for getting/setting an optional
   * name for a nodeset name.
   */
  std::string& nodeset_name(boundary_id_type id);

  /**
   * Returns a the id of the requested boundary by name.  Throws an error
   * if a sideset or nodeset by name is not found
   */
  boundary_id_type get_id_by_name(const std::string& name) const;

  /**
   * Return a writeable reference to the whole sideset name map
   */
  std::map<boundary_id_type, std::string>& set_sideset_name_map ()
  { return _ss_id_to_name; }
  const std::map<boundary_id_type, std::string>& get_sideset_name_map () const
  { return _ss_id_to_name; }

  /**
   * Return a writeable reference to the whole nodeset name map
   */
  std::map<boundary_id_type, std::string>& set_nodeset_name_map ()
  { return _ns_id_to_name; }
  const std::map<boundary_id_type, std::string>& get_nodeset_name_map () const
  { return _ns_id_to_name; }

  /**
   * Number used for internal use. This is the return value
   * if a boundary condition is not specified.
   */
  static const boundary_id_type invalid_id;


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
                boundary_id_type> _boundary_node_id;

  /**
   * Data structure that maps edges of elements
   * to boundary ids. This is only relevant in 3D.
   */
  std::multimap<const Elem*,
                std::pair<unsigned short int, boundary_id_type> >
                                             _boundary_edge_id;

  /**
   * Data structure that maps sides of elements
   * to boundary ids.
   */
  std::multimap<const Elem*,
                std::pair<unsigned short int, boundary_id_type> >
                                             _boundary_side_id;

  /**
   * A collection of user-specified boundary ids for sides, edges and nodes.
   * See _side_boundary_ids, _edge_boundary_ids and _node_boundary_ids
   * for sets containing IDs for only sides and only nodes, respectively.
   */
  std::set<boundary_id_type> _boundary_ids;

  /**
   * Set of user-specified boundary IDs for sides *only*.  Note: _boundary_ids
   * is the union of this set, _edge_boundary_ids and _node_boundary_ids.
   */
  std::set<boundary_id_type> _side_boundary_ids;

  /**
   * Set of user-specified boundary IDs for edges *only*.
   * This is only relevant in 3D. Note: _boundary_ids
   * is the union of this set, _side_boundary_ids
   * and _node_boundary_ids.
   */
  std::set<boundary_id_type> _edge_boundary_ids;

  /**
   * Set of user-specified boundary IDs for nodes *only*.  Note: _boundary_ids
   * is the union of this set, _edge_boundary_ids and _side_boundary_ids.
   */
  std::set<boundary_id_type> _node_boundary_ids;

  /**
   * This structure maintains the mapping of named side sets
   * for file formats that support named blocks.  Currently
   * this is only implemented for ExodusII
   */
  std::map<boundary_id_type, std::string> _ss_id_to_name;

  /**
   * This structure maintains the mapping of named node sets
   * for file formats that support named blocks.  Currently
   * this is only implemented for ExodusII
   */
  std::map<boundary_id_type, std::string> _ns_id_to_name;


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
//       libMesh::out << "  (" << np.first->id()
//		      << ", "  << np.second
//		      << ")"  << std::endl;
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
//       libMesh::out << "  (" << sp.first->id()
//		      << ", "  << sp.second.first
//		      << ", "  << sp.second.second
//		      << ")"   << std::endl;
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
    Fill(std::map<boundary_id_type, dof_id_type>& im) : id_map(im), cnt(0) {}

    ~Fill()
    {
      id_map[invalid_id] = cnt;
    }

    inline
    void operator() (const boundary_id_type& pos)
    {
      id_map[pos] = cnt++;
    }

  private:
    std::map<boundary_id_type, dof_id_type>& id_map;
    dof_id_type cnt;
  };

};






// ------------------------------------------------------------
// BoundaryInfo inline methods
inline
void BoundaryInfo::remove (const Node* node)
{
  libmesh_assert(node);

  // Erase everything associated with node
  _boundary_node_id.erase (node);
}



inline
void BoundaryInfo::remove (const Elem* elem)
{
  libmesh_assert(elem);

  // Erase everything associated with elem
  _boundary_edge_id.erase (elem);
  _boundary_side_id.erase (elem);
}

} // namespace libMesh

#endif // LIBMESH_BOUNDARY_INFO_H
