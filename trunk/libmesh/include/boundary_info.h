// $Id: boundary_info.h,v 1.14 2003-05-14 11:54:36 ddreyer Exp $

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



#ifndef __boundary_info_h__
#define __boundary_info_h__

// C++ includes
#include <map>
#include <set>

// Local includes
#include "elem.h"
#include "boundary_data.h"


// Forward declarations
class MeshBase;
class BoundaryMesh;



/**
 * The \p BoundaryInfo class contains information relevant
 * to boundary conditions.  It can also build a mesh that
 * just includes boundary elements/faces.
 */

//------------------------------------------------------
// BoundaryInfo class definition
class BoundaryInfo
{
 public:

  /**
   * Constructor.  Takes a reference to the mesh.
   */ 
  BoundaryInfo (const MeshBase& m);

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
   * Allows the \p buondary_mesh to be used like any other mesh.
   * Before this is called the \p boundary_mesh data structure is
   * empty.
   */
  void sync (BoundaryMesh& boundary_mesh);
  
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
   * Returns the number of user-specified boundary ids.
   */
  unsigned int n_boundary_ids () const { return boundary_ids.size(); }

  /**
   * Returns the boundary id associated with \p Node \p node.
   * Returns \p invalid_id if the node is not found, so \p invalid_id
   * can be thought of as a "default" boundary id.
   */
  short int boundary_id (const Node* node) const;

  /**
   * Returns the boundary id associated with the \p side side of
   * element \p elem.  Note that only one id per side is allowed,
   * however multiple sides per element are allowed.  Returns \p invalid_id
   * if the \p side does not have an associated boundary id, hence
   * \p invalid_id can be used as the default boundary id.
   */
  short int boundary_id (const Elem* elem,
			 const unsigned short int side) const;

  /**
   * @returns the number of element-based boundary conditions.
   */
  unsigned int n_boundary_conds () const
  { return boundary_side_id.size(); }
  
  /**
   * Creates a list of nodes and ids for those nodes.
   */
  void build_node_list (std::vector<unsigned int>& nl,
			std::vector<short int>&    il) const;

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
  { return boundary_ids; }


  /**
   * Attach a \p BoundaryData object to this \p BoundaryInfo.
   * Do this @e prior to reading the mesh that this \p BoundaryInfo
   * belongs to!  Note that this object does @e not own the
   * \p bd and therefore does not delete \p bd when itself
   * goes out of scope.
   */
  void attach_boundary_data (BoundaryData* bd);

  /**
   * @returns true when this object has a \p BoundaryData attached.
   */
  bool has_boundary_data () const;

  /**
   * @returns a const reference to the \p BoundaryData
   * object that this object owns.
   */
  const BoundaryData& get_boundary_data () const;

  /**
   * @returns a writable reference to the \p BoundaryData
   * object that this object owns.
   */
  BoundaryData& get_boundary_data ();

  /**
   * Print the boundary information data structure.
   */
  void print_info () const;
  
  /**
   * Number used for internal use. This is the return value
   * if a boundary condition is not specified.
   */
  static const short int invalid_id;


 private:


  /**
   * The Mesh this boundary info pertains to.
   */
  const MeshBase& mesh;

  /**
   * Data structure that maps nodes in the mesh
   * to boundary ids.
   */  
  std::map<const Node*,
	   short int> boundary_node_id;

  /**
   * Data structure that maps sides of elements
   * to boundary ids.
   */
  std::multimap<const Elem*,
                std::pair<unsigned short int, short int> >
                                             boundary_side_id;

  /**
   * A collection of user-specified boundary ids.
   */
  std::set<short int> boundary_ids;

  /**
   * The boundary data object that holds values
   * associated with nodes & elements
   */
  BoundaryData* _boundary_data;



  /**
   * Functor class for printing a single node's info
   * To be used with "for_each".
   */
  class PrintNodeInfo 
  {
  public:
    
    inline
    void operator() (const std::pair<const Node*, short int>& np) const
    {
      std::cout << "  (" << np.first->id()
		<< ", "  << np.second
		<< ")"  << std::endl;
    }
    
  };


  /**
   * Functor class for printing a single side's info.
   * To be used with "for_each".
   */
  class PrintSideInfo
  {
  public:
    PrintSideInfo() {}
    
    inline
    void operator() (const std::pair<const Elem*, std::pair<unsigned short int,short int> >& sp) const
    {
      std::cout << "  (" << sp.first->id()
		<< ", "  << sp.second.first
		<< ", "  << sp.second.second 
		<< ")"   << std::endl;
    }

  private:
  };


  
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
// BoundaryData inline methods
inline  
bool BoundaryInfo::has_boundary_data () const
{
  return (_boundary_data != NULL);
}



inline
const BoundaryData& BoundaryInfo::get_boundary_data () const
{
  return *_boundary_data;
}



inline
BoundaryData& BoundaryInfo::get_boundary_data ()
{
  return *_boundary_data;
}


#endif
