// $Id: boundary_info.h,v 1.6 2003-01-31 17:13:35 spetersen Exp $

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
#include "boundary_mesh.h"


// Forward declarations
class Mesh;



/**
 * The \p BoundaryInfo class contains information relevant
 * to boundary conditions.  It can also write a mesh that
 * just includes boundary elements/faces.
 */

//------------------------------------------------------
// BoundaryInfo class definition
class BoundaryInfo
{
 public:

  /**
   * Constructor.  Optionally takes a pointer to the mesh, defaults
   * to \p NULL.
   */ 
  BoundaryInfo(unsigned int d, const Mesh& m);

  /**
   * Destructor.  Not much to do.
   */ 
  ~BoundaryInfo();

  /**
   * Clears the underlying data structures.
   * Returns the object to a pristine state
   * with no data stored.
   */
  void clear();
  
  /**
   * Close the data structures and prepare for use.
   * Synchronizes the \p boundary_mesh
   * data structures with the \p mesh data structures.
   * Allows the \p buondary_mesh to be used like any other mesh.
   * Before this is called the \p boundary_mesh data structure is
   * empty.
   */
  void sync();
  
  /**
   * Add node number \p node with boundary id \p id to the boundary
   * information data structures.
   */ 
  void add_node(const unsigned int node,
		const short int id);

  /**
   * Add side \p side of element number \p elem with boundary id \p id
   * to the boundary information data structure.
   */
  void add_side(const unsigned int elem,
		const unsigned short int side,
		const short int id);

  /**
   * Add side \p side of element \p elem with boundary id \p id
   * to the boundary information data structure.
   */
  void add_side(const Elem* elem,
		const unsigned short int side,
		const short int id);

  /**
   * Returns the number of user-specified boundary ids.
   */
  unsigned int n_boundary_ids() const { return boundary_ids.size(); }

  /**
   * Read boundary data in the shanee format,
   * taking a string for the name of the file.
   * Probably deprecated.
   */
  void read_shanee_boundary(const std::string& name);
  
  /**
   * Read boundary data in the shanee format,
   * taking an istream which represents the file.
   * Probably deprecated.
   */
  void read_shanee_boundary(std::istream& in);

  /**
   * Returns the boundary id associated with node number \p node.
   * Returns \p invalid_id if the node is not found, so \p invalid_id
   * can be thought of as a "default" boundary id.
   */
  short int boundary_id(const unsigned int node) const;

  /**
   * Returns the boundary id associated with the \p side side of
   * element number \p elem.  Note that only one id per side is allowed,
   * however multiple sides per element are allowed.  Returns \p invalid_id
   * if the \p side does not have an associated boundary id, hence
   * \p invalid_id can be used as the default boundary id.
   */
  short int boundary_id(const Elem* elem,
			const unsigned short int side) const;

  /**
   * Returns the boundary id associated with the \p side side of
   * element number \p elem.  Note that only one id per side is allowed,
   * however multiple sides per element are allowed.  Returns \p invalid_id
   * if the \p side does not have an associated boundary id, hence
   * \p invalid_id can be used as the default boundary id.
   */
  short int boundary_id(const unsigned int elem,
			const unsigned short int side) const;

  /**
   * @returns the number of element-based boundary conditions.
   */
  unsigned int n_boundary_conds() const
  { return elem_list.size(); }

  /**
   * @returns a list of nodes that have boundary conditions.
   */
  const std::vector<unsigned int>& get_node_list() const
  { assert(node_list.size() == node_id_list.size()); return node_list; };

  /**
   * @returns a list of elements that have boundary conditions.
   */
   const std::vector<unsigned int>& get_elem_list() const
  { assert(elem_list.size() == elem_id_list.size()); return elem_list; };

  /**
   * @returns the side of each element that has a boundary condition.
   */
  const std::vector<unsigned short int>& get_side_list() const
  { assert(elem_list.size() == side_list.size()); return side_list; };

  /**
   * @returns a list of boundary condition ids for the nodes.  This
   * vector is the same size as \p node_list.
   */
  const std::vector<short int>& get_node_id_list() const
  { assert(node_list.size() == node_id_list.size()); return node_id_list; };

  /**
   * @returns a list of boundary condition ids for the elements.
   * This vector is the same size as \p elem_list.
   */
  const std::vector<short int>& get_elem_id_list() const
  {  assert(elem_list.size() == elem_id_list.size()); return elem_id_list; };

  /**
   * @returns the user-specified boundary ids.
   */
  const std::set<short int>& get_boundary_ids() const
  { return boundary_ids; };

  /**
   * Add boundary values for node \p node with id \p id to the boundary
   * information data structures.
   */
  void add_boundary_values(const unsigned int node,
			   const std::vector<real> values,
			   const short int id);

  /**
   * @returns a reference to the user-specified boundary values.
   */
  const std::vector<std::pair<unsigned int,
                 std::vector<real> > >& get_boundary_values() const
  { return  boundary_values; };

  /**
   * @returns the boundary values specified for node \p node.
   */
  vector<real> get_boundary_values(unsigned int node) const;

  /**
   * Print the boundary information data structure.
   */
  void print_info() const;
  
  /**
   * Number used for internal use. This is the return value
   * if a boundary condition is not specified.
   */
  static const short int invalid_id;
  
  /**
   * A mesh describing just the boundary.
   */
  BoundaryMesh boundary_mesh;


 private:

  const unsigned int dim;
  
  const Mesh& mesh;
  
  
  std::map<unsigned int, short int> boundary_node_id;
  
  std::multimap<const Elem*,
                std::pair<unsigned short int, short int> >
                                             boundary_side_id;
  
  std::set<short int> boundary_ids;


  std::vector<unsigned int>       node_list;
  std::vector<short int>          node_id_list;
  
  std::vector<unsigned int>       elem_list;
  std::vector<unsigned short int> side_list;
  std::vector<short int>          elem_id_list;

  // a vector for boundary values
  std::vector<std::pair<unsigned int,
              std::vector<real> > >   boundary_values;
  

  // New code
  
  /**
   * Functor class for printing a single node's info
   * To be used with "for_each".
   */
  class PrintNodeInfo 
  {
  public:
    
    inline
    void operator() (const std::pair<const unsigned int, short int>& np) const
    {
      std::cout << "  (" << np.first
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
    PrintSideInfo(const std::map<const Elem*, unsigned int>& m) : estn(m) {}
    
    inline
    void operator() (const std::pair<const Elem*, std::pair<unsigned short int,short int> >& sp) const
    {
      std::map<const Elem*, unsigned int>::const_iterator e_num_it =
	estn.find(sp.first);

      assert (e_num_it != estn.end());
      
      std::cout << "  (" << e_num_it->second
		<< ", "  << sp.second.first
		<< ", "  << sp.second.second 
		<< ")"   << std::endl;
    }

  private:
    const std::map<const Elem*, unsigned int>& estn;
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


#endif
