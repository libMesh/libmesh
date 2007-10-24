// $Id: mesh_refinement.h,v 1.14 2003-06-03 16:45:50 benkirk Exp $

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



#ifndef __mesh_refinement_h__
#define __mesh_refinement_h__



#include "mesh_config.h"

#ifdef ENABLE_AMR

// C++ Includes   -----------------------------------
#include <vector>
#include <set>

#if   defined(HAVE_HASH_MAP)
# include <hash_map>
#elif defined(HAVE_EXT_HASH_MAP)
# include <ext/hash_map>
#else
# include <map>
#endif

// Local Includes -----------------------------------
#include "mesh_common.h"

class MeshBase;
class Point;
class Node;
class ErrorVector;



/**
 * This is the \p MeshRefinement class.  This class implements
 * adaptive mesh refinement algorithms for a \p MeshBase.
 *
 * @author Benjamin S. Kirk, 2002-2003.
 */


// ------------------------------------------------------------
// MeshRefinement class definition
class MeshRefinement
{
public:

  /**
   * Constructor.
   */
  MeshRefinement (MeshBase& mesh);

  /**
   * Destructor. Deletes all the elements that are currently stored.
   */
  ~MeshRefinement ();

  /**
   * Deletes all the data that are currently stored.
   */
  void clear ();  

  /**
   * Flags elements for coarsening and refinement based on
   * the computed error passed in \p error_per_cell.  The two
   * fractions \p refine_fraction and \p coarsen_fraction must be in
   * \f$ [0,1] \f$.
   */
  void flag_elements_by_error_fraction (const ErrorVector& error_per_cell,
					const Real refine_fraction  = 0.3,
					const Real coarsen_fraction = 0.0,
					const unsigned int max_level = static_cast<unsigned int>(-1));

  /**
   * Flags elements for coarsening and refinement based on
   * the computed error passed in \p error_per_cell.  This method picks
   * the top \p refine_fraction * \p n_elem elements for refinement and
   * the bottom \p coarsen_fraction * \p n_elem elements for coarsening.
   * The two fractions \p refine_fraction and \p coarsen_fraction must be
   * in \f$ [0,1] \f$.
   */
  void flag_elements_by_elem_fraction (const ErrorVector& error_per_cell,
				       const Real refine_fraction  = 0.3,
				       const Real coarsen_fraction = 0.0,
				       const unsigned int max_level = static_cast<unsigned int>(-1));

  /**
   * Flags elements for coarsening and refinement based on
   * the computed error passed in \p error_per_cell.  This method picks
   * the top \p refine_fraction * \p stddev + \p mean elements for refinement
   * and the bottom \p mean - \p coarsen_fraction * \p stddev elements for
   * coarsening. The two fractions \p refine_fraction and \p coarsen_fraction
   * must be in \f$ [0,1] \f$.
   */
  void flag_elements_by_mean_stddev (const ErrorVector& error_per_cell,
				     const Real refine_fraction  = 1.0,
				     const Real coarsen_fraction = 0.0,
				     const unsigned int max_level = static_cast<unsigned int>(-1));
		      
  /**
   * Refines and coarsens user-requested elements. Will also
   * refine/coarsen additional elements to satisy level-one rule.
   */
  void refine_and_coarsen_elements (const bool maintain_level_one=true);
  
  /**
   * Coarsens user-requested elements.
   */
  void coarsen_elements ();
  
  /**
   * Refines user-requested elements.
   */
  void refine_elements ();
  
  /**
   * Uniformly refines the mesh \p n times.
   */
  void uniformly_refine (unsigned int n=1);
  
  /**
   * Add point \p p to the mesh. The function returns a pointer to
   * the new node.  The \p key tells the method where to look.
   */
  Node* add_point (const Point& p, const unsigned int key);
  
  /**
   * @returns The index of the next unused element number
   * in the \p element vector.  If all the entries are
   * full it returns the size of the vector.
   */
  unsigned int new_element_number ();




private:



  //------------------------------------------------------
  // "Smoothing" algorthms for refined meshes
  
  /**
   * This algorithm restricts the maximim level mismatch
   * at any node in the mesh.  Calling this with \p max_mismatch
   * equal to 1 would transform this mesh:
   \verbatim
   o---o---o---o---o-------o-------o
   |   |   |   |   |       |       |
   |   |   |   |   |       |       |
   o---o---o---o---o       |       |
   |   |   |   |   |       |       |
   |   |   |   |   |       |       |
   o---o---o---o---o-------o-------o
   |   |   |   |   |       |       |
   |   |   |   |   |       |       |
   o---o---o---o---o       |       |
   |   |   |   |   |       |       |
   |   |   |   |   |       |       |
   o---o---o---o---o-------o-------o
   |       |       |               |
   |       |       |               |
   |       |       |               |
   |       |       |               |
   |       |       |               |
   o-------o-------o               |
   |       |       |               |
   |       |       |               |
   |       |       |               |
   |       |       |               |
   |       |       |               |
   o-------o-------o---------------o     
   \endverbatim

   * into this:
   
   \verbatim
   o---o---o---o---o-------o-------o
   |   |   |   |   |       |       |
   |   |   |   |   |       |       |
   o---o---o---o---o       |       |
   |   |   |   |   |       |       |
   |   |   |   |   |       |       |
   o---o---o---o---o-------o-------o
   |   |   |   |   |       |       |
   |   |   |   |   |       |       |
   o---o---o---o---o       |       |
   |   |   |   |   |       |       |
   |   |   |   |   |       |       |
   o---o---o---o---o-------o-------o
   |       |       |       :       |
   |       |       |       :       |
   |       |       |       :       |
   |       |       |       :       |
   |       |       |       :       |
   o-------o-------o.......o.......o
   |       |       |       :       |
   |       |       |       :       |
   |       |       |       :       |
   |       |       |       :       |
   |       |       |       :       |
   o-------o-------o-------o-------o     
   \endverbatim
   by refining the indicated element
   
   */
  bool limit_level_mismatch_at_node (const unsigned int max_mismatch);

  /**
   * This algorithm selects an element for refinement
   * if all of its neighbors are (or will be) refined.
   * This algorithm will transform this mesh:
   \verbatim
   o---o---o---o---o---o---o
   |   |   |   |   |   |   |
   |   |   |   |   |   |   |
   o---o---o---o---o---o---o
   |   |   |   |   |   |   |
   |   |   |   |   |   |   |
   o---o---o---o---o---o---o
   |   |   |       |   |   |
   |   |   |       |   |   |
   o---o---o       o---o---o
   |   |   |       |   |   |
   |   |   |       |   |   |
   o---o---o---o---o---o---o
   |   |   |   |   |   |   |
   |   |   |   |   |   |   |
   o---o---o---o---o---o---o
   |   |   |   |   |   |   |
   |   |   |   |   |   |   |
   o---o---o---o---o---o---o
   \endverbatim

   into this:
   \verbatim
   o---o---o---o---o---o---o
   |   |   |   |   |   |   |
   |   |   |   |   |   |   |
   o---o---o---o---o---o---o
   |   |   |   |   |   |   |
   |   |   |   |   |   |   |
   o---o---o---o---o---o---o
   |   |   |   :   |   |   |
   |   |   |   :   |   |   |
   o---o---o...o...o---o---o
   |   |   |   :   |   |   |
   |   |   |   :   |   |   |
   o---o---o---o---o---o---o
   |   |   |   |   |   |   |
   |   |   |   |   |   |   |
   o---o---o---o---o---o---o
   |   |   |   |   |   |   |
   |   |   |   |   |   |   |
   o---o---o---o---o---o---o
   \endverbatim

   by refining the indicated element
   
   */
  bool eliminate_unrefined_patches ();
				     
  /**
   * @returns The index of the next unused node number
   * in the \p nodes vector.  If all the entries are
   * full it returns the size of the vector.
   */
  unsigned int new_node_number ();



  //---------------------------------------------
  // Utility algorithms
  
  /**
   * Updates the \p unused_elements data structure to
   * be compatible with the current state of the mesh.
   */
  void update_unused_elements ();
  
  /**
   * Updates the \p unused_nodes data structure to
   * be compatible with the current state of the mesh.
   */
  void update_unused_nodes () {}

  /**
   * Sets the refinement flag to \p Elem::DO_NOTHING
   * for each element in the mesh.
   */
  void clean_refinement_flags ();
  
  /**
   * Take user-specified coarsening flags and augment them
   * so that level-one dependency is satisfied.
   */
  bool make_coarsening_compatible (const bool);

  /**
   * Take user-specified refinement flags and augment them
   * so that level-one dependency is satisfied.
   */
  bool make_refinement_compatible (const bool);
  
  /**
   * Data structure that holds the new nodes information.
   * The key is a pointer to the element that created the node.
   * For efficiency we will use a hashed multimap if it is
   * available, otherwise a regular multimap.
   */
#if   defined(HAVE_HASH_MAP)    
    typedef std::hash_multimap<unsigned int, Node*> map_type;    
#elif defined(HAVE_EXT_HASH_MAP)
# if  __GNUC__ >= 3
    typedef __gnu_cxx::hash_multimap<unsigned int, Node*> map_type;
# else
    DIE A HORRIBLE DEATH
# endif
#else
    typedef std::multimap<unsigned int, Node*> map_type;
#endif
  
  map_type _new_nodes_map;

  /**
   * Data structure that holds the indices of elements
   * that have been removed from the mesh but not
   * yet removed from the \p _elements vector.
   */
  std::set<unsigned int> _unused_elements;

  /**
   * Reference to the mesh.
   */
  MeshBase& mesh;
};

#endif // end #ifdef ENABLE_AMR
#endif // end #ifndef __mesh_refinement_h__ 
