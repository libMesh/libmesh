// $Id: mesh_refinement.h,v 1.18 2006-06-12 23:26:04 roystgnr Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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



#include "libmesh_config.h"

#ifdef ENABLE_AMR

// C++ Includes   -----------------------------------
#include <vector>
#include <list>

#if   defined(HAVE_HASH_MAP)
# include <hash_map>
#elif defined(HAVE_EXT_HASH_MAP)
# include <ext/hash_map>
#else
# include <map>
#endif

// Local Includes -----------------------------------
#include "libmesh_common.h"
#include "libmesh.h" // libMesh::invalid_uint

//#include "node.h" // remove this when debugging is done, only need forward dec.

class MeshBase;
class Point;
class Node;
class Elem;
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

private:
  // Both the copy ctor and the assignment operator are 
  // declared private but not implemented.  This is the 
  // standard practice to prevent them from being used.
  MeshRefinement (const MeshRefinement&);
  MeshRefinement& operator=(const MeshRefinement&);

public:

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
					const unsigned int max_level = libMesh::invalid_uint);

  /**
   * Flags elements for coarsening and refinement based on
   * the computed error passed in \p error_per_cell.  This method refines
   * the worst elements with errors greater than
   * \p global_tolerance / n_elem, flagging at most
   * \p refine_fraction * n_elem.  It coarsens elements with errors less than
   * \p coarsen_threshold * \p global_tolerance / n_elem, flagging at most
   * \p coarsen_fraction * n_elem.
   *
   * The three fractions \p refine_fraction \p coarsen_fraction and
   * \p coarsen_threshold should be in \f$ [0,1] \f$.
   */
  void flag_elements_by_error_tolerance (const ErrorVector& error_per_cell,
					 const Real global_tolerance,
					 const Real coarsen_threshold = 0.1,
					 const Real refine_fraction  = 0.3,
					 const Real coarsen_fraction = 0.3);

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
				       const unsigned int max_level = libMesh::invalid_uint);

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
				     const unsigned int max_level = libMesh::invalid_uint);

  /**
   * Takes a mesh whose elements are flagged for h refinement and coarsening,
   * and switches those flags to request p refinement and coarsening instead.
   */
  void switch_h_to_p_refinement();

  /**
   * Takes a mesh whose elements are flagged for h refinement and coarsening,
   * and adds flags to request p refinement and coarsening of the same elements.
   */
  void add_p_to_h_refinement();

  /**
   * Refines and coarsens user-requested elements. Will also
   * refine/coarsen additional elements to satisy level-one rule.
   * It is possible that for a given set of refinement flags there
   * is actually no change upon calling this member function.  Consequently,
   * this function returns \p true if the mesh actually changed (hence
   * data needs to be projected) and \p false otherwise.
   */
  bool refine_and_coarsen_elements (const bool maintain_level_one=true);

  /**
   * Only coarsens the user-requested elements. Some elements
   * will not be coarsened to satisfy the level one rule.
   * It is possible that for a given set of refinement flags there
   * is actually no change upon calling this member function.  Consequently,
   * this function returns \p true if the mesh actually changed (hence
   * data needs to be projected) and \p false otherwise.
   */
  bool coarsen_elements (const bool maintain_level_one=true);

  /**
   * Only refines the user-requested elements. 
   * It is possible that for a given set of refinement flags there
   * is actually no change upon calling this member function.  Consequently,
   * this function returns \p true if the mesh actually changed (hence
   * data needs to be projected) and \p false otherwise.
   */
  bool refine_elements (const bool maintain_level_one=true);
  
  /**
   * Uniformly refines the mesh \p n times.
   */
  void uniformly_refine (unsigned int n=1);
  
  /**
   * Attempts to uniformly coarsen the mesh \p n times.
   */
  void uniformly_coarsen (unsigned int n=1);
  
  /**
   * Uniformly p refines the mesh \p n times.
   */
  void uniformly_p_refine (unsigned int n=1);
  
  /**
   * Attempts to uniformly p coarsen the mesh \p n times.
   */
  void uniformly_p_coarsen (unsigned int n=1);

  /**
   * Returns true if and only if the mesh is level one smooth
   * Returns false otherwise
   * Aborts the program if assert_yes is true and 
   * the mesh is not level one smooth
   */
  bool test_level_one (bool assert_yes = false);
  
  /**
   * Returns true if and only if the mesh has no elements
   * flagged to be coarsened or refined
   * Returns false otherwise
   * Aborts the program if assert_yes is true and 
   * the mesh has flagged elements
   */
  bool test_unflagged (bool assert_yes = false);
  
  /**
   * Add point \p p to the mesh. The function returns a pointer to
   * the new node.  The \p key tells the method where to look.
   */
  Node* add_point (const Point& p, const unsigned int key);

  /**
   * Adds the element \p elem to the mesh.
   */
  Elem* add_elem (Elem* elem);

  /**
   * @returns a constant reference to the \p MeshBase object associated
   * with this object.
   */
  const MeshBase& get_mesh () const { return _mesh; }

  /**
   * @returns a writeable reference to the \p MeshBase object associated
   * with this object.
   */
  MeshBase&       get_mesh ()       { return _mesh; }



private:

  /**
   * Coarsens user-requested elements.  Both coarsen_elements
   * and refine_elements used to be in the public interface for the
   * MeshRefinement object.  Unfortunately, without proper
   * preparation (make_refinement_compatible, make_coarsening_compatible)
   * at least coarsen_elements() did not work alone.  By making them
   * private, we signal to the user that they are not part of the
   * interface.
   *
   * It is possible that for a given set of refinement flags there
   * is actually no change upon calling this member function.  Consequently,
   * this function returns \p true if the mesh actually changed (hence
   * data needs to be projected) and \p false otherwise.
   */
  bool _coarsen_elements ();
  
  /**
   * Refines user-requested elements.
   *
   * It is possible that for a given set of refinement flags there
   * is actually no change upon calling this member function.  Consequently,
   * this function returns \p true if the mesh actually changed (hence
   * data needs to be projected) and \p false otherwise.
   */
  bool _refine_elements ();



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


  //---------------------------------------------
  // Utility algorithms

  /**
   * Updates the \p _new_nodes_map
   */
  void update_nodes_map ();

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
# if   (__GNUC__ == 3) && (__GNUC_MINOR__ == 0) // gcc 3.0   
    typedef std::hash_multimap<unsigned int, Node*> map_type;
# elif (__GNUC__ >= 3)                          // gcc 3.1 & newer
    typedef __gnu_cxx::hash_multimap<unsigned int, Node*> map_type;
# else
    DIE A HORRIBLE DEATH
# endif
#else
    typedef std::multimap<unsigned int, Node*> map_type;
#endif
  
  map_type _new_nodes_map;


  /**
   * Reference to the mesh.
   */
  MeshBase& _mesh;

};

#endif // end #ifdef ENABLE_AMR
#endif // end #ifndef __mesh_refinement_h__ 
