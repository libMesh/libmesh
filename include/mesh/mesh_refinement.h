// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_MESH_REFINEMENT_H
#define LIBMESH_MESH_REFINEMENT_H



#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_ENABLE_AMR

// Local Includes -----------------------------------
#include "libmesh/libmesh_common.h"
#include "libmesh/libmesh.h" // libMesh::invalid_uint
#include "libmesh/topology_map.h"
#include "libmesh/parallel_object.h"

// C++ Includes   -----------------------------------
#include <vector>

namespace libMesh
{

// Forward Declarations -----------------------------
class MeshBase;
class Point;
class Node;
class ErrorVector;
class PeriodicBoundaries;
class Elem;
class PointLocatorBase;

/**
 * This is the \p MeshRefinement class.  This class implements
 * adaptive mesh refinement algorithms for a \p MeshBase.
 *
 * \author Benjamin S. Kirk
 * \date 2002-2007
 */
class MeshRefinement : public ParallelObject
{
public:

  /**
   * Constructor.
   */
  explicit
  MeshRefinement (MeshBase & mesh);

private:
  // Both the copy ctor and the assignment operator are
  // declared private but not implemented.  This is the
  // standard practice to prevent them from being used.
  MeshRefinement (const MeshRefinement &);
  MeshRefinement & operator=(const MeshRefinement &);

public:

  /**
   * Abstract base class to be used for user-specified
   * element flagging.  This can be used instead of or to
   * augment traditional error indicator based refinement.
   * This simply provides a base class that can be derived
   * from and then passed to the
   * \p flag_elements_by () method.
   */
  class ElementFlagging
  {
  public:
    /**
     * Destructor.  Virtual because we will have virtual functions.
     */
    virtual ~ElementFlagging () {}

    /**
     * Callback function to be used for marking elements for refinement.
     */
    virtual void flag_elements () = 0;
  };

  /**
   * Sets the \p PeriodicBoundaries pointer.
   */
  void set_periodic_boundaries_ptr(PeriodicBoundaries * pb_ptr);

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
   *
   * All the function arguments except error_per_cell
   * have been deprecated, and will be removed in
   * future libMesh releases - to control these parameters,
   * set the corresponding member variables.
   */
  void flag_elements_by_error_fraction (const ErrorVector & error_per_cell,
                                        const Real refine_fraction  = 0.3,
                                        const Real coarsen_fraction = 0.0,
                                        const unsigned int max_level = libMesh::invalid_uint);

  /**
   * Flags elements for coarsening and refinement based on
   * the computed error passed in \p error_per_cell.  This method refines
   * the worst elements with errors greater than
   * \p absolute_global_tolerance / n_active_elem, flagging at most
   * \p refine_fraction * n_active_elem
   * It coarsens elements with errors less than
   * \p coarsen_threshold * \p global_tolerance / n_active_elem,
   * flagging at most
   * \p coarsen_fraction * n_active_elem
   *
   * The three fractions \p refine_fraction \p coarsen_fraction and
   * \p coarsen_threshold should be in \f$ [0,1] \f$.
   */
  void flag_elements_by_error_tolerance (const ErrorVector & error_per_cell);

  /**
   * Flags elements for coarsening and refinement based on
   * the computed error passed in \p error_per_cell.  This method attempts to
   * produce a mesh with slightly more than \p nelem_target active elements,
   * trading element refinement for element coarsening when their error
   * ratios exceed \p coarsen_threshold.  It flags no more than
   * \p refine_fraction * n_elem elements for refinement and flags no
   * more than \p coarsen_fraction * n_elem elements for coarsening.
   * This method returns true if it has done all the AMR/C it can do
   * in a single step, or false if further adaptive steps may be required
   * to produce a mesh with a narrow error distribution and the right
   * number of elements.
   */
  bool flag_elements_by_nelem_target (const ErrorVector & error_per_cell);

  /**
   * Flags elements for coarsening and refinement based on
   * the computed error passed in \p error_per_cell.  This method picks
   * the top \p refine_fraction * \p n_elem elements for refinement and
   * the bottom \p coarsen_fraction * \p n_elem elements for coarsening.
   * The two fractions \p refine_fraction and \p coarsen_fraction must be
   * in \f$ [0,1] \f$.
   *
   * All the function arguments except error_per_cell
   * have been deprecated, and will be removed in
   * future libMesh releases - to control these parameters,
   * set the corresponding member variables.
   */
  void flag_elements_by_elem_fraction (const ErrorVector & error_per_cell,
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
   *
   * All the function arguments except error_per_cell
   * have been deprecated, and will be removed in
   * future libMesh releases - to control these parameters,
   * set the corresponding member variables.
   */
  void flag_elements_by_mean_stddev (const ErrorVector & error_per_cell,
                                     const Real refine_fraction  = 1.0,
                                     const Real coarsen_fraction = 0.0,
                                     const unsigned int max_level = libMesh::invalid_uint);

  /**
   * Flag elements based on a function object.  The class \p ElementFlagging
   * defines a mechanism for implementing refinement strategies.
   */
  void flag_elements_by (ElementFlagging & element_flagging);

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
   * refine/coarsen additional elements to satisfy level-one rule.
   * It is possible that for a given set of refinement flags there
   * is actually no change upon calling this member function.  Consequently,
   * this function returns \p true if the mesh actually changed (hence
   * data needs to be projected) and \p false otherwise.
   *
   * This function used to take an argument, \p maintain_level_one -
   * new code should use face_level_mismatch_limit() instead.
   */
  bool refine_and_coarsen_elements ();

  /**
   * Only coarsens the user-requested elements. Some elements
   * will not be coarsened to satisfy the level one rule.
   * It is possible that for a given set of refinement flags there
   * is actually no change upon calling this member function.  Consequently,
   * this function returns \p true if the mesh actually changed (hence
   * data needs to be projected) and \p false otherwise.
   *
   * This function used to take an argument, \p maintain_level_one -
   * new code should use face_level_mismatch_limit() instead.
   */
  bool coarsen_elements ();

  /**
   * Only refines the user-requested elements.
   * It is possible that for a given set of refinement flags there
   * is actually no change upon calling this member function.  Consequently,
   * this function returns \p true if the mesh actually changed (hence
   * data needs to be projected) and \p false otherwise.
   *
   * This function used to take an argument, \p maintain_level_one -
   * new code should use face_level_mismatch_limit() instead.
   */
  bool refine_elements ();

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
   * Sets the refinement flag to \p Elem::DO_NOTHING
   * for each element in the mesh.
   */
  void clean_refinement_flags ();

  /**
   * Returns true if and only if the mesh is level one smooth
   * Returns false otherwise
   * Aborts the program if libmesh_assert_yes is true and
   * the mesh is not level one smooth
   */
  bool test_level_one (bool libmesh_assert_yes = false);

  /**
   * Returns true if and only if the mesh has no elements
   * flagged to be coarsened or refined
   * Returns false otherwise
   * Aborts the program if libmesh_assert_yes is true and
   * the mesh has flagged elements
   */
  bool test_unflagged (bool libmesh_assert_yes = false);

  /**
   * Add a node to the mesh.  The node should be node n of child c of
   * parent Elem parent.  The function returns a pointer to a suitable
   * existing node, or creates a new node and returns a pointer to it
   * if necessary.
   * The processor_id is assigned to any newly created node.
   */
  Node * add_node (const Elem & parent,
                   unsigned int child,
                   unsigned int node,
                   processor_id_type proc_id);

  /**
   * Adds the element \p elem to the mesh.
   */
  Elem * add_elem (Elem * elem);

  /**
   * @returns a constant reference to the \p MeshBase object associated
   * with this object.
   */
  const MeshBase & get_mesh () const { return _mesh; }

  /**
   * @returns a writeable reference to the \p MeshBase object associated
   * with this object.
   */
  MeshBase & get_mesh () { return _mesh; }

  /**
   * If \p coarsen_by_parents is true, complete groups of sibling elements
   * (elements with the same parent) will be flagged for coarsening.
   * This should make the coarsening more likely to occur as requested.
   *
   * \p coarsen_by_parents is true by default.
   */
  bool & coarsen_by_parents();

  /**
   * The \p refine_fraction sets either a desired target or a desired
   * maximum number of elements to flag for refinement, depending on which
   * flag_elements_by method is called.
   *
   * \p refine_fraction must be in \f$ [0,1] \f$, and is 0.3 by default.
   */
  Real & refine_fraction();

  /**
   * The \p coarsen_fraction sets either a desired target or a desired
   * maximum number of elements to flag for coarsening, depending on which
   * flag_elements_by method is called.
   *
   * \p coarsen_fraction must be in \f$ [0,1] \f$, and is 0 by default.
   */
  Real & coarsen_fraction();

  /**
   * The \p max_h_level is the greatest refinement level an element should
   * reach.
   *
   * \p max_h_level is unlimited (libMesh::invalid_uint) by default
   */
  unsigned int & max_h_level();

  /**
   * The \p coarsen_threshold provides hysteresis in AMR/C strategies.
   * Refinement of elements with error estimate E will be done even
   * at the expense of coarsening elements whose children's accumulated
   * error does not exceed \p coarsen_threshold * E.
   *
   * \p coarsen_threshold must be in \f$ [0,1] \f$, and is 0.1 by default.
   */
  Real & coarsen_threshold();

  /**
   * If \p nelem_target is set to a nonzero value, methods like
   * flag_elements_by_nelem_target() will attempt to keep the number
   * of active elements in the mesh close to nelem_target.
   *
   * \p nelem_target is 0 by default.
   */
  dof_id_type & nelem_target();

  /**
   * If \p absolute_global_tolerance is set to a nonzero value, methods
   * like flag_elements_by_global_tolerance() will attempt to reduce
   * the global error of the mesh (defined as the square root of the
   * sum of the squares of the errors on active elements) to below
   * this tolerance.
   *
   * \p absolute_global_tolerance is 0 by default.
   */
  Real & absolute_global_tolerance();

  /**
   * If \p face_level_mismatch_limit is set to a nonzero value, then
   * refinement and coarsening will produce meshes in which the
   * refinement level of two face neighbors will not differ by more than
   * that limit.  If \p face_level_mismatch_limit is 0, then level
   * differences will be unlimited.
   *
   * \p face_level_mismatch_limit is 1 by default.  Currently the only
   * supported options are 0 and 1.
   */
  unsigned char & face_level_mismatch_limit();

  /**
   * If \p edge_level_mismatch_limit is set to a nonzero value, then
   * refinement and coarsening will produce meshes in which the
   * refinement level of two edge neighbors will not differ by more than
   * that limit.  If \p edge_level_mismatch_limit is 0, then level
   * differences will be unlimited.
   *
   * \p edge_level_mismatch_limit is 0 by default.
   */
  unsigned char & edge_level_mismatch_limit();

  /**
   * If \p node_level_mismatch_limit is set to a nonzero value, then
   * refinement and coarsening will produce meshes in which the
   * refinement level of two nodal neighbors will not differ by more than
   * that limit.  If \p node_level_mismatch_limit is 0, then level
   * differences will be unlimited.
   *
   * \p node_level_mismatch_limit is 0 by default.
   */
  unsigned char & node_level_mismatch_limit();

  /**
   * If \p overrefined_boundary_limit is set to a nonnegative value,
   * then refinement and coarsening will produce meshes in which the
   * refinement level of a boundary element is no more than that many
   * levels greater than the level of any of its interior neighbors.
   *
   * This may be counter-intuitive in the 1D-embedded-in-3D case: an
   * edge has *more* interior neighbors than a face containing that
   * edge.
   *
   * If \p overrefined_boundary_limit is negative, then level
   * differences will be unlimited.
   *
   * \p overrefined_boundary_limit is 0 by default.  This implies that
   * adaptive coarsening can only be done on an interior element if
   * any boundary elements on its sides are simultaneously coarsened.
   */
  signed char & overrefined_boundary_limit();

  /**
   * If \p underrefined_boundary_limit is set to a nonnegative value,
   * then refinement and coarsening will produce meshes in which the
   * refinement level of an element is no more than that many
   * levels greater than the level of any boundary elements on its
   * sides.
   *
   * If \p underrefined_boundary_limit is negative, then level
   * differences will be unlimited.
   *
   * \p underrefined_boundary_limit is 0 by default.  This implies that
   * adaptive coarsening can only be done on a boundary element if
   * any interior elements it is on the side of are simultaneously
   * coarsened.
   */
  signed char & underrefined_boundary_limit();


  /**
   * Copy refinement flags on ghost elements from their
   * local processors.  Return true if any flags changed.
   */
  bool make_flags_parallel_consistent ();

  /**
   * Returns the state of the _enforce_mismatch_limit_prior_to_refinement flag.
   * Defaults to false.
   * Deprecated - use enforce_mismatch_limit_prior_to_refinement() instead.
   */
  bool get_enforce_mismatch_limit_prior_to_refinement();

  /**
   * Set _enforce_mismatch_limit_prior_to_refinement option.
   * Defaults to false.
   * Deprecated - use enforce_mismatch_limit_prior_to_refinement() instead.
   */
  void set_enforce_mismatch_limit_prior_to_refinement(bool enforce);

  /**
   * Get/set the _enforce_mismatch_limit_prior_to_refinement flag.
   * The default value for this flag is false.
   */
  bool & enforce_mismatch_limit_prior_to_refinement();

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

  /**
   * Smooths refinement flags according to current settings.
   *
   * It is possible that for a given set of refinement flags there
   * is actually no change upon calling this member function.  Consequently,
   * this function returns \p true if the flags actually changed (hence
   * data needs to be projected) and \p false otherwise.
   */
  void _smooth_flags (bool refining, bool coarsening);

  //------------------------------------------------------
  // "Smoothing" algorthms for refined meshes

  /**
   * This algorithm restricts the maximum level mismatch
   * at any node in the mesh.  Calling this with \p max_mismatch
   * equal to 1 would transform this mesh:
   * \verbatim
   * o---o---o---o---o-------o-------o
   * |   |   |   |   |       |       |
   * |   |   |   |   |       |       |
   * o---o---o---o---o       |       |
   * |   |   |   |   |       |       |
   * |   |   |   |   |       |       |
   * o---o---o---o---o-------o-------o
   * |   |   |   |   |       |       |
   * |   |   |   |   |       |       |
   * o---o---o---o---o       |       |
   * |   |   |   |   |       |       |
   * |   |   |   |   |       |       |
   * o---o---o---o---o-------o-------o
   * |       |       |               |
   * |       |       |               |
   * |       |       |               |
   * |       |       |               |
   * |       |       |               |
   * o-------o-------o               |
   * |       |       |               |
   * |       |       |               |
   * |       |       |               |
   * |       |       |               |
   * |       |       |               |
   * o-------o-------o---------------o
   * \endverbatim
   *
   * into this:
   *
   * \verbatim
   * o---o---o---o---o-------o-------o
   * |   |   |   |   |       |       |
   * |   |   |   |   |       |       |
   * o---o---o---o---o       |       |
   * |   |   |   |   |       |       |
   * |   |   |   |   |       |       |
   * o---o---o---o---o-------o-------o
   * |   |   |   |   |       |       |
   * |   |   |   |   |       |       |
   * o---o---o---o---o       |       |
   * |   |   |   |   |       |       |
   * |   |   |   |   |       |       |
   * o---o---o---o---o-------o-------o
   * |       |       |       :       |
   * |       |       |       :       |
   * |       |       |       :       |
   * |       |       |       :       |
   * |       |       |       :       |
   * o-------o-------o.......o.......o
   * |       |       |       :       |
   * |       |       |       :       |
   * |       |       |       :       |
   * |       |       |       :       |
   * |       |       |       :       |
   * o-------o-------o-------o-------o
   * \endverbatim
   * by refining the indicated element
   */
  bool limit_level_mismatch_at_node (const unsigned int max_mismatch);

  /*
   * This algorithm restricts the maximum level mismatch
   * at any edge in the mesh.  See the ASCII art in the comment of
   * limit_level_mismatch_at_node, and pretend they're hexes.
   */
  bool limit_level_mismatch_at_edge (const unsigned int max_mismatch);

  /*
   * This algorithm flags interior elements for refinement as needed
   * to prevent corresponding boundary element refinement mismatch
   * from exceeding the given limit.
   */
  bool limit_overrefined_boundary (const signed char max_mismatch);

  /*
   * This algorithm flags boundary elements for refinement as needed
   * to prevent corresponding interior element refinement mismatch
   * from exceeding the given limit.
   */
  bool limit_underrefined_boundary (const signed char max_mismatch);

  /**
   * This algorithm selects an element for refinement
   * if all of its neighbors are (or will be) refined.
   * This algorithm will transform this mesh:
   * \verbatim
   * o---o---o---o---o---o---o
   * |   |   |   |   |   |   |
   * |   |   |   |   |   |   |
   * o---o---o---o---o---o---o
   * |   |   |   |   |   |   |
   * |   |   |   |   |   |   |
   * o---o---o---o---o---o---o
   * |   |   |       |   |   |
   * |   |   |       |   |   |
   * o---o---o       o---o---o
   * |   |   |       |   |   |
   * |   |   |       |   |   |
   * o---o---o---o---o---o---o
   * |   |   |   |   |   |   |
   * |   |   |   |   |   |   |
   * o---o---o---o---o---o---o
   * |   |   |   |   |   |   |
   * |   |   |   |   |   |   |
   * o---o---o---o---o---o---o
   * \endverbatim
   *
   * into this:
   * \verbatim
   * o---o---o---o---o---o---o
   * |   |   |   |   |   |   |
   * |   |   |   |   |   |   |
   * o---o---o---o---o---o---o
   * |   |   |   |   |   |   |
   * |   |   |   |   |   |   |
   * o---o---o---o---o---o---o
   * |   |   |   :   |   |   |
   * |   |   |   :   |   |   |
   * o---o---o...o...o---o---o
   * |   |   |   :   |   |   |
   * |   |   |   :   |   |   |
   * o---o---o---o---o---o---o
   * |   |   |   |   |   |   |
   * |   |   |   |   |   |   |
   * o---o---o---o---o---o---o
   * |   |   |   |   |   |   |
   * |   |   |   |   |   |   |
   * o---o---o---o---o---o---o
   * \endverbatim
   *
   * by refining the indicated element
   */
  bool eliminate_unrefined_patches ();


  //---------------------------------------------
  // Utility algorithms

  /**
   * Calculates the error on all coarsenable parents.
   * error_per_parent[parent_id] stores this error if parent_id corresponds
   * to a coarsenable parent, and stores -1 otherwise.
   */
  void create_parent_error_vector (const ErrorVector & error_per_cell,
                                   ErrorVector & error_per_parent,
                                   Real & parent_error_min,
                                   Real & parent_error_max);

  /**
   * Updates the \p _new_nodes_map
   */
  void update_nodes_map ();

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
   * Local dispatch function for getting the correct topological
   * neighbor from the Elem class
   */
  Elem * topological_neighbor (Elem * elem,
                               const PointLocatorBase * point_locator,
                               const unsigned int side);

  /**
   * Local dispatch function for checking the correct has_neighbor
   * function from the Elem class
   */
  bool has_topological_neighbor (Elem * elem,
                                 const PointLocatorBase * point_locator,
                                 Elem * neighbor);

  /**
   * Data structure that holds the new nodes information.
   */
  TopologyMap _new_nodes_map;

  /**
   * Reference to the mesh.
   */
  MeshBase & _mesh;

  /**
   * For backwards compatibility, we initialize this
   * as false and then set it to true if the user uses
   * any of the refinement parameter accessor functions
   */
  bool _use_member_parameters;

  /**
   * Refinement parameter values
   */

  bool _coarsen_by_parents;

  Real _refine_fraction;

  Real _coarsen_fraction;

  unsigned int _max_h_level;

  Real _coarsen_threshold;

  dof_id_type _nelem_target;

  Real _absolute_global_tolerance;

  unsigned char _face_level_mismatch_limit;
  unsigned char _edge_level_mismatch_limit;
  unsigned char _node_level_mismatch_limit;

  signed char _overrefined_boundary_limit;
  signed char _underrefined_boundary_limit;

  /**
   * This option enforces the mismatch level prior to refinement by checking
   * if refining any element marked for refinement \b would cause a mismatch
   * greater than the limit. Applies to all mismatch methods.
   *
   * Calling this with \p node_level_mismatch_limit() = 1
   * would transform this mesh:
   * \verbatim
   * o-------o-------o-------o-------o
   * |       |       |       |       |
   * |       |       |       |       |
   * |       |       |       |       |
   * |       |       |       |       |
   * |       |       |       |       |
   * o-------o---o---o-------o-------o
   * |       |   :   |       |       |
   * |       |   :   |       |       |
   * |       o...o...o       |       |
   * |       |   :   |       |       |
   * |       |   :   |       |       |
   * o-------o---o---o-------o-------o
   * |       |       |               |
   * |       |       |               |
   * |       |       |               |
   * |       |       |               |
   * |       |       |               |
   * o-------o-------o               |
   * |       |       |               |
   * |       |       |               |
   * |       |       |               |
   * |       |       |               |
   * |       |       |               |
   * o-------o-------o---------------o
   * \endverbatim
   *
   * into this:
   *
   * \verbatim
   * o-------o-------o-------o-------o
   * |       |       |       |       |
   * |       |       |       |       |
   * |       |       |       |       |
   * |       |       |       |       |
   * |       |       |       |       |
   * o-------o-------o-------o-------o
   * |       |       |       |       |
   * |       |       |       |       |
   * |       |       |       |       |
   * |       |       |       |       |
   * |       |       |       |       |
   * o-------o-------o-------o-------o
   * |       |       |       :       |
   * |       |       |       :       |
   * |       |       |       :       |
   * |       |       |       :       |
   * |       |       |       :       |
   * o-------o-------o.......o.......o
   * |       |       |       :       |
   * |       |       |       :       |
   * |       |       |       :       |
   * |       |       |       :       |
   * |       |       |       :       |
   * o-------o-------o-------o-------o
   * \endverbatim
   * by moving the refinement flag to the indicated element.
   *
   * Default value is false.
   */
  bool _enforce_mismatch_limit_prior_to_refinement;

  /**
   * This helper function enforces the desired mismatch limits prior
   * to refinement.  It is called from the
   * MeshRefinement::limit_level_mismatch_at_edge() and
   * MeshRefinement::limit_level_mismatch_at_node() functions.
   * Returns true if this enforcement caused the refinement flags for
   * elem to change, false otherwise.
   */
  enum NeighborType {POINT, EDGE};
  bool enforce_mismatch_limit_prior_to_refinement(Elem * elem,
                                                  NeighborType nt,
                                                  unsigned max_mismatch);

#ifdef LIBMESH_ENABLE_PERIODIC
  PeriodicBoundaries * _periodic_boundaries;
#endif
};



// ------------------------------------------------------------
// MeshRefinement class inline members

inline bool & MeshRefinement::coarsen_by_parents()
{
  _use_member_parameters = true;
  return _coarsen_by_parents;
}

inline Real & MeshRefinement::refine_fraction()
{
  _use_member_parameters = true;
  return _refine_fraction;
}

inline Real & MeshRefinement::coarsen_fraction()
{
  _use_member_parameters = true;
  return _coarsen_fraction;
}

inline unsigned int & MeshRefinement::max_h_level()
{
  _use_member_parameters = true;
  return _max_h_level;
}

inline Real & MeshRefinement::coarsen_threshold()
{
  _use_member_parameters = true;
  return _coarsen_threshold;
}

inline dof_id_type & MeshRefinement::nelem_target()
{
  _use_member_parameters = true;
  return _nelem_target;
}

inline Real & MeshRefinement::absolute_global_tolerance()
{
  _use_member_parameters = true;
  return _absolute_global_tolerance;
}

inline unsigned char & MeshRefinement::face_level_mismatch_limit()
{
  return _face_level_mismatch_limit;
}

inline unsigned char & MeshRefinement::edge_level_mismatch_limit()
{
  return _edge_level_mismatch_limit;
}

inline unsigned char & MeshRefinement::node_level_mismatch_limit()
{
  return _node_level_mismatch_limit;
}

inline signed char & MeshRefinement::overrefined_boundary_limit()
{
  return _overrefined_boundary_limit;
}

inline signed char & MeshRefinement::underrefined_boundary_limit()
{
  return _underrefined_boundary_limit;
}

inline bool MeshRefinement::get_enforce_mismatch_limit_prior_to_refinement()
{
  libmesh_deprecated();
  return enforce_mismatch_limit_prior_to_refinement();
}

inline void MeshRefinement::set_enforce_mismatch_limit_prior_to_refinement(bool enforce)
{
  libmesh_deprecated();
  enforce_mismatch_limit_prior_to_refinement() = enforce;
}

inline bool & MeshRefinement::enforce_mismatch_limit_prior_to_refinement()
{
  return _enforce_mismatch_limit_prior_to_refinement;
}



} // namespace libMesh

#endif // end #ifdef LIBMESH_ENABLE_AMR
#endif // LIBMESH_MESH_REFINEMENT_H
