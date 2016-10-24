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



#ifndef LIBMESH_MESH_BASE_H
#define LIBMESH_MESH_BASE_H



// Local Includes -----------------------------------
#include "libmesh/auto_ptr.h"
#include "libmesh/boundary_info.h"
#include "libmesh/dof_object.h" // for invalid_processor_id
#include "libmesh/enum_elem_type.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/multi_predicates.h"
#include "libmesh/partitioner.h" // UniquePtr needs a real declaration
#include "libmesh/point_locator_base.h"
#include "libmesh/variant_filter_iterator.h"
#include "libmesh/parallel_object.h"

// C++ Includes   -----------------------------------
#include <cstddef>
#include <string>

namespace libMesh
{

// forward declarations
class Elem;
class GhostingFunctor;
class Node;
class Point;
class MeshData;

template <class MT>
class MeshInput;


/**
 * This is the \p MeshBase class. This class provides all the data necessary
 * to describe a geometric entity.  It allows for the description of a
 * \p dim dimensional object that lives in \p LIBMESH_DIM-dimensional space.
 * \par
 * A mesh is made of nodes and elements, and this class provides data
 * structures to store and access both.  A mesh may be partitioned into a
 * number of subdomains, and this class provides that functionality.
 * Furthermore, this class provides functions for reading and writing a
 * mesh to disk in various formats.
 *
 * \author  Benjamin S. Kirk
 */
class MeshBase : public ParallelObject
{
public:

  /**
   * Constructor.  Takes \p dim, the dimension of the mesh.
   * The mesh dimension can be changed (and may automatically be
   * changed by mesh generation/loading) later.
   */
  MeshBase (const Parallel::Communicator & comm_in,
            unsigned char dim=1);

#ifndef LIBMESH_DISABLE_COMMWORLD
  /**
   * Deprecated constructor.  Takes \p dim, the dimension of the mesh.
   * The mesh dimension can be changed (and may automatically be
   * changed by mesh generation/loading) later.
   */
  MeshBase (unsigned char dim=1);
#endif

  /**
   * Copy-constructor.
   */
  MeshBase (const MeshBase & other_mesh);

  /**
   * Virtual "copy constructor"
   */
  virtual UniquePtr<MeshBase> clone() const = 0;

  /**
   * Destructor.
   */
  virtual ~MeshBase ();

  /**
   * A partitioner to use at each prepare_for_use()
   */
  virtual UniquePtr<Partitioner> & partitioner() { return _partitioner; }

  /**
   * The information about boundary ids on the mesh
   */
  const BoundaryInfo & get_boundary_info() const { return *boundary_info; }

  /**
   * Writeable information about boundary ids on the mesh
   */
  BoundaryInfo & get_boundary_info() { return *boundary_info; }

  /**
   * Deletes all the data that are currently stored.
   */
  virtual void clear ();

  /**
   * @returns \p true if the mesh has been prepared via a call
   * to \p prepare_for_use, \p false otherwise.
   */
  bool is_prepared () const
  { return _is_prepared; }

  /**
   * @returns \p true if all elements and nodes of the mesh
   * exist on the current processor, \p false otherwise
   */
  virtual bool is_serial () const
  { return true; }

  /**
   * Gathers all elements and nodes of the mesh onto
   * every processor
   */
  virtual void allgather () {}

  /**
   * When supported, deletes all nonlocal elements of the mesh
   * except for "ghosts" which touch a local element, and deletes
   * all nodes which are not part of a local or ghost element
   */
  virtual void delete_remote_elements () {}

  /**
   * @returns the logical dimension of the mesh; i.e. the manifold
   * dimension of the elements in the mesh.  If we ever support
   * multi-dimensional meshes (e.g. hexes and quads in the same mesh)
   * then this will return the largest such dimension.
   */
  unsigned int mesh_dimension () const;

  /**
   * Resets the logical dimension of the mesh. If the mesh has
   * elements of multiple dimensions, this should be set to the largest
   * dimension. E.g. if the mesh has 1D and 2D elements, this should
   * be set to 2. If the mesh has 2D and 3D elements, this should be
   * set to 3.
   */
  void set_mesh_dimension (unsigned char d)
  { _elem_dims.clear(); _elem_dims.insert(d); }

  /**
   * @returns set of dimensions of elements present in the mesh.
   */
  const std::set<unsigned char> & elem_dimensions() const
  { return _elem_dims; }

  /**
   * Returns the "spatial dimension" of the mesh.  The spatial
   * dimension is defined as:
   *
   *   1 - for an exactly x-aligned mesh of 1D elements
   *   2 - for an exactly x-y planar mesh of 2D elements
   *   3 - otherwise
   *
   * No tolerance checks are performed to determine whether the Mesh
   * is x-aligned or x-y planar, only strict equality with zero in the
   * higher dimensions is checked.  Also, x-z and y-z planar meshes are
   * considered to have spatial dimension == 3.
   *
   * The spatial dimension is updated during prepare_for_use() based
   * on the dimensions of the various elements present in the Mesh,
   * but is *never automatically decreased* by this function.
   *
   * For example, if the user calls set_spatial_dimension(2) and then
   * later inserts 3D elements into the mesh,
   * Mesh::spatial_dimension() will return 3 after the next call to
   * prepare_for_use().  On the other hand, if the user calls
   * set_spatial_dimension(3) and then inserts only x-aligned 1D
   * elements into the Mesh, mesh.spatial_dimension() will remain 3.
   */
  unsigned int spatial_dimension () const;

  /**
   * Sets the "spatial dimension" of the Mesh.  See the documentation
   * for Mesh::spatial_dimension() for more information.
   */
  void set_spatial_dimension(unsigned char d);

  /**
   * Returns the number of nodes in the mesh. This function and others must
   * be defined in derived classes since the MeshBase class has no specific
   * storage for nodes or elements.  The standard n_nodes() function
   * may return a cached value on distributed meshes, and so can be
   * called by any processor at any time.
   */
  virtual dof_id_type n_nodes () const = 0;

  /**
   * Returns the number of nodes in the mesh. This function and others must
   * be defined in derived classes since the MeshBase class has no specific
   * storage for nodes or elements.  The parallel_n_nodes() function
   * returns a newly calculated parallel-synchronized value on
   * distributed meshes, and so must be called in parallel only.
   */
  virtual dof_id_type parallel_n_nodes () const = 0;

  /**
   * Returns the number of nodes on processor \p proc.
   */
  dof_id_type n_nodes_on_proc (const processor_id_type proc) const;

  /**
   * Returns the number of nodes on the local processor.
   */
  dof_id_type n_local_nodes () const
  { return this->n_nodes_on_proc (this->processor_id()); }

  /**
   * Returns the number of nodes owned by no processor.
   */
  dof_id_type n_unpartitioned_nodes () const
  { return this->n_nodes_on_proc (DofObject::invalid_processor_id); }

  /**
   * Returns a number greater than or equal to the maximum node id in the
   * mesh.
   */
  virtual dof_id_type max_node_id () const = 0;

#ifdef LIBMESH_ENABLE_UNIQUE_ID
  /**
   * Returns the next unique id to be used.
   */
  unique_id_type next_unique_id() { return _next_unique_id; }

  /**
   * Sets the next unique id to be used.
   */
  void set_next_unique_id(unique_id_type id) { _next_unique_id = id; }
#endif

  /**
   * Reserves space for a known number of nodes.
   * Note that this method may or may not do anything, depending
   * on the actual \p Mesh implementation.  If you know the number
   * of nodes you will add and call this method before repeatedly
   * calling \p add_point() the implementation will be more efficient.
   */
  virtual void reserve_nodes (const dof_id_type nn) = 0;

  /**
   * Returns the number of elements in the mesh.  The standard
   * n_elem() function may return a cached value on distributed
   * meshes, and so can be called by any processor at any time.
   */
  virtual dof_id_type n_elem () const = 0;

  /**
   * Returns the number of elements in the mesh.  The
   * parallel_n_elem() function returns a newly calculated
   * parallel-synchronized value on distributed meshes, and so must be
   * called in parallel only.
   */
  virtual dof_id_type parallel_n_elem () const = 0;

  /**
   * Returns a number greater than or equal to the maximum element id in the
   * mesh.
   */
  virtual dof_id_type max_elem_id () const = 0;

  /**
   * Returns a number greater than or equal to the maximum unique_id in the
   * mesh.
   */
#ifdef LIBMESH_ENABLE_UNIQUE_ID
  virtual unique_id_type parallel_max_unique_id () const = 0;
#endif

  /**
   * Reserves space for a known number of elements.
   * Note that this method may or may not do anything, depending
   * on the actual \p Mesh implementation.  If you know the number
   * of elements you will add and call this method before repeatedly
   * calling \p add_point() the implementation will be more efficient.
   */
  virtual void reserve_elem (const dof_id_type ne) = 0;

  /**
   * Updates parallel caches so that methods like n_elem()
   * accurately reflect changes on other processors
   */
  virtual void update_parallel_id_counts () = 0;

  /**
   * Returns the number of active elements in the mesh.  Implemented
   * in terms of active_element_iterators.
   */
  virtual dof_id_type n_active_elem () const = 0;

  /**
   * Returns the number of elements on processor \p proc.
   */
  dof_id_type n_elem_on_proc (const processor_id_type proc) const;

  /**
   * Returns the number of elements on the local processor.
   */
  dof_id_type n_local_elem () const
  { return this->n_elem_on_proc (this->processor_id()); }

  /**
   * Returns the number of elements owned by no processor.
   */
  dof_id_type n_unpartitioned_elem () const
  { return this->n_elem_on_proc (DofObject::invalid_processor_id); }

  /**
   * Returns the number of active elements on processor \p proc.
   */
  dof_id_type n_active_elem_on_proc (const processor_id_type proc) const;

  /**
   * Returns the number of active elements on the local processor.
   */
  dof_id_type n_active_local_elem () const
  { return this->n_active_elem_on_proc (this->processor_id()); }

  /**
   * This function returns the number of elements that will be written
   * out in the Tecplot format.  For example, a 9-noded quadrilateral will
   * be broken into 4 linear sub-elements for plotting purposes.  Thus, for
   * a mesh of 2 \p QUAD9 elements  \p n_tecplot_elem() will return 8.
   * Implemented in terms of element_iterators.
   */
  dof_id_type n_sub_elem () const;

  /**
   * Same, but only counts active elements.
   */
  dof_id_type n_active_sub_elem () const;

  /**
   * Return a constant reference (for reading only) to the
   * \f$ i^{th} \f$ point, which should be present in this processor's
   * subset of the mesh data structure.
   */
  virtual const Point & point (const dof_id_type i) const = 0;

  /**
   * Return a constant reference (for reading only) to the
   * \f$ i^{th} \f$ node, which should be present in this processor's
   * subset of the mesh data structure.
   */
  virtual const Node & node_ref (const dof_id_type i) const {
    return *this->node_ptr(i);
  }

  /**
   * Return a reference to the \f$ i^{th} \f$ node, which should be
   * present in this processor's subset of the mesh data structure.
   */
  virtual Node & node_ref (const dof_id_type i) {
    return *this->node_ptr(i);
  }

  /**
   * Return a constant reference (for reading only) to the
   * \f$ i^{th} \f$ node, which should be present in this processor's
   * subset of the mesh data structure.
   *
   * This method is deprecated; use the less confusingly-named
   * node_ref() instead.
   */
  virtual const Node & node (const dof_id_type i) const {
    return *this->node_ptr(i);
  }

  /**
   * Return a reference to the \f$ i^{th} \f$ node, which should be
   * present in this processor's subset of the mesh data structure.
   *
   * This method is deprecated; use the less confusingly-named
   * node_ref() instead.
   */
  virtual Node & node (const dof_id_type i) {
    return *this->node_ptr(i);
  }

  /**
   * Return a pointer to the \f$ i^{th} \f$ node, which should be
   * present in this processor's subset of the mesh data structure.
   */
  virtual const Node * node_ptr (const dof_id_type i) const = 0;

  /**
   * Return a writeable pointer to the \f$ i^{th} \f$ node, which
   * should be present in this processor's subset of the mesh data
   * structure.
   */
  virtual Node * node_ptr (const dof_id_type i) = 0;

  /**
   * Return a pointer to the \f$ i^{th} \f$ node, or NULL if no such
   * node exists in this processor's mesh data structure.
   */
  virtual const Node * query_node_ptr (const dof_id_type i) const = 0;

  /**
   * Return a writeable pointer to the \f$ i^{th} \f$ node, or NULL if
   * no such node exists in this processor's mesh data structure.
   */
  virtual Node * query_node_ptr (const dof_id_type i) = 0;

  /**
   * Return a reference to the \f$ i^{th} \f$ element, which should be
   * present in this processor's subset of the mesh data structure.
   */
  virtual const Elem & elem_ref (const dof_id_type i) const {
    return *this->elem_ptr(i);
  }

  /**
   * Return a writeable reference to the \f$ i^{th} \f$ element, which
   * should be present in this processor's subset of the mesh data
   * structure.
   */
  virtual Elem & elem_ref (const dof_id_type i) {
    return *this->elem_ptr(i);
  }

  /**
   * Return a pointer to the \f$ i^{th} \f$ element, which should be
   * present in this processor's subset of the mesh data structure.
   */
  virtual const Elem * elem_ptr (const dof_id_type i) const = 0;

  /**
   * Return a writeable pointer to the \f$ i^{th} \f$ element, which
   * should be present in this processor's subset of the mesh data
   * structure.
   */
  virtual Elem * elem_ptr (const dof_id_type i) = 0;

  /**
   * Return a pointer to the \f$ i^{th} \f$ element, which should be
   * present in this processor's subset of the mesh data structure.
   *
   * This method is deprecated; use the less confusingly-named
   * elem_ptr() instead.
   */
  virtual const Elem * elem (const dof_id_type i) const {
    return this->elem_ptr(i);
  }

  /**
   * Return a writeable pointer to the \f$ i^{th} \f$ element, which
   * should be present in this processor's subset of the mesh data
   * structure.
   *
   * This method is deprecated; use the less confusingly-named
   * elem_ptr() instead.
   */
  virtual Elem * elem (const dof_id_type i) {
    return this->elem_ptr(i);
  }

  /**
   * Return a pointer to the \f$ i^{th} \f$ element, or NULL if no
   * such element exists in this processor's mesh data structure.
   */
  virtual const Elem * query_elem_ptr (const dof_id_type i) const = 0;

  /**
   * Return a writeable pointer to the \f$ i^{th} \f$ element, or NULL
   * if no such element exists in this processor's mesh data structure.
   */
  virtual Elem * query_elem_ptr (const dof_id_type i) = 0;

  /**
   * Return a pointer to the \f$ i^{th} \f$ element, or NULL if no
   * such element exists in this processor's mesh data structure.
   *
   * This method is deprecated; use the less confusingly-named
   * query_elem_ptr() instead.
   */
  virtual const Elem * query_elem (const dof_id_type i) const {
    return this->query_elem_ptr(i);
  }

  /**
   * Return a writeable pointer to the \f$ i^{th} \f$ element, or NULL
   * if no such element exists in this processor's mesh data structure.
   *
   * This method is deprecated; use the less confusingly-named
   * query_elem_ptr() instead.
   */
  virtual Elem * query_elem (const dof_id_type i) {
    return this->query_elem_ptr(i);
  }

  /**
   * Add a new \p Node at \p Point \p p to the end of the vertex array,
   * with processor_id \p procid.
   * Use DofObject::invalid_processor_id (default) to add a node to all
   * processors, or this->processor_id() to add a node to the local
   * processor only.
   * If adding a node locally, passing an \p id other than
   * DofObject::invalid_id will set that specific node id.  Only
   * do this in parallel if you are manually keeping ids consistent.
   */
  virtual Node * add_point (const Point & p,
                            const dof_id_type id = DofObject::invalid_id,
                            const processor_id_type proc_id =
                            DofObject::invalid_processor_id) = 0;

  /**
   * Add \p Node \p n to the end of the vertex array.
   */
  virtual Node * add_node (Node * n) = 0;

  /**
   * Insert \p Node \p n into the Mesh at a location consistent with
   * n->id(), allocating extra storage if necessary.  Will error
   * rather than overwriting an existing Node.  Primarily intended for
   * use with the mesh_inserter_iterator, only use if you know what
   * you are doing...
   */
  virtual Node * insert_node(Node * n) = 0;

  /**
   * Removes the Node n from the mesh.
   */
  virtual void delete_node (Node * n) = 0;

  /**
   * Changes the id of node \p old_id, both by changing node(old_id)->id() and
   * by moving node(old_id) in the mesh's internal container.  No element with
   * the id \p new_id should already exist.
   */
  virtual void renumber_node (dof_id_type old_id, dof_id_type new_id) = 0;

  /**
   * Add elem \p e to the end of the element array.
   * To add an element locally, set e->processor_id() before adding it.
   * To ensure a specific element id, call e->set_id() before adding it;
   * only do this in parallel if you are manually keeping ids consistent.
   *
   * Users should call MeshBase::prepare_for_use() after elements are
   * added to and/or deleted from the mesh.
   */
  virtual Elem * add_elem (Elem * e) = 0;

  /**
   * Insert elem \p e to the element array, preserving its id
   * and replacing/deleting any existing element with the same id.
   *
   * Users should call MeshBase::prepare_for_use() after elements are
   * added to and/or deleted from the mesh.
   */
  virtual Elem * insert_elem (Elem * e) = 0;

  /**
   * Removes element \p e from the mesh. Note that calling this
   * method may produce isolated nodes, i.e. nodes not connected
   * to any element.  This method must be implemented in derived classes
   * in such a way that it does not invalidate element iterators.
   *
   * Users should call MeshBase::prepare_for_use() after elements are
   * added to and/or deleted from the mesh.
   */
  virtual void delete_elem (Elem * e) = 0;

  /**
   * Changes the id of element \p old_id, both by changing elem(old_id)->id()
   * and by moving elem(old_id) in the mesh's internal container.  No element
   * with the id \p new_id should already exist.
   */
  virtual void renumber_elem (dof_id_type old_id, dof_id_type new_id) = 0;

  /**
   * Locate element face (edge in 2D) neighbors.  This is done with the help
   * of a \p std::map that functions like a hash table.
   * After this routine is called all the elements with a \p NULL neighbor
   * pointer are guaranteed to be on the boundary.  Thus this routine is
   * useful for automatically determining the boundaries of the domain.
   * If reset_remote_elements is left to false, remote neighbor links are not
   * reset and searched for in the local mesh.  If reset_current_list is
   * left as true, then any existing links will be reset before initiating
   * the algorithm, while honoring the value of the reset_remote_elements
   * flag.
   */
  virtual void find_neighbors (const bool reset_remote_elements = false,
                               const bool reset_current_list    = true) = 0;

  /**
   * After partitoning a mesh it is useful to renumber the nodes and elements
   * so that they lie in contiguous blocks on the processors.  This method
   * does just that.
   */
  virtual void renumber_nodes_and_elements () = 0;

  /**
   * There is no reason for a user to ever call this function.
   *
   * This function restores a previously broken element/node numbering such that
   * \p mesh.node_ref(n).id() == n.
   */
  virtual void fix_broken_node_and_element_numbering () = 0;


#ifdef LIBMESH_ENABLE_AMR
  /**
   * Delete subactive (i.e. children of coarsened) elements.
   * This removes all elements descended from currently active
   * elements in the mesh.
   */
  virtual bool contract () = 0;
#endif

  /**
   * Prepare a newly created (or read) mesh for use.
   * This involves 4 steps:
   *  1.) call \p find_neighbors()
   *  2.) call \p partition()
   *  3.) call \p renumber_nodes_and_elements()
   *  4.) call \p cache_elem_dims()
   *
   * The argument to skip renumbering is now deprecated - to prevent a
   * mesh from being renumbered, set allow_renumbering(false).
   */
  void prepare_for_use (const bool skip_renumber_nodes_and_elements=false, const bool skip_find_neighbors=false);

  /**
   * Call the default partitioner (currently \p metis_partition()).
   */
  virtual void partition (const unsigned int n_parts);

  void partition ()
  { this->partition(this->n_processors()); }

  /**
   * Redistribute elements between processors.  This gets called
   * automatically by the Partitioner, and is a no-op in the case of a
   * ReplicatedMesh or serialized DistributedMesh
   */
  virtual void redistribute () {}

  /**
   * Recalculate any cached data after elements and nodes have been
   * repartitioned.
   */
  virtual void update_post_partitioning () {}

  /**
   * If false is passed in then this mesh will no longer be renumbered
   * when being prepared for use.  This may slightly adversely affect
   * performance during subsequent element access, particulary when
   * using a distributed mesh.
   */
  void allow_renumbering(bool allow) { _skip_renumber_nodes_and_elements = !allow; }
  bool allow_renumbering() const { return !_skip_renumber_nodes_and_elements; }

  /**
   * If true is passed in then this mesh will no longer be (re)partitioned.
   * It would probably be a bad idea to call this on a Serial Mesh _before_
   * the first partitioning has happened... because no elements would get assigned
   * to your processor pool.
   *
   * Note that turning on skip_partitioning() can have adverse effects on your
   * performance when using AMR... ie you could get large load imbalances.
   *
   * However you might still want to use this if the communication and computation
   * of the rebalance and repartition is too high for your application.
   */
  void skip_partitioning(bool skip) { _skip_partitioning = skip; }
  bool skip_partitioning() const { return _skip_partitioning; }

  /**
   * Adds a functor which can specify ghosting requirements for use on
   * distributed meshes.  Multiple ghosting functors can be added; any
   * element which is required by any functor will be ghosted.
   *
   * GhostingFunctor memory must be managed by the code which calls
   * this function; the GhostingFunctor lifetime is expected to extend
   * until either the functor is removed or the Mesh is destructed.
   */
  void add_ghosting_functor(GhostingFunctor & ghosting_functor)
  { _ghosting_functors.insert(&ghosting_functor); }

  /**
   * Removes a functor which was previously added to the set of
   * ghosting functors.
   */
  void remove_ghosting_functor(GhostingFunctor & ghosting_functor);

  /**
   * Beginning of range of ghosting functors
   */
  std::set<GhostingFunctor *>::const_iterator ghosting_functors_begin() const
  { return _ghosting_functors.begin(); }

  /**
   * End of range of ghosting functors
   */
  std::set<GhostingFunctor *>::const_iterator ghosting_functors_end() const
  { return _ghosting_functors.end(); }

  /**
   * Constructs a list of all subdomain identifiers in the global mesh.
   * Subdomains correspond to separate subsets of the mesh which could correspond
   * e.g. to different materials in a solid mechanics application,
   * or regions where different physical processes are important.  The subdomain
   * mapping is independent from the parallel decomposition.
   */
  void subdomain_ids (std::set<subdomain_id_type> & ids) const;

  /**
   * Returns the number of subdomains in the global mesh. Subdomains correspond
   * to separate subsets of the mesh which could correspond e.g. to different
   * materials in a solid mechanics application, or regions where different
   * physical processes are important.  The subdomain mapping is independent
   * from the parallel decomposition.
   */
  subdomain_id_type n_subdomains () const;

  /**
   * Returns the number of partitions which have been defined via
   * a call to either mesh.partition() or by building a Partitioner
   * object and calling partition.  Note that the partitioner objects
   * are responsible for setting this value.
   */
  unsigned int n_partitions () const
  { return _n_parts; }

  /**
   * @returns a string containing relevant information
   * about the mesh.
   */
  std::string get_info () const;

  /**
   * Prints relevant information about the mesh.
   */
  void print_info (std::ostream & os=libMesh::out) const;

  /**
   * Equivalent to calling print_info() above, but now you can write:
   * Mesh mesh;
   * libMesh::out << mesh << std::endl;
   */
  friend std::ostream & operator << (std::ostream & os, const MeshBase & m);

  /**
   * Interfaces for reading/writing a mesh to/from a file.  Must be
   * implemented in derived classes.
   */
  virtual void read  (const std::string & name,
                      MeshData * mesh_data=libmesh_nullptr,
                      bool skip_renumber_nodes_and_elements=false,
                      bool skip_find_neighbors=false) = 0;
  virtual void write (const std::string & name, MeshData * mesh_data=libmesh_nullptr) = 0;

  /**
   * Converts a mesh with higher-order
   * elements into a mesh with linear elements.  For
   * example, a mesh consisting of \p Tet10 will be converted
   * to a mesh with \p Tet4 etc.
   */
  virtual void all_first_order () = 0;

  /**
   * Converts a (conforming, non-refined) mesh with linear
   * elements into a mesh with second-order elements.  For
   * example, a mesh consisting of \p Tet4 will be converted
   * to a mesh with \p Tet10 etc.  Note that for some elements
   * like \p Hex8 there exist @e two higher order equivalents,
   * \p Hex20 and \p Hex27.  When \p full_ordered is \p true
   * (default), then \p Hex27 is built.  Otherwise, \p Hex20
   * is built.  The same holds obviously for \p Quad4, \p Prism6
   * ...
   */
  virtual void all_second_order (const bool full_ordered=true) = 0;

  /**
   * We need an empty, generic class to act as a predicate for this
   * and derived mesh classes.
   */
  typedef Predicates::multi_predicate Predicate;

  /**
   * structs for the element_iterator's.
   * Note that these iterators were designed so that derived mesh classes could use the
   * _same_ base class iterators interchangeably.  Their definition comes later in the
   * header file.
   */
  struct element_iterator;
  struct const_element_iterator;

  /**
   * structs for the node_iterator's.
   * Note that these iterators were designed so that derived mesh classes could use the
   * _same_ base class iterators interchangeably.  Their definition comes later in the
   * header file.
   */
  struct node_iterator;
  struct const_node_iterator;

  /**
   * In a few (very rare) cases, the user may have manually tagged the
   * elements with specific processor IDs by hand, without using a
   * partitioner.  In this case, the Mesh will not know that the total
   * number of partitions, _n_parts, has changed, unless you call this
   * function.  This is an O(N active elements) calculation.  The return
   * value is the number of partitions, and _n_parts is also set by
   * this function.
   */
  unsigned int recalculate_n_partitions();

  /**
   * \p returns a pointer to a \p PointLocatorBase object for this
   * mesh, constructing a master PointLocator first if necessary.
   * This should never be used in threaded or non-parallel_only code,
   * and so is deprecated.
   */
  const PointLocatorBase & point_locator () const;

  /**
   * \p returns a pointer to a subordinate \p PointLocatorBase object
   * for this mesh, constructing a master PointLocator first if
   * necessary.  This should not be used in threaded or
   * non-parallel_only code unless the master has already been
   * constructed.
   */
  UniquePtr<PointLocatorBase> sub_point_locator () const;

  /**
   * Releases the current \p PointLocator object.
   */
  void clear_point_locator ();

  /**
   * Verify id and processor_id consistency of our elements and
   * nodes containers.
   * Calls libmesh_assert() on each possible failure.
   * Currently only implemented on DistributedMesh; a serial data
   * structure is much harder to get out of sync.
   */
  virtual void libmesh_assert_valid_parallel_ids() const {}

  /**
   * Returns a writable reference for getting/setting an optional
   * name for a subdomain.
   */
  std::string & subdomain_name(subdomain_id_type id);
  const std::string & subdomain_name(subdomain_id_type id) const;

  /**
   * Returns the id of the named subdomain if it exists,
   * Elem::invalid_subdomain_id otherwise.
   */
  subdomain_id_type get_id_by_name(const std::string & name) const;

  //
  // element_iterator accessors
  //

  /**
   * Iterate over all the elements in the Mesh.
   */
  virtual element_iterator elements_begin () = 0;
  virtual element_iterator elements_end () = 0;
  virtual const_element_iterator elements_begin () const = 0;
  virtual const_element_iterator elements_end () const = 0;

  /**
   * Iterate over elements for which elem->ancestor() returns true.
   */
  virtual element_iterator ancestor_elements_begin () = 0;
  virtual element_iterator ancestor_elements_end () = 0;
  virtual const_element_iterator ancestor_elements_begin () const = 0;
  virtual const_element_iterator ancestor_elements_end () const = 0;

  /**
   * Iterate over elements for which elem->subactive() returns true.
   */
  virtual element_iterator subactive_elements_begin () = 0;
  virtual element_iterator subactive_elements_end () = 0;
  virtual const_element_iterator subactive_elements_begin () const = 0;
  virtual const_element_iterator subactive_elements_end () const = 0;

  /**
   * Iterate over elements for which elem->is_semilocal() returns true for the current processor.
   */
  virtual element_iterator semilocal_elements_begin () = 0;
  virtual element_iterator semilocal_elements_end () = 0;
  virtual const_element_iterator semilocal_elements_begin () const = 0;
  virtual const_element_iterator semilocal_elements_end () const = 0;

  /**
   * Iterate over elements which are on or have a neighbor on the current processor.
   */
  virtual element_iterator facelocal_elements_begin () = 0;
  virtual element_iterator facelocal_elements_end () = 0;
  virtual const_element_iterator facelocal_elements_begin () const = 0;
  virtual const_element_iterator facelocal_elements_end () const = 0;

  /**
   * Iterate over elements of a given level.
   */
  virtual element_iterator level_elements_begin (unsigned int level) = 0;
  virtual element_iterator level_elements_end (unsigned int level) = 0;
  virtual const_element_iterator level_elements_begin (unsigned int level) const = 0;
  virtual const_element_iterator level_elements_end (unsigned int level) const = 0;

  /**
   * Iterate over all elements with a specified processor id.
   */
  virtual element_iterator pid_elements_begin (processor_id_type proc_id) = 0;
  virtual element_iterator pid_elements_end (processor_id_type proc_id) = 0;
  virtual const_element_iterator pid_elements_begin (processor_id_type proc_id) const = 0;
  virtual const_element_iterator pid_elements_end (processor_id_type proc_id) const = 0;

  /**
   * Iterate over all elements with a specified geometric type.
   */
  virtual element_iterator type_elements_begin (ElemType type) = 0;
  virtual element_iterator type_elements_end (ElemType type) = 0;
  virtual const_element_iterator type_elements_begin (ElemType type) const = 0;
  virtual const_element_iterator type_elements_end (ElemType type) const = 0;

  /**
   * Iterate over unpartitioned elements in the Mesh.
   */
  virtual element_iterator unpartitioned_elements_begin () = 0;
  virtual element_iterator unpartitioned_elements_end () = 0;
  virtual const_element_iterator unpartitioned_elements_begin () const = 0;
  virtual const_element_iterator unpartitioned_elements_end () const = 0;

  /**
   * Iterate over active unpartitioned elements in the Mesh.
   */
  virtual element_iterator active_unpartitioned_elements_begin () = 0;
  virtual element_iterator active_unpartitioned_elements_end () = 0;
  virtual const_element_iterator active_unpartitioned_elements_begin () const = 0;
  virtual const_element_iterator active_unpartitioned_elements_end () const = 0;

  /**
   * Iterate over "ghost" elements in the Mesh.  A ghost element is
   * one which is *not* local, but *is* semilocal.
   */
  virtual element_iterator ghost_elements_begin () = 0;
  virtual element_iterator ghost_elements_end () = 0;
  virtual const_element_iterator ghost_elements_begin () const = 0;
  virtual const_element_iterator ghost_elements_end () const = 0;

  /**
   * Iterate over elements in the Mesh where the solution (as
   * distributed by the given DofMap) can be evaluated, for the given
   * variable var_num, or for all variables by default.
   */
  virtual element_iterator evaluable_elements_begin
    (const DofMap & dof_map,
     unsigned int var_num = libMesh::invalid_uint) = 0;
  virtual element_iterator evaluable_elements_end
    (const DofMap & dof_map,
     unsigned int var_num = libMesh::invalid_uint) = 0;
  virtual const_element_iterator evaluable_elements_begin
    (const DofMap & dof_map,
     unsigned int var_num = libMesh::invalid_uint) const = 0;
  virtual const_element_iterator evaluable_elements_end
    (const DofMap & dof_map,
     unsigned int var_num = libMesh::invalid_uint) const = 0;

#ifdef LIBMESH_ENABLE_AMR
  /**
   * Iterate over all elements with a specified refinement flag.
   */
  virtual element_iterator flagged_elements_begin (unsigned char rflag) = 0;
  virtual element_iterator flagged_elements_end (unsigned char rflag) = 0;
  virtual const_element_iterator flagged_elements_begin (unsigned char rflag) const = 0;
  virtual const_element_iterator flagged_elements_end (unsigned char rflag) const = 0;

  /**
   * Iterate over all elements with a specified refinement flag on a
   * specified processor.
   */
  virtual element_iterator flagged_pid_elements_begin (unsigned char rflag,
                                                       processor_id_type pid) = 0;
  virtual element_iterator flagged_pid_elements_end (unsigned char rflag,
                                                     processor_id_type pid) = 0;
  virtual const_element_iterator flagged_pid_elements_begin (unsigned char rflag,
                                                             processor_id_type pid) const = 0;
  virtual const_element_iterator flagged_pid_elements_end (unsigned char rflag,
                                                           processor_id_type pid) const = 0;
#endif

  /**
   * Active, local, and negation forms of the element iterators described above.
   * An "active" element is an element without children (i.e. has not been refined).
   * A "local" element is one whose processor_id() matches the current processor.
   */
  virtual element_iterator active_elements_begin () = 0;
  virtual element_iterator active_elements_end () = 0;
  virtual const_element_iterator active_elements_begin () const = 0;
  virtual const_element_iterator active_elements_end () const = 0;

  virtual element_iterator local_elements_begin () = 0;
  virtual element_iterator local_elements_end () = 0;
  virtual const_element_iterator local_elements_begin () const = 0;
  virtual const_element_iterator local_elements_end () const = 0;

  virtual element_iterator active_semilocal_elements_begin () = 0;
  virtual element_iterator active_semilocal_elements_end () = 0;
  virtual const_element_iterator active_semilocal_elements_begin () const = 0;
  virtual const_element_iterator active_semilocal_elements_end () const = 0;

  virtual element_iterator active_type_elements_begin (ElemType type) = 0;
  virtual element_iterator active_type_elements_end (ElemType type) = 0;
  virtual const_element_iterator active_type_elements_begin (ElemType type) const = 0;
  virtual const_element_iterator active_type_elements_end (ElemType type) const = 0;

  virtual element_iterator active_pid_elements_begin (processor_id_type proc_id) = 0;
  virtual element_iterator active_pid_elements_end (processor_id_type proc_id) = 0;
  virtual const_element_iterator active_pid_elements_begin (processor_id_type proc_id) const = 0;
  virtual const_element_iterator active_pid_elements_end (processor_id_type proc_id) const = 0;

  virtual element_iterator active_subdomain_elements_begin (subdomain_id_type subdomain_id) = 0;
  virtual element_iterator active_subdomain_elements_end (subdomain_id_type subdomain_id) = 0;
  virtual const_element_iterator active_subdomain_elements_begin (subdomain_id_type subdomain_id) const = 0;
  virtual const_element_iterator active_subdomain_elements_end (subdomain_id_type subdomain_id) const = 0;

  virtual element_iterator active_local_subdomain_elements_begin (subdomain_id_type subdomain_id) = 0;
  virtual element_iterator active_local_subdomain_elements_end (subdomain_id_type subdomain_id) = 0;
  virtual const_element_iterator active_local_subdomain_elements_begin (subdomain_id_type subdomain_id) const = 0;
  virtual const_element_iterator active_local_subdomain_elements_end (subdomain_id_type subdomain_id) const = 0;

  virtual element_iterator local_level_elements_begin (unsigned int level) = 0;
  virtual element_iterator local_level_elements_end (unsigned int level) = 0;
  virtual const_element_iterator local_level_elements_begin (unsigned int level) const = 0;
  virtual const_element_iterator local_level_elements_end (unsigned int level) const = 0;

  virtual element_iterator local_not_level_elements_begin (unsigned int level) = 0;
  virtual element_iterator local_not_level_elements_end (unsigned int level) = 0;
  virtual const_element_iterator local_not_level_elements_begin (unsigned int level) const = 0;
  virtual const_element_iterator local_not_level_elements_end (unsigned int level) const = 0;

  virtual element_iterator not_level_elements_begin (unsigned int level) = 0;
  virtual element_iterator not_level_elements_end (unsigned int level) = 0;
  virtual const_element_iterator not_level_elements_begin (unsigned int level) const = 0;
  virtual const_element_iterator not_level_elements_end (unsigned int level) const = 0;

  virtual element_iterator active_local_elements_begin () = 0;
  virtual element_iterator active_local_elements_end () = 0;
  virtual const_element_iterator active_local_elements_begin () const = 0;
  virtual const_element_iterator active_local_elements_end () const = 0;

  virtual element_iterator active_not_local_elements_begin () = 0;
  virtual element_iterator active_not_local_elements_end () = 0;
  virtual const_element_iterator active_not_local_elements_begin () const = 0;
  virtual const_element_iterator active_not_local_elements_end () const = 0;

  virtual element_iterator not_local_elements_begin () = 0;
  virtual element_iterator not_local_elements_end () = 0;
  virtual const_element_iterator not_local_elements_begin () const = 0;
  virtual const_element_iterator not_local_elements_end () const = 0;

  virtual element_iterator not_subactive_elements_begin () = 0;
  virtual element_iterator not_subactive_elements_end () = 0;
  virtual const_element_iterator not_subactive_elements_begin () const = 0;
  virtual const_element_iterator not_subactive_elements_end () const = 0;

  virtual element_iterator not_active_elements_begin () = 0;
  virtual element_iterator not_active_elements_end () = 0;
  virtual const_element_iterator not_active_elements_begin () const = 0;
  virtual const_element_iterator not_active_elements_end () const = 0;

  virtual element_iterator not_ancestor_elements_begin () = 0;
  virtual element_iterator not_ancestor_elements_end () = 0;
  virtual const_element_iterator not_ancestor_elements_begin () const = 0;
  virtual const_element_iterator not_ancestor_elements_end () const = 0;

  //
  // node_iterator accessors
  //

  /**
   * Iterate over all the nodes in the Mesh.
   */
  virtual node_iterator nodes_begin () = 0;
  virtual node_iterator nodes_end () = 0;
  virtual const_node_iterator nodes_begin () const = 0;
  virtual const_node_iterator nodes_end () const = 0;

  /**
   * Iterate over only the active nodes in the Mesh.
   */
  virtual node_iterator active_nodes_begin () = 0;
  virtual node_iterator active_nodes_end () = 0;
  virtual const_node_iterator active_nodes_begin () const = 0;
  virtual const_node_iterator active_nodes_end () const = 0;

  /**
   * Iterate over local nodes (nodes whose processor_id() matches the current processor).
   */
  virtual node_iterator local_nodes_begin () = 0;
  virtual node_iterator local_nodes_end () = 0;
  virtual const_node_iterator local_nodes_begin () const = 0;
  virtual const_node_iterator local_nodes_end () const = 0;

  /**
   * Iterate over nodes with processor_id() == proc_id
   */
  virtual node_iterator pid_nodes_begin (processor_id_type proc_id) = 0;
  virtual node_iterator pid_nodes_end (processor_id_type proc_id) = 0;
  virtual const_node_iterator pid_nodes_begin (processor_id_type proc_id) const = 0;
  virtual const_node_iterator pid_nodes_end (processor_id_type proc_id) const = 0;

  /**
   * Iterate over nodes for which BoundaryInfo::has_boundary_id(node, bndry_id) returns true.
   */
  virtual node_iterator bid_nodes_begin (boundary_id_type bndry_id) = 0;
  virtual node_iterator bid_nodes_end (boundary_id_type bndry_id) = 0;
  virtual const_node_iterator bid_nodes_begin (boundary_id_type bndry_id) const = 0;
  virtual const_node_iterator bid_nodes_end (boundary_id_type bndry_id) const = 0;

  /**
   * Iterate over nodes for which BoundaryInfo::n_boundary_ids(node) > 0.
   */
  virtual node_iterator bnd_nodes_begin () = 0;
  virtual node_iterator bnd_nodes_end () = 0;
  virtual const_node_iterator bnd_nodes_begin () const = 0;
  virtual const_node_iterator bnd_nodes_end () const = 0;

  /**
   * Return a writeable reference to the whole subdomain name map
   */
  std::map<subdomain_id_type, std::string> & set_subdomain_name_map ()
  { return _block_id_to_name; }
  const std::map<subdomain_id_type, std::string> & get_subdomain_name_map () const
  { return _block_id_to_name; }


  /**
   * Search the mesh and cache the different dimenions of the elements
   * present in the mesh.  This is done in prepare_for_use(), but can
   * be done manually by other classes after major mesh modifications.
   */
  void cache_elem_dims();

  /**
   * Search the mesh for elements that have a neighboring element
   * of dim+1 and set that element as the interior parent
   */
  void detect_interior_parents();


  /**
   * This class holds the boundary information.  It can store nodes, edges,
   * and faces with a corresponding id that facilitates setting boundary
   * conditions.
   *
   * Direct access to this class will be removed in future libMesh
   * versions.  Use the \p get_boundary_info() accessor instead.
   */
  UniquePtr<BoundaryInfo> boundary_info;


protected:

  /**
   * Returns a writeable reference to the number of partitions.
   */
  unsigned int & set_n_partitions ()
  { return _n_parts; }

  /**
   * The number of partitions the mesh has.  This is set by
   * the partitioners, and may not be changed directly by
   * the user.
   * **NOTE** The number of partitions *need not* equal
   * this->n_processors(), consider for example the case
   * where you simply want to partition a mesh on one
   * processor and view the result in GMV.
   */
  unsigned int _n_parts;

  /**
   * Flag indicating if the mesh has been prepared for use.
   */
  bool _is_prepared;

  /**
   * A \p PointLocator class for this mesh.
   * This will not actually be built unless needed. Further, since we want
   * our \p point_locator() method to be \p const (yet do the dynamic allocating)
   * this needs to be mutable.  Since the PointLocatorBase::build() member is used,
   * and it operates on a constant reference to the mesh, this is OK.
   */
  mutable UniquePtr<PointLocatorBase> _point_locator;

  /**
   * A partitioner to use at each prepare_for_use().
   *
   * This will be built in the constructor of each derived class, but
   * can be replaced by the user through the partitioner() accessor.
   */
  UniquePtr<Partitioner> _partitioner;

#ifdef LIBMESH_ENABLE_UNIQUE_ID
  /**
   * The next available unique id for assigning ids to DOF objects
   */
  unique_id_type _next_unique_id;
#endif

  /**
   * If this is true then no partitioning should be done.
   */
  bool _skip_partitioning;

  /**
   * If this is true then renumbering will be kept to a miniumum.
   *
   * This is set when prepare_for_use() is called.
   */
  bool _skip_renumber_nodes_and_elements;

  /**
   * This structure maintains the mapping of named blocks
   * for file formats that support named blocks.  Currently
   * this is only implemented for ExodusII
   */
  std::map<subdomain_id_type, std::string> _block_id_to_name;

  /**
   * We cache the dimension of the elements present in the mesh.
   * So, if we have a mesh with 1D and 2D elements, this structure
   * will contain 1 and 2.
   */
  std::set<unsigned char> _elem_dims;

  /**
   * The "spatial dimension" of the Mesh.  See the documentation for
   * Mesh::spatial_dimension() for more information.
   */
  unsigned char _spatial_dimension;

  /**
   * The default geometric GhostingFunctor, used to implement standard
   * libMesh element ghosting behavior.  We use a base class pointer
   * here to avoid dragging in more header dependencies.
   */
  UniquePtr<GhostingFunctor> _default_ghosting;

  /**
   * The list of all GhostingFunctor objects to be used when
   * distributing a DistributedMesh.
   *
   * Basically unused by ReplicatedMesh for now, but belongs to
   * MeshBase because the cost is trivial.
   */
  std::set<GhostingFunctor *> _ghosting_functors;

  /**
   * The partitioner class is a friend so that it can set
   * the number of partitions.
   */
  friend class Partitioner;

  /**
   * The MeshInput classes are friends so that they can set the number
   * of partitions.
   */
  friend class MeshInput<MeshBase>;

  /**
   * Make the \p BoundaryInfo class a friend so that
   * it can create and interact with \p BoundaryMesh.
   */
  friend class BoundaryInfo;

private:
  /**
   *  The default shallow assignment operator is a very bad idea, so
   *  we'll make it a compile-time error to try and do it from other
   *  classes and a link-time error to try and do it from this class.
   *  Use clone() if necessary.
   */
  MeshBase & operator= (const MeshBase & other);
};











/**
 * The definition of the element_iterator struct.
 */
struct
MeshBase::element_iterator : variant_filter_iterator<MeshBase::Predicate, Elem *>
{
  // Templated forwarding ctor -- forwards to appropriate variant_filter_iterator ctor
  template <typename PredType, typename IterType>
  element_iterator (const IterType & d,
                    const IterType & e,
                    const PredType & p ) :
    variant_filter_iterator<MeshBase::Predicate, Elem *>(d,e,p) {}
};




/**
 * The definition of the const_element_iterator struct.  It is similar to the regular
 * iterator above, but also provides an additional conversion-to-const ctor.
 */
struct
MeshBase::const_element_iterator : variant_filter_iterator<MeshBase::Predicate,
                                                           Elem * const,
                                                           Elem * const &,
                                                           Elem * const *>
{
  // Templated forwarding ctor -- forwards to appropriate variant_filter_iterator ctor
  template <typename PredType, typename IterType>
  const_element_iterator (const IterType & d,
                          const IterType & e,
                          const PredType & p ) :
    variant_filter_iterator<MeshBase::Predicate, Elem * const, Elem * const &, Elem * const *>(d,e,p)  {}


  // The conversion-to-const ctor.  Takes a regular iterator and calls the appropriate
  // variant_filter_iterator copy constructor.  Note that this one is *not* templated!
  const_element_iterator (const MeshBase::element_iterator & rhs) :
    variant_filter_iterator<Predicate, Elem * const, Elem * const &, Elem * const *>(rhs)
  {
    // libMesh::out << "Called element_iterator conversion-to-const ctor." << std::endl;
  }
};







/**
 * The definition of the node_iterator struct.
 */
struct
MeshBase::node_iterator : variant_filter_iterator<MeshBase::Predicate, Node *>
{
  // Templated forwarding ctor -- forwards to appropriate variant_filter_iterator ctor
  template <typename PredType, typename IterType>
  node_iterator (const IterType & d,
                 const IterType & e,
                 const PredType & p ) :
    variant_filter_iterator<MeshBase::Predicate, Node *>(d,e,p) {}
};




/**
 * The definition of the const_node_iterator struct.  It is similar to the regular
 * iterator above, but also provides an additional conversion-to-const ctor.
 */
struct
MeshBase::const_node_iterator : variant_filter_iterator<MeshBase::Predicate,
                                                        Node * const,
                                                        Node * const &,
                                                        Node * const *>
{
  // Templated forwarding ctor -- forwards to appropriate variant_filter_iterator ctor
  template <typename PredType, typename IterType>
  const_node_iterator (const IterType & d,
                       const IterType & e,
                       const PredType & p ) :
    variant_filter_iterator<MeshBase::Predicate, Node * const, Node * const &, Node * const *>(d,e,p)  {}


  // The conversion-to-const ctor.  Takes a regular iterator and calls the appropriate
  // variant_filter_iterator copy constructor.  Note that this one is *not* templated!
  const_node_iterator (const MeshBase::node_iterator & rhs) :
    variant_filter_iterator<Predicate, Node * const, Node * const &, Node * const *>(rhs)
  {
    // libMesh::out << "Called node_iterator conversion-to-const ctor." << std::endl;
  }
};


} // namespace libMesh

#endif // LIBMESH_MESH_BASE_H
