// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute and/or
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



// library configuration
#include "libmesh/libmesh_config.h"

// C++ includes
#include <algorithm> // for std::min
#include <map>       // for std::multimap
#include <sstream>   // for std::ostringstream


// Local includes
#include "libmesh/boundary_info.h"
#include "libmesh/elem.h"
#include "libmesh/ghost_point_neighbors.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_communication.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/parallel.h"
#include "libmesh/partitioner.h"
#include "libmesh/point_locator_base.h"
#include "libmesh/threads.h"

#include LIBMESH_INCLUDE_UNORDERED_MAP

namespace libMesh
{



// ------------------------------------------------------------
// MeshBase class member functions
MeshBase::MeshBase (const Parallel::Communicator & comm_in,
                    unsigned char d) :
  ParallelObject (comm_in),
  boundary_info  (new BoundaryInfo(*this)),
  _n_parts       (1),
  _is_prepared   (false),
  _point_locator (),
  _count_lower_dim_elems_in_point_locator(true),
  _partitioner   (),
#ifdef LIBMESH_ENABLE_UNIQUE_ID
  _next_unique_id(DofObject::invalid_unique_id),
#endif
  _skip_partitioning(libMesh::on_command_line("--skip-partitioning")),
  _skip_renumber_nodes_and_elements(false),
  _allow_remote_element_removal(true),
  _spatial_dimension(d),
  _default_ghosting(new GhostPointNeighbors(*this))
{
  _elem_dims.insert(d);
  _ghosting_functors.insert(_default_ghosting.get());
  libmesh_assert_less_equal (LIBMESH_DIM, 3);
  libmesh_assert_greater_equal (LIBMESH_DIM, d);
  libmesh_assert (libMesh::initialized());
}


#ifndef LIBMESH_DISABLE_COMMWORLD
MeshBase::MeshBase (unsigned char d) :
  ParallelObject (CommWorld),
  boundary_info  (new BoundaryInfo(*this)),
  _n_parts       (1),
  _is_prepared   (false),
  _point_locator (),
  _partitioner   (),
#ifdef LIBMESH_ENABLE_UNIQUE_ID
  _next_unique_id(DofObject::invalid_unique_id),
#endif
  _skip_partitioning(libMesh::on_command_line("--skip-partitioning")),
  _skip_renumber_nodes_and_elements(false),
  _allow_remote_element_removal(true),
  _spatial_dimension(d),
  _default_ghosting(new GhostPointNeighbors(*this))
{
  libmesh_deprecated();
  _elem_dims.insert(d);
  _ghosting_functors.insert(_default_ghosting.get());
  libmesh_assert_less_equal (LIBMESH_DIM, 3);
  libmesh_assert_greater_equal (LIBMESH_DIM, d);
  libmesh_assert (libMesh::initialized());
}
#endif // !LIBMESH_DISABLE_COMMWORLD



MeshBase::MeshBase (const MeshBase & other_mesh) :
  ParallelObject (other_mesh),
  boundary_info  (new BoundaryInfo(*this)),
  _n_parts       (other_mesh._n_parts),
  _is_prepared   (other_mesh._is_prepared),
  _point_locator (),
  _partitioner   (),
#ifdef LIBMESH_ENABLE_UNIQUE_ID
  _next_unique_id(other_mesh._next_unique_id),
#endif
  _skip_partitioning(libMesh::on_command_line("--skip-partitioning")),
  _skip_renumber_nodes_and_elements(false),
  _allow_remote_element_removal(true),
  _elem_dims(other_mesh._elem_dims),
  _spatial_dimension(other_mesh._spatial_dimension),
  _default_ghosting(new GhostPointNeighbors(*this)),
  _ghosting_functors(other_mesh._ghosting_functors)
{
  // Make sure we don't accidentally delete the other mesh's default
  // ghosting functor; we'll use our own if that's needed.
  if (other_mesh._ghosting_functors.count
        (other_mesh._default_ghosting.get()))
    {
      _ghosting_functors.erase(other_mesh._default_ghosting.get());
      _ghosting_functors.insert(_default_ghosting.get());
    }

  if(other_mesh._partitioner.get())
    {
      _partitioner = other_mesh._partitioner->clone();
    }
}



MeshBase::~MeshBase()
{
  this->clear();

  libmesh_exceptionless_assert (!libMesh::closed());
}



unsigned int MeshBase::mesh_dimension() const
{
  if (!_elem_dims.empty())
    return cast_int<unsigned int>(*_elem_dims.rbegin());
  return 0;
}



unsigned int MeshBase::spatial_dimension () const
{
  return cast_int<unsigned int>(_spatial_dimension);
}



void MeshBase::set_spatial_dimension(unsigned char d)
{
  // The user can set the _spatial_dimension however they wish,
  // libMesh will only *increase* the spatial dimension, however,
  // never decrease it.
  _spatial_dimension = d;
}



void MeshBase::prepare_for_use (const bool skip_renumber_nodes_and_elements, const bool skip_find_neighbors)
{
  parallel_object_only();

  libmesh_assert(this->comm().verify(this->is_serial()));

  // A distributed mesh may have processors with no elements (or
  // processors with no elements of higher dimension, if we ever
  // support mixed-dimension meshes), but we want consistent
  // mesh_dimension anyways.
  //
  // cache_elem_dims() should get the elem_dimensions() and
  // mesh_dimension() correct later, and we don't need it earlier.


  // Renumber the nodes and elements so that they in contiguous
  // blocks.  By default, _skip_renumber_nodes_and_elements is false.
  //
  // We may currently change that by passing
  // skip_renumber_nodes_and_elements==true to this function, but we
  // should use the allow_renumbering() accessor instead.
  //
  // Instances where you if prepare_for_use() should not renumber the nodes
  // and elements include reading in e.g. an xda/r or gmv file. In
  // this case, the ordering of the nodes may depend on an accompanying
  // solution, and the node ordering cannot be changed.

  if (skip_renumber_nodes_and_elements)
    {
      libmesh_deprecated();
      this->allow_renumbering(false);
    }

  // Mesh modification operations might not leave us with consistent
  // id counts, but our partitioner might need that consistency.
  if(!_skip_renumber_nodes_and_elements)
    this->renumber_nodes_and_elements();
  else
    this->update_parallel_id_counts();

  // Let all the elements find their neighbors
  if(!skip_find_neighbors)
    this->find_neighbors();

  // The user may have set boundary conditions.  We require that the
  // boundary conditions were set consistently.  Because we examine
  // neighbors when evaluating non-raw boundary condition IDs, this
  // assert is only valid when our neighbor links are in place.
#ifdef DEBUG
  MeshTools::libmesh_assert_valid_boundary_ids(*this);
#endif

  // Search the mesh for all the dimensions of the elements
  // and cache them.
  this->cache_elem_dims();

  // Search the mesh for elements that have a neighboring element
  // of dim+1 and set that element as the interior parent
  this->detect_interior_parents();

  // Fix up node unique ids in case mesh generation code didn't take
  // exceptional care to do so.
  //  MeshCommunication().make_node_unique_ids_parallel_consistent(*this);

  // We're going to still require that mesh generation code gets
  // element unique ids consistent.
#if defined(DEBUG) && defined(LIBMESH_ENABLE_UNIQUE_ID)
  MeshTools::libmesh_assert_valid_unique_ids(*this);
#endif

  // Reset our PointLocator.  Any old locator is invalidated any time
  // the elements in the underlying elements in the mesh have changed,
  // so we clear it here.
  this->clear_point_locator();

  // Allow our GhostingFunctor objects to reinit if necessary.
  // Do this before partitioning and redistributing, and before
  // deleting remote elements.
  std::set<GhostingFunctor *>::iterator        gf_it = this->ghosting_functors_begin();
  const std::set<GhostingFunctor *>::iterator gf_end = this->ghosting_functors_end();
  for (; gf_it != gf_end; ++gf_it)
    {
      GhostingFunctor *gf = *gf_it;
      libmesh_assert(gf);
      gf->mesh_reinit();
    }

  // Partition the mesh.
  this->partition();

  // If we're using DistributedMesh, we'll probably want it
  // parallelized.
  if (this->_allow_remote_element_removal)
    this->delete_remote_elements();

  if(!_skip_renumber_nodes_and_elements)
    this->renumber_nodes_and_elements();

  // The mesh is now prepared for use.
  _is_prepared = true;

#if defined(DEBUG) && defined(LIBMESH_ENABLE_UNIQUE_ID)
  MeshTools::libmesh_assert_valid_boundary_ids(*this);
  MeshTools::libmesh_assert_valid_unique_ids(*this);
#endif
}



void MeshBase::clear ()
{
  // Reset the number of partitions
  _n_parts = 1;

  // Reset the _is_prepared flag
  _is_prepared = false;

  // Clear boundary information
  this->get_boundary_info().clear();

  // Clear element dimensions
  _elem_dims.clear();

  // Clear our point locator.
  this->clear_point_locator();
}



void MeshBase::remove_ghosting_functor(GhostingFunctor & ghosting_functor)
{
  // We should only be trying to remove ghosting functors we actually
  // have
  libmesh_assert(_ghosting_functors.count(&ghosting_functor));
  _ghosting_functors.erase(&ghosting_functor);
}



void MeshBase::subdomain_ids (std::set<subdomain_id_type> & ids) const
{
  // This requires an inspection on every processor
  parallel_object_only();

  ids.clear();

  const_element_iterator el  = this->active_local_elements_begin();
  const_element_iterator end = this->active_local_elements_end();

  for (; el!=end; ++el)
    ids.insert((*el)->subdomain_id());

  // Some subdomains may only live on other processors
  this->comm().set_union(ids);
}



subdomain_id_type MeshBase::n_subdomains() const
{
  // This requires an inspection on every processor
  parallel_object_only();

  std::set<subdomain_id_type> ids;

  this->subdomain_ids (ids);

  return cast_int<subdomain_id_type>(ids.size());
}




dof_id_type MeshBase::n_nodes_on_proc (const processor_id_type proc_id) const
{
  // We're either counting a processor's nodes or unpartitioned
  // nodes
  libmesh_assert (proc_id < this->n_processors() ||
                  proc_id == DofObject::invalid_processor_id);

  return static_cast<dof_id_type>(std::distance (this->pid_nodes_begin(proc_id),
                                                 this->pid_nodes_end  (proc_id)));
}



dof_id_type MeshBase::n_elem_on_proc (const processor_id_type proc_id) const
{
  // We're either counting a processor's elements or unpartitioned
  // elements
  libmesh_assert (proc_id < this->n_processors() ||
                  proc_id == DofObject::invalid_processor_id);

  return static_cast<dof_id_type>(std::distance (this->pid_elements_begin(proc_id),
                                                 this->pid_elements_end  (proc_id)));
}



dof_id_type MeshBase::n_active_elem_on_proc (const processor_id_type proc_id) const
{
  libmesh_assert_less (proc_id, this->n_processors());
  return static_cast<dof_id_type>(std::distance (this->active_pid_elements_begin(proc_id),
                                                 this->active_pid_elements_end  (proc_id)));
}



dof_id_type MeshBase::n_sub_elem () const
{
  dof_id_type ne=0;

  const_element_iterator       el  = this->elements_begin();
  const const_element_iterator end = this->elements_end();

  for (; el!=end; ++el)
    ne += (*el)->n_sub_elem();

  return ne;
}



dof_id_type MeshBase::n_active_sub_elem () const
{
  dof_id_type ne=0;

  const_element_iterator       el  = this->active_elements_begin();
  const const_element_iterator end = this->active_elements_end();

  for (; el!=end; ++el)
    ne += (*el)->n_sub_elem();

  return ne;
}



std::string MeshBase::get_info() const
{
  std::ostringstream oss;

  oss << " Mesh Information:"                                  << '\n';

  if (!_elem_dims.empty())
    {
      oss << "  elem_dimensions()={";
      std::copy(_elem_dims.begin(),
                --_elem_dims.end(), // --end() is valid if the set is non-empty
                std::ostream_iterator<unsigned int>(oss, ", "));
      oss << cast_int<unsigned int>(*_elem_dims.rbegin());
      oss << "}\n";
    }

  oss << "  spatial_dimension()=" << this->spatial_dimension() << '\n'
      << "  n_nodes()="           << this->n_nodes()           << '\n'
      << "    n_local_nodes()="   << this->n_local_nodes()     << '\n'
      << "  n_elem()="            << this->n_elem()            << '\n'
      << "    n_local_elem()="    << this->n_local_elem()      << '\n'
#ifdef LIBMESH_ENABLE_AMR
      << "    n_active_elem()="   << this->n_active_elem()     << '\n'
#endif
      << "  n_subdomains()="      << static_cast<std::size_t>(this->n_subdomains()) << '\n'
      << "  n_partitions()="      << static_cast<std::size_t>(this->n_partitions()) << '\n'
      << "  n_processors()="      << static_cast<std::size_t>(this->n_processors()) << '\n'
      << "  n_threads()="         << static_cast<std::size_t>(libMesh::n_threads()) << '\n'
      << "  processor_id()="      << static_cast<std::size_t>(this->processor_id()) << '\n';

  return oss.str();
}


void MeshBase::print_info(std::ostream & os) const
{
  os << this->get_info()
     << std::endl;
}


std::ostream & operator << (std::ostream & os, const MeshBase & m)
{
  m.print_info(os);
  return os;
}


void MeshBase::partition (const unsigned int n_parts)
{
  // If we get here and we have unpartitioned elements, we need that
  // fixed.
  if (this->n_unpartitioned_elem() > 0)
    {
      libmesh_assert (partitioner().get());
      libmesh_assert (this->is_serial());
      partitioner()->partition (*this, n_parts);
    }
  // NULL partitioner means don't repartition
  // Non-serial meshes may not be ready for repartitioning here.
  else if(!skip_partitioning() &&
          partitioner().get())
    {
      partitioner()->partition (*this, n_parts);
    }
  else
    {
      // Adaptive coarsening may have "orphaned" nodes on processors
      // whose elements no longer share them.  We need to check for
      // and possibly fix that.
      Partitioner::set_node_processor_ids(*this);

      // Make sure locally cached partition count
      this->recalculate_n_partitions();

      // Make sure any other locally cached data is correct
      this->update_post_partitioning();
    }
}

unsigned int MeshBase::recalculate_n_partitions()
{
  // This requires an inspection on every processor
  parallel_object_only();

  const_element_iterator el  = this->active_local_elements_begin();
  const_element_iterator end = this->active_local_elements_end();

  unsigned int max_proc_id=0;

  for (; el!=end; ++el)
    max_proc_id = std::max(max_proc_id, static_cast<unsigned int>((*el)->processor_id()));

  // The number of partitions is one more than the max processor ID.
  _n_parts = max_proc_id+1;

  this->comm().max(_n_parts);

  return _n_parts;
}



const PointLocatorBase & MeshBase::point_locator () const
{
  libmesh_deprecated();

  if (_point_locator.get() == libmesh_nullptr)
    {
      // PointLocator construction may not be safe within threads
      libmesh_assert(!Threads::in_threads);

      _point_locator.reset (PointLocatorBase::build(TREE_ELEMENTS, *this).release());
    }

  return *_point_locator;
}


UniquePtr<PointLocatorBase> MeshBase::sub_point_locator () const
{
  // If there's no master point locator, then we need one.
  if (_point_locator.get() == libmesh_nullptr)
    {
      // PointLocator construction may not be safe within threads
      libmesh_assert(!Threads::in_threads);

      // And it may require parallel communication
      parallel_object_only();

      _point_locator.reset (PointLocatorBase::build(TREE_ELEMENTS, *this).release());
    }

  // Otherwise there was a master point locator, and we can grab a
  // sub-locator easily.
  return PointLocatorBase::build(TREE_ELEMENTS, *this, _point_locator.get());
}



void MeshBase::clear_point_locator ()
{
  _point_locator.reset(libmesh_nullptr);
}



void MeshBase::set_count_lower_dim_elems_in_point_locator(bool count_lower_dim_elems)
{
  _count_lower_dim_elems_in_point_locator = count_lower_dim_elems;
}



bool MeshBase::get_count_lower_dim_elems_in_point_locator() const
{
  return _count_lower_dim_elems_in_point_locator;
}



std::string & MeshBase::subdomain_name(subdomain_id_type id)
{
  return _block_id_to_name[id];
}

const std::string & MeshBase::subdomain_name(subdomain_id_type id) const
{
  // An empty string to return when no matching subdomain name is found
  static const std::string empty;

  std::map<subdomain_id_type, std::string>::const_iterator iter = _block_id_to_name.find(id);
  if (iter == _block_id_to_name.end())
    return empty;
  else
    return iter->second;
}




subdomain_id_type MeshBase::get_id_by_name(const std::string & name) const
{
  // Linear search over the map values.
  std::map<subdomain_id_type, std::string>::const_iterator
    iter = _block_id_to_name.begin(),
    end_iter = _block_id_to_name.end();

  for ( ; iter != end_iter; ++iter)
    if (iter->second == name)
      return iter->first;

  // If we made it here without returning, we don't have a subdomain
  // with the requested name, so return Elem::invalid_subdomain_id.
  return Elem::invalid_subdomain_id;
}

void MeshBase::cache_elem_dims()
{
  // This requires an inspection on every processor
  parallel_object_only();

  // Need to clear _elem_dims first in case all elements of a
  // particular dimension have been deleted.
  _elem_dims.clear();

  const_element_iterator el  = this->active_elements_begin();
  const_element_iterator end = this->active_elements_end();

  for (; el!=end; ++el)
    _elem_dims.insert((*el)->dim());

  // Some different dimension elements may only live on other processors
  this->comm().set_union(_elem_dims);

  // If the largest element dimension found is larger than the current
  // _spatial_dimension, increase _spatial_dimension.
  unsigned int max_dim = this->mesh_dimension();
  if (max_dim > _spatial_dimension)
    _spatial_dimension = cast_int<unsigned char>(max_dim);

  // _spatial_dimension may need to increase from 1->2 or 2->3 if the
  // mesh is full of 1D elements but they are not x-aligned, or the
  // mesh is full of 2D elements but they are not in the x-y plane.
  // If the mesh is x-aligned or x-y planar, we will end up checking
  // every node's coordinates and not breaking out of the loop
  // early...
  if (_spatial_dimension < 3)
    {
      const_node_iterator node_it  = this->nodes_begin();
      const_node_iterator node_end = this->nodes_end();
      for (; node_it != node_end; ++node_it)
        {
          Node & node = **node_it;

#if LIBMESH_DIM > 1
          // Note: the exact floating point comparison is intentional,
          // we don't want to get tripped up by tolerances.
          if (node(1) != 0.)
            {
              _spatial_dimension = 2;
#if LIBMESH_DIM == 2
              // If libmesh is compiled in 2D mode, this is the
              // largest spatial dimension possible so we can break
              // out.
              break;
#endif
            }
#endif

#if LIBMESH_DIM > 2
          if (node(2) != 0.)
            {
              // Spatial dimension can't get any higher than this, so
              // we can break out.
              _spatial_dimension = 3;
              break;
            }
#endif
        }
    }
}

void MeshBase::detect_interior_parents()
{
  // This requires an inspection on every processor
  parallel_object_only();

  // Check if the mesh contains mixed dimensions. If so, then set interior parents, otherwise return.
  if (this->elem_dimensions().size() == 1)
    return;

  //This map will be used to set interior parents
  LIBMESH_BEST_UNORDERED_MAP<dof_id_type, std::vector<dof_id_type> > node_to_elem;

  const_element_iterator el  = this->active_elements_begin();
  const_element_iterator end = this->active_elements_end();

  for (; el!=end; ++el)
    {
      const Elem * elem = *el;

      // Populating the node_to_elem map, same as MeshTools::build_nodes_to_elem_map
      for (unsigned int n=0; n<elem->n_vertices(); n++)
        {
          libmesh_assert_less (elem->id(), this->max_elem_id());

          node_to_elem[elem->node_id(n)].push_back(elem->id());
        }
    }

  // Automatically set interior parents
  el = this->elements_begin();
  for (; el!=end; ++el)
    {
      Elem * element = *el;

      // Ignore an 3D element or an element that already has an interior parent
      if (element->dim()>=LIBMESH_DIM || element->interior_parent())
        continue;

      // Start by generating a SET of elements that are dim+1 to the current
      // element at each vertex of the current element, thus ignoring interior nodes.
      // If one of the SET of elements is empty, then we will not have an interior parent
      // since an interior parent must be connected to all vertices of the current element
      std::vector< std::set<dof_id_type> > neighbors( element->n_vertices() );

      bool found_interior_parents = false;

      for (dof_id_type n=0; n < element->n_vertices(); n++)
        {
          std::vector<dof_id_type> & element_ids = node_to_elem[element->node_id(n)];
          for (std::vector<dof_id_type>::iterator e_it = element_ids.begin();
               e_it != element_ids.end(); e_it++)
            {
              dof_id_type eid = *e_it;
              if (this->elem_ref(eid).dim() == element->dim()+1)
                neighbors[n].insert(eid);
            }
          if (neighbors[n].size()>0)
            {
              found_interior_parents = true;
            }
          else
            {
              // We have found an empty set, no reason to continue
              // Ensure we set this flag to false before the break since it could have
              // been set to true for previous vertex
              found_interior_parents = false;
              break;
            }
        }

      // If we have successfully generated a set of elements for each vertex, we will compare
      // the set for vertex 0 will the sets for the vertices until we find a id that exists in
      // all sets.  If found, this is our an interior parent id.  The interior parent id found
      // will be the lowest element id if there is potential for multiple interior parents.
      if (found_interior_parents)
        {
          std::set<dof_id_type> & neighbors_0 = neighbors[0];
          for (std::set<dof_id_type>::iterator e_it = neighbors_0.begin();
               e_it != neighbors_0.end(); e_it++)
            {
              found_interior_parents=false;
              dof_id_type interior_parent_id = *e_it;
              for (dof_id_type n=1; n < element->n_vertices(); n++)
                {
                  if (neighbors[n].find(interior_parent_id)!=neighbors[n].end())
                    {
                      found_interior_parents=true;
                    }
                  else
                    {
                      found_interior_parents=false;
                      break;
                    }
                }
              if (found_interior_parents)
                {
                  element->set_interior_parent(this->elem_ptr(interior_parent_id));
                  break;
                }
            }
        }
    }
}

} // namespace libMesh
