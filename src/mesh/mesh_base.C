// The libMesh Finite Element Library.
// Copyright (C) 2002-2021 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include <unordered_map>

// Local includes
#include "libmesh/boundary_info.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/elem.h"
#include "libmesh/ghost_point_neighbors.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_communication.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/parallel.h"
#include "libmesh/parallel_algebra.h"
#include "libmesh/partitioner.h"
#include "libmesh/point_locator_base.h"
#include "libmesh/threads.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/enum_point_locator_type.h"
#include "libmesh/enum_to_string.h"
#include "libmesh/auto_ptr.h" // libmesh_make_unique
#include "libmesh/point_locator_nanoflann.h"

namespace libMesh
{



// ------------------------------------------------------------
// MeshBase class member functions
MeshBase::MeshBase (const Parallel::Communicator & comm_in,
                    unsigned char d) :
  ParallelObject (comm_in),
  boundary_info  (new BoundaryInfo(*this)),
  _n_parts       (1),
  _default_mapping_type(LAGRANGE_MAP),
  _default_mapping_data(0),
  _is_prepared   (false),
  _point_locator (),
  _count_lower_dim_elems_in_point_locator(true),
  _partitioner   (),
#ifdef LIBMESH_ENABLE_UNIQUE_ID
  _next_unique_id(DofObject::invalid_unique_id),
#endif
  _skip_noncritical_partitioning(false),
  _skip_all_partitioning(libMesh::on_command_line("--skip-partitioning")),
  _skip_renumber_nodes_and_elements(false),
  _skip_find_neighbors(false),
  _allow_remote_element_removal(true),
  _spatial_dimension(d),
  _default_ghosting(libmesh_make_unique<GhostPointNeighbors>(*this)),
  _point_locator_close_to_point_tol(0.)
{
  _elem_dims.insert(d);
  _ghosting_functors.insert(_default_ghosting.get());
  libmesh_assert_less_equal (LIBMESH_DIM, 3);
  libmesh_assert_greater_equal (LIBMESH_DIM, d);
  libmesh_assert (libMesh::initialized());
}



MeshBase::MeshBase (const MeshBase & other_mesh) :
  ParallelObject (other_mesh),
  boundary_info  (new BoundaryInfo(*this)),
  _n_parts       (other_mesh._n_parts),
  _default_mapping_type(other_mesh._default_mapping_type),
  _default_mapping_data(other_mesh._default_mapping_data),
  _is_prepared   (other_mesh._is_prepared),
  _point_locator (),
  _count_lower_dim_elems_in_point_locator(other_mesh._count_lower_dim_elems_in_point_locator),
  _partitioner   (),
#ifdef LIBMESH_ENABLE_UNIQUE_ID
  _next_unique_id(other_mesh._next_unique_id),
#endif
  _skip_noncritical_partitioning(false),
  _skip_all_partitioning(libMesh::on_command_line("--skip-partitioning")),
  _skip_renumber_nodes_and_elements(other_mesh._skip_renumber_nodes_and_elements),
  _skip_find_neighbors(other_mesh._skip_find_neighbors),
  _allow_remote_element_removal(other_mesh._allow_remote_element_removal),
  _elem_dims(other_mesh._elem_dims),
  _spatial_dimension(other_mesh._spatial_dimension),
  _default_ghosting(libmesh_make_unique<GhostPointNeighbors>(*this)),
  _point_locator_close_to_point_tol(other_mesh._point_locator_close_to_point_tol)
{
  const GhostingFunctor * const other_default_ghosting = other_mesh._default_ghosting.get();

  for (GhostingFunctor * const gf : other_mesh._ghosting_functors)
    {
      // If the other mesh is using default ghosting, then we will use our own
      // default ghosting
      if (gf == other_default_ghosting)
        {
          _ghosting_functors.insert(_default_ghosting.get());
          continue;
        }

      std::shared_ptr<GhostingFunctor> clone_gf = gf->clone();
      // Some subclasses of GhostingFunctor might not override the
      // clone function yet. If this is the case, GhostingFunctor will
      // return nullptr by default. The clone function should be overridden
      // in all derived classes. This following code ("else") is written
      // for API upgrade. That will allow users gradually to update their code.
      // Once the API upgrade is done, we will come back and delete "else."
      if (clone_gf)
        {
          clone_gf->set_mesh(this);
          add_ghosting_functor(clone_gf);
        }
      else
        {
          libmesh_deprecated();
          add_ghosting_functor(*gf);
        }
    }

  if (other_mesh._partitioner.get())
    _partitioner = other_mesh._partitioner->clone();
}

MeshBase& MeshBase::operator= (MeshBase && other_mesh)
{
  // Move assign as a ParallelObject.
  this->ParallelObject::operator=(other_mesh);

  _n_parts = other_mesh.n_partitions();
  _default_mapping_type = other_mesh.default_mapping_type();
  _default_mapping_data = other_mesh.default_mapping_data();
  _is_prepared = other_mesh.is_prepared();
  _point_locator = std::move(other_mesh._point_locator);
  _count_lower_dim_elems_in_point_locator = other_mesh.get_count_lower_dim_elems_in_point_locator();
  #ifdef LIBMESH_ENABLE_UNIQUE_ID
    _next_unique_id = other_mesh.next_unique_id();
  #endif
  _skip_noncritical_partitioning = other_mesh.skip_noncritical_partitioning();
  _skip_all_partitioning = other_mesh.skip_partitioning();
  _skip_renumber_nodes_and_elements = !(other_mesh.allow_renumbering());
  _skip_find_neighbors = !(other_mesh.allow_find_neighbors());
  _allow_remote_element_removal = other_mesh.allow_remote_element_removal();
  _block_id_to_name = std::move(other_mesh._block_id_to_name);
  _elem_dims = std::move(other_mesh.elem_dimensions());
  _spatial_dimension = other_mesh.spatial_dimension();
  _elem_integer_names = std::move(other_mesh._elem_integer_names);
  _elem_integer_default_values = std::move(other_mesh._elem_integer_default_values);
  _node_integer_names = std::move(other_mesh._node_integer_names);
  _node_integer_default_values = std::move(other_mesh._node_integer_default_values);
  _point_locator_close_to_point_tol = other_mesh.get_point_locator_close_to_point_tol();

  // This relies on our subclasses *not* invalidating pointers when we
  // do their portion of the move assignment later!
  boundary_info = std::move(other_mesh.boundary_info);
  boundary_info->set_mesh(*this);

  // We're *not* really done at this point, but we have the problem
  // that some of our data movement might be expecting subclasses data
  // movement to happen first.  We'll let subclasses handle that by
  // calling our post_dofobject_moves()
  return *this;
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



void MeshBase::set_elem_dimensions(const std::set<unsigned char> & elem_dims)
{
#ifdef DEBUG
  // In debug mode, we call cache_elem_data() and then make sure
  // the result actually agrees with what the user specified.
  parallel_object_only();

  this->cache_elem_data();
  libmesh_assert_msg(_elem_dims == elem_dims, \
                     "Specified element dimensions does not match true element dimensions!");
#endif

  _elem_dims = elem_dims;
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



unsigned int MeshBase::add_elem_integer(const std::string & name,
                                        bool allocate_data,
                                        dof_id_type default_value)
{
  for (auto i : index_range(_elem_integer_names))
    if (_elem_integer_names[i] == name)
      {
        libmesh_assert_less(i, _elem_integer_default_values.size());
        _elem_integer_default_values[i] = default_value;
        return i;
      }

  libmesh_assert_equal_to(_elem_integer_names.size(),
                          _elem_integer_default_values.size());
  _elem_integer_names.push_back(name);
  _elem_integer_default_values.push_back(default_value);
  if (allocate_data)
    this->size_elem_extra_integers();
  return _elem_integer_names.size()-1;
}



std::vector<unsigned int> MeshBase::add_elem_integers(const std::vector<std::string> & names,
                                                      bool allocate_data,
                                                      const std::vector<dof_id_type> * default_values)
{
  libmesh_assert(!default_values || default_values->size() == names.size());
  libmesh_assert_equal_to(_elem_integer_names.size(), _elem_integer_default_values.size());

  std::unordered_map<std::string, std::size_t> name_indices;
  for (auto i : index_range(_elem_integer_names))
    name_indices[_elem_integer_names[i]] = i;

  std::vector<unsigned int> returnval(names.size());

  bool added_an_integer = false;
  for (auto i : index_range(names))
    {
      const std::string & name = names[i];
      auto it = name_indices.find(name);
      if (it != name_indices.end())
        {
          returnval[i] = it->second;
          _elem_integer_default_values[it->second] =
            default_values ? (*default_values)[i] : DofObject::invalid_id;
        }
      else
        {
          returnval[i] = _elem_integer_names.size();
          name_indices[name] = returnval[i];
          _elem_integer_names.push_back(name);
          _elem_integer_default_values.push_back
            (default_values ? (*default_values)[i] : DofObject::invalid_id);
          added_an_integer = true;
        }
    }

  if (allocate_data && added_an_integer)
    this->size_elem_extra_integers();

  return returnval;
}



unsigned int MeshBase::get_elem_integer_index(const std::string & name) const
{
  for (auto i : index_range(_elem_integer_names))
    if (_elem_integer_names[i] == name)
      return i;

  libmesh_error_msg("Unknown elem integer " << name);
  return libMesh::invalid_uint;
}



bool MeshBase::has_elem_integer(const std::string & name) const
{
  for (auto & entry : _elem_integer_names)
    if (entry == name)
      return true;

  return false;
}



unsigned int MeshBase::add_node_integer(const std::string & name,
                                        bool allocate_data,
                                        dof_id_type default_value)
{
  for (auto i : index_range(_node_integer_names))
    if (_node_integer_names[i] == name)
      {
        libmesh_assert_less(i, _node_integer_default_values.size());
        _node_integer_default_values[i] = default_value;
        return i;
      }

  libmesh_assert_equal_to(_node_integer_names.size(),
                          _node_integer_default_values.size());
  _node_integer_names.push_back(name);
  _node_integer_default_values.push_back(default_value);
  if (allocate_data)
    this->size_node_extra_integers();
  return _node_integer_names.size()-1;
}



std::vector<unsigned int> MeshBase::add_node_integers(const std::vector<std::string> & names,
                                                      bool allocate_data,
                                                      const std::vector<dof_id_type> * default_values)
{
  libmesh_assert(!default_values || default_values->size() == names.size());
  libmesh_assert_equal_to(_node_integer_names.size(), _node_integer_default_values.size());

  std::unordered_map<std::string, std::size_t> name_indices;
  for (auto i : index_range(_node_integer_names))
    name_indices[_node_integer_names[i]] = i;

  std::vector<unsigned int> returnval(names.size());

  bool added_an_integer = false;
  for (auto i : index_range(names))
    {
      const std::string & name = names[i];
      auto it = name_indices.find(name);
      if (it != name_indices.end())
        {
          returnval[i] = it->second;
          _node_integer_default_values[it->second] =
            default_values ? (*default_values)[i] : DofObject::invalid_id;
        }
      else
        {
          returnval[i] = _node_integer_names.size();
          name_indices[name] = returnval[i];
          _node_integer_names.push_back(name);
          _node_integer_default_values.push_back
            (default_values ? (*default_values)[i] : DofObject::invalid_id);
          added_an_integer = true;
        }
    }

  if (allocate_data && added_an_integer)
    this->size_node_extra_integers();

  return returnval;
}



unsigned int MeshBase::get_node_integer_index(const std::string & name) const
{
  for (auto i : index_range(_node_integer_names))
    if (_node_integer_names[i] == name)
      return i;

  libmesh_error_msg("Unknown node integer " << name);
  return libMesh::invalid_uint;
}



bool MeshBase::has_node_integer(const std::string & name) const
{
  for (auto & entry : _node_integer_names)
    if (entry == name)
      return true;

  return false;
}



void MeshBase::remove_orphaned_nodes ()
{
  LOG_SCOPE("remove_orphaned_nodes()", "MeshBase");

  // Will hold the set of nodes that are currently connected to elements
  std::unordered_set<Node *> connected_nodes;

  // Loop over the elements.  Find which nodes are connected to at
  // least one of them.
  for (const auto & element : this->element_ptr_range())
    for (auto & n : element->node_ref_range())
      connected_nodes.insert(&n);

  for (const auto & node : this->node_ptr_range())
    if (!connected_nodes.count(node))
      this->delete_node(node);
}



void MeshBase::prepare_for_use (const bool skip_renumber_nodes_and_elements, const bool skip_find_neighbors)
{
  libmesh_deprecated();

  // We only respect the users wish if they tell us to skip renumbering. If they tell us not to
  // skip renumbering but someone previously called allow_renumbering(false), then the latter takes
  // precedence
  if (skip_renumber_nodes_and_elements)
    this->allow_renumbering(false);

  // We always accept the user's value for skip_find_neighbors, in contrast to skip_renumber
  const bool old_allow_find_neighbors = this->allow_find_neighbors();
  this->allow_find_neighbors(!skip_find_neighbors);

  this->prepare_for_use();

  this->allow_find_neighbors(old_allow_find_neighbors);
}

void MeshBase::prepare_for_use (const bool skip_renumber_nodes_and_elements)
{
  libmesh_deprecated();

  // We only respect the users wish if they tell us to skip renumbering. If they tell us not to
  // skip renumbering but someone previously called allow_renumbering(false), then the latter takes
  // precedence
  if (skip_renumber_nodes_and_elements)
    this->allow_renumbering(false);

  this->prepare_for_use();
}

void MeshBase::prepare_for_use ()
{
  LOG_SCOPE("prepare_for_use()", "MeshBase");

  parallel_object_only();

  libmesh_assert(this->comm().verify(this->is_serial()));

  // A distributed mesh may have processors with no elements (or
  // processors with no elements of higher dimension, if we ever
  // support mixed-dimension meshes), but we want consistent
  // mesh_dimension anyways.
  //
  // cache_elem_data() should get the elem_dimensions() and
  // mesh_dimension() correct later, and we don't need it earlier.


  // Renumber the nodes and elements so that they in contiguous
  // blocks.  By default, _skip_renumber_nodes_and_elements is false.
  //
  // Instances where you if prepare_for_use() should not renumber the nodes
  // and elements include reading in e.g. an xda/r or gmv file. In
  // this case, the ordering of the nodes may depend on an accompanying
  // solution, and the node ordering cannot be changed.


  // Mesh modification operations might not leave us with consistent
  // id counts, or might leave us with orphaned nodes we're no longer
  // using, but our partitioner might need that consistency and/or
  // might be confused by orphaned nodes.
  if (!_skip_renumber_nodes_and_elements)
    this->renumber_nodes_and_elements();
  else
    {
      this->remove_orphaned_nodes();
      this->update_parallel_id_counts();
    }

  // Let all the elements find their neighbors
  if (!_skip_find_neighbors)
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
  this->cache_elem_data();

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
  for (auto & gf : _ghosting_functors)
    {
      libmesh_assert(gf);
      gf->mesh_reinit();
    }

  // Partition the mesh unless *all* partitioning is to be skipped.
  // If only noncritical partitioning is to be skipped, the
  // partition() call will still check for orphaned nodes.
  if (!skip_partitioning())
    this->partition();

  // If we're using DistributedMesh, we'll probably want it
  // parallelized.
  if (this->_allow_remote_element_removal)
    this->delete_remote_elements();

  // Much of our boundary info may have been for now-remote parts of the mesh,
  // in which case we don't want to keep local copies of data meant to be
  // local. On the other hand we may have deleted, or the user may have added in
  // a distributed fashion, boundary data that is meant to be global. So we
  // handle both of those scenarios here
  this->get_boundary_info().regenerate_id_sets();

  if (!_skip_renumber_nodes_and_elements)
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
  if (boundary_info)
    boundary_info->clear();

  // Clear element dimensions
  _elem_dims.clear();

  _constraint_rows.clear();

  // Clear our point locator.
  this->clear_point_locator();
}



void MeshBase::remove_ghosting_functor(GhostingFunctor & ghosting_functor)
{
  _ghosting_functors.erase(&ghosting_functor);

  auto it = _shared_functors.find(&ghosting_functor);
  if (it != _shared_functors.end())
    _shared_functors.erase(it);
}



void MeshBase::subdomain_ids (std::set<subdomain_id_type> & ids, const bool global /* = true */) const
{
  // This requires an inspection on every processor
  if (global)
    parallel_object_only();

  ids.clear();

  for (const auto & elem : this->active_local_element_ptr_range())
    ids.insert(elem->subdomain_id());

  if (global)
    {
      // Only include the unpartitioned elements if the user requests the global IDs.
      // In the case of the local subdomain IDs, it doesn't make sense to include the
      // unpartitioned elements because said elements do not have a sense of locality.
      for (const auto & elem : as_range(this->active_unpartitioned_elements_begin(),
                                        this->active_unpartitioned_elements_end()))
        ids.insert(elem->subdomain_id());

      // Some subdomains may only live on other processors
      this->comm().set_union(ids);
    }
}



subdomain_id_type MeshBase::n_subdomains() const
{
  // This requires an inspection on every processor
  parallel_object_only();

  std::set<subdomain_id_type> ids;

  this->subdomain_ids (ids);

  return cast_int<subdomain_id_type>(ids.size());
}



subdomain_id_type MeshBase::n_local_subdomains() const
{
  std::set<subdomain_id_type> ids;

  this->subdomain_ids (ids, /* global = */ false);

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

  for (const auto & elem : this->element_ptr_range())
    ne += elem->n_sub_elem();

  return ne;
}



dof_id_type MeshBase::n_active_sub_elem () const
{
  dof_id_type ne=0;

  for (const auto & elem : this->active_element_ptr_range())
    ne += elem->n_sub_elem();

  return ne;
}



std::string MeshBase::get_info(const unsigned int verbosity /* = 0 */, const bool global /* = true */) const
{
  std::ostringstream oss;

  oss << " Mesh Information:" << '\n';

  if (!_elem_dims.empty())
    {
      oss << "  elem_dimensions()={";
      std::copy(_elem_dims.begin(),
                --_elem_dims.end(), // --end() is valid if the set is non-empty
                std::ostream_iterator<unsigned int>(oss, ", "));
      oss << cast_int<unsigned int>(*_elem_dims.rbegin());
      oss << "}\n";
    }

  oss << "  spatial_dimension()="     << this->spatial_dimension()                            << '\n'
      << "  n_nodes()="               << this->n_nodes()                                      << '\n'
      << "    n_local_nodes()="       << this->n_local_nodes()                                << '\n'
      << "  n_elem()="                << this->n_elem()                                       << '\n'
      << "    n_local_elem()="        << this->n_local_elem()                                 << '\n';
#ifdef LIBMESH_ENABLE_AMR
  oss << "    n_active_elem()="       << this->n_active_elem()                                << '\n';
#endif
  if (global)
    oss << "  n_subdomains()="        << static_cast<std::size_t>(this->n_subdomains())       << '\n';
  else
    oss << "  n_local_subdomains()= " << static_cast<std::size_t>(this->n_local_subdomains()) << '\n';
  oss << "  n_partitions()="          << static_cast<std::size_t>(this->n_partitions())       << '\n'
      << "  n_processors()="          << static_cast<std::size_t>(this->n_processors())       << '\n'
      << "  n_threads()="             << static_cast<std::size_t>(libMesh::n_threads())       << '\n'
      << "  processor_id()="          << static_cast<std::size_t>(this->processor_id())       << '\n'
      << "  is_prepared()="           << (this->is_prepared() ? "true" : "false")             << '\n'
      << "  is_replicated()="         << (this->is_replicated() ? "true" : "false")           << '\n';

  if (verbosity > 0)
    {
      if (global)
        {
          libmesh_parallel_only(this->comm());
          if (this->processor_id() != 0)
            oss << "\n Detailed global get_info() (verbosity > 0) is reduced and output to only rank 0.";
        }

      // Helper for printing element types
      const auto elem_type_helper = [](const std::set<int> & elem_types) {
        std::stringstream ss;
        for (auto it = elem_types.begin(); it != elem_types.end();)
          {
            ss << Utility::enum_to_string((ElemType)*it);
            if (++it != elem_types.end())
              ss << ", ";
          }
        return ss.str();
      };

      // Helper for whether or not the given DofObject is to be included. If we're doing
      // a global reduction, we also count unpartitioned objects on rank 0.
      const auto include_object = [this, &global](const DofObject & dof_object) {
        return this->processor_id() == dof_object.processor_id() ||
               (global &&
                this->processor_id() == 0 &&
                dof_object.processor_id() == DofObject::invalid_processor_id);
      };

      Real volume = 0;

      // Add bounding box information
      const auto bbox = global ? MeshTools::create_bounding_box(*this) : MeshTools::create_local_bounding_box(*this);
      if (!global || this->processor_id() == 0)
        oss << "\n " << (global ? "" : "Local ") << "Mesh Bounding Box:\n"
            << "  Minimum: " << bbox.min()                << "\n"
            << "  Maximum: " << bbox.max()                << "\n"
            << "  Delta:   " << (bbox.max() - bbox.min()) << "\n";

      // Obtain the global or local element types
      std::set<int> elem_types;
      for (const Elem * elem : this->active_local_element_ptr_range())
        elem_types.insert(elem->type());
      if (global)
        {
          // Pick up unpartitioned elems on rank 0
          if (this->processor_id() == 0)
            for (const Elem * elem : as_range(this->active_unpartitioned_elements_begin(),
                                              this->active_unpartitioned_elements_end()))
              elem_types.insert(elem->type());

          this->comm().set_union(elem_types);
        }

      // Add element types
      if (!global || this->processor_id() == 0)
        oss << "\n " << (global ? "" : "Local ") << "Mesh Element Type(s):\n  "
            << elem_type_helper(elem_types) << "\n";

      // Reduce the nodeset ids
      auto nodeset_ids = this->get_boundary_info().get_node_boundary_ids();
      if (global)
        this->comm().set_union(nodeset_ids);

      // Accumulate local information for each nodeset
      struct NodesetInfo
      {
        std::size_t num_nodes = 0;
        BoundingBox bbox;
      };
      std::map<boundary_id_type, NodesetInfo> nodeset_info_map;
      for (const auto & pair : this->get_boundary_info().get_nodeset_map())
        {
          const Node * node = pair.first;
          if (!include_object(*node))
            continue;

          const auto id = pair.second;
          NodesetInfo & info = nodeset_info_map[id];

          ++info.num_nodes;

          if (verbosity > 1)
            info.bbox.union_with(*node);
        }

      // Add nodeset info
      if (!global || this->processor_id() == 0)
        {
          oss << "\n " << (global ? "" : "Local ") << "Mesh Nodesets:\n";
          if (nodeset_ids.empty())
            oss << "  None\n";
        }

      const auto & nodeset_name_map = this->get_boundary_info().get_sideset_name_map();
      for (const auto id : nodeset_ids)
        {
          NodesetInfo & info = nodeset_info_map[id];

          // Reduce the local information for this nodeset if required
          if (global)
            {
              this->comm().sum(info.num_nodes);
              if (verbosity > 1)
                {
                  this->comm().min(info.bbox.min());
                  this->comm().min(info.bbox.max());
                }
            }

          const bool has_name = nodeset_name_map.count(id) && nodeset_name_map.at(id).size();
          const std::string name = has_name ? nodeset_name_map.at(id) : "";
          if (global)
            libmesh_assert(this->comm().verify(name));

          if (global ? this->processor_id() == 0 : info.num_nodes > 0)
            {
              oss << "  Nodeset " << id;
              if (has_name)
                oss << " (" << name << ")";
              oss << ", " << info.num_nodes << " " << (global ? "" : "local ") << "nodes\n";

              if (verbosity > 1)
              {
                oss << "   " << (global ? "Bounding" : "Local bounding") << " box minimum: "
                    << info.bbox.min() << "\n"
                    << "   " << (global ? "Bounding" : "Local bounding") << " box maximum: "
                    << info.bbox.max() << "\n"
                    << "   " << (global ? "Bounding" : "Local bounding") << " box delta: "
                    << (info.bbox.max() - info.bbox.min()) << "\n";
              }
            }
        }

      // Reduce the sideset ids
      auto sideset_ids = this->get_boundary_info().get_side_boundary_ids();
      if (global)
        this->comm().set_union(sideset_ids);

      // Accumulate local information for each sideset
      struct SidesetInfo
      {
        std::size_t num_sides = 0;
        Real volume = 0;
        std::set<int> side_elem_types;
        std::set<int> elem_types;
        std::set<dof_id_type> elem_ids;
        std::set<dof_id_type> node_ids;
        BoundingBox bbox;
      };
      std::map<boundary_id_type, SidesetInfo> sideset_info_map;
      for (const auto & pair : this->get_boundary_info().get_sideset_map())
        {
          const Elem * elem = pair.first;
          if (!include_object(*elem))
            continue;

          const auto id = pair.second.second;
          SidesetInfo & info = sideset_info_map[id];

          const auto s = pair.second.first;
          const auto side = elem->build_side_ptr(s);

          ++info.num_sides;
          info.side_elem_types.insert(side->type());
          info.elem_types.insert(elem->type());
          info.elem_ids.insert(elem->id());

          for (const Node & node : side->node_ref_range())
            if (include_object(node))
              info.node_ids.insert(node.id());

          if (verbosity > 1)
          {
            info.volume += side->volume();
            info.bbox.union_with(side->loose_bounding_box());
          }
        }

      // Add sideset info
      if (!global || this->processor_id() == 0)
        {
          oss << "\n " << (global ? "" : "Local ") << "Mesh Sidesets:\n";
          if (sideset_ids.empty())
            oss << "  None\n";
        }
      const auto & sideset_name_map = this->get_boundary_info().get_sideset_name_map();
      for (const auto id : sideset_ids)
        {
          SidesetInfo & info = sideset_info_map[id];

          auto num_elems = info.elem_ids.size();
          auto num_nodes = info.node_ids.size();

          // Reduce the local information for this sideset if required
          if (global)
            {
              this->comm().sum(info.num_sides);
              this->comm().set_union(info.side_elem_types, 0);
              this->comm().sum(num_elems);
              this->comm().set_union(info.elem_types, 0);
              this->comm().sum(num_nodes);
              if (verbosity > 1)
              {
                this->comm().sum(info.volume);
                this->comm().min(info.bbox.min());
                this->comm().max(info.bbox.max());
              }
            }

          const bool has_name = sideset_name_map.count(id) && sideset_name_map.at(id).size();
          const std::string name = has_name ? sideset_name_map.at(id) : "";
          if (global)
            libmesh_assert(this->comm().verify(name));

          if (global ? this->processor_id() == 0 : info.num_sides > 0)
            {
              oss << "  Sideset " << id;
              if (has_name)
                oss << " (" << name << ")";
              oss << ", " << info.num_sides << " sides (" << elem_type_helper(info.side_elem_types) << ")"
                  << ", " << num_elems << " " << (global ? "" : "local ") << "elems (" << elem_type_helper(info.elem_types) << ")"
                  << ", " << num_nodes << " " << (global ? "" : "local ") << "nodes\n";

              if (verbosity > 1)
              {
                oss << "   " << (global ? "Side" : "Local side") << " volume: " << info.volume << "\n"
                    << "   " << (global ? "Bounding" : "Local bounding") << " box minimum: "
                    << info.bbox.min() << "\n"
                    << "   " << (global ? "Bounding" : "Local bounding") << " box maximum: "
                    << info.bbox.max() << "\n"
                    << "   " << (global ? "Bounding" : "Local bounding") << " box delta: "
                    << (info.bbox.max() - info.bbox.min()) << "\n";
              }
            }
        }

      // Reduce the edgeset ids
      auto edgeset_ids = this->get_boundary_info().get_edge_boundary_ids();
      if (global)
        this->comm().set_union(edgeset_ids);

      // Accumulate local information for each edgeset
      struct EdgesetInfo
      {
        std::size_t num_edges = 0;
        std::set<int> edge_elem_types;
        BoundingBox bbox;
      };
      std::map<boundary_id_type, EdgesetInfo> edgeset_info_map;
      std::unique_ptr<const Elem> edge;

      for (const auto & pair : this->get_boundary_info().get_edgeset_map())
        {
          const Elem * elem = pair.first;
          if (!include_object(*elem))
            continue;

          const auto id = pair.second.second;
          EdgesetInfo & info = edgeset_info_map[id];

          elem->build_edge_ptr(edge, pair.second.first);

          ++info.num_edges;
          info.edge_elem_types.insert(edge->type());

          if (verbosity > 1)
            info.bbox.union_with(edge->loose_bounding_box());
        }

      // Add edgeset info
      if (!global || this->processor_id() == 0)
        {
          oss << "\n " << (global ? "" : "Local ") << "Mesh Edgesets:\n";
          if (edgeset_ids.empty())
            oss << "  None\n";
        }

      const auto & edgeset_name_map = this->get_boundary_info().get_edgeset_name_map();
      for (const auto id : edgeset_ids)
        {
          EdgesetInfo & info = edgeset_info_map[id];

          // Reduce the local information for this edgeset if required
          if (global)
            {
              this->comm().sum(info.num_edges);
              this->comm().set_union(info.edge_elem_types, 0);
              if (verbosity > 1)
                {
                  this->comm().min(info.bbox.min());
                  this->comm().min(info.bbox.max());
                }
            }

          const bool has_name = edgeset_name_map.count(id) && edgeset_name_map.at(id).size();
          const std::string name = has_name ? edgeset_name_map.at(id) : "";
          if (global)
            libmesh_assert(this->comm().verify(name));

          if (global ? this->processor_id() == 0 : info.num_edges > 0)
            {
              oss << "  Edgeset " << id;
              if (has_name)
                oss << " (" << name << ")";
              oss << ", " << info.num_edges << " " << (global ? "" : "local ") << "edges ("
                  << elem_type_helper(info.edge_elem_types) << ")\n";

              if (verbosity > 1)
              {
                oss << "   " << (global ? "Bounding" : "Local bounding") << " box minimum: "
                    << info.bbox.min() << "\n"
                    << "   " << (global ? "Bounding" : "Local bounding") << " box maximum: "
                    << info.bbox.max() << "\n"
                    << "   " << (global ? "Bounding" : "Local bounding") << " box delta: "
                    << (info.bbox.max() - info.bbox.min()) << "\n";
              }
            }
        }

      // Reduce the block IDs and block names
      std::set<subdomain_id_type> subdomains;
      for (const Elem * elem : this->active_element_ptr_range())
        if (include_object(*elem))
          subdomains.insert(elem->subdomain_id());
      if (global)
        this->comm().set_union(subdomains);

      // Accumulate local information for each subdomain
      struct SubdomainInfo
      {
        std::size_t num_elems = 0;
        Real volume = 0;
        std::set<int> elem_types;
        std::set<dof_id_type> active_node_ids;
#ifdef LIBMESH_ENABLE_AMR
        std::size_t num_active_elems = 0;
#endif
        BoundingBox bbox;
      };
      std::map<subdomain_id_type, SubdomainInfo> subdomain_info_map;
      for (const Elem * elem : this->element_ptr_range())
        if (include_object(*elem))
          {
            SubdomainInfo & info = subdomain_info_map[elem->subdomain_id()];

            ++info.num_elems;
            info.elem_types.insert(elem->type());

#ifdef LIBMESH_ENABLE_AMR
            if (elem->active())
              ++info.num_active_elems;
#endif

            for (const Node & node : elem->node_ref_range())
              if (include_object(node) && node.active())
                info.active_node_ids.insert(node.id());

            if (verbosity > 1 && elem->active())
              {
                info.volume += elem->volume();
                info.bbox.union_with(elem->loose_bounding_box());
              }
          }

      // Add subdomain info
      oss << "\n " << (global ? "" : "Local ") << "Mesh Subdomains:\n";
      const auto & subdomain_name_map = this->get_subdomain_name_map();
      for (const auto id : subdomains)
      {
        SubdomainInfo & info = subdomain_info_map[id];

        auto num_active_nodes = info.active_node_ids.size();

        // Reduce the information for this subdomain if needed
        if (global)
          {
            this->comm().sum(info.num_elems);
#ifdef LIBMESH_ENABLE_AMR
            this->comm().sum(info.num_active_elems);
#endif
            this->comm().sum(num_active_nodes);
            this->comm().set_union(info.elem_types, 0);
            if (verbosity > 1)
            {
              this->comm().min(info.bbox.min());
              this->comm().max(info.bbox.max());
              this->comm().sum(info.volume);
            }
          }
        if (verbosity > 1)
          volume += info.volume;

        const bool has_name = subdomain_name_map.count(id);
        const std::string name = has_name ? subdomain_name_map.at(id) : "";
        if (global)
          libmesh_assert(this->comm().verify(name));

        if (!global || this->processor_id() == 0)
          {
            oss << "  Subdomain " << id;
            if (has_name)
              oss << " (" << name << ")";
            oss << ": " << info.num_elems << " " << (global ? "" : "local ") << "elems "
                << "(" << elem_type_helper(info.elem_types);
#ifdef LIBMESH_ENABLE_AMR
            oss << ", " << info.num_active_elems << " active";
#endif
            oss << "), " << num_active_nodes << " " << (global ? "" : "local ") << "active nodes\n";
            if (verbosity > 1)
            {
              oss << "   " << (global ? "Volume" : "Local volume") << ": " << info.volume << "\n";
              oss << "   " << (global ? "Bounding" : "Local bounding") << " box minimum: "
                  << info.bbox.min() << "\n"
                  << "   " << (global ? "Bounding" : "Local bounding") << " box maximum: "
                  << info.bbox.max() << "\n"
                  << "   " << (global ? "Bounding" : "Local bounding") << " box delta: "
                  << (info.bbox.max() - info.bbox.min()) << "\n";
            }
          }
      }

      oss << "  " << (global ? "Global" : "Local") << " mesh volume = " << volume << "\n";

    }

  return oss.str();
}


void MeshBase::print_info(std::ostream & os, const unsigned int verbosity /* = 0 */, const bool global /* = true */) const
{
  os << this->get_info(verbosity, global)
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
  // A nullptr partitioner or a skip_partitioning(true) call or a
  // skip_noncritical_partitioning(true) call means don't repartition;
  // skip_noncritical_partitioning() checks all these.
  else if (!skip_noncritical_partitioning())
    {
      partitioner()->partition (*this, n_parts);
    }
  else
    {
      // Adaptive coarsening may have "orphaned" nodes on processors
      // whose elements no longer share them.  We need to check for
      // and possibly fix that.
      MeshTools::correct_node_proc_ids(*this);

      // Make sure locally cached partition count is correct
      this->recalculate_n_partitions();

      // Make sure any other locally cached data is correct
      this->update_post_partitioning();
    }
}

unsigned int MeshBase::recalculate_n_partitions()
{
  // This requires an inspection on every processor
  parallel_object_only();

  unsigned int max_proc_id=0;

  for (const auto & elem : this->active_local_element_ptr_range())
    max_proc_id = std::max(max_proc_id, static_cast<unsigned int>(elem->processor_id()));

  // The number of partitions is one more than the max processor ID.
  _n_parts = max_proc_id+1;

  this->comm().max(_n_parts);

  return _n_parts;
}



std::unique_ptr<PointLocatorBase> MeshBase::sub_point_locator () const
{
  // If there's no master point locator, then we need one.
  if (_point_locator.get() == nullptr)
    {
      // PointLocator construction may not be safe within threads
      libmesh_assert(!Threads::in_threads);

      // And it may require parallel communication
      parallel_object_only();

#ifdef LIBMESH_ENABLE_NANOFLANN_POINTLOCATOR
      _point_locator = PointLocatorBase::build(NANOFLANN, *this);
#else
      _point_locator = PointLocatorBase::build(TREE_ELEMENTS, *this);
#endif

      if (_point_locator_close_to_point_tol > 0.)
        _point_locator->set_close_to_point_tol(_point_locator_close_to_point_tol);
    }

  // Otherwise there was a master point locator, and we can grab a
  // sub-locator easily.
  return
#ifdef LIBMESH_ENABLE_NANOFLANN_POINTLOCATOR
    PointLocatorBase::build(NANOFLANN, *this, _point_locator.get());
#else
    PointLocatorBase::build(TREE_ELEMENTS, *this, _point_locator.get());
#endif
}



void MeshBase::clear_point_locator ()
{
  _point_locator.reset(nullptr);
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
  libmesh_deprecated();

  this->cache_elem_data();
}

void MeshBase::cache_elem_data()
{
  // This requires an inspection on every processor
  parallel_object_only();

  // Need to clear _elem_dims first in case all elements of a
  // particular dimension have been deleted.
  _elem_dims.clear();
  _mesh_subdomains.clear();

  for (const auto & elem : this->active_element_ptr_range())
  {
    _elem_dims.insert(cast_int<unsigned char>(elem->dim()));
    _mesh_subdomains.insert(elem->subdomain_id());
  }

  if (!this->is_serial())
  {
    // Some different dimension elements may only live on other processors
    this->comm().set_union(_elem_dims);
    this->comm().set_union(_mesh_subdomains);
  }

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
#if LIBMESH_DIM > 1
  if (_spatial_dimension < 3)
    {
      for (const auto & node : this->node_ptr_range())
        {
          // Note: the exact floating point comparison is intentional,
          // we don't want to get tripped up by tolerances.
          if ((*node)(1) != 0.)
            {
              _spatial_dimension = 2;
#if LIBMESH_DIM == 2
              // If libmesh is compiled in 2D mode, this is the
              // largest spatial dimension possible so we can break
              // out.
              break;
#endif
            }

#if LIBMESH_DIM > 2
          if ((*node)(2) != 0.)
            {
              // Spatial dimension can't get any higher than this, so
              // we can break out.
              _spatial_dimension = 3;
              break;
            }
#endif
        }
    }
#endif // LIBMESH_DIM > 1
}

void MeshBase::detect_interior_parents()
{
  // This requires an inspection on every processor
  parallel_object_only();

  // Check if the mesh contains mixed dimensions. If so, then set interior parents, otherwise return.
  if (this->elem_dimensions().size() == 1)
    return;

  //This map will be used to set interior parents
  std::unordered_map<dof_id_type, std::vector<dof_id_type>> node_to_elem;

  for (const auto & elem : this->active_element_ptr_range())
    {
      // Populating the node_to_elem map, same as MeshTools::build_nodes_to_elem_map
      for (auto n : make_range(elem->n_vertices()))
        {
          libmesh_assert_less (elem->id(), this->max_elem_id());

          node_to_elem[elem->node_id(n)].push_back(elem->id());
        }
    }

  // Automatically set interior parents
  for (const auto & element : this->element_ptr_range())
    {
      // Ignore an 3D element or an element that already has an interior parent
      if (element->dim()>=LIBMESH_DIM || element->interior_parent())
        continue;

      // Start by generating a SET of elements that are dim+1 to the current
      // element at each vertex of the current element, thus ignoring interior nodes.
      // If one of the SET of elements is empty, then we will not have an interior parent
      // since an interior parent must be connected to all vertices of the current element
      std::vector<std::set<dof_id_type>> neighbors( element->n_vertices() );

      bool found_interior_parents = false;

      for (auto n : make_range(element->n_vertices()))
        {
          std::vector<dof_id_type> & element_ids = node_to_elem[element->node_id(n)];
          for (const auto & eid : element_ids)
            if (this->elem_ref(eid).dim() == element->dim()+1)
              neighbors[n].insert(eid);

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
          for (const auto & interior_parent_id : neighbors_0)
            {
              found_interior_parents = false;
              for (auto n : make_range(1u, element->n_vertices()))
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



void MeshBase::set_point_locator_close_to_point_tol(Real val)
{
  _point_locator_close_to_point_tol = val;
  if (_point_locator)
    {
      if (val > 0.)
        _point_locator->set_close_to_point_tol(val);
      else
        _point_locator->unset_close_to_point_tol();
    }
}



Real MeshBase::get_point_locator_close_to_point_tol() const
{
  return _point_locator_close_to_point_tol;
}



void MeshBase::size_elem_extra_integers()
{
  const std::size_t new_size = _elem_integer_names.size();
  for (auto elem : this->element_ptr_range())
    elem->add_extra_integers(new_size, _elem_integer_default_values);
}



void MeshBase::size_node_extra_integers()
{
  const std::size_t new_size = _node_integer_names.size();
  for (auto node : this->node_ptr_range())
    node->add_extra_integers(new_size, _node_integer_default_values);
}


std::pair<std::vector<unsigned int>, std::vector<unsigned int>>
MeshBase::merge_extra_integer_names(const MeshBase & other)
{
  std::pair<std::vector<unsigned int>, std::vector<unsigned int>> returnval;
  returnval.first = this->add_elem_integers(other._elem_integer_names, true, &other._elem_integer_default_values);
  returnval.second = this->add_node_integers(other._node_integer_names, true, &other._node_integer_default_values);
  return returnval;
}



void
MeshBase::post_dofobject_moves(MeshBase && other_mesh)
{
  // Now that all the DofObject moving is done, we can move the GhostingFunctor objects
  // which include the _default_ghosting,_ghosting_functors and _shared_functors. We also need
  // to set the mesh object associated with these functors to the assignee mesh.

   // _default_ghosting
  _default_ghosting = std::move(other_mesh._default_ghosting);
  _default_ghosting->set_mesh(this);

  // _ghosting_functors
  _ghosting_functors = std::move(other_mesh._ghosting_functors);

  for (const auto gf : _ghosting_functors )
  {
    gf->set_mesh(this);
  }

  // _shared_functors
  _shared_functors = std::move(other_mesh._shared_functors);

  for (const auto sf : _shared_functors )
  {
    (sf.second)->set_mesh(this);
  }

  // _constraint_rows
  _constraint_rows = std::move(other_mesh._constraint_rows);

  if (other_mesh.partitioner())
    _partitioner = std::move(other_mesh.partitioner());
}


} // namespace libMesh
