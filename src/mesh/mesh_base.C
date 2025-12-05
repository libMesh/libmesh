// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// Local includes
#include "libmesh/boundary_info.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/elem.h"
#include "libmesh/ghost_point_neighbors.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_communication.h"
#include "libmesh/mesh_serializer.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/parallel.h"
#include "libmesh/parallel_algebra.h"
#include "libmesh/parallel_fe_type.h"
#include "libmesh/partitioner.h"
#include "libmesh/point_locator_base.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/threads.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/enum_point_locator_type.h"
#include "libmesh/enum_to_string.h"
#include "libmesh/point_locator_nanoflann.h"
#include "libmesh/elem_side_builder.h"

// C++ includes
#include <algorithm> // for std::min
#include <map>       // for std::multimap
#include <memory>
#include <sstream>   // for std::ostringstream
#include <unordered_map>

#include "libmesh/periodic_boundaries.h"
#include "libmesh/periodic_boundary.h"

namespace libMesh
{



// ------------------------------------------------------------
// MeshBase class member functions
MeshBase::MeshBase (const Parallel::Communicator & comm_in,
                    unsigned char d) :
  ParallelObject (comm_in),
  boundary_info  (new BoundaryInfo(*this)), // BoundaryInfo has protected ctor, can't use std::make_unique
  _n_parts       (1),
  _default_mapping_type(LAGRANGE_MAP),
  _default_mapping_data(0),
  _preparation (),
  _point_locator (),
  _count_lower_dim_elems_in_point_locator(true),
  _partitioner   (),
#ifdef LIBMESH_ENABLE_UNIQUE_ID
  _next_unique_id(DofObject::invalid_unique_id),
#endif
  _interior_mesh(this),
  _skip_noncritical_partitioning(false),
  _skip_all_partitioning(libMesh::on_command_line("--skip-partitioning")),
  _skip_renumber_nodes_and_elements(false),
  _skip_find_neighbors(false),
  _allow_remote_element_removal(true),
  _allow_node_and_elem_unique_id_overlap(false),
  _spatial_dimension(d),
  _default_ghosting(std::make_unique<GhostPointNeighbors>(*this)),
  _point_locator_close_to_point_tol(0.)
{
  _elem_dims.insert(d);
  _ghosting_functors.push_back(_default_ghosting.get());
  libmesh_assert_less_equal (LIBMESH_DIM, 3);
  libmesh_assert_greater_equal (LIBMESH_DIM, d);
  libmesh_assert (libMesh::initialized());
}



MeshBase::MeshBase (const MeshBase & other_mesh) :
  ParallelObject (other_mesh),
  boundary_info  (new BoundaryInfo(*this)), // BoundaryInfo has protected ctor, can't use std::make_unique
  _n_parts       (other_mesh._n_parts),
  _default_mapping_type(other_mesh._default_mapping_type),
  _default_mapping_data(other_mesh._default_mapping_data),
  _preparation (other_mesh._preparation),
  _point_locator (),
  _count_lower_dim_elems_in_point_locator(other_mesh._count_lower_dim_elems_in_point_locator),
  _partitioner   (),
#ifdef LIBMESH_ENABLE_UNIQUE_ID
  _next_unique_id(other_mesh._next_unique_id),
#endif
  // If the other mesh interior_parent pointers just go back to
  // itself, so should we
  _interior_mesh((other_mesh._interior_mesh == &other_mesh) ?
                 this : other_mesh._interior_mesh),
  _skip_noncritical_partitioning(other_mesh._skip_noncritical_partitioning),
  _skip_all_partitioning(other_mesh._skip_all_partitioning),
  _skip_renumber_nodes_and_elements(other_mesh._skip_renumber_nodes_and_elements),
  _skip_find_neighbors(other_mesh._skip_find_neighbors),
  _allow_remote_element_removal(other_mesh._allow_remote_element_removal),
  _allow_node_and_elem_unique_id_overlap(other_mesh._allow_node_and_elem_unique_id_overlap),
  _elem_dims(other_mesh._elem_dims),
  _elem_default_orders(other_mesh._elem_default_orders),
  _supported_nodal_order(other_mesh._supported_nodal_order),
  _mesh_subdomains(other_mesh._mesh_subdomains),
  _elemset_codes_inverse_map(other_mesh._elemset_codes_inverse_map),
  _all_elemset_ids(other_mesh._all_elemset_ids),
  _spatial_dimension(other_mesh._spatial_dimension),
  _default_ghosting(std::make_unique<GhostPointNeighbors>(*this)),
  _point_locator_close_to_point_tol(other_mesh._point_locator_close_to_point_tol)
{
  const GhostingFunctor * const other_default_ghosting = other_mesh._default_ghosting.get();

  for (GhostingFunctor * const gf : other_mesh._ghosting_functors)
    {
      // If the other mesh is using default ghosting, then we will use our own
      // default ghosting
      if (gf == other_default_ghosting)
        {
          _ghosting_functors.push_back(_default_ghosting.get());
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

#ifdef LIBMESH_ENABLE_PERIODIC
  // Deep copy of all periodic boundaries
  if (other_mesh._disjoint_neighbor_boundary_pairs)
    {
      _disjoint_neighbor_boundary_pairs = std::make_unique<PeriodicBoundaries>();

      for (const auto & [id, pb] : *other_mesh._disjoint_neighbor_boundary_pairs)
        if (pb)
          (*_disjoint_neighbor_boundary_pairs)[id] = pb->clone();
    }
#endif

  // _elemset_codes stores pointers to entries in _elemset_codes_inverse_map,
  // so it is not possible to simply copy it directly from other_mesh
  for (const auto & [set, code] : _elemset_codes_inverse_map)
    _elemset_codes.emplace(code, &set);
}

MeshBase& MeshBase::operator= (MeshBase && other_mesh)
{
  LOG_SCOPE("operator=(&&)", "MeshBase");

  // Move assign as a ParallelObject.
  this->ParallelObject::operator=(other_mesh);

  _n_parts = other_mesh.n_partitions();
  _default_mapping_type = other_mesh.default_mapping_type();
  _default_mapping_data = other_mesh.default_mapping_data();
  _preparation = other_mesh._preparation;
  _point_locator = std::move(other_mesh._point_locator);
  _count_lower_dim_elems_in_point_locator = other_mesh.get_count_lower_dim_elems_in_point_locator();
#ifdef LIBMESH_ENABLE_UNIQUE_ID
  _next_unique_id = other_mesh.next_unique_id();
#endif
  // If the other mesh interior_parent pointers just go back to
  // itself, so should we
  _interior_mesh = (other_mesh._interior_mesh == &other_mesh) ?
                   this : other_mesh._interior_mesh;
  _skip_noncritical_partitioning = other_mesh.skip_noncritical_partitioning();
  _skip_all_partitioning = other_mesh.skip_partitioning();
  _skip_renumber_nodes_and_elements = !(other_mesh.allow_renumbering());
  _skip_find_neighbors = !(other_mesh.allow_find_neighbors());
  _allow_remote_element_removal = other_mesh.allow_remote_element_removal();
  _allow_node_and_elem_unique_id_overlap = other_mesh.allow_node_and_elem_unique_id_overlap();
  _block_id_to_name = std::move(other_mesh._block_id_to_name);
  _elem_dims = std::move(other_mesh.elem_dimensions());
  _elem_default_orders = std::move(other_mesh.elem_default_orders());
  _supported_nodal_order = other_mesh.supported_nodal_order();
  _mesh_subdomains = other_mesh._mesh_subdomains;
  _elemset_codes = std::move(other_mesh._elemset_codes);
  _elemset_codes_inverse_map = std::move(other_mesh._elemset_codes_inverse_map);
  _all_elemset_ids = std::move(other_mesh._all_elemset_ids);
  _spatial_dimension = other_mesh.spatial_dimension();
  _elem_integer_names = std::move(other_mesh._elem_integer_names);
  _elem_integer_default_values = std::move(other_mesh._elem_integer_default_values);
  _node_integer_names = std::move(other_mesh._node_integer_names);
  _node_integer_default_values = std::move(other_mesh._node_integer_default_values);
  _point_locator_close_to_point_tol = other_mesh.get_point_locator_close_to_point_tol();

#ifdef LIBMESH_ENABLE_PERIODIC
  // Deep copy of all periodic boundaries:
  // We must clone each PeriodicBoundaryBase in the source map,
  // since unique_ptr cannot be copied and we need independent instances
  if (other_mesh._disjoint_neighbor_boundary_pairs)
    {
      _disjoint_neighbor_boundary_pairs = std::make_unique<PeriodicBoundaries>();

      for (const auto & [id, pb] : *other_mesh._disjoint_neighbor_boundary_pairs)
        if (pb)
          (*_disjoint_neighbor_boundary_pairs)[id] = pb->clone();
    }
#endif

  // This relies on our subclasses *not* invalidating pointers when we
  // do their portion of the move assignment later!
  boundary_info = std::move(other_mesh.boundary_info);
  boundary_info->set_mesh(*this);

#ifdef DEBUG
  // Make sure that move assignment worked for pointers
  for (const auto & [set, code] : _elemset_codes_inverse_map)
    {
      auto it = _elemset_codes.find(code);
      libmesh_assert_msg(it != _elemset_codes.end(),
                         "Elemset code " << code << " not found in _elmset_codes container.");
      libmesh_assert_equal_to(it->second, &set);
    }
#endif

  // We're *not* really done at this point, but we have the problem
  // that some of our data movement might be expecting subclasses data
  // movement to happen first.  We'll let subclasses handle that by
  // calling our post_dofobject_moves()
  return *this;
}


bool MeshBase::operator== (const MeshBase & other_mesh) const
{
  LOG_SCOPE("operator==()", "MeshBase");

  bool is_equal = this->locally_equals(other_mesh);
  this->comm().min(is_equal);
  return is_equal;
}


bool MeshBase::locally_equals (const MeshBase & other_mesh) const
{
  // Check whether (almost) everything in the base is equal
  //
  // We don't check _next_unique_id here, because it's expected to
  // change in a DistributedMesh prepare_for_use(); it's conceptually
  // "mutable".
  //
  // We use separate if statements instead of logical operators here,
  // to make it easy to see the failing condition when using a
  // debugger to figure out why a MeshTools::valid_is_prepared(mesh)
  // is failing.
  if (_n_parts != other_mesh._n_parts)
    return false;
  if (_default_mapping_type != other_mesh._default_mapping_type)
    return false;
  if (_default_mapping_data != other_mesh._default_mapping_data)
    return false;
  if (_preparation != other_mesh._preparation)
    return false;
  if (_count_lower_dim_elems_in_point_locator !=
        other_mesh._count_lower_dim_elems_in_point_locator)
    return false;

  // We should either both have our own interior parents or both not;
  // but if we both don't then we can't really assert anything else
  // because pointing at the same interior mesh is fair but so is
  // pointing at two different copies of "the same" interior mesh.
  if ((_interior_mesh == this) !=
      (other_mesh._interior_mesh == &other_mesh))
    return false;

  if (_skip_noncritical_partitioning != other_mesh._skip_noncritical_partitioning)
    return false;
  if (_skip_all_partitioning != other_mesh._skip_all_partitioning)
    return false;
  if (_skip_renumber_nodes_and_elements != other_mesh._skip_renumber_nodes_and_elements)
    return false;
  if (_skip_find_neighbors != other_mesh._skip_find_neighbors)
    return false;
  if (_allow_remote_element_removal != other_mesh._allow_remote_element_removal)
    return false;
  if (_allow_node_and_elem_unique_id_overlap != other_mesh._allow_node_and_elem_unique_id_overlap)
    return false;
  if (_spatial_dimension != other_mesh._spatial_dimension)
    return false;
  if (_point_locator_close_to_point_tol != other_mesh._point_locator_close_to_point_tol)
    return false;
  if (_block_id_to_name != other_mesh._block_id_to_name)
    return false;
  if (_elem_dims != other_mesh._elem_dims)
    return false;
  if (_elem_default_orders != other_mesh._elem_default_orders)
    return false;
  if (_supported_nodal_order != other_mesh._supported_nodal_order)
    return false;
  if (_mesh_subdomains != other_mesh._mesh_subdomains)
    return false;
  if (_all_elemset_ids != other_mesh._all_elemset_ids)
    return false;
  if (_elem_integer_names != other_mesh._elem_integer_names)
    return false;
  if (_elem_integer_default_values != other_mesh._elem_integer_default_values)
    return false;
  if (_node_integer_names != other_mesh._node_integer_names)
    return false;
  if (_node_integer_default_values != other_mesh._node_integer_default_values)
    return false;
  if (static_cast<bool>(_default_ghosting) != static_cast<bool>(other_mesh._default_ghosting))
    return false;
  if (static_cast<bool>(_partitioner) != static_cast<bool>(other_mesh._partitioner))
    return false;
  if (*boundary_info != *other_mesh.boundary_info)
    return false;

  // First check whether the "existence" of the two pointers differs (one present, one absent)
  if (static_cast<bool>(_disjoint_neighbor_boundary_pairs) !=
      static_cast<bool>(other_mesh._disjoint_neighbor_boundary_pairs))
    return false;
  // If both exist, compare the contents (Weak Test: just compare sizes like `_ghosting_functors`)
  if (_disjoint_neighbor_boundary_pairs &&
      (_disjoint_neighbor_boundary_pairs->size() != other_mesh._disjoint_neighbor_boundary_pairs->size()))
    return false;

  const constraint_rows_type & other_rows =
    other_mesh.get_constraint_rows();
  for (const auto & [node, row] : this->_constraint_rows)
    {
      const dof_id_type node_id = node->id();
      const Node * other_node = other_mesh.query_node_ptr(node_id);
      if (!other_node)
        return false;

      auto it = other_rows.find(other_node);
      if (it == other_rows.end())
        return false;

      const auto & other_row = it->second;
      if (row.size() != other_row.size())
        return false;

      for (auto i : index_range(row))
        {
          const auto & [elem_pair, coef] = row[i];
          const auto & [other_elem_pair, other_coef] = other_row[i];
          libmesh_assert(elem_pair.first);
          libmesh_assert(other_elem_pair.first);
          if (elem_pair.first->id() !=
              other_elem_pair.first->id() ||
              elem_pair.second !=
              other_elem_pair.second ||
              coef != other_coef)
            return false;
        }
    }

  for (const auto & [elemset_code, elemset_ptr] : this->_elemset_codes)
    if (const auto it = other_mesh._elemset_codes.find(elemset_code);
        it == other_mesh._elemset_codes.end() || *elemset_ptr != *it->second)
      return false;

  // FIXME: we have no good way to compare ghosting functors, since
  // they're in a vector of pointers, and we have no way *at all*
  // to compare ghosting functors, since they don't have operator==
  // defined and we encourage users to subclass them.  We can check if
  // we have the same number, is all.
  if (_ghosting_functors.size() !=
      other_mesh._ghosting_functors.size())
    return false;

  // Same deal for partitioners.  We tested that we both have one or
  // both don't, but are they equivalent?  Let's guess "yes".

  // Now let the subclasses decide whether everything else is equal
  return this->subclass_locally_equals(other_mesh);
}


MeshBase::~MeshBase()
{
  this->MeshBase::clear();

  libmesh_exceptionless_assert (!libMesh::closed());
}



unsigned int MeshBase::mesh_dimension() const
{
  if (!_elem_dims.empty())
    return cast_int<unsigned int>(*_elem_dims.rbegin());
  return 0;
}



void MeshBase::set_elem_dimensions(std::set<unsigned char> elem_dims)
{
#ifdef DEBUG
  // In debug mode, we call cache_elem_data() and then make sure
  // the result actually agrees with what the user specified.
  parallel_object_only();

  this->cache_elem_data();
  libmesh_assert_msg(_elem_dims == elem_dims, \
                     "Specified element dimensions does not match true element dimensions!");
#endif

  _elem_dims = std::move(elem_dims);
}



void MeshBase::add_elemset_code(dof_id_type code, MeshBase::elemset_type id_set)
{
  // Populate inverse map, stealing id_set's resources
  auto [it1, inserted1] = _elemset_codes_inverse_map.emplace(std::move(id_set), code);

  // Reference to the newly inserted (or previously existing) id_set
  const auto & inserted_id_set = it1->first;

  // Keep track of all elemset ids ever added for O(1) n_elemsets()
  // performance. Only need to do this if we didn't know about this
  // id_set before...
  if (inserted1)
    _all_elemset_ids.insert(inserted_id_set.begin(), inserted_id_set.end());

  // Take the address of the newly emplaced set to use in
  // _elemset_codes, avoid duplicating std::set storage
  auto [it2, inserted2] = _elemset_codes.emplace(code, &inserted_id_set);

  // Throw an error if this code already exists with a pointer to a
  // different set of ids.
  libmesh_error_msg_if(!inserted2 && it2->second != &inserted_id_set,
                       "The elemset code " << code << " already exists with a different id_set.");
}



unsigned int MeshBase::n_elemsets() const
{
  return _all_elemset_ids.size();
}

void MeshBase::get_elemsets(dof_id_type elemset_code, MeshBase::elemset_type & id_set_to_fill) const
{
  // If we don't recognize this elemset_code, hand back an empty set
  id_set_to_fill.clear();

  if (const auto it = _elemset_codes.find(elemset_code);
      it != _elemset_codes.end())
    id_set_to_fill.insert(it->second->begin(), it->second->end());
}

dof_id_type MeshBase::get_elemset_code(const MeshBase::elemset_type & id_set) const
{
  auto it = _elemset_codes_inverse_map.find(id_set);
  return (it == _elemset_codes_inverse_map.end()) ? DofObject::invalid_id : it->second;
}

std::vector<dof_id_type> MeshBase::get_elemset_codes() const
{
  std::vector<dof_id_type> ret;
  ret.reserve(_elemset_codes.size());
  for (const auto & pr : _elemset_codes)
    ret.push_back(pr.first);
  return ret;
}

void MeshBase::change_elemset_code(dof_id_type old_code, dof_id_type new_code)
{
  // Look up elemset ids for old_code
  auto it = _elemset_codes.find(old_code);

  // If we don't have the old_code, then do nothing. Alternatively, we
  // could throw an error since trying to change an elemset code you
  // don't have could indicate there's a problem...
  if (it == _elemset_codes.end())
    return;

  // Make copy of the set of elemset ids. We are not changing these,
  // only updating the elemset code it corresponds to.
  elemset_type id_set_copy = *(it->second);

  // Look up the corresponding entry in the inverse map. Note: we want
  // the iterator because we are going to remove it.
  auto inverse_it = _elemset_codes_inverse_map.find(id_set_copy);
  libmesh_error_msg_if(inverse_it == _elemset_codes_inverse_map.end(),
                       "Expected _elemset_codes_inverse_map entry for elemset code " << old_code);

  // Erase entry from inverse map
  _elemset_codes_inverse_map.erase(inverse_it);

  // Erase entry from forward map
  _elemset_codes.erase(it);

  // Add new code with original set of ids.
  this->add_elemset_code(new_code, id_set_copy);

  // We can't update any actual elemset codes if there is no extra integer defined for it.
  if (!this->has_elem_integer("elemset_code"))
    return;

  // Get index of elemset_code extra integer
  unsigned int elemset_index = this->get_elem_integer_index("elemset_code");

  // Loop over all elems and update code
  for (auto & elem : this->element_ptr_range())
    {
      dof_id_type elemset_code =
        elem->get_extra_integer(elemset_index);

      if (elemset_code == old_code)
        elem->set_extra_integer(elemset_index, new_code);
    }
}

void MeshBase::change_elemset_id(elemset_id_type old_id, elemset_id_type new_id)
{
  // Early return if we don't have old_id
  if (!_all_elemset_ids.count(old_id))
    return;

  // Throw an error if the new_id is already used
  libmesh_error_msg_if(_all_elemset_ids.count(new_id),
                       "Cannot change elemset id " << old_id <<
                       " to " << new_id << ", " << new_id << " already exists.");

  // We will build up a new version of the inverse map so we can iterate over
  // the current one without invalidating anything.
  std::map<MeshBase::elemset_type, dof_id_type> new_elemset_codes_inverse_map;
  for (const auto & [id_set, elemset_code] : _elemset_codes_inverse_map)
    {
      auto id_set_copy = id_set;
      if (id_set_copy.count(old_id))
        {
          // Remove old_id, insert new_id
          id_set_copy.erase(old_id);
          id_set_copy.insert(new_id);
        }

      // Store in new version of map
      new_elemset_codes_inverse_map.emplace(id_set_copy, elemset_code);
    }

  // Swap existing map with newly-built one
  _elemset_codes_inverse_map.swap(new_elemset_codes_inverse_map);

  // Reconstruct _elemset_codes map
  _elemset_codes.clear();
  for (const auto & [id_set, elemset_code] : _elemset_codes_inverse_map)
    _elemset_codes.emplace(elemset_code, &id_set);

  // Update _all_elemset_ids
  _all_elemset_ids.erase(old_id);
  _all_elemset_ids.insert(new_id);
}

unsigned int MeshBase::spatial_dimension () const
{
  libmesh_assert(_preparation.has_cached_elem_data);

  return cast_int<unsigned int>(_spatial_dimension);
}



void MeshBase::set_spatial_dimension(unsigned char d)
{
  // The user can set the _spatial_dimension however they wish,
  // libMesh will only *increase* the spatial dimension, however,
  // never decrease it.
  _spatial_dimension = d;
}



unsigned int MeshBase::add_elem_integer(std::string name,
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
  _elem_integer_names.push_back(std::move(name));
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
      if (const auto it = name_indices.find(name);
          it != name_indices.end())
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



unsigned int MeshBase::get_elem_integer_index(std::string_view name) const
{
  for (auto i : index_range(_elem_integer_names))
    if (_elem_integer_names[i] == name)
      return i;

  libmesh_error_msg("Unknown elem integer " << name);
  return libMesh::invalid_uint;
}



bool MeshBase::has_elem_integer(std::string_view name) const
{
  for (auto & entry : _elem_integer_names)
    if (entry == name)
      return true;

  return false;
}



unsigned int MeshBase::add_node_integer(std::string name,
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
  _node_integer_names.push_back(std::move(name));
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
      if (const auto it = name_indices.find(name);
          it != name_indices.end())
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



unsigned int MeshBase::get_node_integer_index(std::string_view name) const
{
  for (auto i : index_range(_node_integer_names))
    if (_node_integer_names[i] == name)
      return i;

  libmesh_error_msg("Unknown node integer " << name);
  return libMesh::invalid_uint;
}



bool MeshBase::has_node_integer(std::string_view name) const
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

  _preparation.has_removed_orphaned_nodes = true;
}



#ifdef LIBMESH_ENABLE_DEPRECATED
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
#endif // LIBMESH_ENABLE_DEPRECATED




void MeshBase::prepare_for_use ()
{
  // Mark everything as unprepared, except for those things we've been
  // told we don't need to prepare, for backwards compatibility
  this->clear_point_locator();
  _preparation = false;
  _preparation.has_neighbor_ptrs = _skip_find_neighbors;
  _preparation.has_removed_remote_elements = !_allow_remote_element_removal;

  this->complete_preparation();
}


void MeshBase::complete_preparation()
{
  LOG_SCOPE("complete_preparation()", "MeshBase");

  parallel_object_only();

  libmesh_assert(this->comm().verify(this->is_serial()));

#ifdef DEBUG
  // If we don't go into this method with valid constraint rows, we're
  // only going to be able to make that worse.
  MeshTools::libmesh_assert_valid_constraint_rows(*this);

  // If this mesh thinks it's already  partially prepared, then in
  // optimized builds we'll trust it, but in debug builds we'll check.
  const bool was_partly_prepared = (_preparation == Preparation());
#endif

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
    {
      if (!_preparation.has_removed_orphaned_nodes ||
          !_preparation.has_synched_id_counts)
        this->renumber_nodes_and_elements();
    }
  else
    {
      if (!_preparation.has_removed_orphaned_nodes)
        this->remove_orphaned_nodes();
      if (!_preparation.has_synched_id_counts)
        this->update_parallel_id_counts();
    }

  // Let all the elements find their neighbors
  if (!_skip_find_neighbors && !_preparation.has_neighbor_ptrs)
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
  if (!_preparation.has_cached_elem_data)
    this->cache_elem_data();

  // Search the mesh for elements that have a neighboring element
  // of dim+1 and set that element as the interior parent
  if (!_preparation.has_interior_parent_ptrs)
    this->detect_interior_parents();

  // Fix up node unique ids in case mesh generation code didn't take
  // exceptional care to do so.
  //  MeshCommunication().make_node_unique_ids_parallel_consistent(*this);

  // We're going to still require that mesh generation code gets
  // element unique ids consistent.
#if defined(DEBUG) && defined(LIBMESH_ENABLE_UNIQUE_ID)
  MeshTools::libmesh_assert_valid_unique_ids(*this);
#endif

  // Allow our GhostingFunctor objects to reinit if necessary.
  // Do this before partitioning and redistributing, and before
  // deleting remote elements.
  if (!_preparation.has_reinit_ghosting_functors)
    this->reinit_ghosting_functors();

  // Partition the mesh unless *all* partitioning is to be skipped.
  // If only noncritical partitioning is to be skipped, the
  // partition() call will still check for orphaned nodes.
  if (!skip_partitioning() && !_preparation.is_partitioned)
    this->partition();
  else if (!this->n_unpartitioned_elem() &&
           !this->n_unpartitioned_nodes())
    _preparation.is_partitioned = true;

  // If we're using DistributedMesh, we'll probably want it
  // parallelized.
  if (this->_allow_remote_element_removal &&
      !_preparation.has_removed_remote_elements)
    this->delete_remote_elements();
  else
    _preparation.has_removed_remote_elements = true;

  // Much of our boundary info may have been for now-remote parts of the mesh,
  // in which case we don't want to keep local copies of data meant to be
  // local. On the other hand we may have deleted, or the user may have added in
  // a distributed fashion, boundary data that is meant to be global. So we
  // handle both of those scenarios here
  if (!_preparation.has_boundary_id_sets)
    this->get_boundary_info().regenerate_id_sets();

  if (!_skip_renumber_nodes_and_elements)
    this->renumber_nodes_and_elements();

  // The mesh is now prepared for use, with the possible exception of
  // partitioning that was supposed to be skipped, and it should know
  // it.
#ifndef NDEBUG
  Preparation completed_preparation = _preparation;
  if (skip_partitioning())
    completed_preparation.is_partitioned = true;
  libmesh_assert(completed_preparation);
#endif

#ifdef DEBUG
  // The if() here avoids both unnecessary work *and* stack overflow
  if (was_partly_prepared)
    libmesh_assert(MeshTools::valid_is_prepared(*this));

  MeshTools::libmesh_assert_valid_boundary_ids(*this);
#ifdef LIBMESH_ENABLE_UNIQUE_ID
  MeshTools::libmesh_assert_valid_unique_ids(*this);
#endif
#endif
}

void
MeshBase::reinit_ghosting_functors()
{
  for (auto & gf : _ghosting_functors)
    {
      libmesh_assert(gf);
      gf->mesh_reinit();
    }

  _preparation.has_reinit_ghosting_functors = true;
}

void MeshBase::clear ()
{
  // Reset the number of partitions
  _n_parts = 1;

  // Reset the preparation flags
  _preparation = false;

  // Clear boundary information
  if (boundary_info)
    boundary_info->clear();

  // Clear cached element data
  _elem_dims.clear();
  _elem_default_orders.clear();
  _supported_nodal_order = MAXIMUM;

  _elemset_codes.clear();
  _elemset_codes_inverse_map.clear();

  _constraint_rows.clear();

  // Clear our point locator.
  this->clear_point_locator();
}



void MeshBase::add_ghosting_functor(GhostingFunctor & ghosting_functor)
{
  // We used to implicitly support duplicate inserts to std::set
#ifdef LIBMESH_ENABLE_DEPRECATED
  _ghosting_functors.erase
    (std::remove(_ghosting_functors.begin(),
                 _ghosting_functors.end(),
                 &ghosting_functor),
     _ghosting_functors.end());
#endif

  // We shouldn't have two copies of the same functor
  libmesh_assert(std::find(_ghosting_functors.begin(),
                           _ghosting_functors.end(),
                           &ghosting_functor) ==
                 _ghosting_functors.end());

  _ghosting_functors.push_back(&ghosting_functor);
}



void MeshBase::remove_ghosting_functor(GhostingFunctor & ghosting_functor)
{
  auto raw_it = std::find(_ghosting_functors.begin(),
                          _ghosting_functors.end(), &ghosting_functor);

  // The DofMap has a "to_mesh" parameter that tells it to avoid
  // registering a new functor with the mesh, but it doesn't keep
  // track of which functors weren't added, so we'll support "remove a
  // functor that isn't there" just like we did with set::erase
  // before.
  if (raw_it != _ghosting_functors.end())
    _ghosting_functors.erase(raw_it);

  // We shouldn't have had two copies of the same functor
  libmesh_assert(std::find(_ghosting_functors.begin(),
                           _ghosting_functors.end(),
                           &ghosting_functor) ==
                 _ghosting_functors.end());

  if (const auto it = _shared_functors.find(&ghosting_functor);
      it != _shared_functors.end())
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
      for (const auto & elem : this->active_unpartitioned_element_ptr_range())
        ids.insert(elem->subdomain_id());

      // Some subdomains may only live on other processors
      this->comm().set_union(ids);
    }
}



void MeshBase::redistribute()
{
  // We now have all elements and nodes redistributed; our ghosting
  // functors should be ready to redistribute and/or recompute any
  // cached data they use too.
  for (auto & gf : as_range(this->ghosting_functors_begin(),
                            this->ghosting_functors_end()))
    gf->redistribute();
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

  if (!_elem_default_orders.empty())
    {
      oss << "  elem_default_orders()={";
      std::transform(_elem_default_orders.begin(),
                     --_elem_default_orders.end(),
                     std::ostream_iterator<std::string>(oss, ", "),
                     [](Order o)
                       { return Utility::enum_to_string<Order>(o); });
      oss << Utility::enum_to_string<Order>(*_elem_default_orders.rbegin());
      oss << "}\n";
    }

  oss << "  supported_nodal_order()=" << this->supported_nodal_order()                        << '\n'
      << "  spatial_dimension()="     << this->spatial_dimension()                            << '\n'
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
  oss << "  n_elemsets()="            << static_cast<std::size_t>(this->n_elemsets())         << '\n';
  if (!_elemset_codes.empty())
    oss << "    n_elemset_codes="     << _elemset_codes.size()                                << '\n';
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
            for (const Elem * elem : this->active_unpartitioned_element_ptr_range())
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
      for (const auto & [node, id] : this->get_boundary_info().get_nodeset_map())
        {
          if (!include_object(*node))
            continue;

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

      const auto & nodeset_name_map = this->get_boundary_info().get_nodeset_name_map();
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
                  this->comm().max(info.bbox.max());
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
      ElemSideBuilder side_builder;
      std::map<boundary_id_type, SidesetInfo> sideset_info_map;
      for (const auto & pair : this->get_boundary_info().get_sideset_map())
        {
          const Elem * elem = pair.first;
          if (!include_object(*elem))
            continue;

          const auto id = pair.second.second;
          SidesetInfo & info = sideset_info_map[id];

          const auto s = pair.second.first;
          const Elem & side = side_builder(*elem, s);

          ++info.num_sides;
          info.side_elem_types.insert(side.type());
          info.elem_types.insert(elem->type());
          info.elem_ids.insert(elem->id());

          for (const Node & node : side.node_ref_range())
            if (include_object(node))
              info.node_ids.insert(node.id());

          if (verbosity > 1)
          {
            info.volume += side.volume();
            info.bbox.union_with(side.loose_bounding_box());
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

  _preparation.is_partitioned = true;
}

void MeshBase::all_second_order (const bool full_ordered)
{
  this->all_second_order_range(this->element_ptr_range(), full_ordered);
}

void MeshBase::all_complete_order ()
{
  this->all_complete_order_range(this->element_ptr_range());
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

  if (const auto iter = _block_id_to_name.find(id);
      iter == _block_id_to_name.end())
    return empty;
  else
    return iter->second;
}




subdomain_id_type MeshBase::get_id_by_name(std::string_view name) const
{
  // Linear search over the map values.
  for (const auto & [sbd_id, sbd_name] : _block_id_to_name)
    if (sbd_name == name)
      return sbd_id;

  // If we made it here without returning, we don't have a subdomain
  // with the requested name, so return Elem::invalid_subdomain_id.
  return Elem::invalid_subdomain_id;
}

#ifdef LIBMESH_ENABLE_DEPRECATED
void MeshBase::cache_elem_dims()
{
  libmesh_deprecated();

  this->cache_elem_data();
}
#endif // LIBMESH_ENABLE_DEPRECATED

void MeshBase::cache_elem_data()
{
  // This requires an inspection on every processor
  parallel_object_only();

  // Need to clear containers first in case all elements of a
  // particular dimension/order/subdomain have been deleted.
  _elem_dims.clear();
  _elem_default_orders.clear();
  _mesh_subdomains.clear();
  _supported_nodal_order = MAXIMUM;

  for (const auto & elem : this->active_element_ptr_range())
  {
    _elem_dims.insert(cast_int<unsigned char>(elem->dim()));
    _elem_default_orders.insert(elem->default_order());
    _mesh_subdomains.insert(elem->subdomain_id());
    _supported_nodal_order =
      static_cast<Order>
        (std::min(static_cast<int>(_supported_nodal_order),
                  static_cast<int>(elem->supported_nodal_order())));
  }

  if (!this->is_serial())
  {
    // Some different dimension/order/subdomain elements may only live
    // on other processors
    this->comm().set_union(_elem_dims);
    this->comm().set_union(_elem_default_orders);
    this->comm().min(_supported_nodal_order);
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
  if (_spatial_dimension < LIBMESH_DIM)
    {
      for (const auto & node : this->node_ptr_range())
        {
          // Note: the exact floating point comparison is intentional,
          // we don't want to get tripped up by tolerances.
          if ((*node)(0) != 0. && _spatial_dimension < 1)
            _spatial_dimension = 1;

          if ((*node)(1) != 0. && _spatial_dimension < 2)
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

  _preparation.has_cached_elem_data = true;
}


void MeshBase::sync_subdomain_name_map()
{
  // This requires every processor
  parallel_object_only();

  this->comm().set_union(_block_id_to_name);
}


void MeshBase::detect_interior_parents()
{
  // This requires an inspection on every processor
  parallel_object_only();

  // This requires up-to-date mesh dimensions in cache
  libmesh_assert(_preparation.has_cached_elem_data);

  // Check if the mesh contains mixed dimensions. If so, then we may
  // have interior parents to set.  Otherwise return.
  if (this->elem_dimensions().size() == 1)
    {
      _preparation.has_interior_parent_ptrs = true;
      return;
    }

  // Do we have interior parent pointers going to a different mesh?
  // If so then we'll still check to make sure that's the only place
  // they go, so we can libmesh_not_implemented() if not.
  const bool separate_interior_mesh = (&(this->interior_mesh()) != this);

  // This map will be used to set interior parents
  std::unordered_map<dof_id_type, std::vector<dof_id_type>> node_to_elem;

  for (const auto & elem : this->element_ptr_range())
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
                  if (neighbors[n].count(interior_parent_id))
                    {
                      found_interior_parents = true;
                    }
                  else
                    {
                      found_interior_parents = false;
                      break;
                    }
                }

              if (found_interior_parents)
                {
                  element->set_interior_parent(this->elem_ptr(interior_parent_id));
                  break;
                }
            }

          // Do we have a mixed dimensional mesh that contains some of
          // its own interior parents, but we already expect to have
          // interior parents on a different mesh?  That's going to
          // take some work to support if anyone needs it.
          if (separate_interior_mesh)
            libmesh_not_implemented_msg
              ("interior_parent() values in multiple meshes are unsupported.");
        }
    }

  _preparation.has_interior_parent_ptrs = true;
}



#ifdef LIBMESH_ENABLE_PERIODIC
  /**
   * Register a pair of boundaries as disjoint neighbor boundary pairs.
   */
  void MeshBase::add_disjoint_neighbor_boundary_pairs(const boundary_id_type b1,
                                             const boundary_id_type b2,
                                             const RealVectorValue & translation)
    {
      // Lazily allocate the container the first time it’s needed
      if (!_disjoint_neighbor_boundary_pairs)
        _disjoint_neighbor_boundary_pairs = std::make_unique<PeriodicBoundaries>();

      PeriodicBoundaries & db = *_disjoint_neighbor_boundary_pairs;

      // Create forward and inverse boundary mappings
      PeriodicBoundary forward(translation);
      PeriodicBoundary inverse(translation * -1.0);

      forward.myboundary       = b1;
      forward.pairedboundary   = b2;
      inverse.myboundary       = b2;
      inverse.pairedboundary   = b1;

      // Add both directions into the container
      db.emplace(b1, forward.clone());
      db.emplace(b2, inverse.clone());
    }

  PeriodicBoundaries * MeshBase::get_disjoint_neighbor_boundary_pairs()
    {
      return _disjoint_neighbor_boundary_pairs.get();
    }

  const PeriodicBoundaries * MeshBase::get_disjoint_neighbor_boundary_pairs() const
    {
      return _disjoint_neighbor_boundary_pairs.get();
    }

    void MeshBase::remove_disjoint_boundary_pair(const boundary_id_type b1,
                                                   const boundary_id_type b2)
    {
      // Nothing to remove if not allocated or empty
      if (!_disjoint_neighbor_boundary_pairs || _disjoint_neighbor_boundary_pairs->empty())
        return;

      auto & pairs = *_disjoint_neighbor_boundary_pairs;

      // Helper to check and erase both directions
      auto erase_if_match = [](boundary_id_type key,
                              boundary_id_type pair,
                              PeriodicBoundaries & pb_map)
        {
          auto it = pb_map.find(key);
          if (it != pb_map.end())
            {
              const auto & pb = *(it->second);
              // Check both directions
              if ((pb.myboundary == key && pb.pairedboundary == pair) ||
                  (pb.pairedboundary == key && pb.myboundary == pair))
                pb_map.erase(it);
            }
        };

      erase_if_match(b1, b2, pairs);
      erase_if_match(b2, b1, pairs);
    }


#endif



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

  for (const auto & sf : _shared_functors )
  {
    (sf.second)->set_mesh(this);
  }

  // _constraint_rows
  _constraint_rows = std::move(other_mesh._constraint_rows);

  if (other_mesh.partitioner())
    _partitioner = std::move(other_mesh.partitioner());
}


void
MeshBase::copy_cached_data(const MeshBase & other_mesh)
{
  this->_spatial_dimension = other_mesh._spatial_dimension;
  this->_elem_dims = other_mesh._elem_dims;
  this->_elem_default_orders = other_mesh._elem_default_orders;
  this->_supported_nodal_order = other_mesh._supported_nodal_order;
  this->_mesh_subdomains = other_mesh._mesh_subdomains;
}


bool MeshBase::nodes_and_elements_equal(const MeshBase & other_mesh) const
{
  for (const auto & other_node : other_mesh.node_ptr_range())
    {
      const Node * node = this->query_node_ptr(other_node->id());
      if (!node)
        return false;
      if (*other_node != *node)
        return false;
    }
  for (const auto & node : this->node_ptr_range())
    if (!other_mesh.query_node_ptr(node->id()))
      return false;

  for (const auto & other_elem : other_mesh.element_ptr_range())
    {
      const Elem * elem = this->query_elem_ptr(other_elem->id());
      if (!elem)
        return false;
      if (!other_elem->topologically_equal(*elem))
        return false;
    }
  for (const auto & elem : this->element_ptr_range())
    if (!other_mesh.query_elem_ptr(elem->id()))
      return false;

  return true;
}


dof_id_type MeshBase::n_constraint_rows() const
{
  dof_id_type n_local_rows=0, n_unpartitioned_rows=0;
  for (const auto & [node, node_constraints] : _constraint_rows)
    {
      // Unpartitioned nodes
      if (node->processor_id() == DofObject::invalid_processor_id)
        n_unpartitioned_rows++;
      else if (node->processor_id() == this->processor_id())
        n_local_rows++;
    }

  this->comm().sum(n_local_rows);

  return n_unpartitioned_rows + n_local_rows;
}


void
MeshBase::copy_constraint_rows(const MeshBase & other_mesh)
{
  LOG_SCOPE("copy_constraint_rows(mesh)", "MeshBase");

  _constraint_rows.clear();

  const auto & other_constraint_rows = other_mesh.get_constraint_rows();
  for (const auto & [other_node, other_node_constraints] : other_constraint_rows)
  {
    const Node * const our_node = this->node_ptr(other_node->id());
    constraint_rows_mapped_type our_node_constraints;
    for (const auto & [other_inner_key_pair, constraint_value] : other_node_constraints)
    {
      const auto & [other_elem, local_node_id] = other_inner_key_pair;
      const Elem * const our_elem = this->elem_ptr(other_elem->id());
      our_node_constraints.emplace_back(std::make_pair(our_elem, local_node_id), constraint_value);
    }
    _constraint_rows[our_node] = std::move(our_node_constraints);
  }
}


template <typename T>
void
MeshBase::copy_constraint_rows(const SparseMatrix<T> & constraint_operator,
                               bool precondition_constraint_operator)
{
  LOG_SCOPE("copy_constraint_rows(mat)", "MeshBase");

  this->_constraint_rows.clear();

  // We're not going to support doing this distributed yet; it'd be
  // pointless unless we temporarily had a linear partitioning to
  // better match the constraint operator.
  MeshSerializer serialize(*this);

  // Our current mesh should already reflect the desired assembly space
  libmesh_error_msg_if(this->n_nodes() != constraint_operator.m(),
                       "Constraint operator matrix with " <<
                       constraint_operator.m() <<
                       "rows does not match this mesh with " <<
                       this->n_nodes() << " nodes");

  // First, find what new unconstrained DoFs we need to add.  We can't
  // iterate over columns in a SparseMatrix, so we'll iterate over
  // rows and keep track of columns.

  // If we have nodes that will work unconstrained, keep track of
  // their node ids and corresponding column indices.
  // existing_unconstrained_nodes[column_id] = node_id
  std::map<dof_id_type, dof_id_type> existing_unconstrained_columns;
  std::set<dof_id_type> existing_unconstrained_nodes;

  // In case we need new nodes, keep track of their columns.
  // columns[j][k] will be the kth row index and value of column j
  typedef
    std::unordered_map<dof_id_type,
                       std::vector<std::pair<dof_id_type, Real>>>
    columns_type;
  columns_type columns(constraint_operator.n());

  // If we need to precondition the constraint operator (e.g.  it's an
  // unpreconditioned extraction operator for a Flex IGA matrix),
  // we'll want to keep track of the sum of each column, because we'll
  // be dividing each column by that sum (Jacobi preconditioning on
  // the right, which then leads to symmetric preconditioning on a
  // physics Jacobian).
  std::unordered_map<dof_id_type, Real> column_sums;

  // Work in parallel, though we'll have to sync shortly
  for (auto i : make_range(constraint_operator.row_start(),
                           constraint_operator.row_stop()))
    {
      std::vector<numeric_index_type> indices;
      std::vector<T> values;

      constraint_operator.get_row(i, indices, values);
      libmesh_assert_equal_to(indices.size(), values.size());

      if (indices.size() == 1 &&
          values[0] == T(1))
        {
          // If we have multiple simple Ui=Uj constraints, let the
          // first one be our "unconstrained" node and let the others
          // be constrained to it.
          if (existing_unconstrained_columns.find(indices[0]) !=
              existing_unconstrained_columns.end())
            {
              const auto j = indices[0];
              columns[j].emplace_back(i, 1);
            }
          else
            {
              existing_unconstrained_nodes.insert(i);
              existing_unconstrained_columns.emplace(indices[0],i);
            }
        }
      else
        for (auto jj : index_range(indices))
          {
            const auto j = indices[jj];
            const Real coef = libmesh_real(values[jj]);
            libmesh_assert_equal_to(coef, values[jj]);
            columns[j].emplace_back(i, coef);
          }
    }

  // Merge data from different processors' slabs of the matrix
  this->comm().set_union(existing_unconstrained_nodes);
  this->comm().set_union(existing_unconstrained_columns);

  std::vector<columns_type> all_columns;
  this->comm().allgather(columns, all_columns);

  columns.clear();
  for (auto p : index_range(all_columns))
    for (auto & [j, subcol] : all_columns[p])
      for (auto [i, v] : subcol)
        columns[j].emplace_back(i,v);

  // Keep track of elements on which unconstrained nodes exist, and
  // their local node indices.
  // node_to_elem_ptrs[node] = [elem_id, local_node_num]
  std::unordered_map<const Node *, std::pair<dof_id_type, unsigned int>> node_to_elem_ptrs;

  // Find elements attached to any existing nodes that will stay
  // unconstrained.  We'll also build a subdomain set here so we don't
  // have to assert that the mesh is already prepared before we pick a
  // new subdomain for any NodeElems we need to add.
  std::set<subdomain_id_type> subdomain_ids;
  for (const Elem * elem : this->element_ptr_range())
    {
      subdomain_ids.insert(elem->subdomain_id());
      for (auto n : make_range(elem->n_nodes()))
        {
          const Node * node = elem->node_ptr(n);
          if (existing_unconstrained_nodes.count(node->id()))
            node_to_elem_ptrs.emplace(node, std::make_pair(elem->id(), n));
        }
    }

  const subdomain_id_type new_sbd_id = *subdomain_ids.rbegin() + 1;

  for (auto j : make_range(constraint_operator.n()))
    {
      // If we already have a good node for this then we're done
      if (existing_unconstrained_columns.count(j))
        continue;

      // Get a half-decent spot to place a new NodeElem for
      // unconstrained DoF(s) here.  Getting a *fully*-decent spot
      // would require finding a Moore-Penrose pseudoinverse, and I'm
      // not going to do that, but scaling a transpose will at least
      // get us a little uniqueness to make visualization reasonable.
      Point newpt;
      Real total_scaling = 0;
      unsigned int total_entries = 0;

      // We'll get a decent initial pid choice here too, if only to
      // aid in later repartitioning.
      std::map<processor_id_type, int> pids;

      auto & column = columns[j];
      for (auto [i, r] : column)
        {
          Node & constrained_node = this->node_ref(i);
          const Point constrained_pt = constrained_node;
          newpt += r*constrained_pt;
          total_scaling += r;
          ++total_entries;
          ++pids[constrained_node.processor_id()];
        }

      if (precondition_constraint_operator)
        column_sums[j] = total_scaling;

      libmesh_error_msg_if
        (!total_entries,
         "Empty column " << j <<
         " found in constraint operator matrix");

      // If we have *cancellation* here then we can end up dividing by
      // zero; try just evenly scaling across all constrained node
      // points instead.
      if (total_scaling > TOLERANCE)
        newpt /= total_scaling;
      else
        newpt /= total_entries;

      Node *n = this->add_point(newpt);
      std::unique_ptr<Elem> elem = Elem::build(NODEELEM);
      elem->set_node(0, n);
      elem->subdomain_id() = new_sbd_id;

      Elem * added_elem = this->add_elem(std::move(elem));
      this->_elem_dims.insert(0);
      this->_elem_default_orders.insert(added_elem->default_order());
      this->_supported_nodal_order =
        static_cast<Order>
          (std::min(static_cast<int>(this->_supported_nodal_order),
                    static_cast<int>(added_elem->supported_nodal_order())));
      this->_mesh_subdomains.insert(new_sbd_id);
      node_to_elem_ptrs.emplace(n, std::make_pair(added_elem->id(), 0));
      existing_unconstrained_columns.emplace(j,n->id());

      // Repartition the new objects *after* adding them, so a
      // DistributedMesh doesn't get confused and think you're not
      // adding them on all processors at once.
      int n_pids = 0;
      for (auto [pid, count] : pids)
        if (count >= n_pids)
          {
            n_pids = count;
            added_elem->processor_id() = pid;
            n->processor_id() = pid;
          }
    }

  // Calculate constraint rows in an indexed form that's easy for us
  // to allgather
  std::unordered_map<dof_id_type,
    std::vector<std::pair<std::pair<dof_id_type, unsigned int>,Real>>>
    indexed_constraint_rows;

  for (auto i : make_range(constraint_operator.row_start(),
                           constraint_operator.row_stop()))
    {
      if (existing_unconstrained_nodes.count(i))
        continue;

      std::vector<numeric_index_type> indices;
      std::vector<T> values;

      constraint_operator.get_row(i, indices, values);

      std::vector<std::pair<std::pair<dof_id_type, unsigned int>, Real>> constraint_row;

      for (auto jj : index_range(indices))
        {
          const dof_id_type node_id =
            existing_unconstrained_columns[indices[jj]];

          Node & constraining_node = this->node_ref(node_id);

          libmesh_assert(node_to_elem_ptrs.count(&constraining_node));

          auto p = node_to_elem_ptrs[&constraining_node];

          Real coef = libmesh_real(values[jj]);
          libmesh_assert_equal_to(coef, values[jj]);

          // If we're preconditioning and we created a nodeelem then
          // we can scale the meaning of that nodeelem's value to give
          // us a better-conditioned matrix after the constraints are
          // applied.
          if (precondition_constraint_operator)
            if (auto sum_it = column_sums.find(indices[jj]);
                sum_it != column_sums.end())
              {
                const Real scaling = sum_it->second;

                if (scaling > TOLERANCE)
                  coef /= scaling;
              }

          constraint_row.emplace_back(std::make_pair(p, coef));
        }

      indexed_constraint_rows.emplace(i, std::move(constraint_row));
    }

  this->comm().set_union(indexed_constraint_rows);

  // Add constraint rows as mesh constraint rows
  for (auto & [node_id, indexed_row] : indexed_constraint_rows)
    {
      Node * constrained_node = this->node_ptr(node_id);

      constraint_rows_mapped_type constraint_row;

      for (auto [p, coef] : indexed_row)
        {
          const Elem * elem = this->elem_ptr(p.first);
          constraint_row.emplace_back
            (std::make_pair(std::make_pair(elem, p.second), coef));
        }

      this->_constraint_rows.emplace(constrained_node,
                                     std::move(constraint_row));
    }
}


void MeshBase::print_constraint_rows(std::ostream & os,
                                     bool print_nonlocal) const
{
  parallel_object_only();

  std::string local_constraints =
    this->get_local_constraints(print_nonlocal);

  if (this->processor_id())
    {
      this->comm().send(0, local_constraints);
    }
  else
    {
      os << "Processor 0:\n";
      os << local_constraints;

      for (auto p : IntRange<processor_id_type>(1, this->n_processors()))
        {
          this->comm().receive(p, local_constraints);
          os << "Processor " << p << ":\n";
          os << local_constraints;
        }
    }
}



std::string MeshBase::get_local_constraints(bool print_nonlocal) const
{
  std::ostringstream os;

  if (print_nonlocal)
    os << "All ";
  else
    os << "Local ";

  os << "Mesh Constraint Rows:"
     << std::endl;

  for (const auto & [node, row] : _constraint_rows)
    {
      const bool local = (node->processor_id() == this->processor_id());

      // Skip non-local dofs if requested
      if (!print_nonlocal && !local)
        continue;

      os << "Constraints for " << (local ? "Local" : "Ghost") << " Node " << node->id()
         << ": \t";

      for (const auto & [elem_and_node, coef] : row)
        os << " ((" << elem_and_node.first->id() << ',' << elem_and_node.second << "), " << coef << ")\t";

      os << std::endl;
    }

  return os.str();
}




// Explicit instantiations for our template function
template LIBMESH_EXPORT void
MeshBase::copy_constraint_rows(const SparseMatrix<Real> & constraint_operator,
                               bool precondition_constraint_operator);

#ifdef LIBMESH_USE_COMPLEX_NUMBERS
template LIBMESH_EXPORT void
MeshBase::copy_constraint_rows(const SparseMatrix<Complex> & constraint_operator,
                               bool precondition_constraint_operator);
#endif


} // namespace libMesh
