// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



// Local includes
#include "libmesh/elem.h"
#include "libmesh/elem_range.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_communication.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/node_range.h"
#include "libmesh/parallel.h"
#include "libmesh/parallel_ghost_sync.h"
#include "libmesh/sphere.h"
#include "libmesh/threads.h"
#include "libmesh/string_to_enum.h"

#ifdef DEBUG
#  include "libmesh/remote_elem.h"
#endif

#include LIBMESH_INCLUDE_UNORDERED_MAP
#include LIBMESH_INCLUDE_UNORDERED_SET

// C++ includes
#include <limits>
#include <numeric> // for std::accumulate
#include <set>



// ------------------------------------------------------------
// anonymous namespace for helper classes
namespace {

using namespace libMesh;

/**
 * SumElemWeight(Range) sums the number of nodes per element
 * for each element in the provided range. The join() method
 * defines how to combine the reduction operation from two
 * distinct instances of this class which may be executed on
 * separate threads.
 */
class SumElemWeight
{
public:
  SumElemWeight () :
    _weight(0)
  {}

  SumElemWeight (SumElemWeight &, Threads::split) :
    _weight(0)
  {}

  void operator()(const ConstElemRange & range)
  {
    for (ConstElemRange::const_iterator it = range.begin(); it !=range.end(); ++it)
      _weight += (*it)->n_nodes();
  }

  dof_id_type weight() const
  { return _weight; }

  // If we don't have threads we never need a join, and icpc yells a
  // warning if it sees an anonymous function that's never used
#if LIBMESH_USING_THREADS
  void join (const SumElemWeight & other)
  { _weight += other.weight(); }
#endif

private:
  dof_id_type _weight;
};


/**
 * FindBBox(Range) computes the bounding box for the objects
 * in the specified range.  This class may be split and subranges
 * can be executed on separate threads.  The join() method
 * defines how the results from two separate threads are combined.
 */
class FindBBox
{
public:
  FindBBox () :
    _vmin(LIBMESH_DIM,  std::numeric_limits<Real>::max()),
    _vmax(LIBMESH_DIM, -std::numeric_limits<Real>::max())
  {}

  FindBBox (FindBBox & other, Threads::split) :
    _vmin(other._vmin),
    _vmax(other._vmax)
  {}

  std::vector<Real> & min() { return _vmin; }
  std::vector<Real> & max() { return _vmax; }

  void operator()(const ConstNodeRange & range)
  {
    for (ConstNodeRange::const_iterator it = range.begin(); it != range.end(); ++it)
      {
        const Node * node = *it;
        libmesh_assert(node);

        for (unsigned int i=0; i<LIBMESH_DIM; i++)
          {
            _vmin[i] = std::min(_vmin[i], (*node)(i));
            _vmax[i] = std::max(_vmax[i], (*node)(i));
          }
      }
  }

  void operator()(const ConstElemRange & range)
  {
    for (ConstElemRange::const_iterator it = range.begin(); it != range.end(); ++it)
      {
        const Elem * elem = *it;
        libmesh_assert(elem);

        for (unsigned int n=0; n<elem->n_nodes(); n++)
          {
            const Point & point = elem->point(n);

            for (unsigned int i=0; i<LIBMESH_DIM; i++)
              {
                _vmin[i] = std::min(_vmin[i], point(i));
                _vmax[i] = std::max(_vmax[i], point(i));
              }
          }
      }
  }

  // If we don't have threads we never need a join, and icpc yells a
  // warning if it sees an anonymous function that's never used
#if LIBMESH_USING_THREADS
  void join (const FindBBox & other)
  {
    for (unsigned int i=0; i<LIBMESH_DIM; i++)
      {
        _vmin[i] = std::min(_vmin[i], other._vmin[i]);
        _vmax[i] = std::max(_vmax[i], other._vmax[i]);
      }
  }
#endif

  MeshTools::BoundingBox bbox () const
  {
    Point pmin(_vmin[0]
#if LIBMESH_DIM > 1
               , _vmin[1]
#endif
#if LIBMESH_DIM > 2
               , _vmin[2]
#endif
               );
    Point pmax(_vmax[0]
#if LIBMESH_DIM > 1
               , _vmax[1]
#endif
#if LIBMESH_DIM > 2
               , _vmax[2]
#endif
               );

    const MeshTools::BoundingBox ret_val(pmin, pmax);

    return ret_val;
  }

private:
  std::vector<Real> _vmin;
  std::vector<Real> _vmax;
};

#ifdef DEBUG
void assert_semiverify_dofobj(const Parallel::Communicator & communicator,
                              const DofObject * d)
{
  if (d)
    {
      const unsigned int n_sys = d->n_systems();

      std::vector<unsigned int> n_vars (n_sys, 0);
      for (unsigned int s = 0; s != n_sys; ++s)
        n_vars[s] = d->n_vars(s);

      const unsigned int tot_n_vars =
        std::accumulate(n_vars.begin(), n_vars.end(), 0);

      std::vector<unsigned int> n_comp (tot_n_vars, 0);
      std::vector<dof_id_type> first_dof (tot_n_vars, 0);

      for (unsigned int s = 0, i=0; s != n_sys; ++s)
        for (unsigned int v = 0; v != n_vars[s]; ++v, ++i)
          {
            n_comp[i] = d->n_comp(s,v);
            first_dof[i] = n_comp[i] ? d->dof_number(s,v,0) : DofObject::invalid_id;
          }

      libmesh_assert(communicator.semiverify(&n_sys));
      libmesh_assert(communicator.semiverify(&n_vars));
      libmesh_assert(communicator.semiverify(&n_comp));
      libmesh_assert(communicator.semiverify(&first_dof));
    }
  else
    {
      const unsigned int * p_ui = libmesh_nullptr;
      const std::vector<unsigned int> * p_vui = libmesh_nullptr;
      const std::vector<dof_id_type> * p_vdid = libmesh_nullptr;

      libmesh_assert(communicator.semiverify(p_ui));
      libmesh_assert(communicator.semiverify(p_vui));
      libmesh_assert(communicator.semiverify(p_vui));
      libmesh_assert(communicator.semiverify(p_vdid));
    }
}
#endif // DEBUG

}


namespace libMesh
{
// Small helper function to make intersect more readable.
bool is_between(Real min, Real check, Real max)
{
  return min <= check && check <= max;
}

bool MeshTools::BoundingBox::intersect (const BoundingBox & other_box) const
{
  // Make local variables first to make thiings more clear in a moment
  const Real & my_min_x = this->first(0);
  const Real & my_max_x = this->second(0);
  const Real & other_min_x = other_box.first(0);
  const Real & other_max_x = other_box.second(0);

  const bool x_int = is_between(my_min_x, other_min_x, my_max_x) || is_between(my_min_x, other_max_x, my_max_x) ||
    is_between(other_min_x, my_min_x, other_max_x) || is_between(other_min_x, my_max_x, other_max_x);

  bool intersection_true = x_int;

#if LIBMESH_DIM > 1
  const Real & my_min_y = this->first(1);
  const Real & my_max_y = this->second(1);
  const Real & other_min_y = other_box.first(1);
  const Real & other_max_y = other_box.second(1);

  const bool y_int = is_between(my_min_y, other_min_y, my_max_y) || is_between(my_min_y, other_max_y, my_max_y) ||
    is_between(other_min_y, my_min_y, other_max_y) || is_between(other_min_y, my_max_y, other_max_y);

  intersection_true = intersection_true && y_int;
#endif

#if LIBMESH_DIM > 2
  const Real & my_min_z = this->first(2);
  const Real & my_max_z = this->second(2);
  const Real & other_min_z = other_box.first(2);
  const Real & other_max_z = other_box.second(2);

  const bool z_int = is_between(my_min_z, other_min_z, my_max_z) || is_between(my_min_z, other_max_z, my_max_z) ||
    is_between(other_min_z, my_min_z, other_max_z) || is_between(other_min_z, my_max_z, other_max_z);

  intersection_true = intersection_true && z_int;
#endif

  return intersection_true;
}

bool MeshTools::BoundingBox::contains_point (const Point & p) const
{
  // Make local variables first to make thiings more clear in a moment
  Real my_min_x = this->first(0);
  Real my_max_x = this->second(0);
  bool x_int = is_between(my_min_x, p(0), my_max_x);

  bool intersection_true = x_int;

#if LIBMESH_DIM > 1
  Real my_min_y = this->first(1);
  Real my_max_y = this->second(1);
  bool y_int = is_between(my_min_y, p(1), my_max_y);

  intersection_true = intersection_true && y_int;
#endif


#if LIBMESH_DIM > 2
  Real my_min_z = this->first(2);
  Real my_max_z = this->second(2);
  bool z_int = is_between(my_min_z, p(2), my_max_z);

  intersection_true = intersection_true && z_int;
#endif

  return intersection_true;
}

// ------------------------------------------------------------
// MeshTools functions
dof_id_type MeshTools::total_weight(const MeshBase & mesh)
{
  if (!mesh.is_serial())
    {
      libmesh_parallel_only(mesh.comm());
      dof_id_type weight = MeshTools::weight (mesh, mesh.processor_id());
      mesh.comm().sum(weight);
      dof_id_type unpartitioned_weight =
        MeshTools::weight (mesh, DofObject::invalid_processor_id);
      return weight + unpartitioned_weight;
    }

  SumElemWeight sew;

  Threads::parallel_reduce (ConstElemRange (mesh.elements_begin(),
                                            mesh.elements_end()),
                            sew);
  return sew.weight();

}



dof_id_type MeshTools::weight(const MeshBase & mesh, const processor_id_type pid)
{
  SumElemWeight sew;

  Threads::parallel_reduce (ConstElemRange (mesh.pid_elements_begin(pid),
                                            mesh.pid_elements_end(pid)),
                            sew);
  return sew.weight();
}



void MeshTools::build_nodes_to_elem_map (const MeshBase & mesh,
                                         std::vector<std::vector<dof_id_type> > & nodes_to_elem_map)
{
  nodes_to_elem_map.resize (mesh.n_nodes());

  MeshBase::const_element_iterator       el  = mesh.elements_begin();
  const MeshBase::const_element_iterator end = mesh.elements_end();

  for (; el != end; ++el)
    for (unsigned int n=0; n<(*el)->n_nodes(); n++)
      {
        libmesh_assert_less ((*el)->node_id(n), nodes_to_elem_map.size());
        libmesh_assert_less ((*el)->id(), mesh.n_elem());

        nodes_to_elem_map[(*el)->node_id(n)].push_back((*el)->id());
      }
}



void MeshTools::build_nodes_to_elem_map (const MeshBase & mesh,
                                         std::vector<std::vector<const Elem *> > & nodes_to_elem_map)
{
  nodes_to_elem_map.resize (mesh.n_nodes());

  MeshBase::const_element_iterator       el  = mesh.elements_begin();
  const MeshBase::const_element_iterator end = mesh.elements_end();

  for (; el != end; ++el)
    for (unsigned int n=0; n<(*el)->n_nodes(); n++)
      {
        libmesh_assert_less ((*el)->node_id(n), nodes_to_elem_map.size());

        nodes_to_elem_map[(*el)->node_id(n)].push_back(*el);
      }
}



void MeshTools::find_boundary_nodes (const MeshBase & mesh,
                                     std::vector<bool> & on_boundary)
{
  // Resize the vector which holds boundary nodes and fill with false.
  on_boundary.resize(mesh.max_node_id());
  std::fill(on_boundary.begin(),
            on_boundary.end(),
            false);

  // Loop over elements, find those on boundary, and
  // mark them as true in on_boundary.
  MeshBase::const_element_iterator       el  = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end = mesh.active_elements_end();

  for (; el != end; ++el)
    {
      const Elem * elem = *el;

      for (unsigned int s=0; s<elem->n_neighbors(); s++)
        if (elem->neighbor_ptr(s) == libmesh_nullptr) // on the boundary
          {
            const UniquePtr<const Elem> side = elem->build_side_ptr(s);

            for (unsigned int n=0; n<side->n_nodes(); n++)
              on_boundary[side->node_id(n)] = true;
          }
    }
}



MeshTools::BoundingBox
MeshTools::bounding_box(const MeshBase & mesh)
{
  // This function must be run on all processors at once
  libmesh_parallel_only(mesh.comm());

  FindBBox find_bbox;

  Threads::parallel_reduce (ConstNodeRange (mesh.local_nodes_begin(),
                                            mesh.local_nodes_end()),
                            find_bbox);

  // and the unpartitioned nodes
  Threads::parallel_reduce (ConstNodeRange (mesh.pid_nodes_begin(DofObject::invalid_processor_id),
                                            mesh.pid_nodes_end(DofObject::invalid_processor_id)),
                            find_bbox);

  // Compare the bounding boxes across processors
  mesh.comm().min(find_bbox.min());
  mesh.comm().max(find_bbox.max());

  return find_bbox.bbox();
}



Sphere
MeshTools::bounding_sphere(const MeshBase & mesh)
{
  BoundingBox bbox = bounding_box(mesh);

  const Real  diag = (bbox.second - bbox.first).norm();
  const Point cent = (bbox.second + bbox.first)/2;

  return Sphere (cent, .5*diag);
}



MeshTools::BoundingBox
MeshTools::processor_bounding_box (const MeshBase & mesh,
                                   const processor_id_type pid)
{
  libmesh_assert_less (pid, mesh.n_processors());

  FindBBox find_bbox;

  Threads::parallel_reduce (ConstElemRange (mesh.pid_elements_begin(pid),
                                            mesh.pid_elements_end(pid)),
                            find_bbox);

  // Compare the bounding boxes across processors
  mesh.comm().min(find_bbox.min());
  mesh.comm().max(find_bbox.max());

  return find_bbox.bbox();
}



Sphere
MeshTools::processor_bounding_sphere (const MeshBase & mesh,
                                      const processor_id_type pid)
{
  BoundingBox bbox = processor_bounding_box(mesh,pid);

  const Real  diag = (bbox.second - bbox.first).norm();
  const Point cent = (bbox.second + bbox.first)/2;

  return Sphere (cent, .5*diag);
}



MeshTools::BoundingBox
MeshTools::subdomain_bounding_box (const MeshBase & mesh,
                                   const subdomain_id_type sid)
{
  FindBBox find_bbox;

  Threads::parallel_reduce
    (ConstElemRange (mesh.active_local_subdomain_elements_begin(sid),
                     mesh.active_local_subdomain_elements_end(sid)),
     find_bbox);

  // Compare the bounding boxes across processors
  mesh.comm().min(find_bbox.min());
  mesh.comm().max(find_bbox.max());

  return find_bbox.bbox();
}



Sphere
MeshTools::subdomain_bounding_sphere (const MeshBase & mesh,
                                      const subdomain_id_type sid)
{
  BoundingBox bbox = subdomain_bounding_box(mesh,sid);

  const Real  diag = (bbox.second - bbox.first).norm();
  const Point cent = (bbox.second + bbox.first)/2;

  return Sphere (cent, .5*diag);
}



void MeshTools::elem_types (const MeshBase & mesh,
                            std::vector<ElemType> & et)
{
  MeshBase::const_element_iterator       el  = mesh.elements_begin();
  const MeshBase::const_element_iterator end = mesh.elements_end();

  // Automatically get the first type
  et.push_back((*el)->type());  ++el;

  // Loop over the rest of the elements.
  // If the current element type isn't in the
  // vector, insert it.
  for (; el != end; ++el)
    if (!std::count(et.begin(), et.end(), (*el)->type()))
      et.push_back((*el)->type());
}



dof_id_type MeshTools::n_elem_of_type (const MeshBase & mesh,
                                       const ElemType type)
{
  return static_cast<dof_id_type>(std::distance(mesh.type_elements_begin(type),
                                                mesh.type_elements_end  (type)));
}



dof_id_type MeshTools::n_active_elem_of_type (const MeshBase & mesh,
                                              const ElemType type)
{
  return static_cast<dof_id_type>(std::distance(mesh.active_type_elements_begin(type),
                                                mesh.active_type_elements_end  (type)));
}

dof_id_type MeshTools::n_non_subactive_elem_of_type_at_level(const MeshBase & mesh,
                                                             const ElemType type,
                                                             const unsigned int level)
{
  dof_id_type cnt = 0;
  // iterate over the elements of the specified type
  MeshBase::const_element_iterator el = mesh.type_elements_begin(type);
  const MeshBase::const_element_iterator end = mesh.type_elements_end(type);

  for(; el!=end; ++el)
    if( ((*el)->level() == level) && !(*el)->subactive())
      cnt++;

  return cnt;
}


unsigned int MeshTools::n_active_local_levels(const MeshBase & mesh)
{
  unsigned int nl = 0;

  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for( ; el != end_el; ++el)
    nl = std::max((*el)->level() + 1, nl);

  return nl;
}



unsigned int MeshTools::n_active_levels(const MeshBase & mesh)
{
  libmesh_parallel_only(mesh.comm());

  unsigned int nl = MeshTools::n_active_local_levels(mesh);

  MeshBase::const_element_iterator el =
    mesh.unpartitioned_elements_begin();
  const MeshBase::const_element_iterator end_el =
    mesh.unpartitioned_elements_end();

  for( ; el != end_el; ++el)
    if ((*el)->active())
      nl = std::max((*el)->level() + 1, nl);

  mesh.comm().max(nl);
  return nl;
}



unsigned int MeshTools::n_local_levels(const MeshBase & mesh)
{
  unsigned int nl = 0;

  MeshBase::const_element_iterator el = mesh.local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.local_elements_end();

  for( ; el != end_el; ++el)
    nl = std::max((*el)->level() + 1, nl);

  return nl;
}



unsigned int MeshTools::n_levels(const MeshBase & mesh)
{
  libmesh_parallel_only(mesh.comm());

  unsigned int nl = MeshTools::n_local_levels(mesh);

  MeshBase::const_element_iterator el =
    mesh.unpartitioned_elements_begin();
  const MeshBase::const_element_iterator end_el =
    mesh.unpartitioned_elements_end();

  for( ; el != end_el; ++el)
    nl = std::max((*el)->level() + 1, nl);

  mesh.comm().max(nl);

  // n_levels() is only valid and should only be called in cases where
  // the mesh is validly distributed (or serialized).  Let's run an
  // expensive test in debug mode to make sure this is such a case.
#ifdef DEBUG
  const unsigned int paranoid_nl = MeshTools::paranoid_n_levels(mesh);
  libmesh_assert_equal_to(nl, paranoid_nl);
#endif
  return nl;
}



unsigned int MeshTools::paranoid_n_levels(const MeshBase & mesh)
{
  libmesh_parallel_only(mesh.comm());

  MeshBase::const_element_iterator el =
    mesh.elements_begin();
  const MeshBase::const_element_iterator end_el =
    mesh.elements_end();

  unsigned int nl = 0;
  for( ; el != end_el; ++el)
    nl = std::max((*el)->level() + 1, nl);

  mesh.comm().max(nl);
  return nl;
}



void MeshTools::get_not_subactive_node_ids(const MeshBase & mesh,
                                           std::set<dof_id_type> & not_subactive_node_ids)
{
  MeshBase::const_element_iterator el           = mesh.elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.elements_end();
  for( ; el != end_el; ++el)
    {
      const Elem * elem = (*el);
      if(!elem->subactive())
        for (unsigned int n=0; n<elem->n_nodes(); ++n)
          not_subactive_node_ids.insert(elem->node_id(n));
    }
}



dof_id_type MeshTools::n_elem (const MeshBase::const_element_iterator & begin,
                               const MeshBase::const_element_iterator & end)
{
  return cast_int<dof_id_type>(std::distance(begin, end));
}



dof_id_type MeshTools::n_nodes (const MeshBase::const_node_iterator & begin,
                                const MeshBase::const_node_iterator & end)
{
  return cast_int<dof_id_type>(std::distance(begin, end));
}



unsigned int MeshTools::n_p_levels (const MeshBase & mesh)
{
  libmesh_parallel_only(mesh.comm());

  unsigned int max_p_level = 0;

  // first my local elements
  MeshBase::const_element_iterator
    el     = mesh.local_elements_begin(),
    end_el = mesh.local_elements_end();

  for( ; el != end_el; ++el)
    max_p_level = std::max((*el)->p_level(), max_p_level);

  // then any unpartitioned objects
  el     = mesh.unpartitioned_elements_begin();
  end_el = mesh.unpartitioned_elements_end();

  for( ; el != end_el; ++el)
    max_p_level = std::max((*el)->p_level(), max_p_level);

  mesh.comm().max(max_p_level);
  return max_p_level + 1;
}



void MeshTools::find_nodal_neighbors(const MeshBase &,
                                     const Node & node,
                                     const std::vector<std::vector<const Elem *> > & nodes_to_elem_map,
                                     std::vector<const Node *> & neighbors)
{
  // We'll refer back to the Node ID several times
  dof_id_type global_id = node.id();

  // We'll construct a std::set<const Node *> for more efficient
  // searching while finding the nodal neighbors, and return it to the
  // user in a std::vector.
  std::set<const Node *> neighbor_set;

  // Iterators to iterate through the elements that include this node
  std::vector<const Elem *>::const_iterator
    el = nodes_to_elem_map[global_id].begin(),
    end_el = nodes_to_elem_map[global_id].end();

  // Look through the elements that contain this node
  // find the local node id... then find the side that
  // node lives on in the element
  // next, look for the _other_ node on that side
  // That other node is a "nodal_neighbor"... save it
  for (; el != end_el; ++el)
    {
      // Grab an Elem pointer to use in the subsequent loop
      const Elem * elem = *el;

      // We only care about active elements...
      if (elem->active())
        {
          // Which local node number is global_id?
          unsigned local_node_number = elem->local_node(global_id);

          // Make sure it was found
          libmesh_assert_not_equal_to(local_node_number, libMesh::invalid_uint);

          // If this element has no edges, the edge-based algorithm below doesn't make sense.
          if (elem->n_edges() == 0)
            {
              switch (elem->type())
                {
                case EDGE2:
                  {
                    switch (local_node_number)
                      {
                      case 0:
                        // The other node is a nodal neighbor
                        neighbor_set.insert(elem->node_ptr(1));
                        break;

                      case 1:
                        // The other node is a nodal neighbor
                        neighbor_set.insert(elem->node_ptr(0));
                        break;

                      default:
                        libmesh_error_msg("Invalid local node number: " << local_node_number << " found." << std::endl);
                      }
                    break;
                  }

                case EDGE3:
                  {
                    switch (local_node_number)
                      {
                        // The outside nodes have node 2 as a neighbor
                      case 0:
                      case 1:
                        neighbor_set.insert(elem->node_ptr(2));
                        break;

                        // The middle node has the outer nodes as neighbors
                      case 2:
                        neighbor_set.insert(elem->node_ptr(0));
                        neighbor_set.insert(elem->node_ptr(1));
                        break;

                      default:
                        libmesh_error_msg("Invalid local node number: " << local_node_number << " found." << std::endl);
                      }
                    break;
                  }

                case EDGE4:
                  {
                    switch (local_node_number)
                      {
                      case 0:
                        // The left-middle node is a nodal neighbor
                        neighbor_set.insert(elem->node_ptr(2));
                        break;

                      case 1:
                        // The right-middle node is a nodal neighbor
                        neighbor_set.insert(elem->node_ptr(3));
                        break;

                        // The left-middle node
                      case 2:
                        neighbor_set.insert(elem->node_ptr(0));
                        neighbor_set.insert(elem->node_ptr(3));
                        break;

                        // The right-middle node
                      case 3:
                        neighbor_set.insert(elem->node_ptr(1));
                        neighbor_set.insert(elem->node_ptr(2));
                        break;

                      default:
                        libmesh_error_msg("Invalid local node number: " << local_node_number << " found." << std::endl);
                      }
                    break;
                  }

                default:
                  libmesh_error_msg("Unrecognized ElemType: " << Utility::enum_to_string(elem->type()) << std::endl);
                }
            }

          // Index of the current edge
          unsigned current_edge = 0;

          while (current_edge < elem->n_edges())
            {
              // Find the edge the node is on
              bool found_edge = false;
              for (; current_edge<elem->n_edges(); ++current_edge)
                if ( elem->is_node_on_edge(local_node_number, current_edge) )
                  {
                    found_edge = true;
                    break;
                  }

              // Did we find one?
              if (found_edge)
                {
                  const Node * node_to_save = libmesh_nullptr;

                  // Find another node in this element on this edge
                  for (unsigned other_node_this_edge = 0; other_node_this_edge<elem->n_nodes(); other_node_this_edge++)
                    if ( (elem->is_node_on_edge(other_node_this_edge, current_edge)) && // On the current edge
                         (elem->node_id(other_node_this_edge) != global_id))               // But not the original node
                      {
                        // We've found a nodal neighbor!  Save a pointer to it..
                        node_to_save = elem->node_ptr(other_node_this_edge);
                        break;
                      }

                  // Make sure we found something
                  libmesh_assert(node_to_save != libmesh_nullptr);

                  neighbor_set.insert(node_to_save);
                }

              // Keep looking for edges, node may be on more than one edge
              current_edge++;
            }
        } // if (elem->active())
    } // for

  // Assign the entries from the set to the vector.  Note: this
  // replaces any existing contents in neighbors and modifies its size
  // accordingly.
  neighbors.assign(neighbor_set.begin(), neighbor_set.end());
}



void MeshTools::find_hanging_nodes_and_parents(const MeshBase & mesh,
                                               std::map<dof_id_type, std::vector<dof_id_type> > & hanging_nodes)
{
  MeshBase::const_element_iterator it  = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end = mesh.active_local_elements_end();

  //Loop through all the elements
  for (; it != end; ++it)
    {
      //Save it off for easier access
      const Elem * elem = (*it);

      //Right now this only works for quad4's
      //libmesh_assert_equal_to (elem->type(), QUAD4);
      if(elem->type() == QUAD4)
        {
          //Loop over the sides looking for sides that have hanging nodes
          //This code is inspired by compute_proj_constraints()
          for (unsigned int s=0; s<elem->n_sides(); s++)
            {
              //If not a boundary node
              if (elem->neighbor_ptr(s) != libmesh_nullptr)
                {
                  // Get pointers to the element's neighbor.
                  const Elem * neigh = elem->neighbor_ptr(s);

                  //Is there a coarser element next to this one?
                  if (neigh->level() < elem->level())
                    {
                      const Elem * ancestor = elem;
                      while (neigh->level() < ancestor->level())
                        ancestor = ancestor->parent();
                      unsigned int s_neigh = neigh->which_neighbor_am_i(ancestor);
                      libmesh_assert_less (s_neigh, neigh->n_neighbors());

                      //Couple of helper uints...
                      unsigned int local_node1=0;
                      unsigned int local_node2=0;

                      bool found_in_neighbor = false;

                      //Find the two vertices that make up this side
                      while(!elem->is_node_on_side(local_node1++,s)) { }
                      local_node1--;

                      //Start looking for the second one with the next node
                      local_node2=local_node1+1;

                      //Find the other one
                      while(!elem->is_node_on_side(local_node2++,s)) { }
                      local_node2--;

                      //Pull out their global ids:
                      dof_id_type node1 = elem->node_id(local_node1);
                      dof_id_type node2 = elem->node_id(local_node2);

                      //Now find which node is present in the neighbor
                      //FIXME This assumes a level one rule!
                      //The _other_ one is the hanging node

                      //First look for the first one
                      //FIXME could be streamlined a bit
                      for(unsigned int n=0;n<neigh->n_sides();n++)
                        {
                          if(neigh->node_id(n) == node1)
                            found_in_neighbor=true;
                        }

                      dof_id_type hanging_node=0;

                      if(!found_in_neighbor)
                        hanging_node=node1;
                      else //If it wasn't node1 then it must be node2!
                        hanging_node=node2;

                      //Reset these for reuse
                      local_node1=0;
                      local_node2=0;

                      //Find the first node that makes up the side in the neighbor (these should be the parent nodes)
                      while(!neigh->is_node_on_side(local_node1++,s_neigh)) { }
                      local_node1--;

                      local_node2=local_node1+1;

                      //Find the second node...
                      while(!neigh->is_node_on_side(local_node2++,s_neigh)) { }
                      local_node2--;

                      //Save them if we haven't already found the parents for this one
                      if(hanging_nodes[hanging_node].size()<2)
                        {
                          hanging_nodes[hanging_node].push_back(neigh->node_id(local_node1));
                          hanging_nodes[hanging_node].push_back(neigh->node_id(local_node2));
                        }
                    }
                }
            }
        }
    }
}



#ifdef DEBUG
void MeshTools::libmesh_assert_equal_n_systems (const MeshBase & mesh)
{
  MeshBase::const_element_iterator el =
    mesh.elements_begin();
  const MeshBase::const_element_iterator el_end =
    mesh.elements_end();
  if (el == el_end)
    return;

  const unsigned int n_sys = (*el)->n_systems();

  for (; el != el_end; ++el)
    {
      const Elem * elem = *el;
      libmesh_assert_equal_to (elem->n_systems(), n_sys);
    }

  MeshBase::const_node_iterator node_it =
    mesh.nodes_begin();
  const MeshBase::const_node_iterator node_end =
    mesh.nodes_end();

  if (node_it == node_end)
    return;

  for (; node_it != node_end; ++node_it)
    {
      const Node * node = *node_it;
      libmesh_assert_equal_to (node->n_systems(), n_sys);
    }
}



#ifdef LIBMESH_ENABLE_AMR
void MeshTools::libmesh_assert_old_dof_objects (const MeshBase & mesh)
{
  MeshBase::const_element_iterator el =
    mesh.elements_begin();
  const MeshBase::const_element_iterator el_end =
    mesh.elements_end();

  for (; el != el_end; ++el)
    {
      const Elem * elem = *el;

      if (elem->refinement_flag() == Elem::JUST_REFINED ||
          elem->refinement_flag() == Elem::INACTIVE)
        continue;

      if (elem->has_dofs())
        libmesh_assert(elem->old_dof_object);

      for (unsigned int n=0; n != elem->n_nodes(); ++n)
        {
          const Node & node = elem->node_ref(n);
          if (node.has_dofs())
            libmesh_assert(node.old_dof_object);
        }
    }
}
#else
void MeshTools::libmesh_assert_old_dof_objects (const MeshBase &) {}
#endif // LIBMESH_ENABLE_AMR



void MeshTools::libmesh_assert_valid_node_pointers(const MeshBase & mesh)
{
  const MeshBase::const_element_iterator el_end =
    mesh.elements_end();
  for (MeshBase::const_element_iterator el =
         mesh.elements_begin(); el != el_end; ++el)
    {
      const Elem * elem = *el;
      libmesh_assert (elem);
      while (elem)
        {
          elem->libmesh_assert_valid_node_pointers();
          for (unsigned int n=0; n != elem->n_neighbors(); ++n)
            if (elem->neighbor_ptr(n) &&
                elem->neighbor_ptr(n) != remote_elem)
              elem->neighbor_ptr(n)->libmesh_assert_valid_node_pointers();

          libmesh_assert_not_equal_to (elem->parent(), remote_elem);
          elem = elem->parent();
        }
    }
}


void MeshTools::libmesh_assert_valid_remote_elems(const MeshBase & mesh)
{
  const MeshBase::const_element_iterator el_end =
    mesh.local_elements_end();
  for (MeshBase::const_element_iterator el =
         mesh.local_elements_begin(); el != el_end; ++el)
    {
      const Elem * elem = *el;
      libmesh_assert (elem);

      // We currently don't allow active_local_elements to have
      // remote_elem neighbors
      if (elem->active())
        for (unsigned int n=0; n != elem->n_neighbors(); ++n)
          libmesh_assert_not_equal_to (elem->neighbor_ptr(n), remote_elem);

#ifdef LIBMESH_ENABLE_AMR
      const Elem * parent = elem->parent();
      if (parent)
        libmesh_assert_not_equal_to (parent, remote_elem);

      // We can only be strict about active elements' subactive
      // children
      if (elem->active() && elem->has_children())
        for (unsigned int c=0; c != elem->n_children(); ++c)
          libmesh_assert_not_equal_to (elem->child_ptr(c), remote_elem);
#endif
    }
}


void MeshTools::libmesh_assert_no_links_to_elem(const MeshBase & mesh,
                                                const Elem * bad_elem)
{
  const MeshBase::const_element_iterator el_end =
    mesh.elements_end();
  for (MeshBase::const_element_iterator el =
         mesh.elements_begin(); el != el_end; ++el)
    {
      const Elem * elem = *el;
      libmesh_assert (elem);
      libmesh_assert_not_equal_to (elem->parent(), bad_elem);
      for (unsigned int n=0; n != elem->n_neighbors(); ++n)
        libmesh_assert_not_equal_to (elem->neighbor_ptr(n), bad_elem);
#ifdef LIBMESH_ENABLE_AMR
      if (elem->has_children())
        for (unsigned int c=0; c != elem->n_children(); ++c)
          libmesh_assert_not_equal_to (elem->child_ptr(c), bad_elem);
#endif
    }
}



void MeshTools::libmesh_assert_valid_elem_ids(const MeshBase & mesh)
{
  processor_id_type lastprocid = 0;
  dof_id_type lastelemid = 0;

  const MeshBase::const_element_iterator el_end =
    mesh.active_elements_end();
  for (MeshBase::const_element_iterator el =
         mesh.active_elements_begin(); el != el_end; ++el)
    {
      const Elem * elem = *el;
      libmesh_assert (elem);
      processor_id_type elemprocid = elem->processor_id();
      dof_id_type elemid = elem->id();

      libmesh_assert_greater_equal (elemid, lastelemid);
      libmesh_assert_greater_equal (elemprocid, lastprocid);

      lastelemid = elemid;
      lastprocid = elemprocid;
    }
}



void MeshTools::libmesh_assert_valid_amr_elem_ids(const MeshBase & mesh)
{
  const MeshBase::const_element_iterator el_end =
    mesh.elements_end();
  for (MeshBase::const_element_iterator el =
         mesh.elements_begin(); el != el_end; ++el)
    {
      const Elem * elem = *el;
      libmesh_assert (elem);

      const Elem * parent = elem->parent();

      if (parent)
        {
          libmesh_assert_greater_equal (elem->id(), parent->id());
          libmesh_assert_greater_equal (elem->processor_id(), parent->processor_id());
        }
    }
}



void MeshTools::libmesh_assert_valid_amr_interior_parents(const MeshBase & mesh)
{
  const MeshBase::const_element_iterator el_end =
    mesh.elements_end();
  for (MeshBase::const_element_iterator el =
         mesh.elements_begin(); el != el_end; ++el)
    {
      const Elem * elem = *el;
      libmesh_assert (elem);

      // We can skip to the next element if we're full-dimension
      // and therefore don't have any interior parents
      if (elem->dim() >= LIBMESH_DIM)
        continue;

      const Elem * ip = elem->interior_parent();

      const Elem * parent = elem->parent();

      if (ip && (ip != remote_elem) && parent)
        {
          libmesh_assert_equal_to (ip->top_parent(),
                                   elem->top_parent()->interior_parent());

          if (ip->level() == elem->level())
            libmesh_assert_equal_to (ip->parent(),
                                     parent->interior_parent());
          else
            {
              libmesh_assert_less (ip->level(), elem->level());
              libmesh_assert_equal_to (ip, parent->interior_parent());
            }
        }
    }
}



void MeshTools::libmesh_assert_connected_nodes (const MeshBase & mesh)
{
  std::set<const Node *> used_nodes;

  const MeshBase::const_element_iterator el_end =
    mesh.elements_end();
  for (MeshBase::const_element_iterator el =
         mesh.elements_begin(); el != el_end; ++el)
    {
      const Elem * elem = *el;
      libmesh_assert (elem);

      for (unsigned int n=0; n<elem->n_nodes(); ++n)
        used_nodes.insert(elem->node_ptr(n));
    }

  const MeshBase::const_node_iterator node_end = mesh.nodes_end();

  for (MeshBase::const_node_iterator node_it = mesh.nodes_begin();
       node_it != node_end; ++node_it)
    {
      Node * node = *node_it;
      libmesh_assert(node);
      libmesh_assert(used_nodes.count(node));
    }
}



namespace MeshTools {

void libmesh_assert_valid_boundary_ids(const MeshBase & mesh)
{
  if (mesh.n_processors() == 1)
    return;

  libmesh_parallel_only(mesh.comm());

  const BoundaryInfo & boundary_info = mesh.get_boundary_info();

  dof_id_type pmax_elem_id = mesh.max_elem_id();
  mesh.comm().max(pmax_elem_id);

  for (dof_id_type i=0; i != pmax_elem_id; ++i)
    {
      const Elem *elem = mesh.query_elem_ptr(i);
      unsigned int n_nodes = elem ? elem->n_nodes() : 0;
      unsigned int n_edges = elem ? elem->n_edges() : 0;
      unsigned int n_sides = elem ? elem->n_sides() : 0;
      libmesh_assert(mesh.comm().semiverify
                     (elem ? &n_nodes : libmesh_nullptr));
      libmesh_assert(mesh.comm().semiverify
                     (elem ? &n_edges : libmesh_nullptr));
      libmesh_assert(mesh.comm().semiverify
                     (elem ? &n_sides : libmesh_nullptr));
      mesh.comm().max(n_nodes);
      mesh.comm().max(n_edges);
      mesh.comm().max(n_sides);
      for (unsigned int n=0; n != n_nodes; ++n)
        {
          std::vector<boundary_id_type> bcids;
          if (elem)
            {
              boundary_info.boundary_ids(elem->node_ptr(n), bcids);

              // Ordering of boundary ids shouldn't matter
              std::sort(bcids.begin(), bcids.end());
            }
          libmesh_assert(mesh.comm().semiverify
                         (elem ? &bcids : libmesh_nullptr));
        }

      for (unsigned short e=0; e != n_edges; ++e)
        {
          std::vector<boundary_id_type> bcids;

          if (elem)
            {
              boundary_info.edge_boundary_ids(elem, e, bcids);

              // Ordering of boundary ids shouldn't matter
              std::sort(bcids.begin(), bcids.end());
            }

          libmesh_assert(mesh.comm().semiverify
                         (elem ? &bcids : libmesh_nullptr));

          if (elem)
            {
              boundary_info.raw_edge_boundary_ids(elem, e, bcids);

              // Ordering of boundary ids shouldn't matter
              std::sort(bcids.begin(), bcids.end());
            }

          libmesh_assert(mesh.comm().semiverify
                         (elem ? &bcids : libmesh_nullptr));
        }

      for (unsigned short s=0; s != n_sides; ++s)
        {
          std::vector<boundary_id_type> bcids;

          if (elem)
            {
              boundary_info.boundary_ids(elem, s, bcids);

              // Ordering of boundary ids shouldn't matter
              std::sort(bcids.begin(), bcids.end());
            }

          libmesh_assert(mesh.comm().semiverify
                         (elem ? &bcids : libmesh_nullptr));

          if (elem)
            {
              boundary_info.raw_boundary_ids(elem, s, bcids);

              // Ordering of boundary ids shouldn't matter
              std::sort(bcids.begin(), bcids.end());
            }

          libmesh_assert(mesh.comm().semiverify
                         (elem ? &bcids : libmesh_nullptr));
        }

      for (unsigned short sf=0; sf != 2; ++sf)
        {
          std::vector<boundary_id_type> bcids;

          if (elem)
            {
              boundary_info.shellface_boundary_ids(elem, sf, bcids);

              // Ordering of boundary ids shouldn't matter
              std::sort(bcids.begin(), bcids.end());
            }

          libmesh_assert(mesh.comm().semiverify
                         (elem ? &bcids : libmesh_nullptr));

          if (elem)
            {
              boundary_info.raw_shellface_boundary_ids(elem, sf, bcids);

              // Ordering of boundary ids shouldn't matter
              std::sort(bcids.begin(), bcids.end());
            }

          libmesh_assert(mesh.comm().semiverify
                         (elem ? &bcids : libmesh_nullptr));
        }
    }
}

void libmesh_assert_valid_dof_ids(const MeshBase & mesh)
{
  if (mesh.n_processors() == 1)
    return;

  libmesh_parallel_only(mesh.comm());

  dof_id_type pmax_elem_id = mesh.max_elem_id();
  mesh.comm().max(pmax_elem_id);

  for (dof_id_type i=0; i != pmax_elem_id; ++i)
    assert_semiverify_dofobj(mesh.comm(),
                             mesh.query_elem_ptr(i));

  dof_id_type pmax_node_id = mesh.max_node_id();
  mesh.comm().max(pmax_node_id);

  for (dof_id_type i=0; i != pmax_node_id; ++i)
    assert_semiverify_dofobj(mesh.comm(),
                             mesh.query_node_ptr(i));
}


#ifdef LIBMESH_ENABLE_UNIQUE_ID
void libmesh_assert_valid_unique_ids(const MeshBase &mesh)
{
  libmesh_parallel_only(mesh.comm());

  dof_id_type pmax_elem_id = mesh.max_elem_id();
  mesh.comm().max(pmax_elem_id);

  for (dof_id_type i=0; i != pmax_elem_id; ++i)
    {
      const Elem *elem = mesh.query_elem_ptr(i);
      const unique_id_type unique_id = elem ? elem->unique_id() : 0;
      const unique_id_type * uid_ptr = elem ? &unique_id : libmesh_nullptr;
      libmesh_assert(mesh.comm().semiverify(uid_ptr));
    }

  dof_id_type pmax_node_id = mesh.max_node_id();
  mesh.comm().max(pmax_node_id);

  for (dof_id_type i=0; i != pmax_node_id; ++i)
    {
      const Node *node = mesh.query_node_ptr(i);
      const unique_id_type unique_id = node ? node->unique_id() : 0;
      const unique_id_type * uid_ptr = node ? &unique_id : libmesh_nullptr;
      libmesh_assert(mesh.comm().semiverify(uid_ptr));
    }
}
#endif

template <>
void libmesh_assert_topology_consistent_procids<Elem>(const MeshBase & mesh)
{
  // If we're adaptively refining, check processor ids for consistency
  // between parents and children.
#ifdef LIBMESH_ENABLE_AMR

  // Ancestor elements we won't worry about, but subactive and active
  // elements ought to have parents with consistent processor ids

  const MeshBase::const_element_iterator el_end =
    mesh.elements_end();
  for (MeshBase::const_element_iterator el =
         mesh.elements_begin(); el != el_end; ++el)
    {
      const Elem * elem = *el;
      libmesh_assert(elem);

      if (!elem->active() && !elem->subactive())
        continue;

      const Elem * parent = elem->parent();

      if (parent)
        {
          libmesh_assert(parent->has_children());
          processor_id_type parent_procid = parent->processor_id();
          bool matching_child_id = false;
          for (unsigned int c = 0; c != parent->n_children(); ++c)
            {
              const Elem * child = parent->child_ptr(c);
              libmesh_assert(child);

              // If we've got a remote_elem then we don't know whether
              // it's responsible for the parent's processor id; all
              // we can do is assume it is and let its processor fail
              // an assert if there's something wrong.
              if (child == remote_elem ||
                  child->processor_id() == parent_procid)
                matching_child_id = true;
            }
          libmesh_assert(matching_child_id);
        }
    }
#endif
}



template <>
void libmesh_assert_parallel_consistent_procids<Elem>(const MeshBase & mesh)
{
  if (mesh.n_processors() == 1)
    return;

  libmesh_parallel_only(mesh.comm());

  // We want this test to be valid even when called even after nodes
  // have been added asynchonously but before they're renumbered
  dof_id_type parallel_max_elem_id = mesh.max_elem_id();
  mesh.comm().max(parallel_max_elem_id);

  // Check processor ids for consistency between processors

  for (dof_id_type i=0; i != parallel_max_elem_id; ++i)
    {
      const Elem * elem = mesh.query_elem_ptr(i);

      processor_id_type min_id =
        elem ? elem->processor_id() :
        std::numeric_limits<processor_id_type>::max();
      mesh.comm().min(min_id);

      processor_id_type max_id =
        elem ? elem->processor_id() :
        std::numeric_limits<processor_id_type>::min();
      mesh.comm().max(max_id);

      if (elem)
        {
          libmesh_assert_equal_to (min_id, elem->processor_id());
          libmesh_assert_equal_to (max_id, elem->processor_id());
        }

      if (min_id == mesh.processor_id())
        libmesh_assert(elem);
    }
}



template <>
void libmesh_assert_topology_consistent_procids<Node>(const MeshBase & mesh)
{
  if (mesh.n_processors() == 1)
    return;

  libmesh_parallel_only(mesh.comm());

  // We want this test to be valid even when called even after nodes
  // have been added asynchonously but before they're renumbered
  dof_id_type parallel_max_node_id = mesh.max_node_id();
  mesh.comm().max(parallel_max_node_id);

  std::vector<bool> node_touched_by_me(parallel_max_node_id, false);

  const MeshBase::const_element_iterator el_end =
    mesh.local_elements_end();
  for (MeshBase::const_element_iterator el =
         mesh.local_elements_begin(); el != el_end; ++el)
    {
      const Elem * elem = *el;
      libmesh_assert (elem);

      for (unsigned int i=0; i != elem->n_nodes(); ++i)
        {
          const Node & node = elem->node_ref(i);
          dof_id_type nodeid = node.id();
          node_touched_by_me[nodeid] = true;
        }
    }
  std::vector<bool> node_touched_by_anyone(node_touched_by_me);
  mesh.comm().max(node_touched_by_anyone);

  const MeshBase::const_node_iterator nd_end = mesh.local_nodes_end();
  for (MeshBase::const_node_iterator nd = mesh.local_nodes_begin();
       nd != nd_end; ++nd)
    {
      const Node * node = *nd;
      libmesh_assert(node);

      dof_id_type nodeid = node->id();
      libmesh_assert(!node_touched_by_anyone[nodeid] ||
                     node_touched_by_me[nodeid]);
    }
}



template <>
void libmesh_assert_parallel_consistent_procids<Node>(const MeshBase & mesh)
{
  if (mesh.n_processors() == 1)
    return;

  libmesh_parallel_only(mesh.comm());

  // We want this test to be valid even when called even after nodes
  // have been added asynchonously but before they're renumbered
  dof_id_type parallel_max_node_id = mesh.max_node_id();
  mesh.comm().max(parallel_max_node_id);

  std::vector<bool> node_touched_by_anyone(parallel_max_node_id, false);

  const MeshBase::const_element_iterator el_end =
    mesh.local_elements_end();
  for (MeshBase::const_element_iterator el =
         mesh.local_elements_begin(); el != el_end; ++el)
    {
      const Elem * elem = *el;
      libmesh_assert (elem);

      for (unsigned int i=0; i != elem->n_nodes(); ++i)
        {
          const Node & node = elem->node_ref(i);
          dof_id_type nodeid = node.id();
          node_touched_by_anyone[nodeid] = true;
        }
    }
  mesh.comm().max(node_touched_by_anyone);

  // Check processor ids for consistency between processors
  // on any node an element touches
  for (dof_id_type i=0; i != parallel_max_node_id; ++i)
    {
      if (!node_touched_by_anyone[i])
        continue;

      const Node * node = mesh.query_node_ptr(i);

      processor_id_type min_id =
        node ? node->processor_id() :
        std::numeric_limits<processor_id_type>::max();
      mesh.comm().min(min_id);

      processor_id_type max_id =
        node ? node->processor_id() :
        std::numeric_limits<processor_id_type>::min();
      mesh.comm().max(max_id);

      if (node)
        {
          libmesh_assert_equal_to (min_id, node->processor_id());
          libmesh_assert_equal_to (max_id, node->processor_id());
        }

      if (min_id == mesh.processor_id())
        libmesh_assert(node);
    }
}



} // namespace MeshTools



#ifdef LIBMESH_ENABLE_AMR
void MeshTools::libmesh_assert_valid_refinement_flags(const MeshBase & mesh)
{
  libmesh_parallel_only(mesh.comm());
  if (mesh.n_processors() == 1)
    return;

  std::vector<unsigned char> my_elem_h_state(mesh.max_elem_id(), 255);
  std::vector<unsigned char> my_elem_p_state(mesh.max_elem_id(), 255);

  const MeshBase::const_element_iterator el_end =
    mesh.elements_end();
  for (MeshBase::const_element_iterator el =
         mesh.elements_begin(); el != el_end; ++el)
    {
      const Elem * elem = *el;
      libmesh_assert (elem);
      dof_id_type elemid = elem->id();

      my_elem_h_state[elemid] =
        static_cast<unsigned char>(elem->refinement_flag());

      my_elem_p_state[elemid] =
        static_cast<unsigned char>(elem->p_refinement_flag());
    }
  std::vector<unsigned char> min_elem_h_state(my_elem_h_state);
  mesh.comm().min(min_elem_h_state);

  std::vector<unsigned char> min_elem_p_state(my_elem_p_state);
  mesh.comm().min(min_elem_p_state);

  for (dof_id_type i=0; i!= mesh.max_elem_id(); ++i)
    {
      libmesh_assert(my_elem_h_state[i] == 255 ||
                     my_elem_h_state[i] == min_elem_h_state[i]);
      libmesh_assert(my_elem_p_state[i] == 255 ||
                     my_elem_p_state[i] == min_elem_p_state[i]);
    }
}
#else
void MeshTools::libmesh_assert_valid_refinement_flags(const MeshBase &)
{
}
#endif // LIBMESH_ENABLE_AMR



#ifdef LIBMESH_ENABLE_AMR
void MeshTools::libmesh_assert_valid_refinement_tree(const MeshBase & mesh)
{
  const MeshBase::const_element_iterator el_end =
    mesh.elements_end();
  for (MeshBase::const_element_iterator el =
         mesh.elements_begin(); el != el_end; ++el)
    {
      const Elem * elem = *el;
      libmesh_assert(elem);
      if (elem->has_children())
        for (unsigned int n=0; n != elem->n_children(); ++n)
          {
            libmesh_assert(elem->child_ptr(n));
            if (elem->child_ptr(n) != remote_elem)
              libmesh_assert_equal_to (elem->child_ptr(n)->parent(), elem);
          }
      if (elem->active())
        {
          libmesh_assert(!elem->ancestor());
          libmesh_assert(!elem->subactive());
        }
      else if (elem->ancestor())
        {
          libmesh_assert(!elem->subactive());
        }
      else
        libmesh_assert(elem->subactive());

      if (elem->p_refinement_flag() == Elem::JUST_REFINED)
        libmesh_assert_greater(elem->p_level(), 0);
    }
}
#else
void MeshTools::libmesh_assert_valid_refinement_tree(const MeshBase &)
{
}
#endif // LIBMESH_ENABLE_AMR



void MeshTools::libmesh_assert_valid_neighbors(const MeshBase & mesh,
                                               bool assert_valid_remote_elems)
{
  const MeshBase::const_element_iterator el_end = mesh.elements_end();
  for (MeshBase::const_element_iterator el = mesh.elements_begin();
       el != el_end; ++el)
    {
      const Elem * elem = *el;
      libmesh_assert (elem);
      elem->libmesh_assert_valid_neighbors();
    }

  if (mesh.n_processors() == 1)
    return;

  libmesh_parallel_only(mesh.comm());

  for (dof_id_type i=0; i != mesh.max_elem_id(); ++i)
    {
      const Elem * elem = mesh.query_elem_ptr(i);

      const unsigned int my_n_neigh = elem ? elem->n_neighbors() : 0;
      unsigned int n_neigh = my_n_neigh;
      mesh.comm().max(n_neigh);
      if (elem)
        libmesh_assert_equal_to (my_n_neigh, n_neigh);

      for (unsigned int n = 0; n != n_neigh; ++n)
        {
          dof_id_type my_neighbor = DofObject::invalid_id;
          dof_id_type * p_my_neighbor = libmesh_nullptr;

          // If we have a non-remote_elem neighbor link, then we can
          // verify it.
          if (elem && elem->neighbor_ptr(n) != remote_elem)
            {
              p_my_neighbor = &my_neighbor;
              if (elem->neighbor_ptr(n))
                my_neighbor = elem->neighbor_ptr(n)->id();

              // But wait - if we haven't set remote_elem links yet then
              // some NULL links on ghost elements might be
              // future-remote_elem links, so we can't verify those.
              if (!assert_valid_remote_elems &&
                  !elem->neighbor_ptr(n) &&
                  elem->processor_id() != mesh.processor_id())
                p_my_neighbor = libmesh_nullptr;
            }
          libmesh_assert(mesh.comm().semiverify(p_my_neighbor));
        }
    }
}



#endif // DEBUG



// Functors for correct_node_proc_ids
namespace {

typedef LIBMESH_BEST_UNORDERED_MAP<dof_id_type, processor_id_type> proc_id_map_type;

struct SyncProcIdsFromMap
{
  typedef processor_id_type datum;

  SyncProcIdsFromMap(const proc_id_map_type & _map,
                     MeshBase & _mesh) :
    new_proc_ids(_map), mesh(_mesh) {}

  const proc_id_map_type & new_proc_ids;

  MeshBase & mesh;

  // ------------------------------------------------------------
  void gather_data (const std::vector<dof_id_type> & ids,
                    std::vector<datum> & data)
  {
    // Find the new processor id of each requested node
    data.resize(ids.size());

    for (std::size_t i=0; i != ids.size(); ++i)
      {
        const dof_id_type id = ids[i];
        const proc_id_map_type::const_iterator it = new_proc_ids.find(id);

        // Return the node's new processor id if it has one, or its
        // old processor id if not.
        if (it != new_proc_ids.end())
          data[i] = it->second;
        else
          {
            // We'd better find every node we're asked for
            const Node & node = mesh.node_ref(id);
            data[i] = node.processor_id();
          }
      }
  }

  // ------------------------------------------------------------
  void act_on_data (const std::vector<dof_id_type> & ids,
                    std::vector<datum> proc_ids)
  {
    // Set the node processor ids we've now been informed of
    for (std::size_t i=0; i != ids.size(); ++i)
      {
        Node & node = mesh.node_ref(ids[i]);
        node.processor_id() = proc_ids[i];
      }
  }
};
}



void MeshTools::correct_node_proc_ids (MeshBase & mesh)
{
  // This function must be run on all processors at once
  libmesh_parallel_only(mesh.comm());

  // Fix all nodes' processor ids.  Coarsening may have left us with
  // nodes which are no longer touched by any elements of the same
  // processor id, and for DofMap to work we need to fix that.

  // This is harder now that libMesh no longer requires a distributed
  // mesh to ghost all nodal neighbors: it is possible for two active
  // elements on two different processors to share the same node in
  // such a way that neither processor knows the others' element
  // exists!

  // We require all processors to agree on nodal processor ids before
  // going into this algorithm.
#ifdef DEBUG
  MeshTools::libmesh_assert_parallel_consistent_procids<Node>(mesh);
#endif

  // We build up a set of compatible processor ids for each node
  proc_id_map_type new_proc_ids;

  MeshBase::element_iterator       e_it  = mesh.active_elements_begin();
  const MeshBase::element_iterator e_end = mesh.active_elements_end();
  for (; e_it != e_end; ++e_it)
    {
      Elem * elem = *e_it;
      processor_id_type pid = elem->processor_id();

      for (unsigned int n=0; n != elem->n_nodes(); ++n)
        {
          const Node & node = elem->node_ref(n);
          const dof_id_type id = node.id();
          const proc_id_map_type::iterator it = new_proc_ids.find(id);
          if (it == new_proc_ids.end())
            new_proc_ids.insert(std::make_pair(id,pid));
          else
            it->second = std::min(it->second, pid);
        }
    }

  // Sort the new pids to push to each processor
  std::vector<std::vector<std::pair<dof_id_type, processor_id_type> > >
    ids_to_push(mesh.n_processors());

  for (MeshBase::const_node_iterator n_it = mesh.nodes_begin(),
       n_end = mesh.nodes_end(); n_it != n_end; ++n_it)
    {
      const Node *node = *n_it;
      const dof_id_type id = node->id();
      const proc_id_map_type::iterator it = new_proc_ids.find(id);
      if (it == new_proc_ids.end())
        continue;
      const processor_id_type pid = it->second;
      ids_to_push[node->processor_id()].push_back(std::make_pair(id, pid));
    }

  // Push using non-blocking I/O
  std::vector<Parallel::Request> push_requests(mesh.n_processors());

  for (processor_id_type p=1; p != mesh.n_processors(); ++p)
    {
      const processor_id_type procup =
        cast_int<processor_id_type>
        ((mesh.comm().rank() + p) % mesh.comm().size());

      mesh.comm().send(procup, ids_to_push[procup], push_requests[procup]);
    }

  for (processor_id_type p=0; p != mesh.n_processors(); ++p)
    {
      const processor_id_type procdown =
        cast_int<processor_id_type>
        ((mesh.comm().size() + mesh.comm().rank() - p) %
         mesh.comm().size());

      std::vector<std::pair<dof_id_type, processor_id_type> >
        ids_to_pull;

      if (p)
        mesh.comm().receive(procdown, ids_to_pull);
      else
        ids_to_pull.swap(ids_to_push[procdown]);

      std::vector<std::pair<dof_id_type, processor_id_type> >::iterator
        pulled_ids_it = ids_to_pull.begin(),
        pulled_ids_end = ids_to_pull.end();
      for (; pulled_ids_it != pulled_ids_end; ++pulled_ids_it)
        {
          const dof_id_type id = pulled_ids_it->first;
          const processor_id_type pid = pulled_ids_it->second;
          const proc_id_map_type::iterator it = new_proc_ids.find(id);
          if (it == new_proc_ids.end())
            new_proc_ids.insert(std::make_pair(id,pid));
          else
            it->second = std::min(it->second, pid);
        }
    }

  // Now new_proc_ids is correct for every node we used to own.  Let's
  // ask every other processor about the nodes they used to own.  But
  // first we'll need to keep track of which nodes we used to own,
  // lest we get them confused with nodes we newly own.
  LIBMESH_BEST_UNORDERED_SET<Node *> ex_local_nodes;
  for (MeshBase::node_iterator n_it = mesh.local_nodes_begin(),
       n_end = mesh.local_nodes_end(); n_it != n_end; ++n_it)
    {
      Node *node = *n_it;
      const proc_id_map_type::iterator it = new_proc_ids.find(node->id());
      if (it != new_proc_ids.end() && it->second != mesh.processor_id())
        ex_local_nodes.insert(node);
    }

  // Let's finish with previous I/O before we start more.
  Parallel::wait(push_requests);

  SyncProcIdsFromMap sync(new_proc_ids, mesh);
  Parallel::sync_dofobject_data_by_id
    (mesh.comm(), mesh.nodes_begin(), mesh.nodes_end(), sync);

  // And finally let's update the nodes we used to own.
  for (LIBMESH_BEST_UNORDERED_SET<Node *>::iterator n_it = ex_local_nodes.begin(),
       n_end = ex_local_nodes.end(); n_it != n_end; ++n_it)
    {
      Node *node = *n_it;
      const dof_id_type id = node->id();
      const proc_id_map_type::iterator it = new_proc_ids.find(id);
      libmesh_assert(it != new_proc_ids.end());
      node->processor_id() = it->second;
    }

  // We should still have consistent nodal processor ids coming out of
  // this algorithm.
#ifdef DEBUG
  MeshTools::libmesh_assert_valid_procids<Node>(mesh);
#endif
}



void MeshTools::Private::globally_renumber_nodes_and_elements (MeshBase & mesh)
{
  MeshCommunication().assign_global_indices(mesh);
}

} // namespace libMesh
