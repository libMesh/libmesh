// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_POINT_LOCATOR_TREE_IMPL_H
#define LIBMESH_POINT_LOCATOR_TREE_IMPL_H

// C++ includes

// Local Includes
#include "libmesh/elem.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/point_locator_tree.h"
#include "libmesh/tree.h"

namespace libMesh
{



//------------------------------------------------------------------
// PointLocator methods
template <typename RealType>
PointLocatorTreeTempl<RealType>::PointLocatorTreeTempl (const MeshBase & mesh,
                                    const PointLocatorBase * master) :
  PointLocatorBase (mesh,master),
  _tree            (nullptr),
  _element         (nullptr),
  _out_of_mesh_mode(false),
  _target_bin_size (200),
  _build_type(Trees::NODES)
{
  this->init(_build_type);
}



template <typename RealType>
PointLocatorTreeTempl<RealType>::PointLocatorTreeTempl (const MeshBase & mesh,
                                    const Trees::BuildType build_type,
                                    const PointLocatorBase * master) :
  PointLocatorBase (mesh,master),
  _tree            (nullptr),
  _element         (nullptr),
  _out_of_mesh_mode(false),
  _target_bin_size (200),
  _build_type(build_type)
{
  this->init(_build_type);
}



template <typename RealType>
PointLocatorTreeTempl<RealType>::~PointLocatorTreeTempl ()
{
  this->clear ();
}



template <typename RealType>
void PointLocatorTreeTempl<RealType>::clear ()
{
  // only delete the tree when we are the master
  if (this->_tree != nullptr)
    {
      if (this->_master == nullptr)
        // we own the tree
        delete this->_tree;
      else
        // someone else owns and therefore deletes the tree
        this->_tree = nullptr;

      // make sure operator () throws an assertion
      this->_initialized = false;
    }
}



template <typename RealType>
void PointLocatorTreeTempl<RealType>::init()
{
  this->init(_build_type);
}



template <typename RealType>
void PointLocatorTreeTempl<RealType>::init (Trees::BuildType build_type)
{
  libmesh_assert (!this->_tree);

  if (this->_initialized)
    {
      // Warn that we are already initialized
      libMesh::err << "Warning: PointLocatorTree already initialized!  Will ignore this call..." << std::endl;

      // Further warn if we try to init() again with a different build_type
      if (_build_type != build_type)
        {
          libMesh::err << "Warning: PointLocatorTree is using build_type = " << _build_type << ".\n"
                       << "Your requested build_type, " << build_type << " will not be used!" << std::endl;
        }
    }

  else
    {
      // Let the requested build_type override the _build_type we were
      // constructed with.  This is no big deal since we have not been
      // initialized before.
      _build_type = build_type;

      if (this->_master == nullptr)
        {
          LOG_SCOPE("init(no master)", "PointLocatorTree");

          if (this->_mesh.mesh_dimension() == 3)
            _tree = new Trees::OctTreeTempl<RealType> (this->_mesh, get_target_bin_size(), _build_type);
          else
            {
              // A 1D/2D mesh in 3D space needs special consideration.
              // If the mesh is planar XY, we want to build a QuadTree
              // to search efficiently.  If the mesh is truly a manifold,
              // then we need an octree
#if LIBMESH_DIM > 2
              bool is_planar_xy = false;

              // Build the bounding box for the mesh.  If the delta-z bound is
              // negligibly small then we can use a quadtree.
              {
                BoundingBox bbox = MeshTools::create_bounding_box(this->_mesh);

                const RealType
                  Dx = bbox.second(0) - bbox.first(0),
                  Dz = bbox.second(2) - bbox.first(2);

                // In order to satisfy is_planar_xy the mesh should be planar and should
                // also be in the z=0 plane, since otherwise it is incorrect to use a
                // QuadTree since QuadTrees assume z=0.
                if ( (std::abs(Dz/(Dx + 1.e-20)) < 1e-10) && (std::abs(bbox.second(2)) < 1.e-10) )
                  is_planar_xy = true;
              }

              if (!is_planar_xy)
                _tree = new Trees::OctTreeTempl<RealType> (this->_mesh, get_target_bin_size(), _build_type);
              else
#endif
#if LIBMESH_DIM > 1
                _tree = new Trees::QuadTreeTempl<RealType> (this->_mesh, get_target_bin_size(), _build_type);
#else
              _tree = new Trees::BinaryTreeTempl<RealType> (this->_mesh, get_target_bin_size(), _build_type);
#endif
            }
        }

      else
        {
          // We are _not_ the master.  Let our Tree point to
          // the master's tree.  But for this we first transform
          // the master in a state for which we are friends.
          // And make sure the master has a tree!
          const PointLocatorTree * my_master =
            cast_ptr<const PointLocatorTree *>(this->_master);

          if (my_master->initialized())
            this->_tree = my_master->_tree;
          else
            libmesh_error_msg("ERROR: Initialize master first, then servants!");
        }

      // Not all PointLocators may own a tree, but all of them
      // use their own element pointer.  Let the element pointer
      // be unique for every interpolator.
      // Suppose the interpolators are used concurrently
      // at different locations in the mesh, then it makes quite
      // sense to have unique start elements.
      this->_element = nullptr;
    }

  // ready for take-off
  this->_initialized = true;
}



template <typename RealType>
const ElemTempl<RealType> * PointLocatorTreeTempl<RealType>::operator() (const Point & p,
                                           const std::set<subdomain_id_type> * allowed_subdomains) const
{
  libmesh_assert (this->_initialized);

  LOG_SCOPE("operator()", "PointLocatorTree");

  // If we're provided with an allowed_subdomains list and have a cached element, make sure it complies
  if (allowed_subdomains && this->_element && !allowed_subdomains->count(this->_element->subdomain_id())) this->_element = nullptr;

  // First check the element from last time before asking the tree
  if (this->_element==nullptr || !(this->_element->contains_point(p)))
    {
      // ask the tree
      this->_element = this->_tree->find_element (p,allowed_subdomains);

      if (this->_element == nullptr)
        {
          // If we haven't found the element, we may want to do a linear
          // search using a tolerance.
          if (this->_use_close_to_point_tol)
            {
              if (this->_verbose)
                {
                  libMesh::out << "Performing linear search using close-to-point tolerance "
                               << this->_close_to_point_tol
                               << std::endl;
                }

              this->_element =
                this->perform_linear_search(p,
                                            allowed_subdomains,
                                            /*use_close_to_point*/ true,
                                            this->_close_to_point_tol);

              return this->_element;
            }

          // No element seems to contain this point.  In theory, our
          // tree now correctly handles curved elements.  In
          // out-of-mesh mode this is sometimes expected, and we can
          // just return nullptr without searching further.  Out of
          // out-of-mesh mode, something must have gone wrong.
          libmesh_assert_equal_to (_out_of_mesh_mode, true);

          return this->_element;
        }
    }

  // If we found an element, it should be active
  libmesh_assert (!this->_element || this->_element->active());

  // If we found an element and have a restriction list, they better match
  libmesh_assert (!this->_element || !allowed_subdomains || allowed_subdomains->count(this->_element->subdomain_id()));

  // return the element
  return this->_element;
}


template <typename RealType>
void PointLocatorTreeTempl<RealType>::operator() (const Point & p,
                                   std::set<const Elem *> & candidate_elements,
                                   const std::set<subdomain_id_type> * allowed_subdomains) const
{
  libmesh_assert (this->_initialized);

  LOG_SCOPE("operator() - Version 2", "PointLocatorTree");

  // ask the tree
  this->_tree->find_elements (p, candidate_elements, allowed_subdomains, this->_close_to_point_tol);
}



template <typename RealType>
const ElemTempl<RealType> * PointLocatorTreeTempl<RealType>::perform_linear_search(const Point & p,
                                                     const std::set<subdomain_id_type> * allowed_subdomains,
                                                     bool use_close_to_point,
                                                     Real close_to_point_tolerance) const
{
  LOG_SCOPE("perform_linear_search", "PointLocatorTree");

  // The type of iterator depends on the Trees::BuildType
  // used for this PointLocator.  If it's
  // TREE_LOCAL_ELEMENTS, we only want to double check
  // local elements during this linear search.
  SimpleRange<typename MeshBase::const_element_iterator> r =
    this->_build_type == Trees::LOCAL_ELEMENTS ?
    this->_mesh.active_local_element_ptr_range() :
    this->_mesh.active_element_ptr_range();

  for (const auto & elem : r)
    {
      if (!allowed_subdomains ||
          allowed_subdomains->count(elem->subdomain_id()))
        {
          if (!use_close_to_point)
            {
              if (elem->contains_point(p))
                return elem;
            }
          else
            {
              if (elem->close_to_point(p, close_to_point_tolerance))
                return elem;
            }
        }
    }

  return nullptr;
}


template <typename RealType>
std::set<const ElemTempl<RealType> *> PointLocatorTreeTempl<RealType>::perform_fuzzy_linear_search(const Point & p,
                                                                     const std::set<subdomain_id_type> * allowed_subdomains,
                                                                     Real close_to_point_tolerance) const
{
  LOG_SCOPE("perform_fuzzy_linear_search", "PointLocatorTree");

  std::set<const Elem *> candidate_elements;

  // The type of iterator depends on the Trees::BuildType
  // used for this PointLocator.  If it's
  // TREE_LOCAL_ELEMENTS, we only want to double check
  // local elements during this linear search.
  SimpleRange<typename MeshBase::const_element_iterator> r =
    this->_build_type == Trees::LOCAL_ELEMENTS ?
    this->_mesh.active_local_element_ptr_range() :
    this->_mesh.active_element_ptr_range();

  for (const auto & elem : r)
    if ((!allowed_subdomains || allowed_subdomains->count(elem->subdomain_id())) && elem->close_to_point(p, close_to_point_tolerance))
      candidate_elements.insert(elem);

  return candidate_elements;
}



template <typename RealType>
void PointLocatorTreeTempl<RealType>::enable_out_of_mesh_mode ()
{
  // Out-of-mesh mode should now work properly even on meshes with
  // non-affine elements.
  _out_of_mesh_mode = true;
}


template <typename RealType>
void PointLocatorTreeTempl<RealType>::disable_out_of_mesh_mode ()
{
  _out_of_mesh_mode = false;
}


template <typename RealType>
void PointLocatorTreeTempl<RealType>::set_target_bin_size (unsigned int target_bin_size)
{
  _target_bin_size = target_bin_size;
}


template <typename RealType>
unsigned int PointLocatorTreeTempl<RealType>::get_target_bin_size () const
{
  return _target_bin_size;
}


} // namespace libMesh

#endif // LIBMESH_POINT_LOCATOR_TREE_IMPL_H
