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
PointLocatorTree::PointLocatorTree (const MeshBase & mesh,
                                    const PointLocatorBase * master) :
  PointLocatorBase (mesh,master),
  _tree            (libmesh_nullptr),
  _element         (libmesh_nullptr),
  _out_of_mesh_mode(false),
  _target_bin_size (200),
  _build_type(Trees::NODES)
{
  this->init(_build_type);
}



PointLocatorTree::PointLocatorTree (const MeshBase & mesh,
                                    const Trees::BuildType build_type,
                                    const PointLocatorBase * master) :
  PointLocatorBase (mesh,master),
  _tree            (libmesh_nullptr),
  _element         (libmesh_nullptr),
  _out_of_mesh_mode(false),
  _target_bin_size (200),
  _build_type(build_type)
{
  this->init(_build_type);
}



PointLocatorTree::~PointLocatorTree ()
{
  this->clear ();
}



void PointLocatorTree::clear ()
{
  // only delete the tree when we are the master
  if (this->_tree != libmesh_nullptr)
    {
      if (this->_master == libmesh_nullptr)
        // we own the tree
        delete this->_tree;
      else
        // someone else owns and therefore deletes the tree
        this->_tree = libmesh_nullptr;

      // make sure operator () throws an assertion
      this->_initialized = false;
    }
}



void PointLocatorTree::init()
{
  this->init(_build_type);
}



void PointLocatorTree::init (Trees::BuildType build_type)
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

      if (this->_master == libmesh_nullptr)
        {
          LOG_SCOPE("init(no master)", "PointLocatorTree");

          if (this->_mesh.mesh_dimension() == 3)
            _tree = new Trees::OctTree (this->_mesh, get_target_bin_size(), _build_type);
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
                MeshTools::BoundingBox bbox = MeshTools::bounding_box(this->_mesh);

                const Real
                  Dx = bbox.second(0) - bbox.first(0),
                  Dz = bbox.second(2) - bbox.first(2);

                if (std::abs(Dz/(Dx + 1.e-20)) < 1e-10)
                  is_planar_xy = true;
              }

              if (!is_planar_xy)
                _tree = new Trees::OctTree (this->_mesh, get_target_bin_size(), _build_type);
              else
#endif
#if LIBMESH_DIM > 1
                _tree = new Trees::QuadTree (this->_mesh, get_target_bin_size(), _build_type);
#else
              _tree = new Trees::BinaryTree (this->_mesh, get_target_bin_size(), _build_type);
#endif
            }
        }

      else
        {
          // We are _not_ the master.  Let our Tree point to
          // the master's tree.  But for this we first transform
          // the master in a state for which we are friends.
          // And make sure the master @e has a tree!
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
      this->_element = libmesh_nullptr;
    }

  // ready for take-off
  this->_initialized = true;
}



const Elem * PointLocatorTree::operator() (const Point & p,
                                           const std::set<subdomain_id_type> * allowed_subdomains) const
{
  libmesh_assert (this->_initialized);

  LOG_SCOPE("operator()", "PointLocatorTree");

  // If we're provided with an allowed_subdomains list and have a cached element, make sure it complies
  if (allowed_subdomains && this->_element && !allowed_subdomains->count(this->_element->subdomain_id())) this->_element = libmesh_nullptr;

  // First check the element from last time before asking the tree
  if (this->_element==libmesh_nullptr || !(this->_element->contains_point(p)))
    {
      // ask the tree
      this->_element = this->_tree->find_element (p,allowed_subdomains);

      if (this->_element == libmesh_nullptr)
        {
          // If we haven't found the element, we may want to do a linear
          // search using a tolerance.
          if( _use_close_to_point_tol )
            {
              if(_verbose)
                {
                  libMesh::out << "Performing linear search using close-to-point tolerance "
                               << _close_to_point_tol
                               << std::endl;
                }

              this->_element =
                this->perform_linear_search(p,
                                            allowed_subdomains,
                                            /*use_close_to_point*/ true,
                                            _close_to_point_tol);

              return this->_element;
            }

          // No element seems to contain this point.  In theory, our
          // tree now correctly handles curved elements.  In
          // out-of-mesh mode this is sometimes expected, and we can
          // just return NULL without searching further.  Out of
          // out-of-mesh mode, something must have gone wrong.

          // We'll avoid making this assertion just yet, though,
          // because some users are leaving out_of_mesh_mode disabled
          // as a workaround for old PointLocatorTree behavior.
          // libmesh_assert_equal_to (_out_of_mesh_mode, true);

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


void PointLocatorTree::operator() (const Point & p,
                                   std::set<const Elem *> & candidate_elements,
                                   const std::set<subdomain_id_type> * allowed_subdomains) const
{
  libmesh_assert (this->_initialized);

  LOG_SCOPE("operator() - Version 2", "PointLocatorTree");

  // forward call to perform_linear_search
  candidate_elements = this->perform_fuzzy_linear_search(p, allowed_subdomains, _close_to_point_tol);
}



const Elem * PointLocatorTree::perform_linear_search(const Point & p,
                                                     const std::set<subdomain_id_type> * allowed_subdomains,
                                                     bool use_close_to_point,
                                                     Real close_to_point_tolerance) const
{
  LOG_SCOPE("perform_linear_search", "PointLocatorTree");

  // The type of iterator depends on the Trees::BuildType
  // used for this PointLocator.  If it's
  // TREE_LOCAL_ELEMENTS, we only want to double check
  // local elements during this linear search.
  MeshBase::const_element_iterator pos =
    this->_build_type == Trees::LOCAL_ELEMENTS ?
    this->_mesh.active_local_elements_begin() : this->_mesh.active_elements_begin();

  const MeshBase::const_element_iterator end_pos =
    this->_build_type == Trees::LOCAL_ELEMENTS ?
    this->_mesh.active_local_elements_end() : this->_mesh.active_elements_end();

  for ( ; pos != end_pos; ++pos)
    {
      if (!allowed_subdomains ||
          allowed_subdomains->count((*pos)->subdomain_id()))
        {
          if(!use_close_to_point)
            {
              if ((*pos)->contains_point(p))
                return (*pos);
            }
          else
            {
              if ((*pos)->close_to_point(p, close_to_point_tolerance))
                return (*pos);
            }
        }
    }

  return libmesh_nullptr;
}


std::set<const Elem *> PointLocatorTree::perform_fuzzy_linear_search(const Point & p,
                                                                     const std::set<subdomain_id_type> * allowed_subdomains,
                                                                     Real close_to_point_tolerance) const
{
  LOG_SCOPE("perform_fuzzy_linear_search", "PointLocatorTree");

  std::set<const Elem *> candidate_elements;

  // The type of iterator depends on the Trees::BuildType
  // used for this PointLocator.  If it's
  // TREE_LOCAL_ELEMENTS, we only want to double check
  // local elements during this linear search.
  MeshBase::const_element_iterator pos =
    this->_build_type == Trees::LOCAL_ELEMENTS ?
    this->_mesh.active_local_elements_begin() : this->_mesh.active_elements_begin();

  const MeshBase::const_element_iterator end_pos =
    this->_build_type == Trees::LOCAL_ELEMENTS ?
    this->_mesh.active_local_elements_end() : this->_mesh.active_elements_end();

  for ( ; pos != end_pos; ++pos)
    {
      if ((!allowed_subdomains || allowed_subdomains->count((*pos)->subdomain_id())) &&
          (*pos)->close_to_point(p, close_to_point_tolerance))
        candidate_elements.insert(*pos);
    }

  return candidate_elements;
}



void PointLocatorTree::enable_out_of_mesh_mode ()
{
  // Out-of-mesh mode should now work properly even on meshes with
  // non-affine elements.
  _out_of_mesh_mode = true;
}


void PointLocatorTree::disable_out_of_mesh_mode ()
{
  _out_of_mesh_mode = false;
}


void PointLocatorTree::set_target_bin_size (unsigned int target_bin_size)
{
  _target_bin_size = target_bin_size;
}


unsigned int PointLocatorTree::get_target_bin_size () const
{
  return _target_bin_size;
}


} // namespace libMesh
