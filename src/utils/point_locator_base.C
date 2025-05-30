// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



// Local Includes
#include "libmesh/point_locator_base.h"
#include "libmesh/point_locator_tree.h"
#include "libmesh/elem.h"
#include "libmesh/enum_point_locator_type.h"
#include "libmesh/point_locator_nanoflann.h"

namespace libMesh
{




//------------------------------------------------------------------
// PointLocatorBase methods
PointLocatorBase::PointLocatorBase (const MeshBase & mesh,
                                    const PointLocatorBase * master) :
  _verbose                 (false),
  _master                  (master),
  _mesh                    (mesh),
  _initialized             (false),
  _use_close_to_point_tol  (false),
  _close_to_point_tol      (TOLERANCE),
  _use_contains_point_tol  (false),
  _contains_point_tol      (TOLERANCE)
{
  // If we have a non-nullptr master, inherit its close-to-point tolerances.
  if (_master)
    {
      _use_close_to_point_tol = _master->_use_close_to_point_tol;
      _close_to_point_tol = _master->_close_to_point_tol;
    }
}



PointLocatorBase::~PointLocatorBase () = default;



bool PointLocatorBase::initialized () const
{
  return this->_initialized;
}



std::unique_ptr<PointLocatorBase> PointLocatorBase::build (PointLocatorType t,
                                                           const MeshBase & mesh,
                                                           const PointLocatorBase * master)
{
  switch (t)
    {
    case TREE:
      return std::make_unique<PointLocatorTree>(mesh, /*Trees::NODES,*/ master);

    case TREE_ELEMENTS:
      return std::make_unique<PointLocatorTree>(mesh, Trees::ELEMENTS, master);

    case TREE_LOCAL_ELEMENTS:
      return std::make_unique<PointLocatorTree>(mesh, Trees::LOCAL_ELEMENTS, master);

#ifdef LIBMESH_HAVE_NANOFLANN
    case NANOFLANN:
      return std::make_unique<PointLocatorNanoflann>(mesh, master);
#endif

    default:
      libmesh_error_msg("ERROR: Bad PointLocatorType = " << t);
    }
}

Real PointLocatorBase::get_close_to_point_tol () const
{
  return _close_to_point_tol;
}


void PointLocatorBase::set_close_to_point_tol (Real close_to_point_tol)
{
  _use_close_to_point_tol = true;
  _close_to_point_tol = close_to_point_tol;
}

void PointLocatorBase::unset_close_to_point_tol ()
{
  _use_close_to_point_tol = false;
  _close_to_point_tol = TOLERANCE;
}

void PointLocatorBase::set_contains_point_tol(Real contains_point_tol)
{
  _use_contains_point_tol = true;
  _contains_point_tol = contains_point_tol;
}

void PointLocatorBase::unset_contains_point_tol()
{
  _use_contains_point_tol = false;
  _contains_point_tol = TOLERANCE;
}

Real PointLocatorBase::get_contains_point_tol() const
{
  return _contains_point_tol;
}

const MeshBase & PointLocatorBase::get_mesh () const
{
  return _mesh;
}


const Node *
PointLocatorBase::
locate_node(const Point & p,
            const std::set<subdomain_id_type> * allowed_subdomains,
            Real tol) const
{
  std::set<const Elem *> candidate_elements;
  this->operator()(p, candidate_elements, allowed_subdomains);

  for (const auto & elem : candidate_elements)
    {
      const int elem_n_nodes = elem->n_nodes();
      const Real hmax = elem->hmax();
      const Real dist_tol_sq = (tol * hmax) * (tol * hmax);

      for (int n=0; n != elem_n_nodes; ++n)
        if ((elem->point(n) - p).norm_sq() < dist_tol_sq)
          return elem->node_ptr(n);
    }

  return nullptr;
}

} // namespace libMesh
