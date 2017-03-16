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
#include "libmesh/point_locator_base.h"
#include "libmesh/point_locator_tree.h"

#include "libmesh/elem.h"

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
  _close_to_point_tol      (TOLERANCE)
{
}





PointLocatorBase::~PointLocatorBase ()
{
}



bool PointLocatorBase::initialized () const
{
  return this->_initialized;
}



UniquePtr<PointLocatorBase> PointLocatorBase::build (PointLocatorType t,
                                                     const MeshBase & mesh,
                                                     const PointLocatorBase * master)
{
  switch (t)
    {
    case TREE:
      return UniquePtr<PointLocatorBase>(new PointLocatorTree(mesh, /*Trees::NODES,*/ master));

    case TREE_ELEMENTS:
      return UniquePtr<PointLocatorBase>(new PointLocatorTree(mesh, Trees::ELEMENTS, master));

    case TREE_LOCAL_ELEMENTS:
      return UniquePtr<PointLocatorBase>(new PointLocatorTree(mesh, Trees::LOCAL_ELEMENTS, master));

    default:
      libmesh_error_msg("ERROR: Bad PointLocatorType = " << t);
    }

  libmesh_error_msg("We'll never get here!");
  return UniquePtr<PointLocatorBase>();
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


const Node *
PointLocatorBase::
locate_node(const Point & p,
            const std::set<subdomain_id_type> * allowed_subdomains,
            Real tol) const
{
  std::set<const Elem *> candidate_elements;
  this->operator()(p, candidate_elements, allowed_subdomains);

  for (std::set<const Elem *>::const_iterator
         it = candidate_elements.begin();
       it != candidate_elements.end(); ++it)
    {
      const Elem * elem = *it;
      const int elem_n_nodes = elem->n_nodes();
      const Real hmax = elem->hmax();
      const Real dist_tol_sq = (tol * hmax) * (tol * hmax);

      for (int n=0; n != elem_n_nodes; ++n)
        if ((elem->point(n) - p).norm_sq() < dist_tol_sq)
          return elem->node_ptr(n);
    }

  return libmesh_nullptr;
}

} // namespace libMesh
