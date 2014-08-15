// The libMesh Finite Element Library.
// Copyright (C) 2002-2014 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/point_locator_list.h"

namespace libMesh
{




//------------------------------------------------------------------
// PointLocatorBase methods
PointLocatorBase::PointLocatorBase (const MeshBase& mesh,
                                    const PointLocatorBase* master) :
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



AutoPtr<PointLocatorBase> PointLocatorBase::build (PointLocatorType t,
                                                   const MeshBase& mesh,
                                                   const PointLocatorBase* master)
{
  switch (t)
    {
    case TREE:
      return AutoPtr<PointLocatorBase>(new PointLocatorTree(mesh, /*Trees::NODES,*/ master));

    case TREE_ELEMENTS:
      return AutoPtr<PointLocatorBase>(new PointLocatorTree(mesh, Trees::ELEMENTS, master));

    case TREE_LOCAL_ELEMENTS:
      return AutoPtr<PointLocatorBase>(new PointLocatorTree(mesh, Trees::LOCAL_ELEMENTS, master));

    case LIST:
      return AutoPtr<PointLocatorBase>(new PointLocatorList(mesh, master));

    default:
      libmesh_error_msg("ERROR: Bad PointLocatorType = " << t);
    }

  libmesh_error_msg("We'll never get here!");
  return AutoPtr<PointLocatorBase>();
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

} // namespace libMesh
