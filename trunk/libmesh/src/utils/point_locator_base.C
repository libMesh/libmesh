// $Id: point_locator_base.C,v 1.5 2004-01-03 15:37:44 benkirk Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2004  Benjamin S. Kirk, John W. Peterson
  
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
#include "point_locator_base.h"
#include "mesh.h"
#include "point.h"
#include "point_locator_tree.h"
#include "point_locator_list.h"




//------------------------------------------------------------------
// PointLocatorBase methods
PointLocatorBase::PointLocatorBase (const Mesh& mesh,
				    const PointLocatorBase* master) :
  _master                  (master),
  _mesh                    (mesh),
  _initialized             (false)
{
}





PointLocatorBase::~PointLocatorBase ()
{
}





AutoPtr<PointLocatorBase> PointLocatorBase::build (const PointLocatorType t,
						   const Mesh& mesh,
						   const PointLocatorBase* master)
{
  switch (t)
    {
    case TREE:
      {
	AutoPtr<PointLocatorBase> ap(new PointLocatorTree(mesh,
							  master));
	return ap;
      }

    case LIST:
      {
	AutoPtr<PointLocatorBase> ap(new PointLocatorList(mesh,
							  master));
	return ap;
      }

    default:
      {
	std::cerr << "ERROR: Bad PointLocatorType = " << t << std::endl;
	error();
      }
    }

  error();
  AutoPtr<PointLocatorBase> ap(NULL);
  return ap;
}

