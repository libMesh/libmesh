// $Id: hp_singular.h,v 1.1 2006-10-05 20:50:15 roystgnr Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2006  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __hp_singular_h__
#define __hp_singular_h__

// C++ includes
#include <vector>

// Local Includes
#include "auto_ptr.h"
#include "dense_matrix.h"
#include "dense_vector.h"
#include "fe.h"         // MipsPro requires fe.h and quadrature.h in order to
#include "quadrature.h" //  delete AutoPtrs<> upon destruction
#include "libmesh_common.h"

// Forward Declarations
class Point;
class System;


/**
 * This class uses a user-provided list of singularity locations
 * to choose between h refining and p elevation.
 * Currently we assume that a set of elements has already been flagged
 * for h refinement - any elements which are not adjacent to a
 * user-provided singularity are instead flagged for p refinement.
 *
 * @author Roy H. Stogner, 2006.
 */
class HPSingularity
{
public:

  /**
   * Constructor.
   */
  HPSingularity()
  {
    untested();
  }
  
  /**
   * Destructor.  
   */
  virtual ~HPSingularity () {}


  /**
   * This pure virtual function must be redefined
   * in derived classes to take a mesh flagged for h
   * refinement and potentially change the desired
   * refinement type.
   */
  virtual void select_refinement (System& system);

};


#endif

