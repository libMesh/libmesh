// The libMesh Finite Element Library.
// Copyright (C) 2003  Benjamin S. Kirk

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



// C++ Includes
#include <math.h>

// Mesh library includes
#include "libmesh/libmesh_common.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;





/**
 * This is the exact solution that
 * we are trying to obtain.  We will solve
 *
 * - (u_xx + u_yy) = f
 *
 * and take a finite difference approximation using this
 * function to get f.  This is the well-known "method of
 * manufactured solutions".
 */
Real exact_solution (const Real x,
                     const Real y,
                     const Real z = 0.)
{
  static const Real pi = acos(-1.);

  return cos(.5*pi*x)*sin(.5*pi*y)*cos(.5*pi*z);
}
