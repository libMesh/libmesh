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
 *
 */
Real exact_solution (const Real x,
		     const Real y,
		     const Real t)
{
  const Real xo = 0.2;
  const Real yo = 0.2;
  const Real u  = 0.8;
  const Real v  = 0.8;

  const Real num =
    pow(x - u*t - xo, 2.) +
    pow(y - v*t - yo, 2.);

  const Real den =
    0.01*(4.*t + 1.);

  return exp(-num/den)/(4.*t + 1.);
}
