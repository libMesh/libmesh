// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_RADIAL_BASIS_FUNCTIONS_H
#define LIBMESH_RADIAL_BASIS_FUNCTIONS_H

// C++ includes
#include <limits>

// Local includes
#include "libmesh_common.h"



namespace libMesh
{

  /**
   * Wendland's compactly supported Radial Basis Functions.
   */
  template <unsigned int SpaceDim, unsigned int Continuity>
  class WendlandRBF
  {
  private:
    const Real _rcut;

  public:
    /**
     * Constructor.
     */
    WendlandRBF (const Real r_cut = std::numeric_limits<Real>::max()) :
      _rcut (r_cut)
    {}
  }

  //-------------------------------------------------------
  // Explicit instantiations
} // namespace libMesh


#endif // #define LIBMESH_RADIAL_BASIS_FUNCTIONS_H
