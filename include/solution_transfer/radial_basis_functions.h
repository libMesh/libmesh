// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/libmesh_common.h"
#include "libmesh/utility.h"



namespace libMesh
{

// /**
//  * Simple radial basis function.
//  */
// class SimpleRBF
// {
// private:
//   const Real _rcut;

// public:

//   /**
//    * Constructor.
//    */
//   SimpleRBF (const Real r_cut = 1.) :
//     _rcut (r_cut)
//   {}

//   /**
//    * Evaluate the radial basis function at the reqested location.
//    */
//   Real operator()(Real rad) const
//   {
//     if (rad > _rcut) return 0.;

//     rad /= _rcut;

//     return std::sqrt( 1+ rad*rad );
//   }
// };



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
  WendlandRBF (const Real r_cut = 1.) :
    _rcut (r_cut)
  { libmesh_experimental(); }

  /**
   * Evaluate the radial basis function at the reqested location.
   */
  Real operator()(Real /* rad */) const { libmesh_not_implemented(); return 0.; }
};



//-------------------------------------------------------
// Explicit specializations
template<>
inline
Real WendlandRBF<3,0>::operator()(Real rad) const
{
  if (rad > _rcut) return 0.;

  rad /= _rcut;

  return Utility::pow<2>(1.-rad);
}

template<>
inline
Real WendlandRBF<3,2>::operator()(Real rad) const
{
  if (rad > _rcut) return 0.;

  rad /= _rcut;

  return Utility::pow<4>(1.-rad)*(4.*rad + 1.);
}

template<>
inline
Real WendlandRBF<3,4>::operator()(Real rad) const
{
  if (rad > _rcut) return 0.;

  rad /= _rcut;

  return Utility::pow<6>(1.-rad)*((35.*rad + 18.)*rad + 3.);
}

template<>
inline
Real WendlandRBF<3,8>::operator()(Real rad) const
{
  if (rad > _rcut) return 0.;

  rad /= _rcut;

  return Utility::pow<8>(1.-rad)*(((32.*rad + 25.)*rad + 8.)*rad + 1.);
}


} // namespace libMesh


#endif // #define LIBMESH_RADIAL_BASIS_FUNCTIONS_H
