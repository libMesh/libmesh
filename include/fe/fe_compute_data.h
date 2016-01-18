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


#ifndef LIBMESH_FE_COMPUTE_DATA_H
#define LIBMESH_FE_COMPUTE_DATA_H

// Local includes
#include "libmesh/libmesh.h"

// C++ includes
#include <vector>

namespace libMesh
{

// Forward declarations
class EquationSystems;
class Point;

/**
 * class \p FEComputeData hides arbitrary data to be passed to and from
 * children of \p FEBase through the \p FEInterface::compute_data()
 * method.  This enables the efficient computation of data on
 * the finite element level, while maintaining library integrity.
 * -- With special finite elements @e disabled (like infinite elements),
 * this class wraps the return values of @e all shape functions
 * from \p FEInterface::shape() in a \p std::vector<Number>.
 * -- With enabled infinite elements, this class returns a vector of physically
 * correct shape functions, both for finite and infinite elements.
 */
class FEComputeData
{
public:


  //------------------------------------------------------
  // Conventional FEComputeData
  /**
   * Constructor.  Takes the required input data and clears
   * the output data using \p clear().
   */
  FEComputeData (const EquationSystems & es,
                 const Point & pin) :
    equation_systems(es),
    p(pin)
  {
    this->clear();
  }


  //----------------------------------------------------
  // Input data
  /**
   * Const reference to the \p EquationSystems object
   * that contains simulation-specific data.
   */
  const EquationSystems & equation_systems;
  /**
   * Holds the point where the data are to be computed
   */
  const Point & p;



  //----------------------------------------------------
  // Output data
  /**
   * Storage for the computed shape function values.
   */
  std::vector<Number> shape;


#if defined(LIBMESH_ENABLE_INFINITE_ELEMENTS) && !defined(LIBMESH_USE_COMPLEX_NUMBERS)
  /**
   * Storage for the computed phase lag
   */
  Real phase;

#endif

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
  /**
   * The wave speed.
   */
  Real speed;

#endif

#if defined LIBMESH_ENABLE_INFINITE_ELEMENTS && defined(LIBMESH_USE_COMPLEX_NUMBERS)
  /**
   * The frequency to evaluate shape functions
   * including the wave number depending terms
   */
  Real frequency;

#endif

  /**
   * Clears the @e output data completely.
   */
  void clear () ;


  /**
   * Inits the @e output data to default values, provided
   * the fields are correctly resized.
   */
  void init () ;
};


} // namespace libMesh

#endif // LIBMESH_FE_COMPUTE_DATA_H
