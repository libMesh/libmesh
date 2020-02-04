// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/fe_base.h" // required for the type OutputGradient.

// C++ includes
#include <vector>

namespace libMesh
{

// Forward declarations
class EquationSystems;
template <typename> class PointTempl;
typedef PointTempl<Real> Point;

/**
 * class \p FEComputeData hides arbitrary data to be passed to and from
 * children of \p FEBase through the \p FEInterface::compute_data()
 * method.  This enables the efficient computation of data on
 * the finite element level, while maintaining library integrity.
 * - With special finite elements disabled (like infinite elements),
 *   this class wraps the return values of all shape functions
 *   from \p FEInterface::shape() in a \p std::vector<Number>.
 * - With infinite elements enabled, this class returns a vector of physically
 *   correct shape functions, both for finite and infinite elements.
 *
 * \author Daniel Dreyer
 * \date 2003
 * \brief Helper class used with FEInterface::compute_data().
 */
class FEComputeData
{
public:
  /**
   * Constructor.  Takes the required input data and clears
   * the output data using \p clear().
   */
  FEComputeData (const EquationSystems & es,
                 const Point & pin) :
    equation_systems(es),
    p(pin),
    _need_dshape(false)
  {
    this->clear();
  }

  /**
   * Const reference to the \p EquationSystems object
   * that contains simulation-specific data.
   */
  const EquationSystems & equation_systems;

  /**
   * Holds the point where the data are to be computed
   */
  const Point & p;

  /**
   * Storage for the computed shape function values.
   */
  std::vector<Number> shape;

  /**
   * Storage for the computed shape derivative values.
   */
  std::vector<Gradient> dshape;

  /**
   * Storage for local to global mapping at \p p.
   * This is needed when the gradient in physical
   * coordinates is of interest.
   * FIXME: What kind of type should one use for it?
   *  The matrix-class don't look as if they were made for it
   *  and neither are the TensorTool-members.
   */
  std::vector<std::vector<Real>> local_transform;

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
  /**
   * Storage for the computed phase lag
   */
  Real phase;

  /**
   * The wave speed.
   */
  Real speed;

  /**
   * The frequency to evaluate shape functions
   * including the wave number depending terms.
   * Use imaginary contributions for exponential damping
   */
  Number frequency;
#endif

  /**
   * Clears the output data completely.
   */
  void clear ();

  /**
   * Inits the output data to default values, provided
   * the fields are correctly resized.
   */
  void init ();

  /**
   * Enable the computation of shape gradients (dshape).
   */
  void enable_derivative ();

  /**
   * Disable the computation of shape gradients (dshape).
   * Default is disabled.
   */
  void disable_derivative ()
  {_need_dshape=false; }

  /**
   * Check whether derivatives should be computed or not.
   */
  bool need_derivative ()
  {return _need_dshape; }

private:
  /**
   * variable indicating whether the shape-derivative should be computed or not.
   * Default is false to save time and be compatible with elements where derivatives
   * are not implemented/ cannot be computed.
   */
  bool _need_dshape;

};


} // namespace libMesh

#endif // LIBMESH_FE_COMPUTE_DATA_H
