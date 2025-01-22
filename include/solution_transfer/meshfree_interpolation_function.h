// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef MESHFREE_INTERPOLATION_FUNCTION_H
#define MESHFREE_INTERPOLATION_FUNCTION_H

// libMesh Includes
#include "libmesh/function_base.h"
#include "libmesh/meshfree_interpolation.h"
#include "libmesh/threads.h"

// C++ includes
#include <cstddef>
#include <memory>
#include <vector>

namespace libMesh
{

// Forward Declarations
template <typename T> class DenseVector;

// ------------------------------------------------------------
// MeshfreeInterpolationFunction class definition
class MeshfreeInterpolationFunction : public FunctionBase<Number>
{
public:
  /**
   * Constructor. Requires a MeshfreeInterpolation object.
   */
  MeshfreeInterpolationFunction(const MeshfreeInterpolation & mfi,
                                Threads::spin_mutex & mutex);

  /**
   * Move and copy constructors.
   */
  MeshfreeInterpolationFunction(MeshfreeInterpolationFunction &&) = default;
  MeshfreeInterpolationFunction(const MeshfreeInterpolationFunction &) = default;

  /**
   * Destructor.
   */
  virtual ~MeshfreeInterpolationFunction() = default;

  /**
   * Delete assignment operators.
   */
  MeshfreeInterpolationFunction & operator=(const MeshfreeInterpolationFunction &) = delete;
  MeshfreeInterpolationFunction & operator=(MeshfreeInterpolationFunction &&) = delete;

  /**
   * The actual initialization process.
   */
  void init();

  /**
   * Clears the function.
   */
  void clear();

  /**
   * Returns a new deep copy of the function.
   */
  virtual std::unique_ptr<FunctionBase<Number>> clone() const;

  /**
   * @returns the value at point p and time
   * time, which defaults to zero.
   */
  Number operator()(const Point & p, const Real time = 0.);

  /**
   * Like before, but returns the values in a
   * writable reference.
   */
  void operator()(const Point & p, const Real time, DenseVector<Number> & output);

private:
  const MeshfreeInterpolation & _mfi;
  mutable std::vector<Point> _pts;
  mutable std::vector<Number> _vals;
  Threads::spin_mutex & _mutex;
};

} // namespace libMesh

#endif // MESHFREE_INTERPOLATION_FUNCTION_H
