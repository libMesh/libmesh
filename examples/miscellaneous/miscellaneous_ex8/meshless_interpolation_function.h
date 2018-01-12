// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_MESHLESS_INTERPOLATION_FUNCTION_H
#define LIBMESH_MESHLESS_INTERPOLATION_FUNCTION_H

// libMesh Includes
#include "libmesh/function_base.h"
#include "libmesh/meshfree_interpolation.h"
#include "libmesh/threads.h"
#include "libmesh/auto_ptr.h" // libmesh_make_unique

// C++ includes
#include <cstddef>

namespace libMesh
{



// Forward Declarations
template <typename T>
class DenseVector;


// ------------------------------------------------------------
// MeshlessInterpolationFunction class definition
class MeshlessInterpolationFunction : public FunctionBase<Number>
{
private:
  const MeshfreeInterpolation & _mfi;
  mutable std::vector<Point> _pts;
  mutable std::vector<Number> _vals;
  Threads::spin_mutex & _mutex;

public:

  /**
   * Constructor.  Requires a MeshlessInterpolation object.
   */
  MeshlessInterpolationFunction (const MeshfreeInterpolation & mfi,
                                 Threads::spin_mutex & mutex) :
    _mfi (mfi),
    _mutex(mutex)
  {}


  /**
   * The actual initialization process.
   */
  void init ();

  /**
   * Clears the function.
   */
  void clear ();

  /**
   * Returns a new deep copy of the function.
   */
  virtual std::unique_ptr<FunctionBase<Number>> clone () const;

  /**
   * @returns the value at point p and time
   * time, which defaults to zero.
   */
  Number operator() (const Point & p,
                     const Real time=0.);

  /**
   * Like before, but returns the values in a
   * writable reference.
   */
  void operator() (const Point & p,
                   const Real time,
                   DenseVector<Number> & output);

};



// ------------------------------------------------------------
// MeshlessInterpolationFunction inline methods
inline
Number MeshlessInterpolationFunction::operator() (const Point & p,
                                                  const Real /* time */)
{
  _pts.clear();
  _pts.push_back(p);
  _vals.resize(1);

  Threads::spin_mutex::scoped_lock lock(_mutex);

  _mfi.interpolate_field_data (_mfi.field_variables(),
                               _pts, _vals);

  return _vals.front();
}



inline
void MeshlessInterpolationFunction::operator() (const Point & p,
                                                const Real time,
                                                DenseVector<Number> & output)
{
  output.resize(1);
  output(0) = (*this)(p, time);
}



inline
void MeshlessInterpolationFunction::init ()
{
}



inline
void MeshlessInterpolationFunction::clear ()
{
}



inline
std::unique_ptr<FunctionBase<Number>>
MeshlessInterpolationFunction::clone () const
{
  return libmesh_make_unique<MeshlessInterpolationFunction>(_mfi, _mutex);
}


} // namespace libMesh


#endif // LIBMESH_MESHLESS_INTERPOLATION_FUNCTION_H
