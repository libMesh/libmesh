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

// libMesh Includes
#include "libmesh/meshfree_interpolation_function.h"

// C++ includes
#include <cstddef>
#include <memory>
#include <vector>

namespace libMesh
{

MeshfreeInterpolationFunction::MeshfreeInterpolationFunction(const MeshfreeInterpolation & mfi, Threads::spin_mutex & mutex)
    : _mfi(mfi), _mutex(mutex) {}

Number MeshfreeInterpolationFunction::operator()(const Point & p, const Real /* time */)
{
  _pts.clear();
  _pts.push_back(p);
  _vals.resize(1);

  Threads::spin_mutex::scoped_lock lock(_mutex);

  _mfi.interpolate_field_data(_mfi.field_variables(), _pts, _vals);

  return _vals.front();
}

void MeshfreeInterpolationFunction::operator()(const Point & p, const Real time, DenseVector<Number> & output)
{
  output.resize(1);
  output(0) = (*this)(p, time);
}

void MeshfreeInterpolationFunction::init()
{
  // Initialization logic if any
}

void MeshfreeInterpolationFunction::clear()
{
  _pts.clear();
  _vals.clear();
}

std::unique_ptr<FunctionBase<Number>> MeshfreeInterpolationFunction::clone() const
{
  return std::make_unique<MeshfreeInterpolationFunction>(_mfi, _mutex);
}

} // namespace libMesh
