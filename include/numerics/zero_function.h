// The libMesh Finite Element Library.
// Copyright (C) 2002-2022 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_ZERO_FUNCTION_H
#define LIBMESH_ZERO_FUNCTION_H

// Local includes
#include "libmesh/const_function.h"

// C++ Includes
#include <memory>

namespace libMesh
{

/**
 * ConstFunction that simply returns 0.
 *
 * \author Roy Stogner
 * \date 2012
 * \brief ConstFunction that simply returns 0.
 */
template <typename Output=Number>
class ZeroFunction : public ConstFunction<Output>
{
public:
  ZeroFunction () : ConstFunction<Output>(0) {}

  /**
   * The 5 special functions can be defaulted for this class.
   */
  ZeroFunction (ZeroFunction &&) = default;
  ZeroFunction (const ZeroFunction &) = default;
  ZeroFunction & operator= (const ZeroFunction &) = default;
  ZeroFunction & operator= (ZeroFunction &&) = default;
  virtual ~ZeroFunction () = default;

  virtual std::unique_ptr<FunctionBase<Output>> clone() const override
  {
    return std::make_unique<ZeroFunction<Output>>();
  }
};

} // namespace libMesh

#endif // LIBMESH_ZERO_FUNCTION_H
