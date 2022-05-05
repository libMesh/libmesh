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


#ifndef LIBMESH_CONST_FUNCTION_H
#define LIBMESH_CONST_FUNCTION_H

// Local includes
#include "libmesh/dense_vector.h"
#include "libmesh/function_base.h"
#include "libmesh/point.h"

// C++ includes
#include <memory>
#include <string>

namespace libMesh
{

/**
 * Function that returns a single value that never changes. All
 * overridden virtual functions are documented in function_base.h.
 *
 * \author Roy Stogner
 * \date 2012
 * \brief Function that returns a single value that never changes.
 */
template <typename Output=Number>
class ConstFunction : public FunctionBase<Output>
{
public:
  explicit
  ConstFunction (const Output & c) : _c(c)
  {
    this->_initialized = true;
    this->_is_time_dependent = false;
  }

  /**
   * The 5 special functions can be defaulted for this class.
   */
  ConstFunction (ConstFunction &&) = default;
  ConstFunction (const ConstFunction &) = default;
  ConstFunction & operator= (const ConstFunction &) = default;
  ConstFunction & operator= (ConstFunction &&) = default;
  virtual ~ConstFunction () = default;

  virtual Output operator() (const Point &,
                             const Real = 0) override
  { return _c; }

  virtual void operator() (const Point &,
                           const Real,
                           DenseVector<Output> & output) override
  {
    unsigned int size = output.size();
    for (unsigned int i=0; i != size; ++i)
      output(i) = _c;
  }

  virtual std::unique_ptr<FunctionBase<Output>> clone() const override
  {
    return std::make_unique<ConstFunction<Output>>(_c);
  }

private:
  Output _c;
};

} // namespace libMesh

#endif // LIBMESH_CONST_FUNCTION_H
