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

#ifndef LIBMESH_CONST_FEM_FUNCTION_H
#define LIBMESH_CONST_FEM_FUNCTION_H

// libMesh includes
#include "libmesh/dense_vector.h"
#include "libmesh/fem_function_base.h"
#include "libmesh/int_range.h"
#include "libmesh/auto_ptr.h" // libmesh_make_unique

// C++ includes
#include <string>

namespace libMesh
{

// Forward declarations
template <typename> class PointTempl;
typedef PointTempl<Real> Point;

/**
 * FEMFunction that returns a single value, regardless of the time and
 * location inputs.
 *
 * \author Roy Stogner
 * \date 2012
 * \brief FEMFunction that returns a single value.
 */
template <typename Output=Number>
class ConstFEMFunction : public FEMFunctionBase<Output>
{
public:
  ConstFEMFunction (const Output c) : _c(c) {}

  /**
   * The 5 special functions can be defaulted for this class.
   */
  ConstFEMFunction (ConstFEMFunction &&) = default;
  ConstFEMFunction (const ConstFEMFunction &) = default;
  ConstFEMFunction & operator= (const ConstFEMFunction &) = default;
  ConstFEMFunction & operator= (ConstFEMFunction &&) = default;
  virtual ~ConstFEMFunction () = default;

  virtual std::unique_ptr<FEMFunctionBase<Output>> clone () const
  {return libmesh_make_unique<ConstFEMFunction>(*this); }

  virtual Output operator() (const FEMContext &,
                             const Point &,
                             const Real /* time */ = 0.)
  { return _c; }

  virtual void operator() (const FEMContext &,
                           const Point &,
                           const Real,
                           DenseVector<Output> & output)
  {
    for (auto i : index_range(output))
      output(i) = _c;
  }

private:
  Output _c;
};

} // namespace libMesh;


#endif // LIBMESH_CONST_FEM_FUNCTION_H
