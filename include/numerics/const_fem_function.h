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

#ifndef LIBMESH_CONST_FEM_FUNCTION_H
#define LIBMESH_CONST_FEM_FUNCTION_H

#include <string>

#include "libmesh/dense_vector.h"
#include "libmesh/fem_function_base.h"
#include "libmesh/point.h"

namespace libMesh {

template <typename Output=Number>
class ConstFEMFunction : public FEMFunctionBase<Output>
{
public:
  ConstFEMFunction (const Output c) {_c = c;}

  ~ConstFEMFunction() {}

  virtual AutoPtr<FEMFunctionBase<Output> > clone () const
  {return AutoPtr<FEMFunctionBase<Output> >( new ConstFEMFunction(*this) ); }

  virtual Output operator() (const FEMContext&, const Point&,
			     const Real /* time */ = 0.)
    { return _c; }

  virtual void operator() (const FEMContext&, const Point&,
			   const Real,
			   DenseVector<Output>& output)
  {for(unsigned int i = 0; i < output.size(); i++ )
      output(i) = _c;}

private:
  Output _c;
};

} // namespace libMesh;


#endif // LIBMESH_CONST_FEM_FUNCTION_H
