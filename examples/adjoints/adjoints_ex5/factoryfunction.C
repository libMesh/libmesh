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



// libMesh includes
#include "libmesh/dense_vector.h"
#include "libmesh/factory.h"
#include "libmesh/function_base.h"

#include <memory>

using namespace libMesh;

// An example of a hand-coded function

class ExampleOneFunction : public FunctionBase<Number>
{
  virtual Number operator() (const Point & /*p*/,
                             const Real /*time*/)
  {
    return 1;
  }

  virtual void operator() (const Point & /*p*/,
                           const Real /*time*/,
                           DenseVector<Number> & output)
  {
    for (unsigned int i=0; i != output.size(); ++i)
      output(i) = 1;
  }

  virtual void init() {}
  virtual void clear() {}
  virtual std::unique_ptr<FunctionBase<Number>> clone() const
  {
    return std::make_unique<ExampleOneFunction>();
  }
};

#ifdef LIBMESH_USE_SEPARATE_NAMESPACE
namespace libMesh {
#endif

//-------------------------------------------------
// Full specialization for the Factory<FunctionBase<Number>>
// So we can look up hand-coded functions by name string
template<>
std::map<std::string, Factory<FunctionBase<Number>> *> &
Factory<FunctionBase<Number>>::factory_map()
{
  static std::map<std::string, Factory<FunctionBase<Number>> *> _map;
  return _map;
}

FactoryImp<ExampleOneFunction, FunctionBase<Number>> example_one_factory ("example_one");

#ifdef LIBMESH_USE_SEPARATE_NAMESPACE
} // namespace libMesh
#endif
