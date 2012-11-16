
#ifndef __const_fem_function_h__
#define __const_fem_function_h__

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

  virtual Output operator() (const FEMContext&, const Point& p,
			     const Real time = 0.)
    { return _c; }

  private:
  Output _c;
};

} // namespace libMesh;


#endif // __const_function_h__
