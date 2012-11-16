
#ifndef __zero_function_h__
#define __zero_function_h__

// Local includes
#include "libmesh/const_function.h"

// C++ includes

namespace libMesh {

template <typename Output=Number>
class ZeroFunction : public ConstFunction<Output>
{
public:
  ZeroFunction () : ConstFunction<Output>(0) {}

  virtual AutoPtr<FunctionBase<Output> > clone() const {
    return AutoPtr<FunctionBase<Output> > 
      (new ZeroFunction<Output>());
  }
};

} // namespace libMesh

#endif // __zero_function_h__
