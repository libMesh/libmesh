
#ifndef __zero_function_h__
#define __zero_function_h__

#include "const_function.h"

template <typename Output=Number>
class ZeroFunction : public ConstFunction<Output>
{
public:
  ZeroFunction () : ConstFunction<Output>(0) {}

  virtual AutoPtr<FunctionBase<Output> > clone() {
    return AutoPtr<FunctionBase<Output> > 
      (new ZeroFunction<Output>());
  }
};


#endif // __zero_function_h__
