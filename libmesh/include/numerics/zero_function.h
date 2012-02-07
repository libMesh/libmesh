
#ifndef __zero_function_h__
#define __zero_function_h__

#include <string>

#include "dense_vector.h"
#include "function_base.h"
#include "point.h"

template <typename Output=Number>
class ZeroFunction : public FunctionBase<Output>
{
  ZeroFunction () { this->initialized = true; }

  virtual Output operator() (const Point&,
                             const Real = 0)
    {
      return 0;
    }

  virtual void operator() (const Point&,
                           const Real,
                           DenseVector<Output>& output)
    {
      unsigned int size = output.size();
      for (unsigned int i=0; i != size; ++i)
        output(i) = 0;
    }

  virtual void init() {}

  virtual void clear() {}

  virtual AutoPtr<FunctionBase<Number> > clone() {
    return AutoPtr<FunctionBase<Number> > 
      (new ZeroFunction<Number>());
  }
};


#endif // __zero_function_h__
