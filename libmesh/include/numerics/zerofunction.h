
#ifndef __parsedfunction_h__
#define __parsedfunction_h__

#include <string>

#include "dense_vector.h"
#include "function_base.h"
#include "point.h"

class ZeroFunction : public FunctionBase
{
  virtual Number operator() (const Point&,
                             const Real = 0)
    {
      return 0;
    }

  virtual void operator() (const Point&,
                           const Real,
                           DenseVector<Number>& output)
    {
      unsigned int size = output.size();
      for (unsigned int i=0; i != size; ++i)
        output(i) = 0;
    }

  virtual void init() {}
  virtual void clear() {}
};


#endif // __parsedfunction_h__
