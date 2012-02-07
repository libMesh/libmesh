
#ifndef __parsed_function_h__
#define __parsed_function_h__

#include <cmath>
#include <string>

#include "dense_vector.h"
#include "function_base.h"
#include "point.h"

#ifdef LIBMESH_HAVE_FPARSER

#include "fparser.hh"

template <typename Output=Number>
class ParsedFunction : public FunctionBase<Output>
{
public:
  ParsedFunction (const std::string& expression)
    {
      std::string variables = "x";
#if LIBMESH_DIM > 1
      variables += ",y";
#endif
#if LIBMESH_DIM > 2
      variables += ",z";
#endif
      variables += ",t";

      for (size_t start = 0, end = 0;
           end != std::string::npos;
           end = expression.find( ',', start ))
      {
          // We either want the whole end of the string (end ==
          // npos) or a substring in the middle.
        std::string subexpression =
          expression.substr(start, (end == std::string::npos) ?
                                    std::string::npos : end - start);
    
          // If at end, use start=maxSize.  Else start at next
          // character.
        start = (end == std::string::npos) ?
                std::string::npos : end + 1;

        FunctionParserBase<Output> fp;
        fp.AddConstant("pi", std::acos(Real(-1)));
        fp.AddConstant("e", std::exp(Real(1)));
        fp.Parse(subexpression, variables);
        parsers.push_back(fp);
      }
    }

  virtual Output operator() (const Point& p,
                             const Real time = 0)
    {
      _spacetime[0] = p(0);
#if LIBMESH_DIM > 1
      _spacetime[1] = p(1);
#endif
#if LIBMESH_DIM > 2
      _spacetime[2] = p(2);
#endif
      _spacetime[LIBMESH_DIM] = time;
      return parsers[0].Eval(_spacetime);
    }

  virtual void operator() (const Point& p,
                           const Real time,
                           DenseVector<Output>& output)
    {
      _spacetime[0] = p(0);
#if LIBMESH_DIM > 1
      _spacetime[1] = p(1);
#endif
#if LIBMESH_DIM > 2
      _spacetime[2] = p(2);
#endif
      _spacetime[LIBMESH_DIM] = time;

      unsigned int size = output.size();
      libmesh_assert(size == parsers.size());

      for (unsigned int i=0; i != size; ++i)
        output(i) = parsers[i].Eval(_spacetime);
    }

  /**
   * @returns the vector component \p i at coordinate
   * \p p and time \p time.
   */
  virtual Output component (unsigned int i,
                            const Point& p,
                            Real time)
    {
      _spacetime[0] = p(0);
#if LIBMESH_DIM > 1
      _spacetime[1] = p(1);
#endif
#if LIBMESH_DIM > 2
      _spacetime[2] = p(2);
#endif
      _spacetime[LIBMESH_DIM] = time;

      libmesh_assert(i < parsers.size());

      return parsers[i].Eval(_spacetime);
    }

  virtual void init() {}
  virtual void clear() {}

private:
  Output _spacetime[LIBMESH_DIM+1];
  std::vector<FunctionParserBase<Output> > parsers;
};

#else

template <typename Output>
class ParsedFunction : public FunctionBase<Output>
{
public:
  ParsedFunction (std::string expression)
    {
      libmesh_not_implemented();
    }

  virtual Output operator() (const Point&,
                             const Real time = 0)
    { return 0.; }

  virtual void operator() (const Point&,
                           const Real time,
                           DenseVector<Output>& output) {}

  virtual void init() {}
  virtual void clear() {}
};

#endif // LIBMESH_HAVE_FPARSER

#endif // __parsed_function_h__
