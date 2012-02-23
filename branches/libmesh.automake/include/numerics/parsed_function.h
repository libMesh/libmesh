
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
  ParsedFunction (const std::string& expression, const std::vector<std::string>* additional_vars=NULL,
                  const std::vector<Output>* additional_vals=NULL)
    : _expression(expression),
      // Size the spacetime vector to account for space, time, and any additional
      // variables passed
      //_spacetime(LIBMESH_DIM+1 + (additional_vars ? additional_vars->size() : 0)),
      _additional_vars(additional_vars),
      _additional_vals(additional_vals)
    {
      std::string variables = "x";
#if LIBMESH_DIM > 1
      variables += ",y";
#endif
#if LIBMESH_DIM > 2
      variables += ",z";
#endif
      variables += ",t";

      _spacetime.resize(LIBMESH_DIM+1 + (additional_vars ? additional_vars->size() : 0));

      // If additional vars were passed, append them to the string
      // that we send to the function parser. Also add them to the
      // end of our spacetime vector
      if (additional_vars)
      {
        libmesh_assert(additional_vals != NULL && additional_vars->size() == additional_vals->size());
        for (unsigned int i=0; i < additional_vars->size(); ++i)
        {
          variables += "," + (*additional_vars)[i];
          _spacetime[LIBMESH_DIM+1 + i] = (*additional_vals)[i];
        }
      }

      size_t nextstart = 0, end = 0;

      while (end != std::string::npos)
      {
        // If we're past the end of the string, we can't make any more
        // subparsers
        if (nextstart >= expression.size())
          break;

	// If we're at the start of a brace delimited section, then we
	// parse just that section:
        if (expression[nextstart] == '{')
          {
            nextstart++;
            end = expression.find('}', nextstart);
          }
        // otherwise we parse the whole thing
        else
          end = std::string::npos;

	// We either want the whole end of the string (end == npos) or
	// a substring in the middle.
        std::string subexpression =
          expression.substr(nextstart, (end == std::string::npos) ?
                                        std::string::npos : end - nextstart);

        // fparser can crash on empty expressions
        libmesh_assert(!subexpression.empty());

        // Parse (and optimize if possible) the subexpression.
        // Add some basic constants, to Real precision.
        FunctionParserBase<Output> fp;
        fp.AddConstant("pi", std::acos(Real(-1)));
        fp.AddConstant("e", std::exp(Real(1)));
        fp.Parse(subexpression, variables);
        fp.Optimize();
        parsers.push_back(fp);

	// If at end, use nextstart=maxSize.  Else start at next
	// character.
        nextstart = (end == std::string::npos) ?
                    std::string::npos : end + 1;
      }

      this->_initialized = true;
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

      // The remaining locations in _spacetime are currently fixed at construction
      // but could potentially be made dynamic
      return parsers[0].Eval(&_spacetime[0]);
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

      // The remaining locations in _spacetime are currently fixed at construction
      // but could potentially be made dynamic
      for (unsigned int i=0; i != size; ++i)
        output(i) = parsers[i].Eval(&_spacetime[0]);
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

      // The remaining locations in _spacetime are currently fixed at construction
      // but could potentially be made dynamic
      return parsers[i].Eval(&_spacetime[0]);
    }

  virtual void init() {}
  virtual void clear() {}
  virtual AutoPtr<FunctionBase<Output> > clone() {
    return AutoPtr<FunctionBase<Output> >
      (new ParsedFunction(_expression, _additional_vars, _additional_vals));
  }

private:
  std::string _expression;
  std::vector<FunctionParserBase<Output> > parsers;
  std::vector<Output> _spacetime;

  // Additional variables/values that can be parsed and handled by the function parser
  const std::vector<std::string>* _additional_vars;
  const std::vector<Output>* _additional_vals;
};

#else

template <typename Output>
class ParsedFunction : public FunctionBase<Output>
{
public:
  ParsedFunction (std::string /* expression */)
    {
      libmesh_not_implemented();
    }

  virtual Output operator() (const Point&,
                             const Real /* time */ = 0)
    { return 0.; }

  virtual void operator() (const Point&,
                           const Real /* time */,
                           DenseVector<Output>& /* output */) {}

  virtual void init() {}
  virtual void clear() {}
  virtual AutoPtr<FunctionBase<Output> > clone() {
    return AutoPtr<FunctionBase<Output> >
      (new ParsedFunction<Output>(""));
  }
};

#endif // LIBMESH_HAVE_FPARSER

#endif // __parsed_function_h__
