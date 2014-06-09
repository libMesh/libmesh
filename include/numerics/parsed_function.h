// The libMesh Finite Element Library.
// Copyright (C) 2002-2014 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_PARSED_FUNCTION_H
#define LIBMESH_PARSED_FUNCTION_H

#include "libmesh/libmesh_config.h"
#include "libmesh/function_base.h"

#ifdef LIBMESH_HAVE_FPARSER

// Local includes
#include "libmesh/dense_vector.h"
#include "libmesh/point.h"

// FParser includes
#include "libmesh/fparser.hh"

// C++ includes
#include <algorithm> // std::find
#include <cmath>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>

namespace libMesh {

template <typename Output=Number>
class ParsedFunction : public FunctionBase<Output>
{
public:
  explicit
  ParsedFunction (const std::string& expression, const std::vector<std::string>* additional_vars=NULL,
                  const std::vector<Output>* initial_vals=NULL)
    : _expression(expression)
      // Size the spacetime vector to account for space, time, and any additional
      // variables passed
      //_spacetime(LIBMESH_DIM+1 + (additional_vars ? additional_vars->size() : 0)),
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
        if (initial_vals)
          std::copy(initial_vals->begin(), initial_vals->end(), std::back_inserter(_initial_vals));

        std::copy(additional_vars->begin(), additional_vars->end(), std::back_inserter(_additional_vars));

        for (unsigned int i=0; i < additional_vars->size(); ++i)
          {
            variables += "," + (*additional_vars)[i];
            // Initialize extra variables to the vector passed in or zero
            // Note: The initial_vals vector can be shorter than the additional_vars vector
            _spacetime[LIBMESH_DIM+1 + i] = (initial_vals && i < initial_vals->size()) ? (*initial_vals)[i] : 0;
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
        fp.AddConstant("NaN", std::numeric_limits<Real>::quiet_NaN());
        fp.AddConstant("pi", std::acos(Real(-1)));
        fp.AddConstant("e", std::exp(Real(1)));
        if (fp.Parse(subexpression, variables) != -1) // -1 for success
          libmesh_error_msg("ERROR: FunctionParser is unable to parse expression: " << subexpression << '\n' << fp.ErrorMsg());

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
    return eval(0);
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

    libmesh_assert_equal_to (size, parsers.size());

    // The remaining locations in _spacetime are currently fixed at construction
    // but could potentially be made dynamic
    for (unsigned int i=0; i != size; ++i)
      output(i) = eval(i);
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

    libmesh_assert_less (i, parsers.size());

    // The remaining locations in _spacetime are currently fixed at construction
    // but could potentially be made dynamic
    return eval(i);
  }

  /**
   * @returns the address of a parsed variable so you can supply a parameterized value
   */
  virtual Output & getVarAddress(const std::string & variable_name)
  {
    const std::vector<std::string>::iterator it =
      std::find(_additional_vars.begin(), _additional_vars.end(), variable_name);

    if (it == _additional_vars.end())
      libmesh_error_msg("ERROR: Requested variable not found in parsed function");

    // Iterator Arithmetic (How far from the end of the array is our target address?)
    return _spacetime[_spacetime.size() - (_additional_vars.end() - it)];
  }


  virtual AutoPtr<FunctionBase<Output> > clone() const {
    return AutoPtr<FunctionBase<Output> >
      (new ParsedFunction(_expression, &_additional_vars, &_initial_vals));
  }

private:
  // Evaluate the ith FunctionParser and check the result
  inline Output eval(unsigned int i)
  {
#ifndef DEBUG
    return parsers[i].Eval(&_spacetime[0]);
#else
    libmesh_assert_less(i, parsers.size());

    Output result = parsers[i].Eval(&_spacetime[0]);
    int error_code = parsers[i].EvalError();
    if (error_code)
    {
      libMesh::err << "ERROR: FunctionParser is unable to evaluate the expression at index " << i << " with arguments:\n";
      for (unsigned int j=0; j<_spacetime.size(); ++j)
        libMesh::err << '\t' << _spacetime[j] << '\n';
      libMesh::err << '\n';

      // Currently no API to report error messages, we'll do it manually
      std::string error_message = "Reason: ";

      switch (error_code)
      {
      case 1:
        error_message += "Division by zero";
        break;
      case 2:
        error_message += "Square Root error (negative value)";
        break;
      case 3:
        error_message += "Log error (negative value)";
        break;
      case 4:
        error_message += "Trigonometric error (asin or acos of illegal value)";
        break;
      case 5:
        error_message += "Maximum recursion level reached";
        break;
      default:
        error_message += "Unknown";
        break;
      }
      libmesh_error_msg(error_message);
    }

    return result;
#endif
  }

  std::string _expression;
  std::vector<FunctionParserBase<Output> > parsers;
  std::vector<Output> _spacetime;

  // Additional variables/values that can be parsed and handled by the function parser
  std::vector<std::string> _additional_vars;
  std::vector<Output> _initial_vals;
};


} // namespace libMesh


#else // LIBMESH_HAVE_FPARSER


namespace libMesh {


template <typename Output>
class ParsedFunction : public FunctionBase<Output>
{
public:
  ParsedFunction (std::string /* expression */) : _dummy(0)
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
  virtual Output & getVarAddress(const std::string & /*variable_name*/) { return _dummy; }
  virtual AutoPtr<FunctionBase<Output> > clone() const {
    return AutoPtr<FunctionBase<Output> >
      (new ParsedFunction<Output>(""));
  }
private:
  Output _dummy;
};


} // namespace libMesh


#endif // LIBMESH_HAVE_FPARSER

#endif // LIBMESH_PARSED_FUNCTION_H
