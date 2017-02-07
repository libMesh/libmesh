// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/vector_value.h"
#include "libmesh/point.h"

// FParser includes
#include "libmesh/fparser_ad.hh"

// C++ includes
#include <algorithm> // std::find
#include <cmath>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace libMesh
{

/**
 * FEMFunction that returns a single.
 *
 * \author Roy Stogner
 * \date 2012
 * \brief A Function generated (via FParser) by parsing a mathematical expression.
 */
template <typename Output=Number, typename OutputGradient=Gradient>
class ParsedFunction : public FunctionBase<Output>
{
public:
  explicit
  ParsedFunction (const std::string & expression,
                  const std::vector<std::string> * additional_vars=libmesh_nullptr,
                  const std::vector<Output> * initial_vals=libmesh_nullptr);

  // Re-parse with new expression
  void reparse (const std::string & expression);

  virtual Output operator() (const Point & p,
                             const Real time = 0);

  // Query if the automatic derivative generation was successful
  virtual bool has_derivatives() { return _valid_derivatives; }

  virtual Output dot(const Point & p,
                     const Real time = 0);

  virtual OutputGradient gradient(const Point & p,
                                  const Real time = 0);

  virtual void operator() (const Point & p,
                           const Real time,
                           DenseVector<Output> & output);

  /**
   * @returns the vector component \p i at coordinate
   * \p p and time \p time.
   */
  virtual Output component (unsigned int i,
                            const Point & p,
                            Real time);

  const std::string & expression() { return _expression; }

  /**
   * @returns the address of a parsed variable so you can supply a parameterized value
   */
  virtual Output & getVarAddress(const std::string & variable_name);

  virtual UniquePtr<FunctionBase<Output> > clone() const;

  /**
   * @returns the value of an inline variable.  Will *only* be correct
   * if the inline variable value is independent of input variables,
   * if the inline variable is not redefined within any subexpression,
   * and if the inline variable takes the same value within any
   * subexpressions where it appears.
   */
  Output get_inline_value(const std::string & inline_var_name) const;

  /**
   * Changes the value of an inline variable.  Forever after the
   * variable value will take the given constant, independent of input
   * variables, in every subexpression where it is already defined.
   * Currently only works if the inline variable is not redefined
   * within any one subexpression.
   */
  void set_inline_value(const std::string & inline_var_name,
                        Output newval);

protected:
  // Re-parse with minor changes to expression
  void partial_reparse (const std::string & expression);

  // Helper function for parsing out variable names
  std::size_t find_name (const std::string & varname,
                         const std::string & expr) const;

  // Helper function for determining time-dependence of expression
  bool expression_is_time_dependent( const std::string & expression ) const;

private:
  // Set the _spacetime argument vector
  void set_spacetime(const Point & p,
                     const Real time = 0);

  // Evaluate the ith FunctionParser and check the result
  inline Output eval(FunctionParserADBase<Output> & parser,
                     const std::string & libmesh_dbg_var(function_name),
                     unsigned int libmesh_dbg_var(component_idx)) const;

  std::string _expression;
  std::vector<std::string> _subexpressions;
  std::vector<FunctionParserADBase<Output> > parsers;
  std::vector<Output> _spacetime;

  // derivative functions
  std::vector<FunctionParserADBase<Output> > dx_parsers;
#if LIBMESH_DIM > 1
  std::vector<FunctionParserADBase<Output> > dy_parsers;
#endif
#if LIBMESH_DIM > 2
  std::vector<FunctionParserADBase<Output> > dz_parsers;
#endif
  std::vector<FunctionParserADBase<Output> > dt_parsers;
  bool _valid_derivatives;

  // Variables/values that can be parsed and handled by the function parser
  std::string variables;
  std::vector<std::string> _additional_vars;
  std::vector<Output> _initial_vals;
};


/*----------------------- Inline functions ----------------------------------*/

template <typename Output, typename OutputGradient>
inline
ParsedFunction<Output,OutputGradient>::ParsedFunction (const std::string & expression,
                                                       const std::vector<std::string> * additional_vars,
                                                       const std::vector<Output> * initial_vals) :
  _expression (), // overridden by parse()
  // Size the spacetime vector to account for space, time, and any additional
  // variables passed
  _spacetime (LIBMESH_DIM+1 + (additional_vars ? additional_vars->size() : 0)),
  _valid_derivatives (true),
  _additional_vars (additional_vars ? *additional_vars : std::vector<std::string>()),
  _initial_vals (initial_vals ? *initial_vals : std::vector<Output>())
{
  // time-dependence established in reparse function
  this->reparse(expression);
  this->_initialized = true;
}


template <typename Output, typename OutputGradient>
inline
void
ParsedFunction<Output,OutputGradient>::reparse (const std::string & expression)
{
  variables = "x";
#if LIBMESH_DIM > 1
  variables += ",y";
#endif
#if LIBMESH_DIM > 2
  variables += ",z";
#endif
  variables += ",t";

  // If additional vars were passed, append them to the string
  // that we send to the function parser. Also add them to the
  // end of our spacetime vector
  for (std::size_t i=0; i < _additional_vars.size(); ++i)
    {
      variables += "," + _additional_vars[i];
      // Initialize extra variables to the vector passed in or zero
      // Note: The initial_vals vector can be shorter than the additional_vars vector
      _spacetime[LIBMESH_DIM+1 + i] =
        (i < _initial_vals.size()) ? _initial_vals[i] : 0;
    }

  this->partial_reparse(expression);

  this->_is_time_dependent = this->expression_is_time_dependent(expression);
}

template <typename Output, typename OutputGradient>
inline
Output
ParsedFunction<Output,OutputGradient>::operator() (const Point & p, const Real time)
{
  set_spacetime(p, time);
  return eval(parsers[0], "f", 0);
}

template <typename Output, typename OutputGradient>
inline
Output
ParsedFunction<Output,OutputGradient>::dot (const Point & p, const Real time)
{
  set_spacetime(p, time);
  return eval(dt_parsers[0], "df/dt", 0);
}

template <typename Output, typename OutputGradient>
inline
OutputGradient
ParsedFunction<Output,OutputGradient>::gradient (const Point & p, const Real time)
{
  OutputGradient grad;
  set_spacetime(p, time);

  grad(0) = eval(dx_parsers[0], "df/dx", 0);
#if LIBMESH_DIM > 1
  grad(1) = eval(dy_parsers[0], "df/dy", 0);
#endif
#if LIBMESH_DIM > 2
  grad(2) = eval(dz_parsers[0], "df/dz", 0);
#endif

  return grad;
}

template <typename Output, typename OutputGradient>
inline
void
ParsedFunction<Output,OutputGradient>::operator()
  (const Point & p,
   const Real time,
   DenseVector<Output> & output)
{
  set_spacetime(p, time);

  unsigned int size = output.size();

  libmesh_assert_equal_to (size, parsers.size());

  // The remaining locations in _spacetime are currently fixed at construction
  // but could potentially be made dynamic
  for (unsigned int i=0; i != size; ++i)
    output(i) = eval(parsers[i], "f", i);
}

/**
 * @returns the vector component \p i at coordinate
 * \p p and time \p time.
 */
template <typename Output, typename OutputGradient>
inline
Output
ParsedFunction<Output,OutputGradient>::component (unsigned int i,
                                                  const Point & p,
                                                  Real time)
{
  set_spacetime(p, time);
  libmesh_assert_less (i, parsers.size());

  // The remaining locations in _spacetime are currently fixed at construction
  // but could potentially be made dynamic
  libmesh_assert_less(i, parsers.size());
  return eval(parsers[i], "f", i);
}

/**
 * @returns the address of a parsed variable so you can supply a parameterized value
 */
template <typename Output, typename OutputGradient>
inline
Output &
ParsedFunction<Output,OutputGradient>::getVarAddress (const std::string & variable_name)
{
  const std::vector<std::string>::iterator it =
    std::find(_additional_vars.begin(), _additional_vars.end(), variable_name);

  if (it == _additional_vars.end())
    libmesh_error_msg("ERROR: Requested variable not found in parsed function");

  // Iterator Arithmetic (How far from the end of the array is our target address?)
  return _spacetime[_spacetime.size() - (_additional_vars.end() - it)];
}


template <typename Output, typename OutputGradient>
inline
UniquePtr<FunctionBase<Output> >
ParsedFunction<Output,OutputGradient>::clone() const
{
  return UniquePtr<FunctionBase<Output> >
    (new ParsedFunction(_expression, &_additional_vars, &_initial_vals));
}

template <typename Output, typename OutputGradient>
inline
Output
ParsedFunction<Output,OutputGradient>::get_inline_value (const std::string & inline_var_name) const
{
  libmesh_assert_greater (_subexpressions.size(), 0);

#ifndef NDEBUG
  bool found_var_name = false;
#endif
  Output old_var_value(0.);

  for (std::size_t s=0; s != _subexpressions.size(); ++s)
    {
      const std::string & subexpression = _subexpressions[s];
      const std::size_t varname_i =
        find_name(inline_var_name, subexpression);
      if (varname_i == std::string::npos)
        continue;

      const std::size_t assignment_i =
        subexpression.find(":", varname_i+1);

      libmesh_assert_not_equal_to(assignment_i, std::string::npos);

      libmesh_assert_equal_to(subexpression[assignment_i+1], '=');
      for (unsigned int i = varname_i+1; i != assignment_i; ++i)
        libmesh_assert_equal_to(subexpression[i], ' ');

      std::size_t end_assignment_i =
        subexpression.find(";", assignment_i+1);

      libmesh_assert_not_equal_to(end_assignment_i, std::string::npos);

      std::string new_subexpression =
        subexpression.substr(0, end_assignment_i+1) +
        inline_var_name;

#ifdef LIBMESH_HAVE_FPARSER
      // Parse and evaluate the new subexpression.
      // Add the same constants as we used originally.
      FunctionParserADBase<Output> fp;
      fp.AddConstant("NaN", std::numeric_limits<Real>::quiet_NaN());
      fp.AddConstant("pi", std::acos(Real(-1)));
      fp.AddConstant("e", std::exp(Real(1)));
      if (fp.Parse(new_subexpression, variables) != -1) // -1 for success
        libmesh_error_msg
          ("ERROR: FunctionParser is unable to parse modified expression: "
           << new_subexpression << '\n' << fp.ErrorMsg());

      Output new_var_value = this->eval(fp, new_subexpression, 0);
#ifdef NDEBUG
      return new_var_value;
#else
      if (found_var_name)
        {
          libmesh_assert_equal_to(old_var_value, new_var_value);
        }
      else
        {
          old_var_value = new_var_value;
          found_var_name = true;
        }
#endif

#else
      libmesh_error_msg("ERROR: This functionality requires fparser!");
#endif
    }

  libmesh_assert(found_var_name);
  return old_var_value;
}


template <typename Output, typename OutputGradient>
inline
void
ParsedFunction<Output,OutputGradient>::set_inline_value (const std::string & inline_var_name,
                                                         Output newval)
{
  libmesh_assert_greater (_subexpressions.size(), 0);

#ifndef NDEBUG
  bool found_var_name = false;
#endif
  for (std::size_t s=0; s != _subexpressions.size(); ++s)
    {
      const std::string & subexpression = _subexpressions[s];
      const std::size_t varname_i =
        find_name(inline_var_name, subexpression);
      if (varname_i == std::string::npos)
        continue;

#ifndef NDEBUG
      found_var_name = true;
#endif
      const std::size_t assignment_i =
        subexpression.find(":", varname_i+1);

      libmesh_assert_not_equal_to(assignment_i, std::string::npos);

      libmesh_assert_equal_to(subexpression[assignment_i+1], '=');
      for (unsigned int i = varname_i+1; i != assignment_i; ++i)
        libmesh_assert_equal_to(subexpression[i], ' ');

      std::size_t end_assignment_i =
        subexpression.find(";", assignment_i+1);

      libmesh_assert_not_equal_to(end_assignment_i, std::string::npos);

      std::ostringstream new_subexpression;
      new_subexpression << subexpression.substr(0, assignment_i+2)
                        << std::setprecision(std::numeric_limits<Output>::digits10+2)
#ifdef LIBMESH_USE_COMPLEX_NUMBERS
                        << '(' << newval.real() << '+'
                        << newval.imag() << 'i' << ')'
#else
                        << newval
#endif
                        << subexpression.substr(end_assignment_i,
                                                std::string::npos);
      _subexpressions[s] = new_subexpression.str();
    }

  libmesh_assert(found_var_name);

  std::string new_expression;

  for (std::size_t s=0; s != _subexpressions.size(); ++s)
    {
      new_expression += '{';
      new_expression += _subexpressions[s];
      new_expression += '}';
    }

  this->partial_reparse(new_expression);
}


template <typename Output, typename OutputGradient>
inline
void
ParsedFunction<Output,OutputGradient>::partial_reparse (const std::string & expression)
{
  _expression = expression;
  _subexpressions.clear();
  parsers.clear();

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
      _subexpressions.push_back
        (expression.substr(nextstart, (end == std::string::npos) ?
                           std::string::npos : end - nextstart));

      // fparser can crash on empty expressions
      if (_subexpressions.back().empty())
        libmesh_error_msg("ERROR: FunctionParser is unable to parse empty expression.\n");

      // Parse (and optimize if possible) the subexpression.
      // Add some basic constants, to Real precision.
      FunctionParserADBase<Output> fp;
      fp.AddConstant("NaN", std::numeric_limits<Real>::quiet_NaN());
      fp.AddConstant("pi", std::acos(Real(-1)));
      fp.AddConstant("e", std::exp(Real(1)));
      if (fp.Parse(_subexpressions.back(), variables) != -1) // -1 for success
        libmesh_error_msg
          ("ERROR: FunctionParser is unable to parse expression: "
           << _subexpressions.back() << '\n' << fp.ErrorMsg());

      // use of derivatives is optional. suppress error output on the console
      // use the has_derivatives() method to check if AutoDiff was successful.
      // also enable immediate optimization
      fp.SetADFlags(FunctionParserADBase<Output>::ADSilenceErrors |
                    FunctionParserADBase<Output>::ADAutoOptimize);

      // optimize original function
      fp.Optimize();
      parsers.push_back(fp);

      // generate derivatives through automatic differentiation
      FunctionParserADBase<Output> dx_fp(fp);
      if (dx_fp.AutoDiff("x") != -1) // -1 for success
        _valid_derivatives = false;
      dx_parsers.push_back(dx_fp);
#if LIBMESH_DIM > 1
      FunctionParserADBase<Output> dy_fp(fp);
      if (dy_fp.AutoDiff("y") != -1) // -1 for success
        _valid_derivatives = false;
      dy_parsers.push_back(dy_fp);
#endif
#if LIBMESH_DIM > 2
      FunctionParserADBase<Output> dz_fp(fp);
      if (dz_fp.AutoDiff("z") != -1) // -1 for success
        _valid_derivatives = false;
      dz_parsers.push_back(dz_fp);
#endif
      FunctionParserADBase<Output> dt_fp(fp);
      if (dt_fp.AutoDiff("t") != -1) // -1 for success
        _valid_derivatives = false;
      dt_parsers.push_back(dt_fp);

      // If at end, use nextstart=maxSize.  Else start at next
      // character.
      nextstart = (end == std::string::npos) ?
        std::string::npos : end + 1;
    }
}


template <typename Output, typename OutputGradient>
inline
std::size_t
ParsedFunction<Output,OutputGradient>::find_name (const std::string & varname,
                                                  const std::string & expr) const
{
  const std::size_t namesize = varname.size();
  std::size_t varname_i = expr.find(varname);

  while ((varname_i != std::string::npos) &&
         (((varname_i > 0) &&
           (std::isalnum(expr[varname_i-1]) ||
            (expr[varname_i-1] == '_'))) ||
          ((varname_i+namesize < expr.size()) &&
           (std::isalnum(expr[varname_i+namesize]) ||
            (expr[varname_i+namesize] == '_')))))
    {
      varname_i = expr.find(varname, varname_i+1);
    }

  return varname_i;
}
template <typename Output, typename OutputGradient>
inline
bool
ParsedFunction<Output,OutputGradient>::expression_is_time_dependent( const std::string & expression ) const
{
  bool is_time_dependent = false;

  // By definition, time is "t" for FunctionBase-based objects, so we just need to
  // see if this expression has the variable "t" in it.
  if( this->find_name( std::string("t"), expression ) != std::string::npos )
    is_time_dependent = true;

  return is_time_dependent;
}

// Set the _spacetime argument vector
template <typename Output, typename OutputGradient>
inline
void
ParsedFunction<Output,OutputGradient>::set_spacetime (const Point & p,
                                                      const Real time)
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
}

// Evaluate the ith FunctionParser and check the result
template <typename Output, typename OutputGradient>
inline
Output
ParsedFunction<Output,OutputGradient>::eval (FunctionParserADBase<Output> & parser,
                                             const std::string & libmesh_dbg_var(function_name),
                                             unsigned int libmesh_dbg_var(component_idx)) const
{
#ifndef NDEBUG
  Output result = parser.Eval(&_spacetime[0]);
  int error_code = parser.EvalError();
  if (error_code)
    {
      libMesh::err << "ERROR: FunctionParser is unable to evaluate component "
                   << component_idx
                   << " of expression '"
                   << function_name
                   << "' with arguments:\n";
      for (std::size_t j=0; j<_spacetime.size(); ++j)
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
#else
  return parser.Eval(&_spacetime[0]);
#endif
}

} // namespace libMesh


#else // !LIBMESH_HAVE_FPARSER


namespace libMesh {


template <typename Output=Number>
class ParsedFunction : public FunctionBase<Output>
{
public:
  ParsedFunction (const std::string & /* expression */,
                  const std::vector<std::string> * = libmesh_nullptr,
                  const std::vector<Output> * = libmesh_nullptr) : _dummy(0)
  {
    libmesh_not_implemented();
  }

  virtual Output operator() (const Point &,
                             const Real /* time */ = 0)
  { return 0.; }

  virtual void operator() (const Point &,
                           const Real /* time */,
                           DenseVector<Output> & /* output */) {}

  virtual void init() {}
  virtual void clear() {}
  virtual Output & getVarAddress(const std::string & /*variable_name*/) { return _dummy; }
  virtual UniquePtr<FunctionBase<Output> > clone() const
  {
    return UniquePtr<FunctionBase<Output> > (new ParsedFunction<Output>(""));
  }
private:
  Output _dummy;
};



} // namespace libMesh


#endif // LIBMESH_HAVE_FPARSER

#endif // LIBMESH_PARSED_FUNCTION_H
