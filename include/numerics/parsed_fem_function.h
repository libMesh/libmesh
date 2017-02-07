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



#ifndef LIBMESH_PARSED_FEM_FUNCTION_H
#define LIBMESH_PARSED_FEM_FUNCTION_H

#include "libmesh/libmesh_config.h"

// Local Includes
#include "libmesh/fem_function_base.h"
#include "libmesh/point.h"
#include "libmesh/system.h"

#ifdef LIBMESH_HAVE_FPARSER
// FParser includes
#include "libmesh/fparser.hh"
#endif

// C++ includes
#include <cctype>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

namespace libMesh
{

/**
 * ParsedFEMFunction provides support for FParser-based parsed
 * functions in FEMSystem.
 *
 * \author Roy Stogner
 * \date 2014
 * \brief Support for using parsed functions in FEMSystem.
 */
template <typename Output=Number>
class ParsedFEMFunction : public FEMFunctionBase<Output>
{
public:

  /**
   * Constructor.
   */
  explicit
  ParsedFEMFunction (const System & sys,
                     const std::string & expression,
                     const std::vector<std::string> * additional_vars=libmesh_nullptr,
                     const std::vector<Output> * initial_vals=libmesh_nullptr);

  /**
   * Destructor.
   */
  virtual ~ParsedFEMFunction () {}

  // Re-parse with new expression
  void reparse (const std::string & expression);

  /**
   * Prepares a context object for use.
   */
  virtual void init_context (const FEMContext & c) libmesh_override;

  /**
   * Returns a new copy of the function.  The new copy should be as
   * ``deep'' as necessary to allow independent destruction and
   * simultaneous evaluations of the copies in different threads.
   */
  virtual UniquePtr<FEMFunctionBase<Output> > clone () const libmesh_override;

  // ------------------------------------------------------
  // misc
  /**
   * @returns the scalar value at coordinate
   * \p p and time \p time, which defaults to zero.
   * Purely virtual, so you have to overload it.
   * Note that this cannot be a const method, check \p MeshFunction.
   */
  virtual Output operator() (const FEMContext & c,
                             const Point & p,
                             const Real time = 0.) libmesh_override;

  /**
   * Return function for vectors.
   * Returns in \p output the values of the data at the
   * coordinate \p p and for time \p time.
   */
  void operator() (const FEMContext & c,
                   const Point & p,
                   const Real time,
                   DenseVector<Output> & output) libmesh_override;

  /**
   * @returns the vector component \p i at coordinate
   * \p p and time \p time.
   */
  virtual Output component(const FEMContext & c,
                           unsigned int i,
                           const Point & p,
                           Real time=0.) libmesh_override;

  const std::string & expression() { return _expression; }

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
  // Helper function for reparsing minor changes to expression
  void partial_reparse (const std::string & expression);

  // Helper function for parsing out variable names
  std::size_t find_name (const std::string & varname,
                         const std::string & expr) const;

  // Helper function for evaluating function arguments
  void eval_args(const FEMContext & c,
                 const Point & p,
                 const Real time);

  // Evaluate the ith FunctionParser and check the result
#ifdef LIBMESH_HAVE_FPARSER
  inline Output eval(FunctionParserBase<Output> & parser,
                     const std::string & libmesh_dbg_var(function_name),
                     unsigned int libmesh_dbg_var(component_idx)) const;
#else // LIBMESH_HAVE_FPARSER
  inline Output eval(char & libmesh_dbg_var(parser),
                     const std::string & libmesh_dbg_var(function_name),
                     unsigned int libmesh_dbg_var(component_idx)) const;
#endif

private:
  const System & _sys;
  std::string _expression;
  std::vector<std::string> _subexpressions;
  unsigned int _n_vars,
    _n_requested_vars,
    _n_requested_grad_components,
    _n_requested_hess_components;
  bool _requested_normals;
#ifdef LIBMESH_HAVE_FPARSER
  std::vector<FunctionParserBase<Output> > parsers;
#else
  std::vector<char> parsers;
#endif
  std::vector<Output> _spacetime;

  // Flags for which variables need to be computed

  // _need_var[v] is true iff value(v) is needed
  std::vector<bool> _need_var;

  // _need_var_grad[v*LIBMESH_DIM+i] is true iff grad(v,i) is needed
  std::vector<bool> _need_var_grad;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  // _need_var_hess[v*LIBMESH_DIM*LIBMESH_DIM+i*LIBMESH_DIM+j] is true
  // iff grad(v,i,j) is needed
  std::vector<bool> _need_var_hess;
#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES

  // Variables/values that can be parsed and handled by the function parser
  std::string variables;
  std::vector<std::string> _additional_vars;
  std::vector<Output> _initial_vals;
};


/*----------------------- Inline functions ----------------------------------*/

template <typename Output>
inline
ParsedFEMFunction<Output>::ParsedFEMFunction (const System & sys,
                                              const std::string & expression,
                                              const std::vector<std::string> * additional_vars,
                                              const std::vector<Output> * initial_vals) :
  _sys(sys),
  _expression (), // overridden by parse()
  _n_vars(sys.n_vars()),
  _n_requested_vars(0),
  _n_requested_grad_components(0),
  _n_requested_hess_components(0),
  _requested_normals(false),
  _need_var(_n_vars, false),
  _need_var_grad(_n_vars*LIBMESH_DIM, false),
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  _need_var_hess(_n_vars*LIBMESH_DIM*LIBMESH_DIM, false),
#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES
  _additional_vars (additional_vars ? *additional_vars : std::vector<std::string>()),
  _initial_vals (initial_vals ? *initial_vals : std::vector<Output>())
{
  this->reparse(expression);
}


template <typename Output>
inline
void
ParsedFEMFunction<Output>::reparse (const std::string & expression)
{
  variables = "x";
#if LIBMESH_DIM > 1
  variables += ",y";
#endif
#if LIBMESH_DIM > 2
  variables += ",z";
#endif
  variables += ",t";

  for (unsigned int v=0; v != _n_vars; ++v)
    {
      const std::string & varname = _sys.variable_name(v);
      std::size_t varname_i = find_name(varname, expression);

      // If we didn't find our variable name then let's go to the
      // next.
      if (varname_i == std::string::npos)
        continue;

      _need_var[v] = true;
      variables += ',';
      variables += varname;
      _n_requested_vars++;
    }

  for (unsigned int v=0; v != _n_vars; ++v)
    {
      const std::string & varname = _sys.variable_name(v);

      for (unsigned int d=0; d != LIBMESH_DIM; ++d)
        {
          std::string gradname = std::string("grad_") +
            "xyz"[d] + '_' + varname;
          std::size_t gradname_i = find_name(gradname, expression);

          // If we didn't find that gradient component of our
          // variable name then let's go to the next.
          if (gradname_i == std::string::npos)
            continue;

          _need_var_grad[v*LIBMESH_DIM+d] = true;
          variables += ',';
          variables += gradname;
          _n_requested_grad_components++;
        }
    }

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  for (unsigned int v=0; v != _n_vars; ++v)
    {
      const std::string & varname = _sys.variable_name(v);

      for (unsigned int d1=0; d1 != LIBMESH_DIM; ++d1)
        for (unsigned int d2=0; d2 != LIBMESH_DIM; ++d2)
          {
            std::string hessname = std::string("hess_") +
              "xyz"[d1] + "xyz"[d2] + '_' + varname;
            std::size_t hessname_i = find_name(hessname, expression);

            // If we didn't find that hessian component of our
            // variable name then let's go to the next.
            if (hessname_i == std::string::npos)
              continue;

            _need_var_hess[v*LIBMESH_DIM*LIBMESH_DIM+d1*LIBMESH_DIM+d2]
              = true;
            variables += ',';
            variables += hessname;
            _n_requested_hess_components++;
          }
    }
#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES

  {
    std::size_t nx_i = find_name("n_x", expression);
    std::size_t ny_i = find_name("n_y", expression);
    std::size_t nz_i = find_name("n_z", expression);

    // If we found any requests for normal components then we'll
    // compute normals
    if (nx_i != std::string::npos ||
        ny_i != std::string::npos ||
        nz_i != std::string::npos)
      {
        _requested_normals = true;
        variables += ',';
        variables += "n_x";
        if (LIBMESH_DIM > 1)
          {
            variables += ',';
            variables += "n_y";
          }
        if (LIBMESH_DIM > 2)
          {
            variables += ',';
            variables += "n_z";
          }
      }
  }

  _spacetime.resize
    (LIBMESH_DIM + 1 + _n_requested_vars +
     _n_requested_grad_components + _n_requested_hess_components +
     (_requested_normals ? LIBMESH_DIM : 0) +
     _additional_vars.size());

  // If additional vars were passed, append them to the string
  // that we send to the function parser. Also add them to the
  // end of our spacetime vector
  unsigned int offset = LIBMESH_DIM + 1 + _n_requested_vars +
    _n_requested_grad_components + _n_requested_hess_components;

  for (std::size_t i=0; i < _additional_vars.size(); ++i)
    {
      variables += "," + _additional_vars[i];
      // Initialize extra variables to the vector passed in or zero
      // Note: The initial_vals vector can be shorter than the additional_vars vector
      _spacetime[offset + i] =
        (i < _initial_vals.size()) ? _initial_vals[i] : 0;
    }

  this->partial_reparse(expression);
}

template <typename Output>
inline
void
ParsedFEMFunction<Output>::init_context (const FEMContext & c)
{
  for (unsigned int v=0; v != _n_vars; ++v)
    {
      FEBase * elem_fe;
      c.get_element_fe(v, elem_fe);
      if (_n_requested_vars)
        elem_fe->get_phi();
      if (_n_requested_grad_components)
        elem_fe->get_dphi();
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
      if (_n_requested_hess_components)
        elem_fe->get_d2phi();
#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES
    }

  if (_requested_normals)
    {
      FEBase * side_fe;
      c.get_side_fe(0, side_fe);

      side_fe->get_normals();

      // FIXME: this is a hack to support normals at quadrature
      // points; we don't support normals elsewhere.
      side_fe->get_xyz();
    }
}

template <typename Output>
inline
UniquePtr<FEMFunctionBase<Output> >
ParsedFEMFunction<Output>::clone () const
{
  return UniquePtr<FEMFunctionBase<Output> >
    (new ParsedFEMFunction(_sys, _expression, &_additional_vars, &_initial_vals));
}

template <typename Output>
inline
Output
ParsedFEMFunction<Output>::operator() (const FEMContext & c,
                                       const Point & p,
                                       const Real time)
{
  eval_args(c, p, time);

  return eval(parsers[0], "f", 0);
}



template <typename Output>
inline
void
ParsedFEMFunction<Output>::operator() (const FEMContext & c,
                                       const Point & p,
                                       const Real time,
                                       DenseVector<Output> & output)
{
  eval_args(c, p, time);

  unsigned int size = output.size();

  libmesh_assert_equal_to (size, parsers.size());

  for (unsigned int i=0; i != size; ++i)
    output(i) = eval(parsers[i], "f", i);
}


template <typename Output>
inline
Output
ParsedFEMFunction<Output>::component (const FEMContext & c,
                                      unsigned int i,
                                      const Point & p,
                                      Real time)
{
  eval_args(c, p, time);

  libmesh_assert_less (i, parsers.size());
  return eval(parsers[i], "f", i);
}

template <typename Output>
inline
Output
ParsedFEMFunction<Output>::get_inline_value(const std::string & inline_var_name) const
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
      FunctionParserBase<Output> fp;
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

template <typename Output>
inline
void
ParsedFEMFunction<Output>::set_inline_value (const std::string & inline_var_name,
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


template <typename Output>
inline
void
ParsedFEMFunction<Output>::partial_reparse (const std::string & expression)
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


#ifdef LIBMESH_HAVE_FPARSER
      // Parse (and optimize if possible) the subexpression.
      // Add some basic constants, to Real precision.
      FunctionParserBase<Output> fp;
      fp.AddConstant("NaN", std::numeric_limits<Real>::quiet_NaN());
      fp.AddConstant("pi", std::acos(Real(-1)));
      fp.AddConstant("e", std::exp(Real(1)));
      if (fp.Parse(_subexpressions.back(), variables) != -1) // -1 for success
        libmesh_error_msg
          ("ERROR: FunctionParser is unable to parse expression: "
           << _subexpressions.back() << '\n' << fp.ErrorMsg());
      fp.Optimize();
      parsers.push_back(fp);
#else
      libmesh_error_msg("ERROR: This functionality requires fparser!");
#endif

      // If at end, use nextstart=maxSize.  Else start at next
      // character.
      nextstart = (end == std::string::npos) ?
        std::string::npos : end + 1;
    }
}

template <typename Output>
inline
std::size_t
ParsedFEMFunction<Output>::find_name (const std::string & varname,
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



// Helper function for evaluating function arguments
template <typename Output>
inline
void
ParsedFEMFunction<Output>::eval_args (const FEMContext & c,
                                      const Point & p,
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

  unsigned int request_index = 0;
  for (unsigned int v=0; v != _n_vars; ++v)
    {
      if (!_need_var[v])
        continue;

      c.point_value(v, p, _spacetime[LIBMESH_DIM+1+request_index]);
      request_index++;
    }

  if (_n_requested_grad_components)
    for (unsigned int v=0; v != _n_vars; ++v)
      {
        if (!_need_var_grad[v*LIBMESH_DIM]
#if LIBMESH_DIM > 1
            && !_need_var_grad[v*LIBMESH_DIM+1]
#if LIBMESH_DIM > 2
            && !_need_var_grad[v*LIBMESH_DIM+2]
#endif
#endif
            )
          continue;

        Gradient g;
        c.point_gradient(v, p, g);

        for (unsigned int d=0; d != LIBMESH_DIM; ++d)
          {
            if (!_need_var_grad[v*LIBMESH_DIM+d])
              continue;

            _spacetime[LIBMESH_DIM+1+request_index] = g(d);
            request_index++;
          }
      }

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  if (_n_requested_hess_components)
    for (unsigned int v=0; v != _n_vars; ++v)
      {
        if (!_need_var_hess[v*LIBMESH_DIM*LIBMESH_DIM]
#if LIBMESH_DIM > 1
            && !_need_var_hess[v*LIBMESH_DIM*LIBMESH_DIM+1]
            && !_need_var_hess[v*LIBMESH_DIM*LIBMESH_DIM+2]
            && !_need_var_hess[v*LIBMESH_DIM*LIBMESH_DIM+3]
#if LIBMESH_DIM > 2
            && !_need_var_hess[v*LIBMESH_DIM*LIBMESH_DIM+4]
            && !_need_var_hess[v*LIBMESH_DIM*LIBMESH_DIM+5]
            && !_need_var_hess[v*LIBMESH_DIM*LIBMESH_DIM+6]
            && !_need_var_hess[v*LIBMESH_DIM*LIBMESH_DIM+7]
            && !_need_var_hess[v*LIBMESH_DIM*LIBMESH_DIM+8]
#endif
#endif
            )
          continue;

        Tensor h;
        c.point_hessian(v, p, h);

        for (unsigned int d1=0; d1 != LIBMESH_DIM; ++d1)
          for (unsigned int d2=0; d2 != LIBMESH_DIM; ++d2)
            {
              if (!_need_var_hess[v*LIBMESH_DIM*LIBMESH_DIM+d1*LIBMESH_DIM+d2])
                continue;

              _spacetime[LIBMESH_DIM+1+request_index] = h(d1,d2);
              request_index++;
            }
      }
#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES

  if (_requested_normals)
    {
      FEBase * side_fe;
      c.get_side_fe(0, side_fe);

      const std::vector<Point> & normals = side_fe->get_normals();

      const std::vector<Point> & xyz = side_fe->get_xyz();

      libmesh_assert_equal_to(normals.size(), xyz.size());

      // We currently only support normals at quadrature points!
#ifndef NDEBUG
      bool at_quadrature_point = false;
#endif
      for (std::size_t qp = 0; qp != normals.size(); ++qp)
        {
          if (p == xyz[qp])
            {
              const Point & n = normals[qp];
              for (unsigned int d=0; d != LIBMESH_DIM; ++d)
                {
                  _spacetime[LIBMESH_DIM+1+request_index] = n(d);
                  request_index++;
                }
#ifndef NDEBUG
              at_quadrature_point = true;
#endif
              break;
            }
        }

      libmesh_assert(at_quadrature_point);
    }

  // The remaining locations in _spacetime are currently set at construction
  // but could potentially be made dynamic
}


// Evaluate the ith FunctionParser and check the result
#ifdef LIBMESH_HAVE_FPARSER
template <typename Output>
inline
Output
ParsedFEMFunction<Output>::eval (FunctionParserBase<Output> & parser,
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
#else // LIBMESH_HAVE_FPARSER
template <typename Output>
inline
Output
ParsedFEMFunction<Output>::eval (char & libmesh_dbg_var(parser),
                                 const std::string & libmesh_dbg_var(function_name),
                                 unsigned int libmesh_dbg_var(component_idx)) const
{
  libmesh_error_msg("ERROR: This functionality requires fparser!");
  return Output(0);
}
#endif


} // namespace libMesh

#endif // LIBMESH_PARSED_FEM_FUNCTION_H
