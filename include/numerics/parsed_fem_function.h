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

namespace libMesh
{

// ------------------------------------------------------------
// ParsedFEMFunction class definition
template <typename Output=Number>
class ParsedFEMFunction : public FEMFunctionBase<Output>
{
protected:
  // Helper function for initial parsing
  std::size_t find_name (const std::string & varname,
                         const std::string & expr)
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
  void eval_args(const FEMContext& c, const Point& p, const Real time)
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


    // The remaining locations in _spacetime are currently set at construction
    // but could potentially be made dynamic
  }

public:

  /**
   * Constructor.
   */
  explicit
  ParsedFEMFunction (const System& sys,
                     const std::string& expression,
                     const std::vector<std::string>* additional_vars=NULL,
                     const std::vector<Output>* initial_vals=NULL)
    : _sys(sys),
      _expression(expression),
      _n_vars(sys.n_vars()),
      _n_requested_vars(0),
      _n_requested_grad_components(0),
      _n_requested_hess_components(0),
      _need_var(_n_vars, false),
      _need_var_grad(_n_vars*LIBMESH_DIM, false)
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
    , _need_var_hess(_n_vars*LIBMESH_DIM*LIBMESH_DIM, false)
#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES

  {
    std::string variables = "x";
#if LIBMESH_DIM > 1
    variables += ",y";
#endif
#if LIBMESH_DIM > 2
    variables += ",z";
#endif
    variables += ",t";

    for (unsigned int v=0; v != _n_vars; ++v)
      {
        const std::string & varname = sys.variable_name(v);
        std::size_t varname_i = find_name(varname, _expression);

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
        const std::string & varname = sys.variable_name(v);

        for (unsigned int d=0; d != LIBMESH_DIM; ++d)
          {
            std::string gradname = std::string("grad_") +
              "xyz"[d] + '_' + varname;
            std::size_t gradname_i = find_name(gradname, _expression);

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
        const std::string & varname = sys.variable_name(v);

        for (unsigned int d1=0; d1 != LIBMESH_DIM; ++d1)
          for (unsigned int d2=0; d2 != LIBMESH_DIM; ++d2)
            {
              std::string hessname = std::string("hess_") +
                "xyz"[d1] + "xyz"[d2] + '_' + varname;
              std::size_t hessname_i = find_name(hessname, _expression);

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

    _spacetime.resize
      (LIBMESH_DIM + 1 + _n_requested_vars +
       _n_requested_grad_components + _n_requested_hess_components +
       (additional_vars ? additional_vars->size() : 0));

    // If additional vars were passed, append them to the string
    // that we send to the function parser. Also add them to the
    // end of our spacetime vector
    if (additional_vars)
      {
        if (initial_vals)
          std::copy(initial_vals->begin(), initial_vals->end(), std::back_inserter(_initial_vals));

        std::copy(additional_vars->begin(), additional_vars->end(), std::back_inserter(_additional_vars));

        unsigned int offset = LIBMESH_DIM + 1 + _n_requested_vars +
          _n_requested_grad_components + _n_requested_hess_components;

        for (unsigned int i=0; i < additional_vars->size(); ++i)
          {
            variables += "," + (*additional_vars)[i];
            // Initialize extra variables to the vector passed in or zero
            // Note: The initial_vals vector can be shorter than the additional_vars vector
            _spacetime[offset + i] =
              (initial_vals && i < initial_vals->size()) ?
              (*initial_vals)[i] : 0;
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
        if (subexpression.empty())
          libmesh_error_msg("ERROR: FunctionParser is unable to parse empty expression.\n");


#ifdef LIBMESH_HAVE_FPARSER
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
#else
        libmesh_error_msg("ERROR: This functionality requires fparser!");
#endif

        // If at end, use nextstart=maxSize.  Else start at next
        // character.
        nextstart = (end == std::string::npos) ?
          std::string::npos : end + 1;
      }
  }

  /**
   * Destructor.
   */
  virtual ~ParsedFEMFunction () {}

  /**
   * Prepares a context object for use.
   */
  virtual void init_context (const FEMContext &c) {
    for (unsigned int v=0; v != _n_vars; ++v)
      {
        FEBase* elem_fe;
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
  }

  /**
   * Returns a new copy of the function.  The new copy should be as
   * ``deep'' as necessary to allow independent destruction and
   * simultaneous evaluations of the copies in different threads.
   */
  virtual UniquePtr<FEMFunctionBase<Output> > clone () const {
    return UniquePtr<FEMFunctionBase<Output> >
      (new ParsedFEMFunction
       (_sys, _expression, &_additional_vars, &_initial_vals));
  }

  // ------------------------------------------------------
  // misc
  /**
   * @returns the scalar value at coordinate
   * \p p and time \p time, which defaults to zero.
   * Purely virtual, so you have to overload it.
   * Note that this cannot be a const method, check \p MeshFunction.
   */
  virtual Output operator() (const FEMContext& c, const Point& p,
                             const Real time = 0.)
  {
    eval_args(c, p, time);

#if LIBMESH_HAVE_FPARSER
    return parsers[0].Eval(&_spacetime[0]);
#else
    libmesh_error_msg("ERROR: This functionality requires fparser!");
    return Output(0);
#endif
  }



  /**
   * Return function for vectors.
   * Returns in \p output the values of the data at the
   * coordinate \p p and for time \p time.
   */
  void operator() (const FEMContext& c, const Point& p,
                   const Real time,
                   DenseVector<Output>& output)
  {
    eval_args(c, p, time);

    unsigned int size = output.size();

#ifdef LIBMESH_HAVE_FPARSER
    libmesh_assert_equal_to (size, parsers.size());

    // The remaining locations in _spacetime are currently fixed at construction
    // but could potentially be made dynamic
    for (unsigned int i=0; i != size; ++i)
      output(i) = parsers[i].Eval(&_spacetime[0]);
#else
    libmesh_error_msg("ERROR: This functionality requires fparser!");
#endif
  }


  /**
   * @returns the vector component \p i at coordinate
   * \p p and time \p time.
   */
  virtual Output component(const FEMContext& c, unsigned int i,
                           const Point& p,
                           Real time=0.)
  {
    eval_args(c, p, time);

#ifdef LIBMESH_HAVE_FPARSER
    libmesh_assert_less (i, parsers.size());

    // The remaining locations in _spacetime are currently fixed at construction
    // but could potentially be made dynamic
    return parsers[i].Eval(&_spacetime[0]);
#else
    libmesh_error_msg("ERROR: This functionality requires fparser!");
    return Output(0);
#endif
  }

private:
  const System& _sys;
  std::string _expression;
  unsigned int _n_vars,
    _n_requested_vars,
    _n_requested_grad_components,
    _n_requested_hess_components;
#ifdef LIBMESH_HAVE_FPARSER
  std::vector<FunctionParserBase<Output> > parsers;
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

  // Additional variables/values that can be parsed and handled by the function parser
  std::vector<std::string> _additional_vars;
  std::vector<Output> _initial_vals;
};

} // namespace libMesh

#endif // LIBMESH_PARSED_FEM_FUNCTION_H
