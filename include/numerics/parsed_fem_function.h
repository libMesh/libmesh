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

namespace libMesh
{

// ------------------------------------------------------------
// ParsedFEMFunction class definition
template <typename Output=Number>
class ParsedFEMFunction : public FEMFunctionBase<Output>
{
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
      _n_vars(sys.n_vars())
  {
    std::string variables = "x";
#if LIBMESH_DIM > 1
    variables += ",y";
#endif
#if LIBMESH_DIM > 2
    variables += ",z";
#endif
    variables += ",t";

    _spacetime.resize(LIBMESH_DIM+1 + _n_vars + (additional_vars ? additional_vars->size() : 0));

    for (unsigned int v=0; v != _n_vars; ++v)
      {
        variables += ',';
        variables += sys.variable_name(v);
      }

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
            _spacetime[LIBMESH_DIM+1+_n_vars + i] = (initial_vals && i < initial_vals->size()) ? (*initial_vals)[i] : 0;
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
        elem_fe->get_phi();
      }
  }

  /**
   * Returns a new copy of the function.  The new copy should be as
   * ``deep'' as necessary to allow independent destruction and
   * simultaneous evaluations of the copies in different threads.
   */
  virtual AutoPtr<FEMFunctionBase<Output> > clone () const {
    return AutoPtr<FEMFunctionBase<Output> >
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
    _spacetime[0] = p(0);
#if LIBMESH_DIM > 1
    _spacetime[1] = p(1);
#endif
#if LIBMESH_DIM > 2
    _spacetime[2] = p(2);
#endif
    _spacetime[LIBMESH_DIM] = time;

    for (unsigned int v=0; v != _n_vars; ++v)
      {
        c.point_value(v, p, _spacetime[LIBMESH_DIM+1+v]);
      }

    // The remaining locations in _spacetime are currently fixed at construction
    // but could potentially be made dynamic
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
    _spacetime[0] = p(0);
#if LIBMESH_DIM > 1
    _spacetime[1] = p(1);
#endif
#if LIBMESH_DIM > 2
    _spacetime[2] = p(2);
#endif
    _spacetime[LIBMESH_DIM] = time;

    for (unsigned int v=0; v != _n_vars; ++v)
      {
        c.point_value(v, p, _spacetime[LIBMESH_DIM+1+v]);
      }

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
    _spacetime[0] = p(0);
#if LIBMESH_DIM > 1
    _spacetime[1] = p(1);
#endif
#if LIBMESH_DIM > 2
    _spacetime[2] = p(2);
#endif
    _spacetime[LIBMESH_DIM] = time;

    for (unsigned int v=0; v != _n_vars; ++v)
      {
        c.point_value(v, p, _spacetime[LIBMESH_DIM+1+v]);
      }
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
  unsigned int _n_vars;
#ifdef LIBMESH_HAVE_FPARSER
  std::vector<FunctionParserBase<Output> > parsers;
#endif
  std::vector<Output> _spacetime;

  // Additional variables/values that can be parsed and handled by the function parser
  std::vector<std::string> _additional_vars;
  std::vector<Output> _initial_vals;
};

} // namespace libMesh

#endif // LIBMESH_PARSED_FEM_FUNCTION_H
