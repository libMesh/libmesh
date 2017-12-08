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

// libMesh includes
#include "libmesh/function_base.h"
#include "libmesh/auto_ptr.h" // libmesh_make_unique

// Example includes
#include "laplace_exact_solution.h"

using namespace libMesh;

#ifndef SOLUTION_FUNCTION_H
#define SOLUTION_FUNCTION_H

class SolutionFunction : public FunctionBase<Number>
{
public:

  SolutionFunction (const unsigned int u_var) :
    _u_var(u_var) {}

  ~SolutionFunction() {}

  virtual Number operator() (const Point &,
                             const Real = 0)
  { libmesh_not_implemented(); }

  virtual void operator() (const Point & p,
                           const Real,
                           DenseVector<Number> & output)
  {
    output.zero();
    const Real x=p(0), y=p(1), z=p(2);
    // libMesh assumes each component of the vector-valued variable is stored
    // contiguously.
    output(_u_var)   = soln(0, x, y, z);
    output(_u_var+1) = soln(1, x, y, z);
    output(_u_var+2) = soln(2, x, y, z);
  }

  virtual Number component(unsigned int component_in,
                           const Point & p,
                           const Real)
  {
    const Real x=p(0), y=p(1), z=p(2);
    return soln(component_in, x, y, z);
  }

  virtual std::unique_ptr<FunctionBase<Number>> clone() const
  { return libmesh_make_unique<SolutionFunction>(_u_var); }

private:

  const unsigned int _u_var;
  LaplaceExactSolution soln;
};

//FIXME: PB: We ought to be able to merge the above class with this one
//           through templating, but I'm being lazy.
class SolutionGradient : public FunctionBase<Gradient>
{
public:

  SolutionGradient(const unsigned int u_var)
    : _u_var(u_var) {}
  ~SolutionGradient(){}

  virtual Gradient operator() (const Point &, const Real = 0)
  { libmesh_not_implemented(); }

  virtual void operator() (const Point & p,
                           const Real,
                           DenseVector<Gradient>& output)
  {
    output.zero();
    const Real x=p(0), y=p(1), z=p(2);
    output(_u_var)   = soln(0, x, y, z);
    output(_u_var+1) = soln(1, x, y, z);
    output(_u_var+2) = soln(2, x, y, z);
  }

  virtual Gradient component(unsigned int component_in,
                             const Point & p,
                             const Real)
  {
    const Real x=p(0), y=p(1), z=p(2);
    return soln(component_in, x, y, z);
  }

  virtual std::unique_ptr<FunctionBase<Gradient>> clone() const
  { return libmesh_make_unique<SolutionGradient>(_u_var); }

private:

  const unsigned int _u_var;
  LaplaceExactGradient soln;
};

#endif // SOLUTION_FUNCTION_H
