// The libMesh Finite Element Library.
// Copyright (C) 2002-2023 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef SOLUTION_FUNCTION_H
#define SOLUTION_FUNCTION_H

// libMesh includes
#include "libmesh/function_base.h"

// Example includes
#include "div_grad_exact_solution.h"

// C++ includes
#include <memory>

using namespace libMesh;

template<unsigned int dim>
class SolutionFunction : public FunctionBase<Number>
{
public:

  SolutionFunction() = default;
  ~SolutionFunction() = default;

  virtual Number operator() (const Point &,
                             const Real = 0)
  { libmesh_not_implemented(); }

  virtual void operator() (const Point & p,
                           const Real,
                           DenseVector<Number> & output);

  virtual Number component(unsigned int component_in,
                           const Point & p,
                           const Real);

  virtual std::unique_ptr<FunctionBase<Number>> clone() const
  { return std::make_unique<SolutionFunction>(); }

private:

  DivGradExactSolution soln;
};

template<>
void SolutionFunction<2>::operator() (const Point & p,
                                      const Real,
                                      DenseVector<Number> & output)
{
  output.zero();
  const Real x=p(0), y=p(1);
  // libMesh assumes each component of the vector-valued variable is stored
  // contiguously.
  output(0) = soln(x, y)(0);
  output(1) = soln(x, y)(1);
  // scalar solution
  output(2) = soln.scalar(x, y);
}

template<>
void SolutionFunction<3>::operator() (const Point & p,
                                      const Real,
                                      DenseVector<Number> & output)
{
  output.zero();
  const Real x=p(0), y=p(1), z=p(2);
  // libMesh assumes each component of the vector-valued variable is stored
  // contiguously.
  output(0) = soln(x, y, z)(0);
  output(1) = soln(x, y, z)(1);
  output(2) = soln(x, y, z)(2);
  // scalar solution
  output(3) = soln.scalar(x, y);
}

template<>
Number SolutionFunction<2>::component(unsigned int component_in,
                                      const Point & p,
                                      const Real)
{
  const Real x=p(0), y=p(1);
  if (component_in == 2)
    return soln.scalar(x, y);
  else
    return soln(x, y)(component_in);
}

template<>
Number SolutionFunction<3>::component(unsigned int component_in,
                                      const Point & p,
                                      const Real)
{
  const Real x=p(0), y=p(1), z=p(2);
  if (component_in == 3)
    return soln.scalar(x, y);
  else
    return soln(x, y, z)(component_in);
}

template<unsigned int dim>
class SolutionGradient : public FunctionBase<Gradient>
{
public:

  SolutionGradient() = default;
  ~SolutionGradient() = default;

  virtual Gradient operator() (const Point &, const Real = 0)
  { libmesh_not_implemented(); }

  virtual void operator() (const Point & p,
                           const Real,
                           DenseVector<Gradient> & output);

  virtual Gradient component(unsigned int component_in,
                             const Point & p,
                             const Real);

  virtual std::unique_ptr<FunctionBase<Gradient>> clone() const
  { return std::make_unique<SolutionGradient>(); }

private:

  DivGradExactSolution soln;
};

template<>
void SolutionGradient<2>::operator() (const Point & p,
                                      const Real,
                                      DenseVector<Gradient> & output)
{
  output.zero();
  const Real x=p(0), y=p(1);
  output(0) = soln.grad(x, y).row(0);
  output(1) = soln.grad(x, y).row(1);
}

template<>
void SolutionGradient<3>::operator() (const Point & p,
                                      const Real,
                                      DenseVector<Gradient> & output)
{
  output.zero();
  const Real x=p(0), y=p(1), z=p(2);
  output(0) = soln.grad(x, y, z).row(0);
  output(1) = soln.grad(x, y, z).row(1);
  output(2) = soln.grad(x, y, z).row(2);
}

template<>
Gradient SolutionGradient<2>::component(unsigned int component_in,
                                        const Point & p,
                                        const Real)
{
  const Real x=p(0), y=p(1);
  if (component_in == 2)
    // Only interested in L2 norms for the scalar field
    return {};
  else
    return soln.grad(x, y).row(component_in);
}

template<>
Gradient SolutionGradient<3>::component(unsigned int component_in,
                                        const Point & p,
                                        const Real)
{
  const Real x=p(0), y=p(1), z=p(2);
  if (component_in == 3)
    // Only interested in L2 norms for the scalar field
    return {};
  return soln.grad(x, y, z).row(component_in);
}

#endif // SOLUTION_FUNCTION_H
