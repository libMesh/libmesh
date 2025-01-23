// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

template <unsigned int dim>
class SolutionFunction : public FunctionBase<Number>
{
public:
  SolutionFunction() = default;
  ~SolutionFunction() = default;

  virtual Number operator()(const Point &, const Real = 0) override { libmesh_not_implemented(); }

  virtual void operator()(const Point & p, const Real, DenseVector<Number> & output) override;

  virtual std::unique_ptr<FunctionBase<Number>> clone() const override
  {
    return std::make_unique<SolutionFunction>();
  }

  virtual Number component(unsigned int i, const Point & p, Real time = 0.) override
  {
    const auto size = dim == 2 ? 4 : 5;
    DenseVector<Number> outvec(size);
    (*this)(p, time, outvec);
    return outvec(i);
  }

private:
  DivGradExactSolution soln;
};

template <>
void
SolutionFunction<2>::operator()(const Point & p, const Real, DenseVector<Number> & output)
{
  // We should have our vector variable, our scalar variable, and our higher order scalar variable
  libmesh_assert(output.size() <= 4);
  output.zero();
  const Real x = p(0), y = p(1);
  // libMesh assumes each component of a vector-valued variable is stored
  // contiguously.
  const auto vector = soln(x, y);
  output(0) = vector(0);
  output(1) = vector(1);
  const auto scalar = soln.scalar(x, y);
  output(2) = scalar;
  output(3) = scalar;
}

template <>
void
SolutionFunction<3>::operator()(const Point & p, const Real, DenseVector<Number> & output)
{
  // We should have our vector variable, our scalar variable, and our higher order scalar variable
  libmesh_assert(output.size() <= 5);
  output.zero();
  const Real x = p(0), y = p(1), z = p(2);
  // libMesh assumes each component of the vector-valued variable is stored
  // contiguously.
  const auto vector = soln(x, y, z);
  output(0) = vector(0);
  output(1) = vector(1);
  output(2) = vector(2);
  const auto scalar = soln.scalar(x, y, z);
  output(3) = scalar;
  output(4) = scalar;
}

template <unsigned int dim>
class SolutionGradient : public FunctionBase<Gradient>
{
public:
  SolutionGradient() = default;
  ~SolutionGradient() = default;

  virtual Gradient operator()(const Point &, const Real = 0) override { libmesh_not_implemented(); }

  virtual void operator()(const Point & p, const Real, DenseVector<Gradient> & output) override;

  virtual std::unique_ptr<FunctionBase<Gradient>> clone() const override
  {
    return std::make_unique<SolutionGradient>();
  }

  virtual Gradient component(unsigned int i, const Point & p, Real time = 0.) override
  {
    const auto size = dim == 2 ? 4 : 5;
    DenseVector<Gradient> outvec(size);
    (*this)(p, time, outvec);
    return outvec(i);
  }

private:
  DivGradExactSolution soln;
};

template <>
void
SolutionGradient<2>::operator()(const Point & p, const Real, DenseVector<Gradient> & output)
{
  output.zero();
  const Real x = p(0), y = p(1);
  output(0) = soln.grad(x, y).row(0);
  output(1) = soln.grad(x, y).row(1);
}

template <>
void
SolutionGradient<3>::operator()(const Point & p, const Real, DenseVector<Gradient> & output)
{
  output.zero();
  const Real x = p(0), y = p(1), z = p(2);
  output(0) = soln.grad(x, y, z).row(0);
  output(1) = soln.grad(x, y, z).row(1);
  output(2) = soln.grad(x, y, z).row(2);
}

#endif // SOLUTION_FUNCTION_H
