// The libMesh Finite Element Library.
// Copyright (C) 2002-2024 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef EXACT_SOLN_H
#define EXACT_SOLN_H

#include "libmesh/parameters.h"
#include "libmesh/libmesh_common.h"

namespace libMesh
{
class Point;

class ExactSoln
{
public:
  virtual Real operator()(const Point & p) const = 0;
  virtual Real forcing(const Point & p) const = 0;
};

inline Number
compute_error(const Point & p,
              const Parameters & params,
              const std::string & /*sys_name*/,
              const std::string & unknown_name)
{
  const auto * const error_obj = params.get<const ExactSoln *>(unknown_name + "_exact_sol");
  return (*error_obj)(p);
}

class USoln : public ExactSoln
{
public:
  USoln(const Real nu_in, const bool cavity_in) : nu(nu_in), cavity(cavity_in) {}

  Real operator()(const Point & p) const override
  {
    const auto x = p(0);
    const auto y = p(1);
    if (cavity)
      return sin(y) * cos((1. / 2) * x * pi);
    else
      return sin(1. / 2 * y * pi) * cos(1. / 2 * x * pi);
  }

  Real forcing(const Point & p) const override
  {
    const auto x = p(0);
    const auto y = p(1);
    if (cavity)
    {
      return nu * sin(y) * cos((1. / 2) * x * pi) +
             (1. / 4) * pi * pi * nu * sin(y) * cos((1. / 2) * x * pi) -
             1. / 2 * pi * sin(x) * sin(y) * sin((1. / 2) * y * pi) * cos((1. / 2) * x * pi) +
             sin(x) * cos(y) * cos((1. / 2) * x * pi) * cos((1. / 2) * y * pi) -
             pi * sin(y) * sin(y) * sin((1. / 2) * x * pi) * cos((1. / 2) * x * pi) +
             sin(y) * cos(x);
    }
    else
    {
      const auto quant1 = sin((1. / 2) * y * pi);
      const auto quant2 = cos((1. / 2) * y * pi);
      return (1. / 2) * pi * pi * nu * sin((1. / 2) * y * pi) * cos((1. / 2) * x * pi) -
             1. / 2 * pi * sin((1. / 4) * x * pi) * quant1 * quant1 * cos((1. / 2) * x * pi) -
             1. / 4 * pi * sin((1. / 4) * x * pi) * sin((3. / 2) * y * pi) +
             (1. / 2) * pi * sin((1. / 4) * x * pi) * cos((1. / 2) * x * pi) * quant2 * quant2 -
             pi * sin((1. / 2) * x * pi) * quant1 * quant1 * cos((1. / 2) * x * pi);
    }
  }

private:
  const Real nu;
  const bool cavity;
};

class VSoln : public ExactSoln
{
public:
  VSoln(const Real nu_in, const bool cavity_in) : nu(nu_in), cavity(cavity_in) {}

  Real operator()(const Point & p) const override
  {
    const auto x = p(0);
    const auto y = p(1);
    if (cavity)
      return sin(x) * cos((1. / 2) * y * pi);
    else
      return sin((1. / 4) * x * pi) * cos((1. / 2) * y * pi);
  }

  Real forcing(const Point & p) const override
  {
    const auto x = p(0);
    const auto y = p(1);
    if (cavity)
      return nu * sin(x) * cos((1. / 2) * y * pi) +
             (1. / 4) * pi * pi * nu * sin(x) * cos((1. / 2) * y * pi) -
             pi * sin(x) * sin(x) * sin((1. / 2) * y * pi) * cos((1. / 2) * y * pi) -
             1. / 2 * pi * sin(x) * sin(y) * sin((1. / 2) * x * pi) * cos((1. / 2) * y * pi) +
             sin(y) * cos(x) * cos((1. / 2) * x * pi) * cos((1. / 2) * y * pi) + sin(x) * cos(y);
    else
    {
      const auto quant1 = sin((1. / 4) * x * pi);
      return (5. / 16) * pi * pi * nu * sin((1. / 4) * x * pi) * cos((1. / 2) * y * pi) -
             pi * quant1 * quant1 * sin((1. / 2) * y * pi) * cos((1. / 2) * y * pi) -
             1. / 2 * pi * sin((1. / 4) * x * pi) * sin((1. / 2) * x * pi) *
                 sin((1. / 2) * y * pi) * cos((1. / 2) * y * pi) +
             (1. / 4) * pi * sin((1. / 2) * y * pi) * cos((1. / 4) * x * pi) *
                 cos((1. / 2) * x * pi) * cos((1. / 2) * y * pi) +
             (3. / 2) * pi * cos((1. / 4) * x * pi) * cos((3. / 2) * y * pi);
    }
  }

private:
  const Real nu;
  const bool cavity;
};

class PSoln : public ExactSoln
{
public:
  PSoln(const bool cavity_in) : cavity(cavity_in) {}

  Real operator()(const Point & p) const override
  {
    const auto x = p(0);
    const auto y = p(1);
    if (cavity)
      return sin(x) * sin(y);
    else
      return sin((3. / 2) * y * pi) * cos((1. / 4) * x * pi);
  }

  Real forcing(const Point & p) const override
  {
    const auto x = p(0);
    const auto y = p(1);
    if (cavity)
      return -1. / 2 * pi * sin(x) * sin((1. / 2) * y * pi) -
             1. / 2 * pi * sin(y) * sin((1. / 2) * x * pi);
    else
      return -1. / 2 * pi * sin((1. / 4) * x * pi) * sin((1. / 2) * y * pi) -
             1. / 2 * pi * sin((1. / 2) * x * pi) * sin((1. / 2) * y * pi);
  }

private:
  const bool cavity;
};
}

#endif // EXACT_SOLN_H
