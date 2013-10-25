// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_DIRICHLET_BOUNDARIES_H
#define LIBMESH_DIRICHLET_BOUNDARIES_H

#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_ENABLE_DIRICHLET

// Local Includes -----------------------------------
#include "libmesh/function_base.h"
#include "libmesh/id_types.h"
#include "libmesh/vector_value.h"

// C++ Includes   -----------------------------------
#include <algorithm>
#include <cstddef>
#include <iterator>
#include <set>
#include <string>
#include <vector>

namespace libMesh
{

class DirichletBoundary
{
public:
  DirichletBoundary(const std::set<boundary_id_type> &b_in,
                    const std::vector<unsigned int>& variables_in,
                    const FunctionBase<Number> *f_in,
                    const FunctionBase<Gradient> *g_in = NULL) :
    b(b_in),
    variables(variables_in),
    f(f_in ? f_in->clone() : AutoPtr<FunctionBase<Number> >(NULL)),
    g(g_in ? g_in->clone() : AutoPtr<FunctionBase<Gradient> >(NULL))
    {
      libmesh_assert(f.get());
      f->init();
      if (g.get())
        g->init();
    }

  DirichletBoundary(const std::set<boundary_id_type> &b_in,
                    const std::vector<unsigned int>& variables_in,
                    const FunctionBase<Number> &f_in) :
    b(b_in),
    variables(variables_in),
    f(f_in.clone()),
    g(AutoPtr<FunctionBase<Gradient> >(NULL))
    {
      f->init();
    }


  DirichletBoundary(const std::set<boundary_id_type> &b_in,
                    const std::vector<unsigned int>& variables_in,
                    const FunctionBase<Number> &f_in,
                    const FunctionBase<Gradient> &g_in) :
    b(b_in),
    variables(variables_in),
    f(f_in.clone()),
    g(g_in.clone())
    {
      f->init();
      g->init();
    }


  DirichletBoundary (const DirichletBoundary &dirichlet_in) :
    b(dirichlet_in.b),
    variables(dirichlet_in.variables),
    f(dirichlet_in.f.get() ? dirichlet_in.f->clone() : AutoPtr<FunctionBase<Number> >(NULL)),
    g(dirichlet_in.g.get() ? dirichlet_in.g->clone() : AutoPtr<FunctionBase<Gradient> >(NULL))
    {
      libmesh_assert(f.get());
      f->init();
      if (g.get())
        g->init();
    }

  std::set<boundary_id_type> b;
  std::vector<unsigned int> variables;
  AutoPtr<FunctionBase<Number> > f;
  AutoPtr<FunctionBase<Gradient> > g;
};


/**
 * We're using a class instead of a typedef to allow forward
 * declarations and future flexibility.  Note that std::map has no
 * virtual destructor, so downcasting here would be dangerous.
 */
class DirichletBoundaries : public std::vector<DirichletBoundary *>
{
public:
  DirichletBoundaries() {}

  ~DirichletBoundaries();
};

} // namespace libMesh

#endif // LIBMESH_ENABLE_DIRICHLET

#endif // LIBMESH_DIRICHLET_BOUNDARIES_H
