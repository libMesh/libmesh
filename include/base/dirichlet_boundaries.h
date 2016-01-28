// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/fem_function_base.h"
#include "libmesh/function_base.h"
#include "libmesh/id_types.h"
#include "libmesh/system.h"
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

/**
 * This class allows one to associate dirichlet boundary values with
 * a given set of mesh boundary ids and system variable ids.
 *
 * Dirichlet values must be supplied as the input function "f"; when
 * using some specialized elements, gradient values must be supplied
 * via the input function "g".
 *
 * Dirichlet functions are allowed to return NaN; if this is
 * encountered, then the degree of freedom values in a patch around
 * the location of the returned NaN will be left unconstrained.  E.g.
 * a NaN on a boundary edge in 3D would leave that edge and the two
 * adjoining face interiors unconstrained, but would still permit the
 * other edge and node DoFs around those faces to be constrained.
 */
class DirichletBoundary
{
public:
  DirichletBoundary(const std::set<boundary_id_type> & b_in,
                    const std::vector<unsigned int> & variables_in,
                    const FunctionBase<Number> * f_in,
                    const FunctionBase<Gradient> * g_in = libmesh_nullptr) :
    b(b_in),
    variables(variables_in),
    f(f_in ? f_in->clone() : UniquePtr<FunctionBase<Number> >()),
    g(g_in ? g_in->clone() : UniquePtr<FunctionBase<Gradient> >()),
    f_fem(UniquePtr<FEMFunctionBase<Number> >()),
    g_fem(UniquePtr<FEMFunctionBase<Gradient> >()),
    f_system(libmesh_nullptr)
  {
    libmesh_assert(f.get());
    f->init();
    if (g.get())
      g->init();
  }

  DirichletBoundary(const std::set<boundary_id_type> & b_in,
                    const std::vector<unsigned int> & variables_in,
                    const FunctionBase<Number> & f_in) :
    b(b_in),
    variables(variables_in),
    f(f_in.clone()),
    g(UniquePtr<FunctionBase<Gradient> >()),
    f_fem(UniquePtr<FEMFunctionBase<Number> >()),
    g_fem(UniquePtr<FEMFunctionBase<Gradient> >()),
    f_system(libmesh_nullptr)
  {
    f->init();
  }


  DirichletBoundary(const std::set<boundary_id_type> & b_in,
                    const std::vector<unsigned int> & variables_in,
                    const FunctionBase<Number> & f_in,
                    const FunctionBase<Gradient> & g_in) :
    b(b_in),
    variables(variables_in),
    f(f_in.clone()),
    g(g_in.clone()),
    f_fem(UniquePtr<FEMFunctionBase<Number> >()),
    g_fem(UniquePtr<FEMFunctionBase<Gradient> >()),
    f_system(libmesh_nullptr)
  {
    f->init();
    g->init();
  }


  DirichletBoundary(const std::set<boundary_id_type> & b_in,
                    const std::vector<unsigned int> & variables_in,
                    const System & f_sys_in,
                    const FEMFunctionBase<Number> * f_in,
                    const FEMFunctionBase<Gradient> * g_in = libmesh_nullptr) :
    b(b_in),
    variables(variables_in),
    f(UniquePtr<FunctionBase<Number> >()),
    g(UniquePtr<FunctionBase<Gradient> >()),
    f_fem(f_in ? f_in->clone() : UniquePtr<FEMFunctionBase<Number> >()),
    g_fem(g_in ? g_in->clone() : UniquePtr<FEMFunctionBase<Gradient> >()),
    f_system(&f_sys_in)
  {
    libmesh_assert(f_fem.get());
  }

  DirichletBoundary(const std::set<boundary_id_type> & b_in,
                    const std::vector<unsigned int> & variables_in,
                    const System & f_sys_in,
                    const FEMFunctionBase<Number> & f_in) :
    b(b_in),
    variables(variables_in),
    f(UniquePtr<FunctionBase<Number> >()),
    g(UniquePtr<FunctionBase<Gradient> >()),
    f_fem(f_in.clone()),
    g_fem(UniquePtr<FEMFunctionBase<Gradient> >()),
    f_system(&f_sys_in)
  {
  }


  DirichletBoundary(const std::set<boundary_id_type> & b_in,
                    const std::vector<unsigned int> & variables_in,
                    const System & f_sys_in,
                    const FEMFunctionBase<Number> & f_in,
                    const FEMFunctionBase<Gradient> & g_in) :
    b(b_in),
    variables(variables_in),
    f(UniquePtr<FunctionBase<Number> >()),
    g(UniquePtr<FunctionBase<Gradient> >()),
    f_fem(f_in.clone()),
    g_fem(g_in.clone()),
    f_system(&f_sys_in)
  {
  }




  DirichletBoundary (const DirichletBoundary & dirichlet_in) :
    b(dirichlet_in.b),
    variables(dirichlet_in.variables),
    f(dirichlet_in.f.get() ?
      dirichlet_in.f->clone() : UniquePtr<FunctionBase<Number> >()),
    g(dirichlet_in.g.get() ?
      dirichlet_in.g->clone() : UniquePtr<FunctionBase<Gradient> >()),
    f_fem(dirichlet_in.f_fem.get() ?
          dirichlet_in.f_fem->clone() : UniquePtr<FEMFunctionBase<Number> >()),
    g_fem(dirichlet_in.g_fem.get() ?
          dirichlet_in.g_fem->clone() : UniquePtr<FEMFunctionBase<Gradient> >()),
    f_system(dirichlet_in.f_system)
  {
    libmesh_assert(f.get() || f_fem.get());
    libmesh_assert(!(f.get() && f_fem.get()));
    libmesh_assert(!(f.get() && g_fem.get()));
    libmesh_assert(!(f_fem.get() && g.get()));
    libmesh_assert(!(f_fem.get() && !f_system));
    if (f.get())
      f->init();
    if (g.get())
      g->init();
  }

  std::set<boundary_id_type> b;
  std::vector<unsigned int> variables;

  UniquePtr<FunctionBase<Number> > f;
  UniquePtr<FunctionBase<Gradient> > g;

  UniquePtr<FEMFunctionBase<Number> > f_fem;
  UniquePtr<FEMFunctionBase<Gradient> > g_fem;

  const System * f_system;
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
