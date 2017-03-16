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



#include "libmesh/dirichlet_boundaries.h"

#ifdef LIBMESH_ENABLE_DIRICHLET

// Local Includes
#include "libmesh/composite_fem_function.h"
#include "libmesh/composite_function.h"
#include "libmesh/vector_value.h"


namespace libMesh
{

DirichletBoundary::
DirichletBoundary(const std::set<boundary_id_type> & b_in,
                  const std::vector<unsigned int> & variables_in,
                  const FunctionBase<Number> * f_in,
                  const FunctionBase<Gradient> * g_in) :
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


DirichletBoundary::
DirichletBoundary(const std::set<boundary_id_type> & b_in,
                  const std::vector<unsigned int> & variables_in,
                  const FunctionBase<Number> & f_in,
                  VariableIndexing type) :
  b(b_in),
  variables(variables_in),
  f(UniquePtr<FunctionBase<Number> >()),
  g(UniquePtr<FunctionBase<Gradient> >()),
  f_fem(UniquePtr<FEMFunctionBase<Number> >()),
  g_fem(UniquePtr<FEMFunctionBase<Gradient> >()),
  f_system(libmesh_nullptr)
{
  if (type == LOCAL_VARIABLE_ORDER)
    {
      CompositeFunction<Number> * c =
        new CompositeFunction<Number>();
      c->attach_subfunction(f_in, variables_in);
      f.reset(c);
    }
  else
    {
      f.reset(f_in.clone().release());
    }
  f->init();
}


DirichletBoundary::
DirichletBoundary(const std::set<boundary_id_type> & b_in,
                  const std::vector<unsigned int> & variables_in,
                  const FunctionBase<Number> & f_in,
                  const FunctionBase<Gradient> & g_in,
                  VariableIndexing type) :
  b(b_in),
  variables(variables_in),
  f(UniquePtr<FunctionBase<Number> >()),
  g(UniquePtr<FunctionBase<Gradient> >()),
  f_fem(UniquePtr<FEMFunctionBase<Number> >()),
  g_fem(UniquePtr<FEMFunctionBase<Gradient> >()),
  f_system(libmesh_nullptr)
{
  if (type == LOCAL_VARIABLE_ORDER)
    {
      CompositeFunction<Number> * cf =
        new CompositeFunction<Number>();
      cf->attach_subfunction(f_in, variables_in);
      f.reset(cf);
      CompositeFunction<Gradient> * cg =
        new CompositeFunction<Gradient>();
      cg->attach_subfunction(g_in, variables_in);
      g.reset(cg);
    }
  else
    {
      f.reset(f_in.clone().release());
      g.reset(g_in.clone().release());
    }

  f->init();
  g->init();
}


DirichletBoundary::
DirichletBoundary(const std::set<boundary_id_type> & b_in,
                  const std::vector<unsigned int> & variables_in,
                  const System & f_sys_in,
                  const FEMFunctionBase<Number> * f_in,
                  const FEMFunctionBase<Gradient> * g_in) :
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


DirichletBoundary::
DirichletBoundary(const std::set<boundary_id_type> & b_in,
                  const std::vector<unsigned int> & variables_in,
                  const System & f_sys_in,
                  const FEMFunctionBase<Number> & f_in,
                  VariableIndexing type) :
  b(b_in),
  variables(variables_in),
  f(UniquePtr<FunctionBase<Number> >()),
  g(UniquePtr<FunctionBase<Gradient> >()),
  f_fem(UniquePtr<FEMFunctionBase<Number> >()),
  g_fem(UniquePtr<FEMFunctionBase<Gradient> >()),
  f_system(&f_sys_in)
{
  if (type == LOCAL_VARIABLE_ORDER)
    {
      CompositeFEMFunction<Number> * c = new CompositeFEMFunction<Number>();
      c->attach_subfunction(f_in, variables_in);
      f_fem.reset(c);
    }
  else
    {
      f_fem.reset(f_in.clone().release());
    }
}


DirichletBoundary::
DirichletBoundary(const std::set<boundary_id_type> & b_in,
                  const std::vector<unsigned int> & variables_in,
                  const System & f_sys_in,
                  const FEMFunctionBase<Number> & f_in,
                  const FEMFunctionBase<Gradient> & g_in,
                  VariableIndexing type) :
  b(b_in),
  variables(variables_in),
  f(UniquePtr<FunctionBase<Number> >()),
  g(UniquePtr<FunctionBase<Gradient> >()),
  f_fem(UniquePtr<FEMFunctionBase<Number> >()),
  g_fem(UniquePtr<FEMFunctionBase<Gradient> >()),
  f_system(&f_sys_in)
{
  if (type == LOCAL_VARIABLE_ORDER)
    {
      CompositeFEMFunction<Number> * cf =
        new CompositeFEMFunction<Number>();
      cf->attach_subfunction(f_in, variables_in);
      f_fem.reset(cf);
      CompositeFEMFunction<Gradient> * cg =
        new CompositeFEMFunction<Gradient>();
      cg->attach_subfunction(g_in, variables_in);
      g_fem.reset(cg);
    }
  else
    {
      f_fem.reset(f_in.clone().release());
      g_fem.reset(g_in.clone().release());
    }
}


DirichletBoundary::
DirichletBoundary(const DirichletBoundary & dirichlet_in) :
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


DirichletBoundary::~DirichletBoundary () {}

} // namespace libMesh

#endif // LIBMESH_ENABLE_DIRICHLET
