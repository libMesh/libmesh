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



#include "libmesh/dirichlet_boundaries.h"

#ifdef LIBMESH_ENABLE_DIRICHLET

// Local Includes
#include "libmesh/composite_fem_function.h"
#include "libmesh/composite_function.h"
#include "libmesh/vector_value.h"

// C++ Includes
#include <memory>

namespace libMesh
{

DirichletBoundary::
DirichletBoundary(std::set<boundary_id_type> b_in,
                  std::vector<unsigned int> variables_in,
                  const FunctionBase<Number> * f_in,
                  const FunctionBase<Gradient> * g_in) :
  b(std::move(b_in)),
  variables(std::move(variables_in)),
  f(f_in ? f_in->clone() : nullptr),
  g(g_in ? g_in->clone() : nullptr),
  f_system(nullptr),
  jacobian_tolerance(0.)
{
  libmesh_assert(f);
  f->init();
  if (g)
    g->init();
}


DirichletBoundary::
DirichletBoundary(std::set<boundary_id_type> b_in,
                  std::vector<unsigned int> variables_in,
                  const FunctionBase<Number> & f_in,
                  VariableIndexing type) :
  b(std::move(b_in)),
  variables(std::move(variables_in)),
  f_system(nullptr),
  jacobian_tolerance(0.)
{
  if (type == LOCAL_VARIABLE_ORDER)
    {
      auto c = std::make_unique<CompositeFunction<Number>>();
      c->attach_subfunction(f_in, variables);
      f = std::move(c);
    }
  else
    f = f_in.clone();

  f->init();
}


DirichletBoundary::
DirichletBoundary(std::set<boundary_id_type> b_in,
                  std::vector<unsigned int> variables_in,
                  const FunctionBase<Number> & f_in,
                  const FunctionBase<Gradient> & g_in,
                  VariableIndexing type) :
  b(std::move(b_in)),
  variables(std::move(variables_in)),
  f_system(nullptr),
  jacobian_tolerance(0.)
{
  if (type == LOCAL_VARIABLE_ORDER)
    {
      auto cf = std::make_unique<CompositeFunction<Number>>();
      cf->attach_subfunction(f_in, variables);
      f = std::move(cf);

      auto cg = std::make_unique<CompositeFunction<Gradient>>();
      cg->attach_subfunction(g_in, variables);
      g = std::move(cg);
    }
  else
    {
      f = f_in.clone();
      g = g_in.clone();
    }

  f->init();
  g->init();
}


DirichletBoundary::
DirichletBoundary(std::set<boundary_id_type> b_in,
                  std::vector<unsigned int> variables_in,
                  const System & f_sys_in,
                  const FEMFunctionBase<Number> * f_in,
                  const FEMFunctionBase<Gradient> * g_in) :
  b(std::move(b_in)),
  variables(std::move(variables_in)),
  f_fem(f_in ? f_in->clone() : nullptr),
  g_fem(g_in ? g_in->clone() : nullptr),
  f_system(&f_sys_in),
  jacobian_tolerance(0.)
{
  libmesh_assert(f_fem);
}


DirichletBoundary::
DirichletBoundary(std::set<boundary_id_type> b_in,
                  std::vector<unsigned int> variables_in,
                  const System & f_sys_in,
                  const FEMFunctionBase<Number> & f_in,
                  VariableIndexing type) :
  b(std::move(b_in)),
  variables(std::move(variables_in)),
  f_system(&f_sys_in),
  jacobian_tolerance(0.)
{
  if (type == LOCAL_VARIABLE_ORDER)
    {
      auto c = std::make_unique<CompositeFEMFunction<Number>>();
      c->attach_subfunction(f_in, variables);
      f_fem = std::move(c);
    }
  else
    f_fem = f_in.clone();
}


DirichletBoundary::
DirichletBoundary(std::set<boundary_id_type> b_in,
                  std::vector<unsigned int> variables_in,
                  const System & f_sys_in,
                  const FEMFunctionBase<Number> & f_in,
                  const FEMFunctionBase<Gradient> & g_in,
                  VariableIndexing type) :
  b(std::move(b_in)),
  variables(std::move(variables_in)),
  f_system(&f_sys_in),
  jacobian_tolerance(0.)
{
  if (type == LOCAL_VARIABLE_ORDER)
    {
      auto cf = std::make_unique<CompositeFEMFunction<Number>>();
      cf->attach_subfunction(f_in, variables);
      f_fem = std::move(cf);

      auto cg = std::make_unique<CompositeFEMFunction<Gradient>>();
      cg->attach_subfunction(g_in, variables);
      g_fem = std::move(cg);
    }
  else
    {
      f_fem = f_in.clone();
      g_fem = g_in.clone();
    }
}


DirichletBoundary::
DirichletBoundary(const DirichletBoundary & d_in) :
  b(d_in.b),
  variables(d_in.variables),
  f(d_in.f ? d_in.f->clone() : nullptr),
  g(d_in.g ? d_in.g->clone() : nullptr),
  f_fem(d_in.f_fem ? d_in.f_fem->clone() : nullptr),
  g_fem(d_in.g_fem ? d_in.g_fem->clone() : nullptr),
  f_system(d_in.f_system),
  jacobian_tolerance(d_in.jacobian_tolerance)
{
  libmesh_assert(f || f_fem);
  libmesh_assert(!(f && f_fem));
  libmesh_assert(!(f && g_fem));
  libmesh_assert(!(f_fem && g));
  libmesh_assert(!(f_fem && !f_system));
  if (f)
    f->init();
  if (g)
    g->init();
}


DirichletBoundary & DirichletBoundary::operator= (const DirichletBoundary & rhs)
{
  // Implementation in terms of the copy constructor to avoid code duplication.
  DirichletBoundary tmp(rhs);
  std::swap(tmp, *this); // class must be "MoveAssignable" and "MoveConstructible" for std::swap to work.
  return *this;
}


DirichletBoundary::~DirichletBoundary () = default;

} // namespace libMesh

#endif // LIBMESH_ENABLE_DIRICHLET
