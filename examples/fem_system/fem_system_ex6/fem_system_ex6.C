// The libMesh Finite Element Library.
// Copyright (C) 2002-2026 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// <h1>FEMSystem Example 6 - Nonlinear Diffusion with Automatic Differentiation</h1>
// \author Alexander Lindsay
// \date 2026
//
// This example solves the steady nonlinear diffusion equation
//
//   -div((1 + u^2) grad u) = f
//
// on the unit square with homogeneous Dirichlet boundary conditions,
// with f chosen so that u = sin(pi x) sin(pi y).
//
// The element residual is evaluated with MetaPhysicL dual numbers in
// place of real numbers. Each element degree of freedom is seeded with
// a unit derivative with respect to itself, and the derivatives then
// propagate through the solution value, its gradient, and the
// nonlinear diffusivity. The Jacobian is read off the derivatives of
// the residual, so it is exact without any hand-coded linearization.
// libMesh's numeric types, such as VectorValue, accept dual-number
// components, so the residual is written with the same expressions
// used for real-valued assembly.
//
// Exactness of the Jacobian is visible in the quadratic convergence of
// Newton's method. Setting verify_analytic_jacobians to a positive
// tolerance additionally compares each element Jacobian against a
// finite-difference approximation.

// C++ includes
#include <cmath>
#include <memory>

// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fem_system.h"
#include "libmesh/fem_context.h"
#include "libmesh/steady_solver.h"
#include "libmesh/diff_solver.h"
#include "libmesh/dof_map.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/zero_function.h"
#include "libmesh/exact_solution.h"
#include "libmesh/getpot.h"
#include "libmesh/enum_solver_package.h"
#include "libmesh/int_range.h"
#include "libmesh/vector_value.h"

// MetaPhysicL includes
#ifdef LIBMESH_HAVE_METAPHYSICL
#include "metaphysicl/dualnumber.h"
#include "metaphysicl/dynamicsparsenumberarray.h"
#include "metaphysicl/raw_type.h"
#endif

using namespace libMesh;

#if defined(LIBMESH_HAVE_METAPHYSICL) && !defined(LIBMESH_USE_COMPLEX_NUMBERS)

// A real number together with its derivatives with respect to the
// element degrees of freedom, stored sparsely by local dof index.
typedef MetaPhysicL::DualNumber<Real, MetaPhysicL::DynamicSparseNumberArray<Real, unsigned int>>
    ADReal;

// The manufactured solution
Number
exact_value(const Point & p, const Parameters &, const std::string &, const std::string &)
{
  return std::sin(libMesh::pi * p(0)) * std::sin(libMesh::pi * p(1));
}

Gradient
exact_gradient(const Point & p, const Parameters &, const std::string &, const std::string &)
{
  const Real pi = libMesh::pi;
  return Gradient(pi * std::cos(pi * p(0)) * std::sin(pi * p(1)),
                  pi * std::sin(pi * p(0)) * std::cos(pi * p(1)));
}

// The source term f = -div((1 + u^2) grad u) for the manufactured solution,
// expanded as -(1 + u^2) lap(u) - 2 u |grad u|^2 with lap(u) = -2 pi^2 u.
Real
forcing(const Point & p)
{
  const Real u = libmesh_real(exact_value(p, Parameters(), "", ""));
  const Gradient grad_u = exact_gradient(p, Parameters(), "", "");
  return 2 * libMesh::pi * libMesh::pi * (1 + u * u) * u - 2 * u * grad_u.norm_sq();
}

/**
 * An FEMSystem for -div((1 + u^2) grad u) = f whose element Jacobian
 * is computed by forward-mode automatic differentiation of the
 * element residual.
 */
class NonlinearDiffusionSystem : public FEMSystem
{
public:
  NonlinearDiffusionSystem(EquationSystems & es, const std::string & name, const unsigned int number)
    : FEMSystem(es, name, number)
  {
  }

protected:
  virtual void init_data() override;

  virtual void init_context(DiffContext & context) override;

  virtual bool element_time_derivative(bool request_jacobian, DiffContext & context) override;

private:
  unsigned int _u_var;
};

void
NonlinearDiffusionSystem::init_data()
{
  _u_var = this->add_variable("u", SECOND, LAGRANGE);

  // u = 0 on all four sides of the square
  ZeroFunction<Number> zero;
  this->get_dof_map().add_dirichlet_boundary(
      DirichletBoundary({0, 1, 2, 3}, {_u_var}, zero, LOCAL_VARIABLE_ORDER));

  FEMSystem::init_data();
}

void
NonlinearDiffusionSystem::init_context(DiffContext & context)
{
  FEMContext & c = cast_ref<FEMContext &>(context);

  FEBase * fe = nullptr;
  c.get_element_fe(_u_var, fe);
  fe->get_JxW();
  fe->get_phi();
  fe->get_dphi();
  fe->get_xyz();

  // There are no side integrals.
  FEBase * side_fe = nullptr;
  c.get_side_fe(_u_var, side_fe);
  side_fe->get_nothing();

  FEMSystem::init_context(context);
}

bool
NonlinearDiffusionSystem::element_time_derivative(bool request_jacobian, DiffContext & context)
{
  FEMContext & c = cast_ref<FEMContext &>(context);

  FEBase * fe = nullptr;
  c.get_element_fe(_u_var, fe);
  const std::vector<Real> & JxW = fe->get_JxW();
  const std::vector<std::vector<Real>> & phi = fe->get_phi();
  const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();
  const std::vector<Point> & xyz = fe->get_xyz();

  const unsigned int n_dofs = c.n_dof_indices(_u_var);
  DenseSubVector<Number> & F = c.get_elem_residual(_u_var);
  DenseSubMatrix<Number> & K = c.get_elem_jacobian(_u_var, _u_var);

  // Seed each element degree of freedom with a unit derivative with
  // respect to itself.
  const DenseSubVector<Number> & u_coefs = c.get_elem_solution(_u_var);
  std::vector<ADReal> u_dofs(n_dofs);
  for (const auto j : make_range(n_dofs))
    {
      u_dofs[j] = u_coefs(j);
      u_dofs[j].derivatives().insert(j) = 1;
    }

  // FEMSystem's steady residual convention is F(u) = 0 with F the
  // weak form of div((1 + u^2) grad u) + f.
  std::vector<ADReal> residual(n_dofs, 0);
  for (const auto qp : index_range(JxW))
    {
      ADReal u = 0;
      VectorValue<ADReal> grad_u;
      for (const auto j : make_range(n_dofs))
        {
          u += u_dofs[j] * phi[j][qp];
          grad_u += u_dofs[j] * dphi[j][qp];
        }

      const ADReal diffusivity = 1 + u * u;
      const Real f = forcing(xyz[qp]);

      for (const auto i : make_range(n_dofs))
        residual[i] += JxW[qp] * (-diffusivity * (grad_u * dphi[i][qp]) + f * phi[i][qp]);
    }

  // The residual is the value part of each dual number, and its
  // derivatives are the corresponding Jacobian row.
  for (const auto i : make_range(n_dofs))
    {
      F(i) += MetaPhysicL::raw_value(residual[i]);

      if (request_jacobian)
        {
          const auto & derivatives = residual[i].derivatives();
          for (const auto k : index_range(derivatives.nude_indices()))
            K(i, derivatives.nude_indices()[k]) +=
                derivatives.nude_data()[k] * c.get_elem_solution_derivative();
        }
    }

  return request_jacobian;
}

#endif // LIBMESH_HAVE_METAPHYSICL && !LIBMESH_USE_COMPLEX_NUMBERS

int
main(int argc, char ** argv)
{
  // Initialize libMesh.
  LibMeshInit init(argc, argv);

  // This example requires a linear solver package.
  libmesh_example_requires(libMesh::default_solver_package() != INVALID_SOLVER_PACKAGE,
                           "--enable-petsc, --enable-trilinos, or --enable-eigen");

  // Dual numbers are provided by MetaPhysicL.
#ifndef LIBMESH_HAVE_METAPHYSICL
  libmesh_example_requires(false, "--enable-metaphysicl");
#elif defined(LIBMESH_USE_COMPLEX_NUMBERS)
  // The residual is written for real-valued dual numbers.
  libmesh_example_requires(false, "--disable-complex");
#else

  // The mesh is two-dimensional.
  libmesh_example_requires(2 <= LIBMESH_DIM, "2D support");

  // We use Dirichlet boundary conditions here.
#ifndef LIBMESH_ENABLE_DIRICHLET
  libmesh_example_requires(false, "--enable-dirichlet");
#endif

  // Parse the input file, allowing the command line to override it.
  GetPot infile("fem_system_ex6.in");
  infile.parse_command_line(argc, argv);

  const unsigned int grid_size = infile("grid_size", 16);

  Mesh mesh(init.comm());
  MeshTools::Generation::build_square(mesh, grid_size, grid_size, 0., 1., 0., 1., QUAD9);
  mesh.print_info();

  EquationSystems equation_systems(mesh);
  auto & system = equation_systems.add_system<NonlinearDiffusionSystem>("NonlinearDiffusion");
  system.time_solver = std::make_unique<SteadySolver>(system);

  // A positive tolerance checks each element Jacobian against a
  // finite-difference approximation.
  system.verify_analytic_jacobians = infile("verify_analytic_jacobians", 0.);

  equation_systems.init();

  DiffSolver & solver = *(system.time_solver->diff_solver());
  solver.quiet = false;
  solver.verbose = true;
  solver.max_nonlinear_iterations = infile("max_nonlinear_iterations", 10);
  solver.relative_residual_tolerance = infile("relative_residual_tolerance", 1.e-12);
  solver.relative_step_tolerance = infile("relative_step_tolerance", 0.);
  solver.absolute_residual_tolerance = infile("absolute_residual_tolerance", 0.);

  equation_systems.print_info();

  system.solve();

  ExactSolution exact_sol(equation_systems);
  exact_sol.attach_exact_value(exact_value);
  exact_sol.attach_exact_deriv(exact_gradient);
  exact_sol.compute_error("NonlinearDiffusion", "u");
  libMesh::out << "L2 error: " << exact_sol.l2_error("NonlinearDiffusion", "u") << std::endl
               << "H1 error: " << exact_sol.h1_error("NonlinearDiffusion", "u") << std::endl;

#endif // LIBMESH_HAVE_METAPHYSICL && !LIBMESH_USE_COMPLEX_NUMBERS

  return 0;
}
