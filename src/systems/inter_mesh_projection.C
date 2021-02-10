// The libMesh Finite Element Library.
// Copyright (C) 2002-2020 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// Local includes
#include "libmesh/inter_mesh_projection.h"

namespace libMesh
{
GradientMeshFunction::GradientMeshFunction(MeshFunction * _mesh_function):
  mesh_function(libmesh_make_unique<MeshFunction>(*_mesh_function))
{
  libmesh_experimental();
}

void GradientMeshFunction::operator() (const Point & p, const Real, DenseVector<Gradient> & output)
{
  Real time = 0.0;
  mesh_function->gradient(p, time, output.get_values());
  return;
}

InterMeshProjection::InterMeshProjection(System & _from_system, System & _to_system) :
  from_system(_from_system),
  to_system(_to_system)
{
  libmesh_experimental();
}

void InterMeshProjection::project_system_vectors()
{
  // Number of vectors to be projected
  libmesh_assert_equal_to (to_system.n_vectors(), from_system.n_vectors());

  // Number of variable components in each vector
  unsigned int n_vars = from_system.n_vars();
  libmesh_assert_equal_to (to_system.n_vars(), n_vars);

  // We are going to use the multi-variable MeshFunction, so we can pass
  // a single vector of variables rather than have a MeshFunction for each variable
  std::vector<unsigned int> variables_vector;

  for (unsigned int j = 0; j != n_vars; ++j)
    {
      libmesh_assert_equal_to (to_system.variable_name(j), from_system.variable_name(j));
      variables_vector.push_back(j);
    }

  // Any system holds the solution along with the other vectors system.vectors
  // We will first project the solution and then move to the system.vectors

  // Construct local version of the current system vector
  // This has to be a serial vector
  // Roy's FIXME: Technically it just has to be a ghosted vector whose algebraically
  // ghosted values cover a domain which is a superset of the to-system's domain ...
  // that's hard to do and we can skip it until the poor scalability bites someone.
  std::unique_ptr<NumericVector<Number>> solution_vector_serial = NumericVector<Number>::build(from_system.comm());
  solution_vector_serial->init(from_system.solution->size(), true, SERIAL);

  std::vector<Number> solution_vector;
  from_system.update_global_solution(solution_vector);
  (*solution_vector_serial) = solution_vector;

  // Construct a MeshFunction for the solution
  MeshFunction mesh_func_solution(from_system.get_equation_systems(), *solution_vector_serial, from_system.get_dof_map(), variables_vector);

  mesh_func_solution.init();

  // For some element types (say C1) we also need to pass a gradient evaluation MeshFunction
  // To do this evaluate, a new shim class GradientMeshFunction has been added which redirects
  // gptr::operator evaluations inside projection methods into MeshFunction::gradient calls.
  GradientMeshFunction gptr_solution(&mesh_func_solution);
  gptr_solution.init();

  to_system.project_vector(*to_system.solution, &mesh_func_solution, &gptr_solution);

  // Now loop over the vectors in system.vectors (includes old_nonlin_sol, rhs, adjoints, adjoint_rhs, sensitivity_rhs)
  for (System::vectors_iterator vec = from_system.vectors_begin(), vec_end = from_system.vectors_end(); vec != vec_end; ++vec)
    {
      // The name of this vector
      const std::string & vec_name = vec->first;

      // Construct local version of the current system vector
      // This has to be a serial vector
      // Roy's FIXME: Technically it just has to be a ghosted vector whose algebraically
      // ghosted values cover a domain which is a superset of the to-system's domain ...
      // that's hard to do and we can skip it until the poor scalability bites someone.
      std::unique_ptr<NumericVector<Number>> current_vector_proxy = NumericVector<Number>::build(from_system.comm());
      current_vector_proxy->init(from_system.get_vector(vec_name).size(), true, SERIAL);

      from_system.get_vector(vec_name).localize(*current_vector_proxy);

      // Construct a MeshFunction for the current component
      MeshFunction mesh_func(from_system.get_equation_systems(), *current_vector_proxy, from_system.get_dof_map(), variables_vector);
      mesh_func.init();

      // Project the current system vector, you need to check if this vector is an adjoint to pass
      // the right options to project_vector
      GradientMeshFunction gptr(&mesh_func);
      gptr.init();

      // The fourth argument here is whether the vector is an adjoint solution or not, we will be getting that information
      // via the from_system instead of the to_system in case the user has not set the System::_vector_is_adjoint map to true.
      to_system.project_vector(to_system.get_vector(vec_name), &mesh_func, &gptr, from_system.vector_is_adjoint(vec_name));

    }
  // End loop over the vectors in the system

}
// End InterMeshProjection::project_system_vectors
}
// End namespace libMesh
