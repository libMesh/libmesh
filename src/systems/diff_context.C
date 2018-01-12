// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


#include "libmesh/diff_context.h"
#include "libmesh/diff_system.h"
#include "libmesh/diff_system.h"
#include "libmesh/unsteady_solver.h"
#include "libmesh/auto_ptr.h" // libmesh_make_unique

namespace libMesh
{



DiffContext::DiffContext (const System & sys) :
  time(sys.time),
  system_time(sys.time),
  elem_solution_derivative(1.),
  elem_solution_rate_derivative(1.),
  elem_solution_accel_derivative(1.),
  fixed_solution_derivative(0.),
  _dof_indices_var(sys.n_vars()),
  _deltat(libmesh_nullptr),
  _system(sys),
  _is_adjoint(false)
{
  // Finally initialize solution/residual/jacobian data structures
  unsigned int nv = sys.n_vars();

  _elem_subsolutions.reserve(nv);
  _elem_subresiduals.reserve(nv);
  _elem_subjacobians.resize(nv);
  _elem_subsolution_rates.reserve(nv);
  _elem_subsolution_accels.reserve(nv);
  if (sys.use_fixed_solution)
    _elem_fixed_subsolutions.reserve(nv);

  // If the user resizes sys.qoi, it will invalidate us

  std::size_t n_qoi = sys.qoi.size();
  _elem_qoi.resize(n_qoi);
  _elem_qoi_derivative.resize(n_qoi);
  _elem_qoi_subderivatives.resize(n_qoi);
  for (std::size_t q=0; q != n_qoi; ++q)
    _elem_qoi_subderivatives[q].reserve(nv);

  for (unsigned int i=0; i != nv; ++i)
    {
      _elem_subsolutions.emplace_back(libmesh_make_unique<DenseSubVector<Number>>(_elem_solution));
      _elem_subresiduals.emplace_back(libmesh_make_unique<DenseSubVector<Number>>(_elem_residual));
      for (std::size_t q=0; q != n_qoi; ++q)
        _elem_qoi_subderivatives[q].emplace_back(libmesh_make_unique<DenseSubVector<Number>>(_elem_qoi_derivative[q]));
      _elem_subjacobians[i].reserve(nv);

      // Only make space for these if we're using DiffSystem
      // This is assuming *only* DiffSystem is using elem_solution_rate/accel
      const DifferentiableSystem * diff_system = dynamic_cast<const DifferentiableSystem *>(&sys);
      if (diff_system)
        {
          // Now, we only need these if the solver is unsteady
          if (!diff_system->get_time_solver().is_steady())
            {
              _elem_subsolution_rates.emplace_back(libmesh_make_unique<DenseSubVector<Number>>(_elem_solution_rate));

              // We only need accel space if the TimeSolver is second order
              const UnsteadySolver & time_solver = cast_ref<const UnsteadySolver &>(diff_system->get_time_solver());

              if (time_solver.time_order() >= 2 || !diff_system->get_second_order_vars().empty())
                _elem_subsolution_accels.emplace_back(libmesh_make_unique<DenseSubVector<Number>>(_elem_solution_accel));
            }
        }

      if (sys.use_fixed_solution)
        _elem_fixed_subsolutions.emplace_back(libmesh_make_unique<DenseSubVector<Number>>(_elem_fixed_solution));

      for (unsigned int j=0; j != nv; ++j)
        _elem_subjacobians[i].emplace_back(libmesh_make_unique<DenseSubMatrix<Number>>(_elem_jacobian));
    }
}



DiffContext::~DiffContext ()
{
}


void DiffContext::set_deltat_pointer(Real * dt)
{
  // We may actually want to be able to set this pointer to NULL, so
  // don't report an error for that.
  _deltat = dt;
}


Real DiffContext::get_deltat_value()
{
  libmesh_assert(_deltat);

  return *_deltat;
}


void DiffContext::add_localized_vector (NumericVector<Number> & localized_vector, const System & sys)
{
  // Make an empty pair keyed with a reference to this _localized_vector
  _localized_vectors[&localized_vector] = std::make_pair(DenseVector<Number>(), std::vector<std::unique_ptr<DenseSubVector<Number>>>());

  unsigned int nv = sys.n_vars();

  _localized_vectors[&localized_vector].second.reserve(nv);

  // Fill the DenseSubVector with nv copies of DenseVector
  for (unsigned int i=0; i != nv; ++i)
    _localized_vectors[&localized_vector].second.emplace_back(libmesh_make_unique<DenseSubVector<Number>>(_localized_vectors[&localized_vector].first));
}


DenseVector<Number> & DiffContext::get_localized_vector (const NumericVector<Number> & localized_vector)
{
  return _localized_vectors[&localized_vector].first;
}


const DenseVector<Number> & DiffContext::get_localized_vector (const NumericVector<Number> & localized_vector) const
{
  auto localized_vectors_it = _localized_vectors.find(&localized_vector);
  libmesh_assert(localized_vectors_it != _localized_vectors.end());
  return localized_vectors_it->second.first;
}


DenseSubVector<Number> & DiffContext::get_localized_subvector (const NumericVector<Number> & localized_vector, unsigned int var)
{
  return *_localized_vectors[&localized_vector].second[var];
}


const DenseSubVector<Number> & DiffContext::get_localized_subvector (const NumericVector<Number> & localized_vector, unsigned int var) const
{
  auto localized_vectors_it = _localized_vectors.find(&localized_vector);
  libmesh_assert(localized_vectors_it != _localized_vectors.end());
  return *localized_vectors_it->second.second[var];
}

} // namespace libMesh
