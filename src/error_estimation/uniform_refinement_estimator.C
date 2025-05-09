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


// Local Includes
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/equation_systems.h"
#include "libmesh/error_vector.h"
#include "libmesh/fe.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/quadrature.h"
#include "libmesh/system.h"
#include "libmesh/uniform_refinement_estimator.h"
#include "libmesh/partitioner.h"
#include "libmesh/tensor_tools.h"
#include "libmesh/enum_error_estimator_type.h"
#include "libmesh/enum_norm_type.h"
#include "libmesh/int_range.h"

// C++ includes
#include <algorithm> // for std::fill
#include <sstream>
#include <cstdlib> // *must* precede <cmath> for proper std:abs() on PGI, Sun Studio CC
#include <cmath>    // for sqrt
#include <memory>

#ifdef LIBMESH_ENABLE_AMR

namespace libMesh
{

//-----------------------------------------------------------------
// ErrorEstimator implementations

UniformRefinementEstimator::UniformRefinementEstimator() :
    ErrorEstimator(),
    number_h_refinements(1),
    number_p_refinements(0),
    _extra_order(1)
{
  error_norm = H1;
}



ErrorEstimatorType UniformRefinementEstimator::type() const
{
  return UNIFORM_REFINEMENT;
}


void UniformRefinementEstimator::estimate_error (const System & _system,
                                                 ErrorVector & error_per_cell,
                                                 const NumericVector<Number> * solution_vector,
                                                 bool estimate_parent_error)
{
  LOG_SCOPE("estimate_error()", "UniformRefinementEstimator");
  std::map<const System *, const NumericVector<Number> *> solution_vectors;
  solution_vectors[&_system] = solution_vector;
  this->_estimate_error (nullptr,
                         &_system,
                         &error_per_cell,
                         nullptr,
                         nullptr,
                         &solution_vectors,
                         estimate_parent_error);
}

void UniformRefinementEstimator::estimate_errors (const EquationSystems & _es,
                                                  ErrorVector & error_per_cell,
                                                  const std::map<const System *, SystemNorm> & error_norms,
                                                  const std::map<const System *, const NumericVector<Number> *> * solution_vectors,
                                                  bool estimate_parent_error)
{
  LOG_SCOPE("estimate_errors()", "UniformRefinementEstimator");
  this->_estimate_error (&_es,
                         nullptr,
                         &error_per_cell,
                         nullptr,
                         &error_norms,
                         solution_vectors,
                         estimate_parent_error);
}

void UniformRefinementEstimator::estimate_errors (const EquationSystems & _es,
                                                  ErrorMap & errors_per_cell,
                                                  const std::map<const System *, const NumericVector<Number> *> * solution_vectors,
                                                  bool estimate_parent_error)
{
  LOG_SCOPE("estimate_errors()", "UniformRefinementEstimator");
  this->_estimate_error (&_es,
                         nullptr,
                         nullptr,
                         &errors_per_cell,
                         nullptr,
                         solution_vectors,
                         estimate_parent_error);
}

void UniformRefinementEstimator::_estimate_error (const EquationSystems * _es,
                                                  const System * _system,
                                                  ErrorVector * error_per_cell,
                                                  ErrorMap * errors_per_cell,
                                                  const std::map<const System *, SystemNorm> * _error_norms,
                                                  const std::map<const System *, const NumericVector<Number> *> * solution_vectors,
                                                  bool)
{
  // Get a vector of the Systems we're going to work on,
  // and set up a error_norms map if necessary
  std::vector<System *> system_list;
  auto error_norms = std::make_unique<std::map<const System *, SystemNorm>>();

  if (_es)
    {
      libmesh_assert(!_system);
      libmesh_assert(_es->n_systems());
      _system = &(_es->get_system(0));
      libmesh_assert_equal_to (&(_system->get_equation_systems()), _es);

      libmesh_assert(_es->n_systems());
      for (auto i : make_range(_es->n_systems()))
        // We have to break the rules here, because we can't refine a const System
        system_list.push_back(const_cast<System *>(&(_es->get_system(i))));

      // If we're computing one vector, we need to know how to scale
      // each variable's contributions to it.
      if (_error_norms)
        {
          libmesh_assert(!errors_per_cell);
        }
      else
        // If we're computing many vectors, we just need to know which
        // variables to skip
        {
          libmesh_assert (errors_per_cell);

          _error_norms = error_norms.get();

          for (auto i : make_range(_es->n_systems()))
            {
              const System & sys = _es->get_system(i);
              unsigned int n_vars = sys.n_vars();

              std::vector<Real> weights(n_vars, 0.0);
              for (unsigned int v = 0; v != n_vars; ++v)
                {
                  if (!errors_per_cell->count(std::make_pair(&sys, v)))
                    continue;

                  weights[v] = 1.0;
                }
              (*error_norms)[&sys] =
                SystemNorm(std::vector<FEMNormType>(n_vars, error_norm.type(0)),
                           weights);
            }
        }
    }
  else
    {
      libmesh_assert(_system);
      // We have to break the rules here, because we can't refine a const System
      system_list.push_back(const_cast<System *>(_system));

      libmesh_assert(!_error_norms);
      (*error_norms)[_system] = error_norm;
      _error_norms = error_norms.get();
    }

  // An EquationSystems reference will be convenient.
  // We have to break the rules here, because we can't refine a const System
  EquationSystems & es =
    const_cast<EquationSystems &>(_system->get_equation_systems());

  // The current mesh
  MeshBase & mesh = es.get_mesh();

  // The dimensionality of the mesh
  const unsigned int dim = mesh.mesh_dimension();

  // Resize the error_per_cell vectors to be
  // the number of elements, initialize them to 0.
  if (error_per_cell)
    {
      error_per_cell->clear();
      error_per_cell->resize (mesh.max_elem_id(), 0.);
    }
  else
    {
      libmesh_assert(errors_per_cell);
      for (const auto & pr : *errors_per_cell)
        {
          ErrorVector * e = pr.second.get();
          e->clear();
          e->resize(mesh.max_elem_id(), 0.);
        }
    }

  // We'll want to back up all coarse grid vectors
  std::vector<std::map<std::string, std::unique_ptr<NumericVector<Number>>>> coarse_vectors(system_list.size());
  std::vector<std::unique_ptr<NumericVector<Number>>> coarse_solutions(system_list.size());
  std::vector<std::unique_ptr<NumericVector<Number>>> coarse_local_solutions(system_list.size());
  // And make copies of projected solutions
  std::vector<std::unique_ptr<NumericVector<Number>>> projected_solutions(system_list.size());

  // And we'll need to temporarily change solution projection settings
  std::vector<bool> old_projection_settings(system_list.size());

  // And it'll be best to avoid any repartitioning
  std::unique_ptr<Partitioner> old_partitioner(mesh.partitioner().release());

  // And we can't allow any renumbering
  const bool old_renumbering_setting = mesh.allow_renumbering();
  mesh.allow_renumbering(false);

  for (auto i : index_range(system_list))
    {
      System & system = *system_list[i];

      // Check for valid error_norms
      libmesh_assert (_error_norms->find(&system) !=
                      _error_norms->end());

      // Back up the solution vector
      coarse_solutions[i] = system.solution->clone();
      coarse_local_solutions[i] = system.current_local_solution->clone();

      // Back up all other coarse grid vectors
      for (System::vectors_iterator vec = system.vectors_begin(); vec !=
             system.vectors_end(); ++vec)
        {
          // The (string) name of this vector
          const std::string & var_name = vec->first;

          coarse_vectors[i][var_name] = vec->second->clone();
        }

      // Use a non-standard solution vector if necessary
      if (solution_vectors &&
          solution_vectors->find(&system) != solution_vectors->end() &&
          solution_vectors->find(&system)->second &&
          solution_vectors->find(&system)->second != system.solution.get())
        {
          NumericVector<Number> * newsol =
            const_cast<NumericVector<Number> *>
            (solution_vectors->find(&system)->second);
          newsol->swap(*system.solution);
          system.update();
        }

      // Make sure the solution is projected when we refine the mesh
      old_projection_settings[i] = system.project_solution_on_reinit();
      system.project_solution_on_reinit() = true;
    }

  // Find the number of coarse mesh elements, to make it possible
  // to find correct coarse elem ids later
  const dof_id_type max_coarse_elem_id = mesh.max_elem_id();
#ifndef NDEBUG
  // n_coarse_elem is only used in an assertion later so
  // avoid declaring it unless asserts are active.
  const dof_id_type n_coarse_elem = mesh.n_elem();
#endif

  // Uniformly refine the mesh
  MeshRefinement mesh_refinement(mesh);

  libmesh_assert (number_h_refinements > 0 || number_p_refinements > 0);

  // FIXME: this may break if there is more than one System
  // on this mesh but estimate_error was still called instead of
  // estimate_errors
  for (unsigned int i = 0; i != number_h_refinements; ++i)
    {
      mesh_refinement.uniformly_refine(1);
      es.reinit();
    }

  for (unsigned int i = 0; i != number_p_refinements; ++i)
    {
      mesh_refinement.uniformly_p_refine(1);
      es.reinit();
    }

  for (auto i : index_range(system_list))
    {
      System & system = *system_list[i];

      // Copy the projected coarse grid solutions, which will be
      // overwritten by solve()
      projected_solutions[i] = NumericVector<Number>::build(system.comm());
      projected_solutions[i]->init(system.solution->size(),
                                   system.solution->local_size(),
                                   system.get_dof_map().get_send_list(),
                                   true, GHOSTED);
      system.solution->localize(*projected_solutions[i],
                                system.get_dof_map().get_send_list());
    }

  // Are we doing a forward or an adjoint solve?
  bool solve_adjoint = false;
  if (solution_vectors)
    {
      System * sys = system_list[0];
      libmesh_assert (solution_vectors->find(sys) !=
                      solution_vectors->end());
      const NumericVector<Number> * vec = solution_vectors->find(sys)->second;
      for (auto j : make_range(sys->n_qois()))
        {
          std::ostringstream adjoint_name;
          adjoint_name << "adjoint_solution" << j;

          if (vec == sys->request_vector(adjoint_name.str()))
            {
              solve_adjoint = true;
              break;
            }
        }
    }

  // Get the uniformly refined solution.

  if (_es)
    {
      // Even if we had a decent preconditioner, valid matrix etc. before
      // refinement, we don't any more.
      for (auto i : make_range(_es->n_systems()))
        es.get_system(i).disable_cache();

      // No specified vectors == forward solve
      if (!solution_vectors)
        es.solve();
      else
        {
          libmesh_assert_equal_to (solution_vectors->size(), es.n_systems());
          libmesh_assert (solution_vectors->find(system_list[0]) !=
                          solution_vectors->end());
          libmesh_assert(solve_adjoint ||
                         (solution_vectors->find(system_list[0])->second ==
                          system_list[0]->solution.get()) ||
                         !solution_vectors->find(system_list[0])->second);

#ifdef DEBUG
          for (const auto & sys : system_list)
            {
              libmesh_assert (solution_vectors->find(sys) !=
                              solution_vectors->end());
              const NumericVector<Number> * vec = solution_vectors->find(sys)->second;
              if (solve_adjoint)
                {
                  bool found_vec = false;
                  for (auto j : make_range(sys->n_qois()))
                    {
                      std::ostringstream adjoint_name;
                      adjoint_name << "adjoint_solution" << j;

                      if (vec == sys->request_vector(adjoint_name.str()))
                        {
                          found_vec = true;
                          break;
                        }
                    }
                  libmesh_assert(found_vec);
                }
              else
                libmesh_assert(vec == sys->solution.get() || !vec);
            }
#endif

          if (solve_adjoint)
            {
              std::vector<unsigned int> adjs(system_list.size(),
                                             libMesh::invalid_uint);
              // Set up proper initial guesses
              for (auto i : index_range(system_list))
                {
                  System * sys = system_list[i];
                  libmesh_assert (solution_vectors->find(sys) !=
                                  solution_vectors->end());
                  const NumericVector<Number> * vec = solution_vectors->find(sys)->second;
                  for (auto j : make_range(sys->n_qois()))
                    {
                      std::ostringstream adjoint_name;
                      adjoint_name << "adjoint_solution" << j;

                      if (vec == sys->request_vector(adjoint_name.str()))
                        {
                          adjs[i] = j;
                          break;
                        }
                    }
                  libmesh_assert_not_equal_to (adjs[i], libMesh::invalid_uint);
                  sys->get_adjoint_solution(adjs[i]) = *sys->solution;
                }

              es.adjoint_solve();

              // Put the adjoint_solution into solution for
              // comparisons
              for (auto i : index_range(system_list))
                {
                  system_list[i]->get_adjoint_solution(adjs[i]).swap(*system_list[i]->solution);
                  system_list[i]->update();
                }
            }
          else
            es.solve();
        }
    }
  else
    {
      System * sys = system_list[0];

      // Even if we had a decent preconditioner, valid matrix etc. before
      // refinement, we don't any more.
      sys->disable_cache();

      // No specified vectors == forward solve
      if (!solution_vectors)
        sys->solve();
      else
        {
          libmesh_assert (solution_vectors->find(sys) !=
                          solution_vectors->end());

          const NumericVector<Number> * vec = solution_vectors->find(sys)->second;

          libmesh_assert(solve_adjoint ||
                         (solution_vectors->find(sys)->second ==
                          sys->solution.get()) ||
                         !solution_vectors->find(sys)->second);

          if (solve_adjoint)
            {
              unsigned int adj = libMesh::invalid_uint;
              for (unsigned int j=0, n_qois = sys->n_qois();
                   j != n_qois; ++j)
                {
                  std::ostringstream adjoint_name;
                  adjoint_name << "adjoint_solution" << j;

                  if (vec == sys->request_vector(adjoint_name.str()))
                    {
                      adj = j;
                      break;
                    }
                }
              libmesh_assert_not_equal_to (adj, libMesh::invalid_uint);

              // Set up proper initial guess
              sys->get_adjoint_solution(adj) = *sys->solution;
              sys->adjoint_solve();
              // Put the adjoint_solution into solution for
              // comparisons
              sys->get_adjoint_solution(adj).swap(*sys->solution);
              sys->update();
            }
          else
            sys->solve();
        }
    }

  // Get the error in the uniformly refined solution(s).
  for (auto sysnum : index_range(system_list))
    {
      System & system = *system_list[sysnum];

      unsigned int n_vars = system.n_vars();

      DofMap & dof_map = system.get_dof_map();

      const SystemNorm & system_i_norm =
        _error_norms->find(&system)->second;

      NumericVector<Number> * projected_solution = projected_solutions[sysnum].get();

      // Loop over all the variables in the system
      for (unsigned int var=0; var<n_vars; var++)
        {
          // Get the error vector to fill for this system and variable
          ErrorVector * err_vec = error_per_cell;
          if (!err_vec)
            {
              libmesh_assert(errors_per_cell);
              err_vec =
                (*errors_per_cell)[std::make_pair(&system,var)].get();
            }

          // The type of finite element to use for this variable
          const FEType & fe_type = dof_map.variable_type (var);

          // Finite element object for each fine element
          std::unique_ptr<FEBase> fe (FEBase::build (dim, fe_type));

          // Build and attach an appropriate quadrature rule
          std::unique_ptr<QBase> qrule = fe_type.default_quadrature_rule(dim, _extra_order);
          fe->attach_quadrature_rule (qrule.get());

          const std::vector<Real> &  JxW = fe->get_JxW();
          const std::vector<std::vector<Real>> & phi = fe->get_phi();
          const std::vector<std::vector<RealGradient>> & dphi =
            fe->get_dphi();
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
          const std::vector<std::vector<RealTensor>> & d2phi =
            fe->get_d2phi();
#endif

          // The global DOF indices for the fine element
          std::vector<dof_id_type> dof_indices;

          // Iterate over all the active elements in the fine mesh
          // that live on this processor.
          for (const auto & elem : mesh.active_local_element_ptr_range())
            {
              // Find the element id for the corresponding coarse grid element
              const Elem * coarse = elem;
              dof_id_type e_id = coarse->id();
              while (e_id >= max_coarse_elem_id)
                {
                  libmesh_assert (coarse->parent());
                  coarse = coarse->parent();
                  e_id = coarse->id();
                }

              Real L2normsq = 0., H1seminormsq = 0.;
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
              Real H2seminormsq = 0.;
#endif

              // reinitialize the element-specific data
              // for the current element
              fe->reinit (elem);

              // Get the local to global degree of freedom maps
              dof_map.dof_indices (elem, dof_indices, var);

              // The number of quadrature points
              const unsigned int n_qp = qrule->n_points();

              // The number of shape functions
              const unsigned int n_sf =
                cast_int<unsigned int>(dof_indices.size());

              //
              // Begin the loop over the Quadrature points.
              //
              for (unsigned int qp=0; qp<n_qp; qp++)
                {
                  Number u_fine = 0., u_coarse = 0.;

                  Gradient grad_u_fine, grad_u_coarse;
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
                  Tensor grad2_u_fine, grad2_u_coarse;
#endif

                  // Compute solution values at the current
                  // quadrature point.  This requires a sum
                  // over all the shape functions evaluated
                  // at the quadrature point.
                  for (unsigned int i=0; i<n_sf; i++)
                    {
                      u_fine            += phi[i][qp]*system.current_solution (dof_indices[i]);
                      u_coarse          += phi[i][qp]*(*projected_solution) (dof_indices[i]);
                      grad_u_fine       += dphi[i][qp]*system.current_solution (dof_indices[i]);
                      grad_u_coarse     += dphi[i][qp]*(*projected_solution) (dof_indices[i]);
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
                      grad2_u_fine      += d2phi[i][qp]*system.current_solution (dof_indices[i]);
                      grad2_u_coarse    += d2phi[i][qp]*(*projected_solution) (dof_indices[i]);
#endif
                    }

                  // Compute the value of the error at this quadrature point
                  const Number val_error = u_fine - u_coarse;

                  // Add the squares of the error to each contribution
                  if (system_i_norm.type(var) == L2 ||
                      system_i_norm.type(var) == H1 ||
                      system_i_norm.type(var) == H2)
                    {
                      L2normsq += JxW[qp] * system_i_norm.weight_sq(var) *
                        TensorTools::norm_sq(val_error);
                      libmesh_assert_greater_equal (L2normsq, 0.);
                    }


                  // Compute the value of the error in the gradient at this
                  // quadrature point
                  if (system_i_norm.type(var) == H1 ||
                      system_i_norm.type(var) == H2 ||
                      system_i_norm.type(var) == H1_SEMINORM)
                    {
                      Gradient grad_error = grad_u_fine - grad_u_coarse;

                      H1seminormsq += JxW[qp] * system_i_norm.weight_sq(var) *
                        grad_error.norm_sq();
                      libmesh_assert_greater_equal (H1seminormsq, 0.);
                    }

                  // Compute the value of the error in the hessian at this
                  // quadrature point
                  if (system_i_norm.type(var) == H2 ||
                      system_i_norm.type(var) == H2_SEMINORM)
                    {
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
                      Tensor grad2_error = grad2_u_fine - grad2_u_coarse;

                      H2seminormsq += JxW[qp] * system_i_norm.weight_sq(var) *
                        grad2_error.norm_sq();
                      libmesh_assert_greater_equal (H2seminormsq, 0.);
#else
                      libmesh_error_msg
                        ("libMesh was not configured with --enable-second");
#endif
                    }
                } // end qp loop

              if (system_i_norm.type(var) == L2 ||
                  system_i_norm.type(var) == H1 ||
                  system_i_norm.type(var) == H2)
                (*err_vec)[e_id] +=
                  static_cast<ErrorVectorReal>(L2normsq);
              if (system_i_norm.type(var) == H1 ||
                  system_i_norm.type(var) == H2 ||
                  system_i_norm.type(var) == H1_SEMINORM)
                (*err_vec)[e_id] +=
                  static_cast<ErrorVectorReal>(H1seminormsq);

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
              if (system_i_norm.type(var) == H2 ||
                  system_i_norm.type(var) == H2_SEMINORM)
                (*err_vec)[e_id] +=
                  static_cast<ErrorVectorReal>(H2seminormsq);
#endif
            } // End loop over active local elements
        } // End loop over variables

      // Don't bother projecting the solution; we'll restore from backup
      // after coarsening
      system.project_solution_on_reinit() = false;
    }


  // Uniformly coarsen the mesh, without projecting the solution
  libmesh_assert (number_h_refinements > 0 || number_p_refinements > 0);

  for (unsigned int i = 0; i != number_h_refinements; ++i)
    {
      mesh_refinement.uniformly_coarsen(1);
      // FIXME - should the reinits here be necessary? - RHS
      es.reinit();
    }

  for (unsigned int i = 0; i != number_p_refinements; ++i)
    {
      mesh_refinement.uniformly_p_coarsen(1);
      es.reinit();
    }

  // We should be back where we started
  libmesh_assert_equal_to (n_coarse_elem, mesh.n_elem());

  // Each processor has now computed the error contributions
  // for its local elements.  We need to sum the vector
  // and then take the square-root of each component.  Note
  // that we only need to sum if we are running on multiple
  // processors, and we only need to take the square-root
  // if the value is nonzero.  There will in general be many
  // zeros for the inactive elements.

  if (error_per_cell)
    {
      // First sum the vector of estimated error values
      this->reduce_error(*error_per_cell, es.comm());

      // Compute the square-root of each component.
      LOG_SCOPE("std::sqrt()", "UniformRefinementEstimator");
      for (auto & val : *error_per_cell)
        if (val != 0.)
          val = std::sqrt(val);
    }
  else
    {
      for (const auto & pr : *errors_per_cell)
        {
          ErrorVector & e = *(pr.second);
          // First sum the vector of estimated error values
          this->reduce_error(e, es.comm());

          // Compute the square-root of each component.
          LOG_SCOPE("std::sqrt()", "UniformRefinementEstimator");
          for (auto & val : e)
            if (val != 0.)
              val = std::sqrt(val);
        }
    }

  // Restore old solutions and clean up the heap
  for (auto i : index_range(system_list))
    {
      System & system = *system_list[i];

      system.project_solution_on_reinit() = old_projection_settings[i];

      // Restore the coarse solution vectors and delete their copies
      *system.solution = *coarse_solutions[i];
      *system.current_local_solution = *coarse_local_solutions[i];

      for (System::vectors_iterator vec = system.vectors_begin(); vec !=
             system.vectors_end(); ++vec)
        {
          // The (string) name of this vector
          const std::string & var_name = vec->first;

          system.get_vector(var_name) = *coarse_vectors[i][var_name];

          coarse_vectors[i][var_name]->clear();
        }
    }

  // Restore old partitioner and renumbering settings
  mesh.partitioner().reset(old_partitioner.release());
  mesh.allow_renumbering(old_renumbering_setting);
}

} // namespace libMesh

#endif // #ifdef LIBMESH_ENABLE_AMR
