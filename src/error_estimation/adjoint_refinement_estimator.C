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

// C++ includes
#include <algorithm> // for std::fill
#include <cstdlib> // *must* precede <cmath> for proper std:abs() on PGI, Sun Studio CC
#include <cmath>    // for sqrt
#include <set>
#include <sstream> // for ostringstream

// Local Includes
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/equation_systems.h"
#include "libmesh/error_vector.h"
#include "libmesh/fe.h"
#include "libmesh/fe_interface.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/quadrature.h"
#include "libmesh/system.h"
#include "libmesh/implicit_system.h"
#include "libmesh/partitioner.h"
#include "libmesh/adjoint_refinement_estimator.h"

#include LIBMESH_INCLUDE_UNORDERED_MAP
#include LIBMESH_INCLUDE_UNORDERED_SET

#ifdef LIBMESH_ENABLE_AMR

namespace libMesh
{

//-----------------------------------------------------------------
// AdjointRefinementErrorEstimator

// As of 10/2/2012, this function implements a 'brute-force' adjoint based QoI
// error estimator, using the following relationship:
// Q(u) - Q(u_h) \aprrox - R( (u_h)_(h/2), z_(h/2) ) .
// where: Q(u) is the true QoI
// u_h is the approximate primal solution on the current FE space
// Q(u_h) is the approximate QoI
// z_(h/2) is the adjoint corresponding to Q, on a richer FE space
// (u_h)_(h/2) is a projection of the primal solution on the richer FE space
// By richer FE space, we mean a grid that has been refined once and a polynomial order
// that has been increased once, i.e. one h and one p refinement

// Both a global QoI error estimate and element wise error indicators are included
// Note that the element wise error indicators slightly over estimate the error in
// each element

void AdjointRefinementEstimator::estimate_error (const System & _system,
                                                 ErrorVector & error_per_cell,
                                                 const NumericVector<Number> * solution_vector,
                                                 bool /*estimate_parent_error*/)
{
  // We have to break the rules here, because we can't refine a const System
  System & system = const_cast<System &>(_system);

  // An EquationSystems reference will be convenient.
  EquationSystems & es = system.get_equation_systems();

  // The current mesh
  MeshBase & mesh = es.get_mesh();

  // Get coarse grid adjoint solutions.  This should be a relatively
  // quick (especially with preconditioner reuse) way to get a good
  // initial guess for the fine grid adjoint solutions.  More
  // importantly, subtracting off a coarse adjoint approximation gives
  // us better local error indication, and subtracting off *some* lift
  // function is necessary for correctness if we have heterogeneous
  // adjoint Dirichlet conditions.

  // Solve the adjoint problem(s) on the coarse FE space
  // Only if the user didn't already solve it for us
  if (!system.is_adjoint_already_solved())
    system.adjoint_solve(_qoi_set);

  // Loop over all the adjoint problems and, if any have heterogenous
  // Dirichlet conditions, get the corresponding coarse lift
  // function(s)
  for (std::size_t j=0; j != system.qoi.size(); j++)
    {
      // Skip this QoI if it is not in the QoI Set or if there are no
      // heterogeneous Dirichlet boundaries for it
      if (_qoi_set.has_index(j) &&
          system.get_dof_map().has_adjoint_dirichlet_boundaries(j))
        {
          std::ostringstream liftfunc_name;
          liftfunc_name << "adjoint_lift_function" << j;
          NumericVector<Number> & liftvec =
            system.add_vector(liftfunc_name.str());

          system.get_dof_map().enforce_constraints_exactly
            (system, &liftvec, true);
        }
    }

  // We'll want to back up all coarse grid vectors
  std::map<std::string, NumericVector<Number> *> coarse_vectors;
  for (System::vectors_iterator vec = system.vectors_begin(); vec !=
         system.vectors_end(); ++vec)
    {
      // The (string) name of this vector
      const std::string & var_name = vec->first;

      coarse_vectors[var_name] = vec->second->clone().release();
    }
  // Back up the coarse solution and coarse local solution
  NumericVector<Number> * coarse_solution =
    system.solution->clone().release();
  NumericVector<Number> * coarse_local_solution =
    system.current_local_solution->clone().release();

  // And we'll need to temporarily change solution projection settings
  bool old_projection_setting;
  old_projection_setting = system.project_solution_on_reinit();

  // Make sure the solution is projected when we refine the mesh
  system.project_solution_on_reinit() = true;

  // And it'll be best to avoid any repartitioning
  UniquePtr<Partitioner> old_partitioner(mesh.partitioner().release());

  // And we can't allow any renumbering
  const bool old_renumbering_setting = mesh.allow_renumbering();
  mesh.allow_renumbering(false);

  // Use a non-standard solution vector if necessary
  if (solution_vector && solution_vector != system.solution.get())
    {
      NumericVector<Number> * newsol =
        const_cast<NumericVector<Number> *> (solution_vector);
      newsol->swap(*system.solution);
      system.update();
    }

  // Resize the error_per_cell vector to be
  // the number of elements, initialized to 0.
  error_per_cell.clear();
  error_per_cell.resize (mesh.max_elem_id(), 0.);

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

  // Copy the projected coarse grid solutions, which will be
  // overwritten by solve()
  std::vector<NumericVector<Number> *> coarse_adjoints;
  for (std::size_t j=0; j != system.qoi.size(); j++)
    {
      if (_qoi_set.has_index(j))
        {
          NumericVector<Number> * coarse_adjoint =
            NumericVector<Number>::build(mesh.comm()).release();

          // Can do "fast" init since we're overwriting this in a sec
          coarse_adjoint->init(system.get_adjoint_solution(j),
                               /* fast = */ true);

          *coarse_adjoint = system.get_adjoint_solution(j);

          coarse_adjoints.push_back(coarse_adjoint);
        }
      else
        coarse_adjoints.push_back(static_cast<NumericVector<Number> *>(libmesh_nullptr));
    }

  // Rebuild the rhs with the projected primal solution
  (dynamic_cast<ImplicitSystem &>(system)).assembly(true, false);
  NumericVector<Number> & projected_residual = (dynamic_cast<ExplicitSystem &>(system)).get_vector("RHS Vector");
  projected_residual.close();

  // Solve the adjoint problem(s) on the refined FE space
  system.adjoint_solve(_qoi_set);

  // Now that we have the refined adjoint solution and the projected primal solution,
  // we first compute the global QoI error estimate

  // Resize the computed_global_QoI_errors vector to hold the error estimates for each QoI
  computed_global_QoI_errors.resize(system.qoi.size());

  // Loop over all the adjoint solutions and get the QoI error
  // contributions from all of them.  While we're looping anyway we'll
  // pull off the coarse adjoints
  for (std::size_t j=0; j != system.qoi.size(); j++)
    {
      // Skip this QoI if not in the QoI Set
      if (_qoi_set.has_index(j))
        {
          // If the adjoint solution has heterogeneous dirichlet
          // values, then to get a proper error estimate here we need
          // to subtract off a coarse grid lift function.  In any case
          // we can get a better error estimate by separating off a
          // coarse representation of the adjoint solution, so we'll
          // use that for our lift function.
          system.get_adjoint_solution(j) -= *coarse_adjoints[j];

          computed_global_QoI_errors[j] = projected_residual.dot(system.get_adjoint_solution(j));
        }
    }

  // Done with the global error estimates, now construct the element wise error indicators

  // We ought to account for 'spill-over' effects while computing the
  // element error indicators This happens because the same dof is
  // shared by multiple elements, one way of mitigating this is to
  // scale the contribution from each dof by the number of elements it
  // belongs to We first obtain the number of elements each node
  // belongs to

  // A map that relates a node id to an int that will tell us how many elements it is a node of
  LIBMESH_BEST_UNORDERED_MAP<dof_id_type, unsigned int>shared_element_count;

  // To fill this map, we will loop over elements, and then in each element, we will loop
  // over the nodes each element contains, and then query it for the number of coarse
  // grid elements it was a node of

  // Keep track of which nodes we have already dealt with
  LIBMESH_BEST_UNORDERED_SET<dof_id_type> processed_node_ids;

  // We will be iterating over all the active elements in the fine mesh that live on
  // this processor
  {
    MeshBase::const_element_iterator elem_it = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator elem_end = mesh.active_local_elements_end();

    // Start loop over elems
    for(; elem_it != elem_end; ++elem_it)
      {
        // Pointer to this element
        const Elem * elem = *elem_it;

        // Loop over the nodes in the element
        for(unsigned int n=0; n != elem->n_nodes(); ++n)
          {
            // Get a reference to the current node
            const Node & node = elem->node_ref(n);

            // Get the id of this node
            dof_id_type node_id = node.id();

            // If we havent already processed this node, do so now
            if(processed_node_ids.find(node_id) == processed_node_ids.end())
              {
                // Declare a neighbor_set to be filled by the find_point_neighbors
                std::set<const Elem *> fine_grid_neighbor_set;

                // Call find_point_neighbors to fill the neighbor_set
                elem->find_point_neighbors(node, fine_grid_neighbor_set);

                // A vector to hold the coarse grid parents neighbors
                std::vector<dof_id_type> coarse_grid_neighbors;

                // Iterators over the fine grid neighbors set
                std::set<const Elem *>::iterator fine_neighbor_it = fine_grid_neighbor_set.begin();
                const std::set<const Elem *>::iterator fine_neighbor_end = fine_grid_neighbor_set.end();

                // Loop over all the fine neighbors of this node
                for(; fine_neighbor_it != fine_neighbor_end ; ++fine_neighbor_it)
                  {
                    // Pointer to the current fine neighbor element
                    const Elem * fine_elem = *fine_neighbor_it;

                    // Find the element id for the corresponding coarse grid element
                    const Elem * coarse_elem = fine_elem;
                    for (unsigned int j = 0; j != number_h_refinements; ++j)
                      {
                        libmesh_assert (coarse_elem->parent());

                        coarse_elem = coarse_elem->parent();
                      }

                    // Loop over the existing coarse neighbors and check if this one is
                    // already in there
                    const dof_id_type coarse_id = coarse_elem->id();
                    std::size_t j = 0;
                    for (; j != coarse_grid_neighbors.size(); j++)
                      {
                        // If the set already contains this element break out of the loop
                        if(coarse_grid_neighbors[j] == coarse_id)
                          {
                            break;
                          }
                      }

                    // If we didn't leave the loop even at the last element,
                    // this is a new neighbour, put in the coarse_grid_neighbor_set
                    if(j == coarse_grid_neighbors.size())
                      {
                        coarse_grid_neighbors.push_back(coarse_id);
                      }

                  } // End loop over fine neighbors

                // Set the shared_neighbour index for this node to the
                // size of the coarse grid neighbor set
                shared_element_count[node_id] =
                  cast_int<unsigned int>(coarse_grid_neighbors.size());

                // Add this node to processed_node_ids vector
                processed_node_ids.insert(node_id);

              } // End if not processed node

          } // End loop over nodes

      }  // End loop over elems
  }

  // Get a DoF map, will be used to get the nodal dof_indices for each element
  DofMap & dof_map = system.get_dof_map();

  // The global DOF indices, we will use these later on when we compute the element wise indicators
  std::vector<dof_id_type> dof_indices;

  // Localize the global rhs and adjoint solution vectors (which might be shared on multiple processsors) onto a
  // local ghosted vector, this ensures each processor has all the dof_indices to compute an error indicator for
  // an element it owns
  UniquePtr<NumericVector<Number> > localized_projected_residual = NumericVector<Number>::build(system.comm());
  localized_projected_residual->init(system.n_dofs(), system.n_local_dofs(), system.get_dof_map().get_send_list(), false, GHOSTED);
  projected_residual.localize(*localized_projected_residual, system.get_dof_map().get_send_list());

  // Each adjoint solution will also require ghosting; for efficiency we'll reuse the same memory
  UniquePtr<NumericVector<Number> > localized_adjoint_solution = NumericVector<Number>::build(system.comm());
  localized_adjoint_solution->init(system.n_dofs(), system.n_local_dofs(), system.get_dof_map().get_send_list(), false, GHOSTED);

  // We will loop over each adjoint solution, localize that adjoint
  // solution and then loop over local elements
  for (std::size_t i=0; i != system.qoi.size(); i++)
    {
      // Skip this QoI if not in the QoI Set
      if (_qoi_set.has_index(i))
        {
          // Get the weight for the current QoI
          Real error_weight = _qoi_set.weight(i);

          (system.get_adjoint_solution(i)).localize(*localized_adjoint_solution, system.get_dof_map().get_send_list());

          // Loop over elements
          MeshBase::const_element_iterator elem_it = mesh.active_local_elements_begin();
          const MeshBase::const_element_iterator elem_end = mesh.active_local_elements_end();

          for(; elem_it != elem_end; ++elem_it)
            {
              // Pointer to the element
              const Elem * elem = *elem_it;

              // Go up number_h_refinements levels up to find the coarse parent
              const Elem * coarse = elem;

              for (unsigned int j = 0; j != number_h_refinements; ++j)
                {
                  libmesh_assert (coarse->parent());

                  coarse = coarse->parent();
                }

              const dof_id_type e_id = coarse->id();

              // Get the local to global degree of freedom maps for this element
              dof_map.dof_indices (elem, dof_indices);

              // We will have to manually do the dot products.
              Number local_contribution = 0.;

              for (std::size_t j=0; j != dof_indices.size(); j++)
                {
                  // The contribution to the error indicator for this element from the current QoI
                  local_contribution += (*localized_projected_residual)(dof_indices[j]) * (*localized_adjoint_solution)(dof_indices[j]);
                }

              // Multiply by the error weight for this QoI
              local_contribution *= error_weight;

              // FIXME: we're throwing away information in the
              // --enable-complex case
              error_per_cell[e_id] += static_cast<ErrorVectorReal>
                (std::abs(local_contribution));

            } // End loop over elements

        } // End if belong to QoI set

    } // End loop over QoIs

  for (std::size_t j=0; j != system.qoi.size(); j++)
    {
      if (_qoi_set.has_index(j))
        {
          delete coarse_adjoints[j];
        }
    }

  // Don't bother projecting the solution; we'll restore from backup
  // after coarsening
  system.project_solution_on_reinit() = false;

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

  // Restore old solutions and clean up the heap
  system.project_solution_on_reinit() = old_projection_setting;

  // Restore the coarse solution vectors and delete their copies
  *system.solution = *coarse_solution;
  delete coarse_solution;
  *system.current_local_solution = *coarse_local_solution;
  delete coarse_local_solution;

  for (System::vectors_iterator vec = system.vectors_begin(); vec !=
         system.vectors_end(); ++vec)
    {
      // The (string) name of this vector
      const std::string & var_name = vec->first;

      // If it's a vector we already had (and not a newly created
      // vector like an adjoint rhs), we need to restore it.
      std::map<std::string, NumericVector<Number> *>::iterator it =
        coarse_vectors.find(var_name);
      if (it != coarse_vectors.end())
        {
          NumericVector<Number> * coarsevec = it->second;
          system.get_vector(var_name) = *coarsevec;

          coarsevec->clear();
          delete coarsevec;
        }
    }

  // Restore old partitioner and renumbering settings
  mesh.partitioner().reset(old_partitioner.release());
  mesh.allow_renumbering(old_renumbering_setting);

  // Fiinally sum the vector of estimated error values.
  this->reduce_error(error_per_cell, system.comm());

  // We don't take a square root here; this is a goal-oriented
  // estimate not a Hilbert norm estimate.
} // end estimate_error function

}// namespace libMesh

#endif // #ifdef LIBMESH_ENABLE_AMR
