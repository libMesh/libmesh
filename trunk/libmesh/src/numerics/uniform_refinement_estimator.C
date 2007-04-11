// $Id: uniform_refinement_estimator.C,v 1.11 2007-04-11 23:06:42 roystgnr Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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
#include <cmath>    // for sqrt


// Local Includes
#include "dof_map.h"
#include "elem.h"
#include "equation_systems.h"
#include "error_vector.h"
#include "fe.h"
#include "fe_interface.h"
#include "libmesh_common.h"
#include "libmesh_logging.h"
#include "mesh.h"
#include "mesh_refinement.h"
#include "numeric_vector.h"
#include "quadrature.h"
#include "system.h"
#include "uniform_refinement_estimator.h"

//-----------------------------------------------------------------
// ErrorEstimator implementations
void UniformRefinementEstimator::estimate_error (const System& _system,
					         ErrorVector& error_per_cell,
					         bool estimate_parent_error)
{
  START_LOG("estimate_error()", "UniformRefinementEstimator");
  this->_estimate_error
    (NULL, &_system, &error_per_cell, NULL, NULL, estimate_parent_error);
  STOP_LOG("estimate_error()", "UniformRefinementEstimator");
}

void UniformRefinementEstimator::estimate_errors (const EquationSystems& _es,
					         ErrorVector& error_per_cell,
                                                 std::map<const System*, std::vector<float> >& component_scales,
					         bool estimate_parent_error)
{
  START_LOG("estimate_errors()", "UniformRefinementEstimator");
  this->_estimate_error
    (&_es, NULL, &error_per_cell, NULL, &component_scales, estimate_parent_error);
  STOP_LOG("estimate_errors()", "UniformRefinementEstimator");
}

void UniformRefinementEstimator::estimate_errors (const EquationSystems& _es,
                                                  ErrorMap& errors_per_cell,
                                                  bool estimate_parent_error)
{
  START_LOG("estimate_errors()", "UniformRefinementEstimator");
  this->_estimate_error
    (&_es, NULL, NULL, &errors_per_cell, NULL, estimate_parent_error);
  STOP_LOG("estimate_errors()", "UniformRefinementEstimator");
}

void UniformRefinementEstimator::_estimate_error (const EquationSystems* _es,
                                                 const System* _system,
					         ErrorVector* error_per_cell,
                                                 ErrorMap* errors_per_cell,
                                                 std::map<const System*, std::vector<float> > *_component_scales,
					         bool)
{
  // Get a vector of the Systems we're going to work on,
  // and set up a component_scales map if necessary
  std::vector<System *> system_list;
  std::map<const System*, std::vector<float> > *component_scales;

  if (_es)
    {
      assert(!_system);
      assert(_es->n_systems());
      _system = &(_es->get_system(0));
      assert(&(_system->get_equation_systems()) == _es);

      assert(_es->n_systems());
      for (unsigned int i=0; i != _es->n_systems(); ++i)
      // We have to break the rules here, because we can't refine a const System
        system_list.push_back(const_cast<System *>(&(_es->get_system(i))));

      // If we're computing one vector, we need to know how to scale
      // each variable's contributions to it.
      if (_component_scales)
        {
          assert(!errors_per_cell);
          component_scales = _component_scales;
        }
      else
      // If we're computing many vectors, we just need to know which
      // variables to skip
        {
          assert (errors_per_cell);

          component_scales = new std::map<const System*, std::vector<float> >;
          for (unsigned int i=0; i!= _es->n_systems(); ++i)
            {
              const System &sys = _es->get_system(i);
              unsigned int n_vars = sys.n_vars();

              std::vector<float> cs(sys.n_vars(), 0.0);
              for (unsigned int v = 0; v != n_vars; ++v)
                {
                  if (errors_per_cell->find(std::make_pair(&sys, v)) ==
                      errors_per_cell->end())
                    continue;

                  cs[v] = 1.0;
                }
              (*component_scales)[_system] = cs;
            }
        }
    }
  else
    {
      assert(_system);
      // We have to break the rules here, because we can't refine a const System
      system_list.push_back(const_cast<System *>(_system));

      assert(!_component_scales);
      component_scales = new std::map<const System*, std::vector<float> >;
      (*component_scales)[_system] = component_scale;
    }

  // An EquationSystems reference will be convenient.
  // We have to break the rules here, because we can't refine a const System
  EquationSystems& es =
    const_cast<EquationSystems &>(_system->get_equation_systems());

  // The current mesh
  Mesh& mesh = es.get_mesh();

  // The dimensionality of the mesh
  const unsigned int dim = mesh.mesh_dimension();
  
  // Resize the error_per_cell vectors to be
  // the number of elements, initialize them to 0.
  if (error_per_cell)
    {
      error_per_cell->clear();
      error_per_cell->resize (mesh.n_elem(), 0.);
    }
  else
    {
      assert(errors_per_cell);
      for (ErrorMap::iterator i = errors_per_cell->begin();
           i != errors_per_cell->end(); ++i)
        {
          ErrorVector *e = i->second;
          e->clear();
          e->resize(mesh.n_elem(), 0.);
        }
    }

  // component_mask has long since been deprecated
  if (!component_mask.empty())
    deprecated();

  // We'll want to back up all coarse grid vectors
  std::vector<std::map<std::string, NumericVector<Number> *> >
    coarse_vectors(system_list.size());
  std::vector<NumericVector<Number> *>
    coarse_solutions(system_list.size());
  std::vector<NumericVector<Number> *>
    coarse_local_solutions(system_list.size());
  // And make copies of projected solutions
  std::vector<NumericVector<Number> *>
    projected_solutions(system_list.size());

  // And we'll need to temporarily change solution projection settings
  std::vector<bool> old_projection_settings(system_list.size());

  for (unsigned int i=0; i != system_list.size(); ++i)
    {
      System &system = *system_list[i];


      // Check for valid component_scales
      if (!(*component_scales)[&system].empty())
        {
          if ((*component_scales)[&system].size() != system.n_vars())
	    {
	      std::cerr << "ERROR: component_scale is the wrong size:"
		        << std::endl
		        << " component_scales[" << i << "].scale()=" 
                        << (*component_scales)[&system].size()
		        << std::endl
		        << " n_vars=" << system.n_vars()
		        << std::endl;
	      error();
	    }
        }
      else
        {
          // No specified scaling.  Scale all variables by one.
          (*component_scales)[&system].resize (system.n_vars());
          std::fill ((*component_scales)[&system].begin(),
                     (*component_scales)[&system].end(), 1.0);
        }
  
      // Back up the solution vector
      coarse_solutions[i] = system.solution->clone().release();
      coarse_local_solutions[i] =
        system.current_local_solution->clone().release();

      // Back up all other coarse grid vectors
      for (System::vectors_iterator vec = system.vectors_begin(); vec !=
           system.vectors_end(); ++vec)
        {
          // The (string) name of this vector
          const std::string& var_name = vec->first;

          coarse_vectors[i][var_name] = vec->second->clone().release();
        }

      // Make sure the solution is projected when we refine the mesh
      old_projection_settings[i] = system.project_solution_on_reinit();
      system.project_solution_on_reinit() = true;
    }

  // Find the number of coarse mesh elements, to make it possible
  // to find correct coarse elem ids later
  const unsigned int n_coarse_elem = mesh.n_elem();

  // Uniformly refine the mesh
  MeshRefinement mesh_refinement(mesh);

  assert (number_h_refinements > 0 || number_p_refinements > 0);

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

  for (unsigned int i=0; i != system_list.size(); ++i)
    {
      System &system = *system_list[i];
      
      // Copy the projected coarse grid solutions, which will be
      // overwritten by solve()
//      projected_solutions[i] = system.solution->clone().release();
      projected_solutions[i] = NumericVector<Number>::build().release();
      projected_solutions[i]->init(system.solution->size(), system.solution->size());
      system.solution->localize(*projected_solutions[i],
                                system.get_dof_map().get_send_list());
    }

  // Get the uniformly refined solution.

  if (_es)
    es.solve();
  else
    system_list[0]->solve();
  
  // Get the error in the uniformly refined solution(s).

  for (unsigned int i=0; i != system_list.size(); ++i)
    {
      System &system = *system_list[i];

      unsigned int n_vars = system.n_vars();

      DofMap &dof_map = system.get_dof_map();

      std::vector<float> &_component_scale =
        (*component_scales)[&system];

      NumericVector<Number> *projected_solution = projected_solutions[i];

      // Loop over all the variables in the system
      for (unsigned int var=0; var<n_vars; var++)
        {
          // Possibly skip this variable
          if (!_component_scale.empty())
	    if (_component_scale[var] == 0.0) continue;

          // The type of finite element to use for this variable
          const FEType& fe_type = dof_map.variable_type (var);
      
          // Finite element object for each fine element
          AutoPtr<FEBase> fe (FEBase::build (dim, fe_type));

          // Build and attach an appropriate quadrature rule
          AutoPtr<QBase> qrule = fe_type.default_quadrature_rule(dim);
          fe->attach_quadrature_rule (qrule.get());
      
          const std::vector<Real>&  JxW = fe->get_JxW();
          const std::vector<std::vector<Real> >& phi = fe->get_phi();
          const std::vector<std::vector<RealGradient> >& dphi =
            fe->get_dphi();
#ifdef ENABLE_SECOND_DERIVATIVES
          const std::vector<std::vector<RealTensor> >& d2phi =
            fe->get_d2phi();
#endif

          // The global DOF indices for the fine element
          std::vector<unsigned int> dof_indices;
      
          // Iterate over all the active elements in the fine mesh
          // that live on this processor.
          MeshBase::const_element_iterator       elem_it  = mesh.active_local_elements_begin();
          const MeshBase::const_element_iterator elem_end = mesh.active_local_elements_end(); 

          for (; elem_it != elem_end; ++elem_it)
	    {
	      // e is necessarily an active element on the local processor
	      const Elem* elem = *elem_it;

              // Find the element id for the corresponding coarse grid element
              const Elem* coarse = elem;
              unsigned int e_id = coarse->id();
              while (e_id >= n_coarse_elem)
                {
                  assert (coarse->parent());
                  coarse = coarse->parent();
                  e_id = coarse->id();
                }
          
              double L2normsq = 0., H1seminormsq = 0., H2seminormsq = 0.;

              // reinitialize the element-specific data
              // for the current element
              fe->reinit (elem);

              // Get the local to global degree of freedom maps
              dof_map.dof_indices (elem, dof_indices, var);

              // The number of quadrature points
              const unsigned int n_qp = qrule->n_points();

              // The number of shape functions
              const unsigned int n_sf = dof_indices.size();

              //
              // Begin the loop over the Quadrature points.
              //
              for (unsigned int qp=0; qp<n_qp; qp++)
                {
                  Number u_fine = 0., u_coarse = 0.;

#ifndef USE_COMPLEX_NUMBERS
                  RealGradient grad_u_fine, grad_u_coarse;
#else
                  // Gradient     grad_u_fine;
                  RealGradient grad_u_fine_re, grad_u_coarse_re;
                  RealGradient grad_u_fine_im, grad_u_coarse_im;
#endif
#ifdef ENABLE_SECOND_DERIVATIVES
  #ifndef USE_COMPLEX_NUMBERS
                  RealTensor grad2_u_fine, grad2_u_coarse;
  #else
                  RealTensor grad2_u_fine_re, grad2_u_coarse_re;
                  RealTensor grad2_u_fine_im, grad2_u_coarse_im;
  #endif
#endif

                  // Compute solution values at the current
                  // quadrature point.  This reqiures a sum
                  // over all the shape functions evaluated
                  // at the quadrature point.
                  for (unsigned int i=0; i<n_sf; i++)
                    {
                      u_fine            += phi[i][qp]*system.current_solution (dof_indices[i]);
                      u_coarse          += phi[i][qp]*(*projected_solution) (dof_indices[i]);
#ifndef USE_COMPLEX_NUMBERS
                      grad_u_fine       += dphi[i][qp]*system.current_solution (dof_indices[i]);
                      grad_u_coarse     += dphi[i][qp]*(*projected_solution) (dof_indices[i]);
#else
                      grad_u_fine_re    += dphi[i][qp]*system.current_solution (dof_indices[i]).real();
                      grad_u_fine_im    += dphi[i][qp]*system.current_solution (dof_indices[i]).imag();
                      grad_u_coarse_re  += dphi[i][qp]*(*projected_solution) (dof_indices[i]).real();
                      grad_u_coarse_im  += dphi[i][qp]*(*projected_solution) (dof_indices[i]).imag();
#endif
#ifdef ENABLE_SECOND_DERIVATIVES
  #ifndef USE_COMPLEX_NUMBERS
                      grad2_u_fine      += d2phi[i][qp]*system.current_solution (dof_indices[i]);
                      grad2_u_coarse    += d2phi[i][qp]*(*projected_solution) (dof_indices[i]);
  #else
                      grad2_u_fine_re   += d2phi[i][qp]*system.current_solution (dof_indices[i]).real();
                      grad2_u_fine_im   += d2phi[i][qp]*system.current_solution (dof_indices[i]).imag();
                      grad2_u_coarse_re += d2phi[i][qp]*(*projected_solution) (dof_indices[i]).real();
                      grad2_u_coarse_im += d2phi[i][qp]*(*projected_solution) (dof_indices[i]).imag();
  #endif
#endif
                    }

#ifdef USE_COMPLEX_NUMBERS
                  Gradient grad_u_fine (grad_u_fine_re, grad_u_fine_im);
                  Gradient grad_u_coarse (grad_u_coarse_re, grad_u_coarse_im);
  #ifdef ENABLE_SECOND_DERIVATIVES
                  Tensor grad2_u_fine (grad2_u_fine_re, grad2_u_fine_im);
                  Tensor grad2_u_coarse (grad2_u_coarse_re, grad2_u_coarse_im);
  #endif
#endif

                  // Compute the value of the error at this quadrature point
                  const Number val_error = u_fine - u_coarse;

                  // Add the squares of the error to each contribution
#ifndef USE_COMPLEX_NUMBERS
                  L2normsq += JxW[qp]*(val_error*val_error);
#else
                  L2normsq += JxW[qp]*std::norm(val_error);
#endif

                  // Compute the value of the error in the gradient at this
                  // quadrature point
                  if (_sobolev_order > 0)
                    {
                      Gradient grad_error = grad_u_fine - grad_u_coarse;

#ifndef USE_COMPLEX_NUMBERS
                      H1seminormsq += JxW[qp]*(grad_error*grad_error);
#else
                      H1seminormsq += JxW[qp]*std::abs(grad_error*grad_error);
#endif
                    }


#ifdef ENABLE_SECOND_DERIVATIVES
                  // Compute the value of the error in the hessian at this
                  // quadrature point
                  if (_sobolev_order > 1)
                    {
                      Tensor grad2_error = grad2_u_fine - grad2_u_coarse;

                      H2seminormsq += JxW[qp]*std::abs(grad2_error.contract(grad2_error));
                    }
#endif

                } // end qp loop

              assert (L2normsq     >= 0.);
              assert (H1seminormsq >= 0.);

              if (error_per_cell)
                {
                  (*error_per_cell)[e_id] += L2normsq;
                  if (_sobolev_order > 0)
                    (*error_per_cell)[e_id] += H1seminormsq;
                  if (_sobolev_order > 1)
                    (*error_per_cell)[e_id] += H2seminormsq;
                }
              else
                {
                  assert(errors_per_cell);
                  ErrorVector &e =
                    *((*errors_per_cell)[std::make_pair(&system,var)]);

                  e[e_id] += L2normsq;
                  if (_sobolev_order > 0)
                    e[e_id] += H1seminormsq;
                  if (_sobolev_order > 1)
                    e[e_id] += H2seminormsq;
                }
            } // End loop over active local elements
        } // End loop over variables

      // Don't bother projecting the solution; we'll restore from backup
      // after coarsening
      system.project_solution_on_reinit() = false;
    }


  // Uniformly coarsen the mesh, without projecting the solution
  assert (number_h_refinements > 0 || number_p_refinements > 0);

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
  assert(n_coarse_elem == mesh.n_elem());

  // Each processor has now computed the error contribuions
  // for its local elements.  We need to sum the vector
  // and then take the square-root of each component.  Note
  // that we only need to sum if we are running on multiple
  // processors, and we only need to take the square-root
  // if the value is nonzero.  There will in general be many
  // zeros for the inactive elements.

  // First sum the vector of estimated error values
  this->reduce_error(*error_per_cell);

  // Compute the square-root of each component.
  START_LOG("std::sqrt()", "UniformRefinementEstimator");
  for (unsigned int i=0; i<error_per_cell->size(); i++)
    if ((*error_per_cell)[i] != 0.)
      (*error_per_cell)[i] = std::sqrt((*error_per_cell)[i]);
  STOP_LOG("std::sqrt()", "UniformRefinementEstimator");

  // Restore old solutions and clean up the heap
  for (unsigned int i=0; i != system_list.size(); ++i)
    {
      System &system = *system_list[i];

      system.project_solution_on_reinit() = old_projection_settings[i];
  
      // Restore the coarse solution vectors and delete their copies
      *system.solution = *coarse_solutions[i];
      delete coarse_solutions[i];
      *system.current_local_solution = *coarse_local_solutions[i];
      delete coarse_local_solutions[i];
      delete projected_solutions[i];

      for (System::vectors_iterator vec = system.vectors_begin(); vec !=
           system.vectors_end(); ++vec)
        {
          // The (string) name of this vector
          const std::string& var_name = vec->first;

          system.get_vector(var_name) = *coarse_vectors[i][var_name];

          coarse_vectors[i][var_name]->clear();
          delete coarse_vectors[i][var_name];
        }
    }

  if (!_component_scales)
    delete component_scales;
}
