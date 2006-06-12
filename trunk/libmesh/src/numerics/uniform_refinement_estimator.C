// $Id: uniform_refinement_estimator.C,v 1.2 2006-06-12 21:51:14 roystgnr Exp $

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
					  std::vector<float>& error_per_cell)
{
  START_LOG("estimate_error()", "UniformRefinementEstimator");

  // We have to break the rules here, because we can't refine a const System
  System& system = const_cast<System &>(_system);

  // The current mesh
  Mesh& mesh = system.get_mesh();

  // The dimensionality of the mesh
  const unsigned int dim = mesh.mesh_dimension();
  
  // The number of variables in the system
  const unsigned int n_vars = system.n_vars();
  
  // The DofMap for this system
  const DofMap& dof_map = system.get_dof_map();

  // Resize the error_per_cell vector to be
  // the number of elements, initialize it to 0.
  error_per_cell.resize (mesh.n_elem());
  std::fill (error_per_cell.begin(), error_per_cell.end(), 0.);

  // Check for a valid component_scale
  if (!component_scale.empty())
    {
      if (component_scale.size() != n_vars)
	{
	  std::cerr << "ERROR: component_scale is the wrong size:"
		    << std::endl
		    << " component_scale.scale()=" << component_scale.size()
		    << std::endl
		    << " n_vars=" << n_vars
		    << std::endl;
	  error();
	}
    }
  else
    {
      // No specified scaling.  Scale all variables by one.
      component_scale.resize (n_vars);
      std::fill (component_scale.begin(), component_scale.end(), 1.0);
    }
  
  NumericVector<Number> * debugging;

  // Back up the coarse grid vectors
  std::map<std::string, NumericVector<Number> *> coarse_vectors;
  for (System::vectors_iterator vec = system.vectors_begin(); vec !=
       system.vectors_end(); ++vec)
    {
      // The (string) name of this vector
      const std::string& var_name = vec->first;

      coarse_vectors[var_name] = vec->second->clone().release();
      debugging = coarse_vectors[var_name];
    }

  // Find the number of coarse mesh elements, to make it possible
  // to find correct coarse elem ids later
  const unsigned int n_coarse_elem = mesh.n_elem();

  // Uniformly refine the mesh, making sure the solution is projected
  bool old_projection_setting = system.project_solution_on_reinit();
  system.project_solution_on_reinit() = true;

  MeshRefinement mesh_refinement(mesh);

  assert (number_h_refinements > 0 || number_p_refinements > 0);

  // FIXME: this will probably break if there is more than one System
  // on this mesh
  for (unsigned int i = 0; i != number_h_refinements; ++i)
    {
      mesh_refinement.uniformly_refine(1);
      system.get_equation_systems().reinit();
    }
      
  for (unsigned int i = 0; i != number_p_refinements; ++i)
    {
      mesh_refinement.uniformly_p_refine(1);
      system.get_equation_systems().reinit();
    }

  system.project_solution_on_reinit() = old_projection_setting;
      
  // Copy the projected coarse grid solution, which will be
  // overwritten by solve()
  AutoPtr<NumericVector<Number> > projected_solution = system.solution->clone();

  // Get the uniformly refined solution.

  system.solve();
  
  // Get the error in the uniformly refined solution.

  // Loop over all the variables in the system
  for (unsigned int var=0; var<n_vars; var++)
    {
      // Possibly skip this variable
      if (!component_scale.empty())
	if (component_scale[var] == 0.0) continue;

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
  #ifdef ENABLE_SECOND_DERIVATIVES
              Tensor grad2_u_fine (grad2_u_fine_re, grad2_u_fine_im);
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

                  H2seminormsq += JxW[qp]*(grad2_error.contract(grad2_error));
                }
#endif

            } // end qp loop

          assert (L2normsq     >= 0.);
          assert (H1seminormsq >= 0.);

          error_per_cell[e_id] = L2normsq;
          if (_sobolev_order > 0)
            error_per_cell[e_id] += H1seminormsq;
          if (_sobolev_order > 1)
            error_per_cell[e_id] += H2seminormsq;

        } // End loop over active local elements
    } // End loop over variables


  // Uniformly coarsen the mesh, without projecting the solution
  old_projection_setting = system.project_solution_on_reinit();
  system.project_solution_on_reinit() = false;

  assert (number_h_refinements > 0 || number_p_refinements > 0);

  for (unsigned int i = 0; i != number_h_refinements; ++i)
    {
      mesh_refinement.uniformly_coarsen(1);
      // FIXME - should the reinits here be necessary? - RHS
      system.get_equation_systems().reinit();
    }
      
  for (unsigned int i = 0; i != number_p_refinements; ++i)
    {
      mesh_refinement.uniformly_p_coarsen(1);
      system.get_equation_systems().reinit();
    }

  system.project_solution_on_reinit() = old_projection_setting;

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
  this->reduce_error(error_per_cell);

  // Compute the square-root of each component.
  START_LOG("std::sqrt()", "UniformRefinementEstimator");
  for (unsigned int i=0; i<error_per_cell.size(); i++)
    if (error_per_cell[i] != 0.)
      error_per_cell[i] = std::sqrt(error_per_cell[i]);
  STOP_LOG("std::sqrt()", "UniformRefinementEstimator");

  // Restore the coarse solution vectors and delete their copies
  for (System::vectors_iterator vec = system.vectors_begin(); vec !=
       system.vectors_end(); ++vec)
    {
      // The (string) name of this vector
      const std::string& var_name = vec->first;

      system.get_vector(var_name) = *coarse_vectors[var_name];

      coarse_vectors[var_name]->clear();
      delete coarse_vectors[var_name];
    }

  STOP_LOG("estimate_error()", "UniformRefinementEstimator");
}
