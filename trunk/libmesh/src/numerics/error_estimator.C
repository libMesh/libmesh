// $Id: error_estimator.C,v 1.4 2003-05-19 21:21:13 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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


// libMesh Configuration
#include "mesh_config.h"

// C++ includes
#include <algorithm>
#include <math.h>

#ifdef HAVE_MPI
namespace Mpi
{
  extern "C"{
#include <mpi.h>
  }
}
#endif


// Local Includes
#include "error_estimator.h"
#include "equation_systems.h"
#include "steady_system.h"
#include "dof_map.h"
#include "mesh.h"
#include "fe.h"
#include "fe_interface.h"
#include "quadrature_gauss.h"
#include "mesh_logging.h"




//-----------------------------------------------------------------
// ErrorEstimator implementations
void ErrorEstimator::flux_jump (const EquationSystems& es,
				const std::string& name,
				std::vector<float>& error_per_cell)
{
  START_LOG("flux_jump()", "ErrorEstimator");
  
  /*

  - e & f are global element ids
  
  Case (1.) Elements are at the same level, e<f
            Compute the flux jump on the face and
	    add it as a contribution to error_per_cell[e]
	    and error_per_cell[f]
  
                   ----------------------
		  |           |          |
		  |           |    f     |
		  |           |          |
		  |    e      |---> n    | 
		  |           |          |
		  |           |          |
                   ----------------------


   Case (2.) The neighbor is at a higher level.
             Compute the flux jump on e's face and
	     add it as a contribution to error_per_cell[e]
	     and error_per_cell[f]

                   ----------------------
		  |     |     |          |
		  |     |  e  |---> n    |
		  |     |     |          |
		  |-----------|    f     | 
		  |     |     |          |
		  |     |     |          |
		  |     |     |          |
                   ----------------------
*/
   
  // The current mesh
  const Mesh& mesh = es.get_mesh();

  // The dimensionality of the mesh
  const unsigned int dim = mesh.mesh_dimension();
  
  // The System object to estimate the error for
  const SteadySystem& system = es.get_system<SteadySystem>(name);

  // The DofMap for this system
  const DofMap& dof_map = system.get_dof_map();

  // Resize the error_per_cell vector to be
  // the number of elements, initialize it to 0.
  error_per_cell.resize (mesh.n_elem());
  std::fill (error_per_cell.begin(), error_per_cell.end(), 0.);

  

  
  // Loop over all the variables in the system
  for (unsigned int var=0; var<system.n_vars(); var++)
    {
      // The type of finite element to use for this variable
      const FEType& fe_type = dof_map.variable_type (var);

      // Finite element objects for the same face from
      // different sides
      AutoPtr<FEBase> fe_e (FEBase::build (dim, fe_type));
      AutoPtr<FEBase> fe_f (FEBase::build (dim, fe_type));

      // Build an appropriate Gaussian quadrature rule
      QGauss qrule (dim-1, fe_type.default_quadrature_order());

      // Tell the finite element for element e about the quadrature
      // rule.  The finite element for element f need not know about it
      fe_e->attach_quadrature_rule (&qrule);
      
      // By convention we will always do the integration
      // on the face of element e.  Get its Jacobian values, etc..
      const std::vector<Real>&  JxW_face     = fe_e->get_JxW();
      const std::vector<Point>& qface_point  = fe_e->get_xyz();
      const std::vector<Point>& face_normals = fe_e->get_normals();

      // The quadrature points on element f.  These will be computed
      // from the quadrature points on element e.
      std::vector<Point> qp_f;
      
      // The shape function gradients on elements e & f
      const std::vector<std::vector<Point> > & dphi_e = fe_e->get_dphi();
      const std::vector<std::vector<Point> > & dphi_f = fe_f->get_dphi();
      
      // The global DOF indices for elements e & f
      std::vector<unsigned int> dof_indices_e;
      std::vector<unsigned int> dof_indices_f;


      
      // Iterate over all the active elements in the mesh
      // that live on this processor.
      const_active_local_elem_iterator       elem_it (mesh.elements_begin());
      const const_active_local_elem_iterator elem_end(mesh.elements_end());

      for (; elem_it != elem_end; ++elem_it)
	{
	  // e is necessarily an active element on the local processor
	  const Elem* e = *elem_it;
	  const unsigned int e_id = e->id();
	  
	  // Loop over the neighbors of element e
	  for (unsigned int s_e=0; s_e<e->n_neighbors(); s_e++)
	    if (e->neighbor(s_e) != NULL) // e is not on the boundary
	      {
		const Elem* f           = e->neighbor(s_e);
		const unsigned int f_id = f->id();
		
		if (   //-------------------------------------
		    ((f->active()) &&
		     (f->level() == e->level()) &&
		     (e_id < f_id))                 // Case 1.
		    
		    || //-------------------------------------
		    
		    (f->level() < e->level())       // Case 2.
		    
		    )  //-------------------------------------
		  {		    
		    // Update the shape functions on side s_e of
		    // element e
		    fe_e->reinit (e, s_e);

		    // Build the side
		    AutoPtr<Elem> side (e->side(s_e));

		    // Get the maximum h for this side
		    const Real h = side->hmax();
		    
		    // Get the DOF indices for the two elements
		    dof_map.dof_indices (e, dof_indices_e, var);
		    dof_map.dof_indices (f, dof_indices_f, var);

		    // The number of DOFS on each element
		    const unsigned int n_dofs_e = dof_indices_e.size();
		    const unsigned int n_dofs_f = dof_indices_f.size();

		    // The number of quadrature points
		    const unsigned int n_qp = qrule.n_points();

		    // Find the location of the quadrature points
		    // on element f
		    FEInterface::inverse_map (dim, fe_type, f, qface_point, qp_f);

		    // Compute the shape functions on element f
		    // at the quadrature points of element e
		    fe_f->reinit (f, &qp_f);

		    // The error contribution from this face
		    Real error = 0.;

		    
		    
		    // loop over the integration points on the face
		    for (unsigned int qp=0; qp<n_qp; qp++)
		      {
			// The solution gradient from each element
			Point grad_e, grad_f;
			
			// Compute the solution gradient on element e
			for (unsigned int i=0; i<n_dofs_e; i++)
			  grad_e.add_scaled (dphi_e[i][qp],
					     system.current_solution(dof_indices_e[i]));
			
			// Compute the solution gradient on element f
			for (unsigned int i=0; i<n_dofs_f; i++)
			  grad_f.add_scaled (dphi_f[i][qp],
					     system.current_solution(dof_indices_f[i]));


			// The flux jump at the face 
			const Real jump = (grad_e - grad_f)*face_normals[qp];

			// Integrate the error on the face
			error += JxW_face[qp]*h*jump*jump;			
			
		      } // End quadrature point loop

		    // Add the error contribution to elements e & f
		    error_per_cell[e_id] += error;
		    error_per_cell[f_id] += error;
		    
		    
		  }
	      } // if (e->neigbor(s_e) != NULL)
	  
	} // End loop over active local elements
      
    } // End loop over variables


  // Each processor has now computed the error contribuions
  // for its local elements.  We need to sum the vector
  // and then take the square-root of each component.  Note
  // that we only need to sum if we are running on multiple
  // processors, and we only need to take the square-root
  // if the value is nonzero.  There will in general be many
  // zeros for the inactive elements.
#ifdef HAVE_MPI
  if (libMesh::n_processors() > 1)
    {
      using namespace Mpi;
      
      MPI_Allreduce (&error_per_cell[0], &error_per_cell[0],
		     error_per_cell.size(),
		     MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    }  
#endif

  // Compute the square-root of each component.
  for (unsigned int i=0; i<error_per_cell.size(); i++)
    if (error_per_cell[i] != 0.)
      error_per_cell[i] = sqrt(error_per_cell[i]);
  
  STOP_LOG("flux_jump()", "ErrorEstimator");
}
