// $Id: discontinuity_measure.C,v 1.4 2006-09-19 17:50:52 roystgnr Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2006  Benjamin S. Kirk, John W. Peterson
  
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
#include "libmesh_common.h"
#include "discontinuity_measure.h"
#include "dof_map.h"
#include "error_vector.h"
#include "fe.h"
#include "fe_interface.h"
#include "quadrature_gauss.h"
#include "libmesh_logging.h"
#include "elem.h"
#include "mesh.h"
#include "system.h"

//-----------------------------------------------------------------
// ErrorEstimator implementations
void DiscontinuityMeasure::estimate_error (const System& system,
					   ErrorVector& error_per_cell,
					   bool estimate_parent_error)
{
  //  START_LOG("jumps()", "DiscontinuityMeasure");
  
  /*

  Case (1.) Elements are at the same level, e<f
            Compute the function jump on the face and
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
             Compute the function jump on e's face and
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
  const Mesh& mesh = system.get_mesh();

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

  // Declare a vector of floats which is as long as
  // error_per_cell above, and fill with zeros.  This vector will be
  // used to keep track of the number of sides on each active
  // element which are either:
  // 1) an internal side
  // 2) an side on a boundary for which a boundary condition
  //    function has been specified.
  // The error estimator can be scaled by the number of internal sides
  // which the element actually has to obtain a more uniform measure
  // of the error.  Use floats instead of ints since in case 2 (above)
  // f gets 1/2 of a flux face contribution from each of his
  // neighbors
  std::vector<float> n_internal_sides (error_per_cell.size());

  // Check for the use of component_mask
  this->convert_component_mask_to_scale();
  
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
  


  // Implement 1D estimator separately
  if(dim == 1)
  {
    // Declare temporary vectors which represent the locations
    // (on the reference element) of the left and right nodes.
    // These are used to trick the FE reinit() functions into
    // recomputing finite element data at the nodes.
    std::vector<Point> left_edge (1);
    std::vector<Point> right_edge(1);
    left_edge[0]  = Point(-1., 0., 0.);
    right_edge[0] = Point( 1., 0., 0.);
    
    for (unsigned int var=0; var<n_vars; var++)
    {
      // Possibly skip this variable
      if (!component_scale.empty())
        if (component_scale[var] == 0.0) continue;

      // The type of finite element to use for this variable
      const FEType& fe_type = dof_map.variable_type (var);

      // Finite element objects for the same face from
      // different sides
      AutoPtr<FEBase> fe_e (FEBase::build (dim, fe_type));
      AutoPtr<FEBase> fe_f (FEBase::build (dim, fe_type));

      const std::vector<std::vector<Real> > & phi_e =
        fe_e->get_phi();
      const std::vector<std::vector<Real> > & phi_f = 
        fe_f->get_phi();
      std::vector<unsigned int> dof_indices_e;
      std::vector<unsigned int> dof_indices_f;


      // Iterate over all the active elements in the mesh
      // that live on this processor.
      MeshBase::const_element_iterator       elem_it  = 
        mesh.active_local_elements_begin();
      const MeshBase::const_element_iterator elem_end = 
        mesh.active_local_elements_end(); 

      for ( ; elem_it != elem_end; ++elem_it)
      {
        // e is necessarily an active element on the local processor
        const Elem* e = *elem_it;
        const unsigned int e_id = e->id();

        for (unsigned int n_e=0; n_e<e->n_neighbors(); n_e++)
        {
          if (e->neighbor(n_e) != NULL) // e is not on the boundary
          {

            const Elem* f = e->neighbor(n_e);
            const unsigned int f_id = f->id();
	    bool case1 = false, case2 = false;

	    // Case 1.	    
            if ( f->active() && (f->level() == e->level()) && (e_id < f_id))
	      case1 = true;

	    // Case 2.
            else if (f->level() < e->level())
	      case2 = true;

	    // Compute jumps if we are in case 1 or case 2.
	    if (case1 || case2)
	      {
		Real error = 0.;

		const Real h_e = e->hmax();
		const Real h_f = f->hmax();

		dof_map.dof_indices (e, dof_indices_e, var);
		dof_map.dof_indices (f, dof_indices_f, var);

		// The number of DOFS on each element
		const unsigned int n_dofs_e = dof_indices_e.size();
		const unsigned int n_dofs_f = dof_indices_f.size();

		if (n_e == 0) // e is not the left edge of mesh
		  {
		    fe_e->reinit(e, &left_edge);
		    fe_f->reinit(f, &right_edge);
		  }
	      
		else if (n_e == 1) // e is not the right edge of the mesh
		  {
		    fe_e->reinit(e, &right_edge);
		    fe_f->reinit(f, &left_edge);
		  }

		else
		  {
		    std::cerr << "A 1D element cannot have more than 2 neighbors!"
			      << std::endl;
		    error();
		  }

		// The solution gradient from each element
		Number val_e = 0., val_f = 0.;

		// Compute the solution gradient at left edge of element e
		for (unsigned int i=0; i<n_dofs_e; i++)
		  val_e += phi_e[i][0] *
			   system.current_solution(dof_indices_e[i]);

		// Compute the solution gradient at right edge of element f
		for (unsigned int i=0; i<n_dofs_f; i++)
		  val_f += phi_f[i][0] *
			   system.current_solution(dof_indices_f[i]);

		Number jump = val_e - val_f;
#ifndef USE_COMPLEX_NUMBERS
		error += jump*jump;
#else
		error += std::norm(jump);
#endif
		// Add the error contribution to element e
		assert(e_id < error_per_cell.size());
		assert(f_id < error_per_cell.size());
		error_per_cell[e_id] += h_e*error*component_scale[var];
		error_per_cell[f_id] += h_f*error*component_scale[var];
	      } // end if (case1 || case2)
          } // end e->neighbor(n_e) != NULL
        } // end loop over neighbors
      } // end loop over active elements
    } // end loop over variables in the system

    // First sum the vector
    this->reduce_error(error_per_cell);

    // Compute the square-root of each component.
    START_LOG("std::sqrt()", "DiscontinuityMeasure");
    for (unsigned int i=0; i<error_per_cell.size(); i++)
      if (error_per_cell[i] != 0.)
        error_per_cell[i] = std::sqrt(error_per_cell[i]);
    STOP_LOG("std::sqrt()", "DiscontinuityMeasure");

    return;
  } // end dim == 1

      
  //--------- 2D and 3D code ---------//
      
  
  // Loop over all the variables in the system
  for (unsigned int var=0; var<n_vars; var++)
    {
      // Possibly skip this variable
      if (!component_scale.empty())
	if (component_scale[var] == 0.0) continue;

      // The (string) name of this variable
      const std::string& var_name = system.variable_name(var);
      
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

      // The quadrature points on element f.  These will be computed
      // from the quadrature points on element e.
      std::vector<Point> qp_f;
      
      // The shape function gradients on elements e & f
      const std::vector<std::vector<Real> > & phi_e = fe_e->get_phi();
      const std::vector<std::vector<Real> > & phi_f = fe_f->get_phi();
      
      // The global DOF indices for elements e & f
      std::vector<unsigned int> dof_indices_e;
      std::vector<unsigned int> dof_indices_f;


      
      // Iterate over all the active elements in the mesh
      // that live on this processor.
      MeshBase::const_element_iterator       elem_it  = mesh.active_local_elements_begin();
      const MeshBase::const_element_iterator elem_end = mesh.active_local_elements_end(); 

      for (; elem_it != elem_end; ++elem_it)
	{
	  // e is necessarily an active element on the local processor
	  const Elem* e = *elem_it;
	  const unsigned int e_id = e->id();

	  // Loop over the neighbors of element e
	  for (unsigned int n_e=0; n_e<e->n_neighbors(); n_e++)
	    {
	      if (e->neighbor(n_e) != NULL) // e is not on the boundary
		{
		  const Elem* f           = e->neighbor(n_e);
		  const unsigned int f_id = f->id();
		  bool case1 = false, case2 = false;

		  // Case 1.	    
		  if ( f->active() && (f->level() == e->level()) && (e_id < f_id))
		    case1 = true;

		  // Case 2.
		  else if (f->level() < e->level())
		    case2 = true;
		  
		  // Compute flux jumps if we are in case 1 or case 2.
		  if (case1 || case2)
		    {		    
                      // Update the shape functions on side s_e of
                      // element e
                      START_LOG("fe_e->reinit()", "DiscontinuityMeasure");
                      fe_e->reinit (e, n_e);
                      STOP_LOG("fe_e->reinit()", "DiscontinuityMeasure");

                      // Build the side
                      START_LOG("construct side", "DiscontinuityMeasure");
                      AutoPtr<Elem> side (e->build_side(n_e));
                      STOP_LOG("construct side", "DiscontinuityMeasure");

                      // Get the maximum h for this side
                      START_LOG("side->hmax()", "DiscontinuityMeasure");
                      const Real h = side->hmax();
                      STOP_LOG("side->hmax()", "DiscontinuityMeasure");
                      
		      // Get the DOF indices for the two elements
		      START_LOG("dof_indices()", "DiscontinuityMeasure");
		      dof_map.dof_indices (e, dof_indices_e, var);
		      dof_map.dof_indices (f, dof_indices_f, var);
		      STOP_LOG("dof_indices()", "DiscontinuityMeasure");

		      // The number of DOFS on each element
		      const unsigned int n_dofs_e = dof_indices_e.size();
		      const unsigned int n_dofs_f = dof_indices_f.size();

		      // The number of quadrature points
		      const unsigned int n_qp = qrule.n_points();

		      // Find the location of the quadrature points
		      // on element f
		      START_LOG("inverse_map()", "DiscontinuityMeasure");
		      FEInterface::inverse_map (dim, fe_type, f, qface_point, qp_f);
		      STOP_LOG("inverse_map()", "DiscontinuityMeasure");
		      
		      // Compute the shape functions on element f
		      // at the quadrature points of element e
		      START_LOG("fe_f->reinit()", "DiscontinuityMeasure");
		      fe_f->reinit (f, &qp_f);
		      STOP_LOG("fe_f->reinit()", "DiscontinuityMeasure");
		      
		      // The error contribution from this face
		      Real error = 0.;

		    
		      START_LOG("jump integral", "DiscontinuityMeasure");
		      // loop over the integration points on the face
		      for (unsigned int qp=0; qp<n_qp; qp++)
			{
			  // The solution gradient from each element
			  Number val_e = 0., val_f = 0.;
			
			  // Compute the solution gradient on element e
			  for (unsigned int i=0; i<n_dofs_e; i++)
			    val_e += phi_e[i][qp] *
				     system.current_solution(dof_indices_e[i]);
			
			  // Compute the solution gradient on element f
			  for (unsigned int i=0; i<n_dofs_f; i++)
			    val_f += phi_f[i][qp] *
				     system.current_solution(dof_indices_f[i]);
			
			  // The flux jump at the face 
			  const Number jump = (val_e - val_f);

			  // The flux jump squared.  If using complex numbers,
			  // std::norm(z) returns |z|^2, where |z| is the modulus of z.
#ifndef USE_COMPLEX_NUMBERS
			  const Real jump2 = jump*jump;
#else
			  const Real jump2 = std::norm(jump);
#endif

			  // Integrate the error on the face.  The error is
			  // scaled by an additional power of h, where h is
			  // the maximum side length for the element.  This
			  // arises in the definition of the indicator.
			  error += JxW_face[qp]*h*jump2;			
			
			} // End quadrature point loop
		      STOP_LOG("jump integral", "DiscontinuityMeasure");
		      
		      // Add the error contribution to elements e & f
                      assert(e_id < error_per_cell.size());
                      assert(f_id < error_per_cell.size());
		      error_per_cell[e_id] += error*component_scale[var];
		      error_per_cell[f_id] += error*component_scale[var];

		      // Increment the number of flux faces for e
		      n_internal_sides[e_id]++;

		      // In case 1, e and f are at the same level, so
		      // increment the number of flux faces for f as well.
		      if (case1)
			n_internal_sides[f_id]++;

		      // In case 2, the number of additional flux faces for
		      // f depends on the difference in levels between e and
		      // f and the physical dimension of the mesh.  
		      else if (case2)
			{
			  // With a difference of n levels between elements e and f,
			  // We compute a fractional flux face for element f by adding:
			  // 1/2^n in 2D
			  // 1/4^n in 3D
			  // each time.  This code will get hit 2^n times in 2D and 4^n
			  // times in 3D so that the final flux face count for element f
			  // will be an integer value.
			  const unsigned int divisor = 1 << (dim-1)*(e->level() - f->level());
			  
			  // Add a fractional flux face to element f.
			  n_internal_sides[f_id] += 1.0 / static_cast<Real>(divisor);
			}
		    } // end if (case1 || case2)
		} // if (e->neigbor(n_e) != NULL)

	      // Otherwise, e is on the boundary.  If it happens to
	      // be on a Dirichlet boundary, we need not do anything.
	      // On the other hand, if e is on a Neumann (flux) boundary
	      // with grad(u).n = g, we need to compute the additional residual
	      // (h * \int |g - grad(u_h).n|^2 dS)^(1/2).
	      // We can only do this with some knowledge of the boundary
	      // conditions, i.e. the user must have attached an appropriate
	      // BC function.
	      else
		{
		  if (this->_bc_function != NULL)
		    {
		      START_LOG("boundary integrals", "DiscontinuityMeasure");
		      // here();
		  
		      // Update the shape functions on side s_e of element e
		      fe_e->reinit (e, n_e);

		      // The reinitialization also recomputes the locations of
		      // the quadrature points on the side.  By checking if the
		      // first quadrature point on the side is on a flux boundary
		      // for a particular variable, we will determine if the whole
		      // element is on a flux boundary (assuming quadrature points
		      // are strictly contained in the side).
		      if (this->_bc_function(system, qface_point[0], var_name).first)
			{
			  // Build the side
			  AutoPtr<Elem> side (e->build_side(n_e));

			  // Get the maximum h for this side
			  const Real h = side->hmax();
		    
			  // Get the DOF indices 
			  dof_map.dof_indices (e, dof_indices_e, var);

			  // The number of DOFS on each element
			  const unsigned int n_dofs_e = dof_indices_e.size();

			  // The number of quadrature points
			  const unsigned int n_qp = qrule.n_points();

			  // The error contribution from this face
			  Real error = 0.;
		    
			  // loop over the integration points on the face.
			  for (unsigned int qp=0; qp<n_qp; qp++)
			    {
			      // Value of the imposed flux BC at this quadrature point.
			      const std::pair<bool,Real> flux_bc =
				this->_bc_function(system, qface_point[qp], var_name);

			      // Be sure the BC function still thinks we're on the 
			      // flux boundary.
			      assert (flux_bc.first == true);
			      
			      // The solution from each element
			      Number val_e = 0.;
			
			      // Compute the solution gradient on element e
			      for (unsigned int i=0; i<n_dofs_e; i++)
				val_e += phi_e[i][qp] *
					 system.current_solution(dof_indices_e[i]);

			      // The difference between the desired BC and the approximate solution. 
			      const Number jump = flux_bc.second - val_e;

			      // The flux jump squared.  If using complex numbers,
			      // std::norm(z) returns |z|^2, where |z| is the modulus of z.
#ifndef USE_COMPLEX_NUMBERS
			      const Real jump2 = jump*jump;
#else
			      const Real jump2 = std::norm(jump);
#endif

			      
// 			      std::cout << "Error contribution from "
// 					<< var_name
// 					<< " flux BC: "
// 					<< JxW_face[qp]*h*jump2
// 					<< std::endl;

			      
			      // Integrate the error on the face.  The error is
			      // scaled by an additional power of h, where h is
			      // the maximum side length for the element.  This
			      // arises in the definition of the indicator.
			      error += JxW_face[qp]*h*jump2;			
			
			    } // End quadrature point loop

			  // Add the error contribution to element e 
                          assert(e_id < error_per_cell.size());
			  error_per_cell[e_id] += error*component_scale[var];

			  // Increment the number of flux faces for element e
			  n_internal_sides[e_id]++;
			  
			} // end if side on flux boundary
		      
		      STOP_LOG("boundary integrals", "DiscontinuityMeasure");
		    } // end if _bc_function != NULL
		} // end if (e->neighbor(n_e) == NULL)
	    } // end loop over neighbors
	} // End loop over active local elements
    } // End loop over variables


  
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
  START_LOG("std::sqrt()", "DiscontinuityMeasure");
  for (unsigned int i=0; i<error_per_cell.size(); i++)
    if (error_per_cell[i] != 0.)
      error_per_cell[i] = std::sqrt(error_per_cell[i]);
  STOP_LOG("std::sqrt()", "DiscontinuityMeasure");


  if (this->scale_by_n_internal_sides)
    {
      // Sum the vector of flux face counts
      this->reduce_error(n_internal_sides);

      // Sanity check: Make sure the number of flux faces is
      // always an integer value
#ifdef DEBUG
      for (unsigned int i=0; i<n_internal_sides.size(); ++i)
	assert (n_internal_sides[i] == static_cast<float>(static_cast<unsigned int>(n_internal_sides[i])) );
#endif
  
      // Scale the error by the number of flux faces for each element
      for (unsigned int i=0; i<n_internal_sides.size(); ++i)
	{
	  if (n_internal_sides[i] == 0.0) // inactive or non-local element
	    continue;
      
	  //std::cout << "Element " << i << " has " << n_internal_sides[i] << " flux faces." << std::endl;
	  error_per_cell[i] /= static_cast<Real>(n_internal_sides[i]); 
	}
    }
  
  //  STOP_LOG("flux_jumps()", "DiscontinuityMeasure");
}







void
DiscontinuityMeasure::attach_flux_bc_function (std::pair<bool,Real> fptr(const System& system,
									const Point& p,
									const std::string& var_name))
{
  assert (fptr != NULL);
  
  _bc_function = fptr;
}
