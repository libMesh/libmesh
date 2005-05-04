// $Id: system_projection.C,v 1.12 2005-05-04 21:24:42 roystgnr Exp $

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
#include <vector>

// Local includes
#include "system.h"
#include "mesh.h"
#include "elem.h"
#include "libmesh.h"
#include "dof_map.h"
#include "fe_type.h"
#include "fe_interface.h"
#include "numeric_vector.h"
#include "libmesh_logging.h"

#include "dense_matrix.h"
#include "mesh_base.h"
#include "fe_base.h"
#include "quadrature_gauss.h"



// ------------------------------------------------------------
// System implementation
void System::project_vector (NumericVector<Number>& vector) const
{
  // Create a copy of the vector, which currently
  // contains the old data.
  AutoPtr<NumericVector<Number> >
    old_vector (vector.clone());

  // Project the old vector to the new vector
  this->project_vector (*old_vector, vector);
}


/**
 * This method projects the vector
 * via L2 projections or nodal
 * interpolations on each element.
 */
void System::project_vector (const NumericVector<Number>& old_vector,
				 NumericVector<Number>& new_vector) const
{
  START_LOG ("project_vector()", "System");

  /**
   * This method projects a solution from an old mesh to a current, refined
   * mesh.  The input vector \p old_vector gives the solution on the
   * old mesh, while the \p new_vector gives the solution (to be computed)
   * on the new mesh.
   */
  new_vector.clear();

  // Resize the new vector.

  // If the old vector was uniprocessor, make the new
  // vector uniprocessor
  if (old_vector.size() == old_vector.local_size())
    new_vector.init (this->n_dofs(),
		     this->n_dofs());

  // Otherwise it is a parallel, distributed vector.
  else
    new_vector.init (this->n_dofs(),
		     this->n_local_dofs());
  // Note that the above will have zeroed the new_vector
   
#ifdef ENABLE_AMR

  // A vector for Lagrange element interpolation, indicating if we
  // have visited a DOF yet
  std::vector<bool> already_done (this->n_dofs(), false);

  // The number of variables in this system
  const unsigned int n_variables = this->n_vars();

  // The dimensionality of the current mesh
  const unsigned int dim = _mesh.mesh_dimension();

  // The DofMap for this system
  const DofMap& dof_map = this->get_dof_map();

  // The element matrix and RHS for non-Lagrange element projection
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;
  DenseVector<Number> Ue;


  // Loop over all the variables in the system
  for (unsigned int var=0; var<n_variables; var++)
    {
      // Get a FE object of the appropriate type
      const FEType& fe_type = dof_map.variable_type(var);     
      AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));      

      // Prepare variables for non-Lagrange projection
      int order = fe->get_order();
      QGauss qrule       (dim, libMeshEnums::Order(order*2));
      fe->attach_quadrature_rule (&qrule);

      // The values of the shape functions at the quadrature
      // points on the child element.
      const std::vector<std::vector<Real> >& phi_values =
	fe->get_phi();

      // The Jacobian * quadrature weight at the quadrature points
      const std::vector<Real>& JxW =
	fe->get_JxW();
     
      // The XYZ locations of the quadrature points on the
      // child element
      const std::vector<Point>& xyz_values =
	fe->get_xyz();


      // The global DOF indices
      std::vector<unsigned int> new_dof_indices, old_dof_indices;
   
      // Iterators for the active elements on local processor
      MeshBase::element_iterator       elem_it =
		      _mesh.active_local_elements_begin();
      const MeshBase::element_iterator elem_end = 
		      _mesh.active_local_elements_end();
   
      for ( ; elem_it != elem_end; ++elem_it)
	{
	  const Elem* elem = *elem_it;
	  const Elem* parent = elem->parent();

	  // Update the DOF indices for this element based on
          // the new mesh
	  dof_map.dof_indices (elem, new_dof_indices, var);

	  // The number of DOFs on this child
	  const unsigned int new_n_dofs = new_dof_indices.size();

          // Zero the interpolated values
          Ue.resize (new_n_dofs); Ue.zero();

	  // Update the DOF indices based on the old mesh.
	  // This is done in one of three ways:
	  // 1.) If the child was just refined then it was not
	  //     active in the previous mesh & hence has no solution
	  //     values on it.  In this case simply project or
	  //     interpolate the solution from the parent, who was
          //     active in the previous mesh
	  // 2.) If the child was just coarsened, obtaining a
	  //     well-defined solution may require doing independent
	  //     projections on nodes, edges, faces, and interiors
	  // 3.) If the child was active in the previous
	  //     mesh, we can just copy coefficients directly
	  if (elem->refinement_flag() == Elem::JUST_REFINED)
	    {
	      assert (parent != NULL);
	 
	      dof_map.old_dof_indices (parent, old_dof_indices, var);
            }
	  else if (elem->refinement_flag() == Elem::JUST_COARSENED)
	    {
	      assert (elem->has_children());
            }
	  else
	    {
	      dof_map.old_dof_indices (elem, old_dof_indices, var);

	      assert (old_dof_indices.size() == new_n_dofs);
	    }

	  const unsigned int old_n_dofs = old_dof_indices.size();

          if (fe_type.family != LAGRANGE) {

	    // For refined elements, we do an L2 projection
	    if (elem->refinement_flag() == Elem::JUST_REFINED)
	      {
	        // Update the fe object based on the current child
	        fe->reinit (elem);
	  
	        // The number of quadrature points on the child
	        const unsigned int n_qp = qrule.n_points();

	        // Reinitialize the element matrix and vector for
	        // the current element.  Note that this will zero them
	        // before they are summed.
	        Ke.resize (new_n_dofs, new_n_dofs); Ke.zero();
	        Fe.resize (new_n_dofs); Fe.zero();
	  
	  
	        // Loop over the quadrature points
	        for (unsigned int qp=0; qp<n_qp; qp++)
	          {
	            // The solution value at the quadrature point	      
	            Number val = libMesh::zero;

		    // The location of the quadrature point
		    // on the parent element
		    const Point q_point =
		    FEInterface::inverse_map (dim, fe_type,
					      parent, xyz_values[qp]);

		    // Sum the function values * the DOF values
		    // at the point of interest to get the function value
		    // (Note that the # of DOFs on the parent need not be the
		    //  same as on the child!)
		    for (unsigned int i=0; i<old_n_dofs; i++)
		      {
		        val += (old_vector(old_dof_indices[i])*
			        FEInterface::shape(dim, fe_type, parent,
						   i, q_point));
		      }

	            // Now \p val contains the solution value of variable
	            // \p var at the qp'th quadrature point on the child
	            // element.  It has been interpolated from the parent
	            // in case the child was just refined.  Next we will
	            // construct the L2-projection matrix for the element.

	            // Construct the Mass Matrix
	            for (unsigned int i=0; i<new_n_dofs; i++)
		      for (unsigned int j=0; j<new_n_dofs; j++)
		        Ke(i,j) += JxW[qp]*phi_values[i][qp]*phi_values[j][qp];

	            // Construct the RHS
	            for (unsigned int i=0; i<new_n_dofs; i++)
		      Fe(i) += JxW[qp]*phi_values[i][qp]*val;
	      
	          } // end qp loop

                Ke.cholesky_solve(Fe, Ue);
	      }
            else if (elem->refinement_flag() == Elem::JUST_COARSENED)
	      {
		// FIXME: proper non-Lagrange coarsening will take
                // some work
                error();
              }
	    // For unrefined uncoarsened elements, we just copy DoFs
	    else
	      {
		for (unsigned int i=0; i<new_n_dofs; i++)
		  Ue(i) = old_vector(old_dof_indices[i]);
	      }
          } else { // fe type is Lagrange
	    // Loop over the DOFs on the element
	    for (unsigned int new_local_dof=0;
	         new_local_dof<new_n_dofs; new_local_dof++)
	      {
	        // The global DOF number for this local DOF
	        const unsigned int new_global_dof =
		  new_dof_indices[new_local_dof];

	        // The global DOF might lie outside of the bounds of a
	        // distributed vector.  Check for that and possibly continue
	        // the loop
	        if ((new_global_dof <  new_vector.first_local_index()) ||
		    (new_global_dof >= new_vector.last_local_index()))
		  continue;

	        // We might have already computed the solution for this DOF.
	        // This is likely in the case of a shared node, particularly
	        // at the corners of an element.  Check to see if that is the
	        // case
	        if (already_done[new_global_dof] == true)
		  continue;

		already_done[new_global_dof] = true;
	      
	        if (elem->refinement_flag() == Elem::JUST_REFINED)
		  {
		    // The location of the child's node on the parent element
		    const Point point =
		      FEInterface::inverse_map (dim, fe_type, parent,
					        elem->point(new_local_dof));
		  
		    // Sum the function values * the DOF values
		    // at the point of interest to get the function value
		    // (Note that the # of DOFs on the parent need not be the
		    //  same as on the child!)
		    for (unsigned int old_local_dof=0;
		         old_local_dof<old_n_dofs; old_local_dof++)
		      {
		        const unsigned int old_global_dof =
			  old_dof_indices[old_local_dof];
		      
		        Ue(new_local_dof) += 
                          (old_vector(old_global_dof)*
			  FEInterface::shape(dim, fe_type, parent,
					     old_local_dof, point));
		      }
		  }
	        else
		  {
		    // Get the old global DOF index
		    const unsigned int old_global_dof =
		      old_dof_indices[new_local_dof];

		    Ue(new_local_dof) = old_vector(old_global_dof);
		  }
             } // end local DOF loop

          }  // end fe_type if()

	  for (unsigned int i = 0; i < new_n_dofs; i++) 
	    if (Ue(i) != 0.)
              new_vector.set(new_dof_indices[i], Ue(i));
        }  // end elem loop
    } // end variables loop

  new_vector.close();

#else

  // AMR is disabled: simply copy the vector
  new_vector = old_vector;
  
#endif // #ifdef ENABLE_AMR

  STOP_LOG("project_vector()", "System");
}
