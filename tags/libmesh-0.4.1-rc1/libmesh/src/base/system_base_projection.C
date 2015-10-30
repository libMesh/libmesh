// $Id: system_base_projection.C,v 1.5 2003-09-02 18:02:40 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002-2003  Benjamin S. Kirk, John W. Peterson
  
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
#include "system_base.h"
#include "mesh.h"
#include "libmesh.h"
#include "dof_map.h"
#include "fe_type.h"
#include "fe_interface.h"
#include "numeric_vector.h"
#include "mesh_logging.h"



// ------------------------------------------------------------
// SystemBase implementation
void SystemBase::project_vector (NumericVector<Number>& vector) const
{
  // Create a copy of the vector, which currently
  // contains the old data.
  AutoPtr<NumericVector<Number> >
    old_vector (vector.clone());

  // Project the old vector to the new vector
  this->project_vector (*old_vector, vector);
}



void SystemBase::project_vector (const NumericVector<Number>& old_vector,
				 NumericVector<Number>& new_vector) const
{
#ifdef ENABLE_AMR
  
  START_LOG ("project_vector()", "SystemBase");
  
  /**
   * This method projects a solution from an old mesh to a current, refined
   * mesh.  The input vector \p old_vector gives the solution on the
   * old mesh, while the \p new_vector gives the solution (to be computed)
   * on the new mesh.
   */
  new_vector.clear();

  // If the old vector was uniprocessor, make the new
  // vector uniprocessor
  if (old_vector.size() == old_vector.local_size())
    new_vector.init (this->n_dofs(),
		     this->n_dofs());

  // Otherwise it is a parallel, distributed vector.
  else
    new_vector.init (this->n_dofs(),
		     this->n_local_dofs());
    
  // Note that the init above will have zeroed the new_vector

  // A vector indicating if we have visited a DOF yet
  std::vector<bool> already_done (this->n_dofs(), false);
  

  // The number of variables in this system
  const unsigned int n_variables = this->n_vars();

  // The dimensionality of the current mesh
  const unsigned int dim = _mesh.mesh_dimension();

  // The DofMap for this system
  const DofMap& dof_map = this->get_dof_map();


  
  // Loop over all the variables in the system
  for (unsigned int var=0; var<n_variables; var++)
    {
      // Get a FE object of the appropriate type
      const FEType& fe_type = dof_map.variable_type(var);

      // This method is only good for Lagrange elements
      assert (fe_type.family == LAGRANGE);

      // The global DOF indices on the
      // old and new meshes
      std::vector<unsigned int> new_dof_indices, old_dof_indices;

      
      
      // Iterators for the active elements on local processor
      const_active_local_elem_iterator       elem_it (_mesh.elements_begin());
      const const_active_local_elem_iterator elem_end(_mesh.elements_end());
      
      for ( ; elem_it != elem_end; ++elem_it)
	{
	  const Elem* elem = *elem_it;
	  const Elem* parent = elem->parent();

	  // Update the DOF indices for this element based on
	  // the new mesh
	  dof_map.dof_indices (elem, new_dof_indices, var);
	  
	  // The number of DOFs on this element
	  const unsigned int new_n_dofs = new_dof_indices.size();

	  // This only works for Lagrange elements, where the
	  // number of dofs <= the number of nodes. Check it.
	  assert (new_n_dofs <= elem->n_nodes());


	  // Update the DOF indices based on
	  // the old mesh.  There are two possibilities here:
	  // 1.) If the child was just refined then it was not
	  //     active in the previous mesh & hence has no DOFs
	  //     on the old mesh. In this case simply use the DOFs
	  //     from the parent, who was active in the previous mesh
	  // 2.) Otherwise the child was active in the previous
	  //     mesh, and we can just use its DOFs from the old
	  //     mesh
	  if (elem->refinement_flag() == Elem::JUST_REFINED)
	    {
	      // Case 1.  This is the first time we might need the parent
	      assert (parent != NULL);
	      dof_map.old_dof_indices (parent, old_dof_indices, var);

	      assert (old_dof_indices.size() <= parent->n_nodes());
	    }
	  else
	    {
	      // Case 2.
	      dof_map.old_dof_indices (elem, old_dof_indices, var);
	      assert (old_dof_indices.size() <= elem->n_nodes());
	    }

	  // The number of DOFs on the old mesh
	  const unsigned int old_n_dofs = old_dof_indices.size();



	  
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
	      else
		already_done[new_global_dof] = true;
	      
	      
	      // The solution value at the node, computed from the old mesh
	      Number val = libMesh::zero;

	      // Interpolate the old solution at the node
	      // on the child.  This is done in one of two ways.
	      // 1.) If the child was just refined then it was not
	      //     active in the previous mesh & hence has no solution
	      //     values on it.  In this case simply interpolate the
	      //     solution from the parent, who was active in the
	      //     previous mesh
	      // 2.) Otherwise the child was active in the previous
	      //     mesh, and we can just interpolate directly
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
		      
		      val += (old_vector(old_global_dof)*
			      FEInterface::shape(dim, fe_type, parent,
						 old_local_dof, point));
		    }
		}
	      else
		{
		  // Get the old global DOF index
		  const unsigned int old_global_dof =
		    old_dof_indices[new_local_dof];

		  val = old_vector(old_global_dof);
		}

	      // Assign the value to the new vector
	      new_vector.set (new_global_dof, val);
	      
	    } // end new dof loop	  
	} // end elem loop
    } // end variables loop

  // Close the new vector so that it is ready for use.
  new_vector.close();
  
  STOP_LOG ("project_vector()", "SystemBase");


#else

  // AMR is disabled...  Simply copy the vector
  new_vector = old_vector;
  
#endif // #ifdef ENABLE_AMR
}



/**
 * This method projects the vector
 * via a full L2 projection. In this function
 * the L2 projection system is constructed.
 * It must be solved subsequently.
 */
// void SystemBase::project_vector (const NumericVector<Number>& old_vector,
// 				 NumericVector<Number>& new_vector) const
// {
// #ifdef ENABLE_AMR

//   /**
//    * This method projects a solution from an old mesh to a current, refined
//    * mesh.  The input vector \p old_vector gives the solution on the
//    * old mesh, while the \p new_vector gives the solution (to be computed)
//    * on the new mesh.
//    */
//   new_vector.clear();

//   // If the old vector was uniprocessor, make the new
//   // vector uniprocessor
//   if (old_vector.size() == old_vector.local_size())
//     new_vector.init (this->n_dofs(),
// 		     this->n_dofs());

//   // Otherwise it is a parallel, distributed vector.
//   else
//     new_vector.init (this->n_dofs(),
// 		     this->n_local_dofs());
    



//   // The number of variables in this system
//   const unsigned int n_variables = this->n_vars();
  
//   // The dimensionality of the current mesh
//   const unsigned int dim = _mesh.mesh_dimension();

//   // The DofMap for this system
//   const DofMap& dof_map = this->get_dof_map();

//   // The element matrix and RHS
//   DenseMatrix<Number> Ke;
//   std::vector<Number> Fe;


  
//   // Loop over all the variables in the system
//   for (unsigned int var=0; var<n_variables; var++)
//     {
//       // Get a FE object of the appropriate type
//       const FEType& fe_type = dof_map.variable_type(var);     
//       AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));      
//       QGauss qrule       (dim, THIRD);

//       fe->attach_quadrature_rule       (&qrule);

//       // The values of the shape functions at the quadrature
//       // points on the child element.
//       const std::vector<std::vector<Real> >& phi_values =
// 	fe->get_phi();

//       // The Jacobian * quadrature weight at the quadrature points
//       const std::vector<Real>& JxW =
// 	fe->get_JxW();
      
//       // The XYZ locations of the quadrature points on the
//       // child element
//       const std::vector<Point>& xyz_values =
// 	fe->get_xyz();

//       // The global DOF indices
//       std::vector<unsigned int> dof_indices, old_dof_indices;

      
      
//       // Iterators for the active elements on local processor
//       const_active_local_elem_iterator       elem_it (_mesh.elements_begin());
//       const const_active_local_elem_iterator elem_end(_mesh.elements_end());
      
//       for ( ; elem_it != elem_end; ++elem_it)
// 	{
// 	  const Elem* elem = *elem_it;
// 	  const Elem* parent = elem->parent();

// 	  // Update the fe object based on the current child
// 	  fe->reinit (elem);

// 	  // Update the DOF indices for this element
// 	  dof_map.dof_indices (elem, dof_indices, var);
	  
// 	  // The number of quadrature points on the child
// 	  const unsigned int n_qp = qrule.n_points();

// 	  // The number of DOFs on this child
// 	  const unsigned int n_dofs = dof_indices.size();

// 	  // Reinitialize the element matrix and vector for
// 	  // the current element.  Note that this will zero them
// 	  // before they are summed.
// 	  Ke.resize (n_dofs, n_dofs);
// 	  Fe.resize (n_dofs);  std::fill (Fe.begin(), Fe.end(), 0.);

	  
	  
// 	  // Loop over the quadrature points
// 	  for (unsigned int qp=0; qp<n_qp; qp++)
// 	    {
// 	      // The solution value at the quadrature point	      
// 	      Number val = libMesh::zero;

// 	      // Interpolate the old solution at the quadrature point
// 	      // on the child.  This is done in one of two ways.
// 	      // 1.) If the child was just refined then it was not
// 	      //     active in the previous mesh & hence has no solution
// 	      //     values on it.  In this case simply interpolate the
// 	      //     solution from the parent, who was active in the
// 	      //     previous mesh
// 	      // 2.) Otherwise the child was active in the previous
// 	      //     mesh, and we can just interpolate directly
// 	      if (elem->refinement_flag() ==
// 		  Elem::JUST_REFINED)
// 		{
// 		  // Sanity check
// 		  assert (parent != NULL);
	 
// 		  // The location of the quadrature point
// 		  // on the parent element
// 		  const Point q_point =
// 		    FEInterface::inverse_map (dim, fe_type,
// 					      parent, xyz_values[qp]);

// 		  dof_map.old_dof_indices (parent, old_dof_indices, var);
		  
// 		  // Sum the function values * the DOF values
// 		  // at the point of interest to get the function value
// 		  // (Note that the # of DOFs on the parent need not be the
// 		  //  same as on the child!)
// 		  for (unsigned int i=0; i<old_dof_indices.size(); i++)
// 		    {
// 		      val += (old_vector(old_dof_indices[i])*
// 			      FEInterface::shape(dim, fe_type, parent,
// 						 i, q_point));
// 		    }
// 		}
// 	      else
// 		{
// 		  dof_map.old_dof_indices (elem, old_dof_indices, var);

// 		  assert (old_dof_indices.size() == n_dofs);
		  
// 		  // Sum all the function values * the DOF values
// 		  // at the quadrature point on the child to get the
// 		  // function value
// 		  for (unsigned int i=0; i<n_dofs; i++)
// 		    {
// 		      val += (old_vector(old_dof_indices[i])*
// 			      phi_values[i][qp]);
// 		    }
// 		}

// 	      // Now \p val contains the solution value of variable
// 	      // \p var at the qp'th quadrature point on the child
// 	      // element.  It has been interpolated from the parent
// 	      // in case the child was just refined.  Next we will
// 	      // construct the L2-projection matrix for the element.

// 	      // Construct the Mass Matrix
// 	      for (unsigned int i=0; i<n_dofs; i++)
// 		for (unsigned int j=0; j<n_dofs; j++)
// 		  Ke(i,j) += JxW[qp]*phi_values[i][qp]*phi_values[j][qp];

// 	      // Construct the RHS
// 	      for (unsigned int i=0; i<n_dofs; i++)
// 		Fe[i] += JxW[qp]*phi_values[i][qp]*val;
	      
	      
// 	    } // end qp loop

// 	  // The mass matrix & RHS are now constructed.  All that
// 	  // remains is to constrain any hanging nodes & insert
// 	  // the element contributions into the global matrix & RHS
// 	  dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

// 	  rhs->add_vector    (Fe, dof_indices);
// 	  matrix->add_matrix (Ke, dof_indices);
	  
// 	} // end elem loop
//     } // end variables loop
  
// #endif // #ifdef ENABLE_AMR
// }
