// $Id: system_base_projection.C,v 1.1 2003-04-30 21:09:29 benkirk Exp $

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



// C++ includes
// Local includes
#include "system_base.h"
#include "libmesh.h"
#include "dof_map.h"
#include "dense_matrix.h"
#include "fe.h"
#include "fe_interface.h"
#include "quadrature_gauss.h"



// ------------------------------------------------------------
// SystemBase implementation
void SystemBase::project_vector (const NumericVector<Number>* old_vector,
				 NumericVector<Number>* new_vector) const
{
#ifdef ENABLE_AMR

  /**
   * This method projects a solution from an old mesh to a current, refined
   * mesh.  The input vector \p old_vector gives the solution on the
   * old mesh, while the \p new_vector gives the solution (to be computed)
   * on the new mesh.
   */
  new_vector->clear();
  new_vector->init (this->n_dofs(),
		    this->n_local_dofs());



  // The number of variables in this system
  const unsigned int n_variables = this->n_vars();
  
  // The dimensionality of the current mesh
  const unsigned int dim = _mesh.mesh_dimension();

  // The DofMap for this system
  const DofMap& dof_map = this->get_dof_map();

  // The element matrix and RHS
  DenseMatrix<Number> Ke;
  std::vector<Number> Fe;


  
  // Loop over all the variables in the system
  for (unsigned int var=0; var<n_variables; var++)
    {
      // Get a FE object of the appropriate type
      const FEType& fe_type = this->get_dof_map().variable_type(var);     
      AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));      
      QGauss qrule       (dim, THIRD);

      fe->attach_quadrature_rule       (&qrule);

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
      std::vector<unsigned int> dof_indices, old_dof_indices;

      
      
      // Iterators for the active elements on local processor
      active_local_elem_iterator       elem_it (_mesh.elements_begin());
      const active_local_elem_iterator elem_end(_mesh.elements_end());
      
      for ( ; elem_it != elem_end; ++elem_it)
	{
	  const Elem* elem = *elem_it;
	  const Elem* parent = elem->parent();

	  // Update the fe object based on the current child
	  fe->reinit (elem);

	  // Update the DOF indices for this element
	  dof_map.dof_indices (elem, dof_indices, var);
	  
	  // The number of quadrature points on the child
	  const unsigned int n_qp = qrule.n_points();

	  // The number of DOFs on this child
	  const unsigned int n_dofs = dof_indices.size();

	  // Reinitialize the element matrix and vector for
	  // the current element.  Note that this will zero them
	  // before they are summed.
	  Ke.resize (n_dofs, n_dofs);
	  Fe.resize (n_dofs);  std::fill (Fe.begin(), Fe.end(), 0.);

	  
	  
	  // Loop over the quadrature points
	  for (unsigned int qp=0; qp<n_qp; qp++)
	    {
	      // The solution value at the quadrature point	      
	      Number val = libMesh::zero;

	      // Interpolate the old solution at the quadrature point
	      // on the child.  This is done in one of two ways.
	      // 1.) If the child was just refined then it was not
	      //     active in the previous mesh & hence has no solution
	      //     values on it.  In this case simply interpolate the
	      //     solution from the parent, who was active in the
	      //     previous mesh
	      // 2.) Otherwise the child was active in the previous
	      //     mesh, and we can just interpolate directly
	      if (elem->refinement_flag() ==
		  Elem::JUST_REFINED)
		{
		  // Sanity check
		  assert (parent != NULL);
	 
		  // The location of the quadrature point
		  // on the parent element
		  const Point q_point =
		    FEInterface::inverse_map (dim, fe_type,
					      parent, xyz_values[qp]);

		  dof_map.old_dof_indices (parent, old_dof_indices, var);
		  
		  // Sum the function values * the DOF values
		  // at the point of interest to get the function value
		  // (Note that the # of DOFs on the parent need not be the
		  //  same as on the child!)
		  for (unsigned int i=0; i<old_dof_indices.size(); i++)
		    {
		      val += ((*old_vector)(old_dof_indices[i])*
			      FEInterface::shape(dim, fe_type, parent,
						 i, q_point));
		    }
		}
	      else
		{
		  dof_map.old_dof_indices (elem, old_dof_indices, var);

		  assert (old_dof_indices.size() == n_dofs);
		  
		  // Sum all the function values * the DOF values
		  // at the quadrature point on the child to get the
		  // function value
		  for (unsigned int i=0; i<n_dofs; i++)
		    {
		      val += ((*old_vector)(old_dof_indices[i])*
			      phi_values[i][qp]);
		    }
		}

	      // Now \p val contains the solution value of variable
	      // \p var at the qp'th quadrature point on the child
	      // element.  It has been interpolated from the parent
	      // in case the child was just refined.  Next we will
	      // construct the L2-projection matrix for the element.

	      // Construct the Mass Matrix
	      for (unsigned int i=0; i<n_dofs; i++)
		for (unsigned int j=0; j<n_dofs; j++)
		  Ke(i,j) += JxW[qp]*phi_values[i][qp]*phi_values[j][qp];

	      // Construct the RHS
	      for (unsigned int i=0; i<n_dofs; i++)
		Fe[i] += JxW[qp]*phi_values[i][qp]*val;
	      
	      
	    } // end qp loop

	  // The mass matrix & RHS are now constructed.  All that
	  // remains is to constrain any hanging nodes & insert
	  // the element contributions into the global matrix & RHS
	  dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

	  rhs->add_vector    (Fe, dof_indices);
	  matrix->add_matrix (Ke, dof_indices);
	  
	} // end elem loop
    } // end variables loop
  
#endif // #ifdef ENABLE_AMR
}
