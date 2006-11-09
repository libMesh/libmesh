// $Id: system_projection.C,v 1.36 2006-11-09 08:17:31 roystgnr Exp $

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
#include "dof_map.h"
#include "fe_interface.h"
#include "numeric_vector.h"
#include "libmesh_logging.h"

#include "dense_matrix.h"
#include "fe_base.h"
#include "quadrature_gauss.h"
#include "dense_vector.h"



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
void System::project_vector (const NumericVector<Number>& old_v,
			     NumericVector<Number>& new_v) const
{
  START_LOG ("project_vector()", "System");

  /**
   * This method projects a solution from an old mesh to a current, refined
   * mesh.  The input vector \p old_v gives the solution on the
   * old mesh, while the \p new_v gives the solution (to be computed)
   * on the new mesh.
   */
  new_v.clear();

  // Resize the new vector and get a serial version.

  NumericVector<Number> *new_vector_ptr;
  AutoPtr<NumericVector<Number> > new_vector_built;
  NumericVector<Number> *local_old_vector;
  AutoPtr<NumericVector<Number> > local_old_vector_built;
  const NumericVector<Number> *old_vector_ptr;

  // If the old vector was uniprocessor, make the new
  // vector uniprocessor
  if (old_v.size() == old_v.local_size())
    {
      new_v.init (this->n_dofs(),
		  this->n_dofs());
      new_vector_ptr = &new_v;
      old_vector_ptr = &old_v;
    }

  // Otherwise it is a parallel, distributed vector, which
  // we need to localize.
  else
    {
      new_v.init (this->n_dofs(),
		  this->n_local_dofs());
      new_vector_built = NumericVector<Number>::build();
      local_old_vector_built = NumericVector<Number>::build();
      new_vector_ptr = new_vector_built.get();
      local_old_vector = local_old_vector_built.get();
      new_vector_ptr->init(this->n_dofs(), this->n_dofs());
      local_old_vector->init(old_v.size(), old_v.size());
      old_v.localize(*local_old_vector);
      local_old_vector->close();
      old_vector_ptr = local_old_vector;
    }
  // Note that the above will have zeroed the new_vector

  NumericVector<Number> &new_vector = *new_vector_ptr;
  const NumericVector<Number> &old_vector = *old_vector_ptr;
   
#ifdef ENABLE_AMR

  // A vector for Lagrange element interpolation, indicating if we
  // have visited a DOF yet
  // FIXME: we should use this for non-Lagrange coarsening, too
  std::vector<bool> already_done (this->n_dofs(), false);

  // The number of variables in this system
  const unsigned int n_variables = this->n_vars();

  // The dimensionality of the current mesh
  const unsigned int dim = _mesh.mesh_dimension();

  // The DofMap for this system
  const DofMap& dof_map = this->get_dof_map();

  // The element matrix and RHS for projections.
  // Note that Ke is always real-valued, whereas
  // Fe may be complex valued if complex number
  // support is enabled
  DenseMatrix<Real> Ke;
  DenseVector<Number> Fe;
  // The new element coefficients
  DenseVector<Number> Ue;


  // Loop over all the variables in the system
  for (unsigned int var=0; var<n_variables; var++)
    {
      // Get FE objects of the appropriate type
      const FEType& base_fe_type = dof_map.variable_type(var);     
      AutoPtr<FEBase> fe (FEBase::build(dim, base_fe_type));      
      AutoPtr<FEBase> fe_coarse (FEBase::build(dim, base_fe_type));      

      // Create FE objects with potentially different p_level
      FEType fe_type, temp_fe_type;

      // Prepare variables for non-Lagrange projection
      AutoPtr<QBase> qrule     (base_fe_type.default_quadrature_rule(dim));
      AutoPtr<QBase> qedgerule (base_fe_type.default_quadrature_rule(1));
      AutoPtr<QBase> qsiderule (base_fe_type.default_quadrature_rule(dim-1));
      std::vector<Point> coarse_qpoints;

      // The values of the shape functions at the quadrature
      // points
      const std::vector<std::vector<Real> >& phi_values =
	fe->get_phi();
      const std::vector<std::vector<Real> >& phi_coarse =
	fe_coarse->get_phi();

      // The gradients of the shape functions at the quadrature
      // points on the child element.
      const std::vector<std::vector<RealGradient> > *dphi_values =
        NULL;
      const std::vector<std::vector<RealGradient> > *dphi_coarse =
        NULL;

      const FEContinuity cont = fe->get_continuity();

      if (cont == C_ONE)
        {
          const std::vector<std::vector<RealGradient> >&
            ref_dphi_values = fe->get_dphi();
          dphi_values = &ref_dphi_values;
          const std::vector<std::vector<RealGradient> >&
            ref_dphi_coarse = fe_coarse->get_dphi();
          dphi_coarse = &ref_dphi_coarse;
        }

      // The Jacobian * quadrature weight at the quadrature points
      const std::vector<Real>& JxW =
	fe->get_JxW();
     
      // The XYZ locations of the quadrature points on the
      // child element
      const std::vector<Point>& xyz_values =
	fe->get_xyz();


      // The global DOF indices
      std::vector<unsigned int> new_dof_indices, old_dof_indices;
      // Side/edge DOF indices
      std::vector<unsigned int> new_side_dofs, old_side_dofs;
   
      // Iterators for the active elements on local processor
      MeshBase::element_iterator       elem_it =
		      _mesh.active_local_elements_begin();
      const MeshBase::element_iterator elem_end = 
		      _mesh.active_local_elements_end();
   
      for ( ; elem_it != elem_end; ++elem_it)
	{
	  const Elem* elem = *elem_it;
	  const Elem* parent = elem->parent();

          // Adjust the FE type for p-refined elements
          fe_type = base_fe_type;
          fe_type.order = static_cast<Order>(fe_type.order +
                                             elem->p_level());

          // We may need to remember the parent's p_level
          unsigned int old_parent_level = 0;

	  // Update the DOF indices for this element based on
          // the new mesh
	  dof_map.dof_indices (elem, new_dof_indices, var);

	  // The number of DOFs on the new element
	  const unsigned int new_n_dofs = new_dof_indices.size();

          // Fixed vs. free DoFs on edge/face projections
          std::vector<char> dof_is_fixed(new_n_dofs, false); // bools
          std::vector<int> free_dof(new_n_dofs, 0);

	  // The element type
	  const ElemType elem_type = elem->type();

	  // The number of nodes on the new element
	  const unsigned int n_nodes = elem->n_nodes();

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
              old_parent_level = parent->p_level();

              // We may have done p refinement or coarsening as well;
              // if so then we need to reset the parent's p level
              // so we can get the right DoFs from it
              if (elem->p_refinement_flag() == Elem::JUST_REFINED)
                {
                  assert(elem->p_level() > 0);
                  (const_cast<Elem *>(parent))->hack_p_level(elem->p_level() - 1);
                }
              else if (elem->p_refinement_flag() == Elem::JUST_COARSENED)
                {
                  (const_cast<Elem *>(parent))->hack_p_level(elem->p_level() + 1);
                }
	 
	      dof_map.old_dof_indices (parent, old_dof_indices, var);
            }
	  else
	    {
	      dof_map.old_dof_indices (elem, old_dof_indices, var);

              if (elem->p_refinement_flag() == Elem::DO_NOTHING)
	        assert (old_dof_indices.size() == new_n_dofs);
              if (elem->p_refinement_flag() == Elem::JUST_COARSENED)
	        assert (elem->has_children());
	    }

	  unsigned int old_n_dofs = old_dof_indices.size();

          if (fe_type.family != LAGRANGE) {

	    // For refined non-Lagrange elements, we do an L2
            // projection
            // FIXME: this will be a suboptimal and ill-defined
            // result if we're using non-nested finite element
            // spaces or if we're on a p-coarsened element!
	    if (elem->refinement_flag() == Elem::JUST_REFINED)
	      {
	        // Update the fe object based on the current child
                fe->attach_quadrature_rule (qrule.get());
	        fe->reinit (elem);
	  
	        // The number of quadrature points on the child
	        const unsigned int n_qp = qrule->n_points();

                FEInterface::inverse_map (dim, fe_type, parent,
					  xyz_values, coarse_qpoints);

                fe_coarse->reinit(parent, &coarse_qpoints);

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

		    // Sum the function values * the DOF values
		    // at the point of interest to get the function value
		    // (Note that the # of DOFs on the parent need not be the
		    //  same as on the child!)
		    for (unsigned int i=0; i<old_n_dofs; i++)
		      {
		        val += (old_vector(old_dof_indices[i])*
			        phi_coarse[i][qp]);
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

                // Fix up the parent's p level in case we changed it
                (const_cast<Elem *>(parent))->hack_p_level(old_parent_level);
	      }
            else if (elem->refinement_flag() == Elem::JUST_COARSENED)
	      {
                FEBase::coarsened_dof_values(old_vector, dof_map,
					     elem, Ue, var, true);
              }
	    // For unrefined uncoarsened elements, we just copy DoFs
	    else
	      {
                // FIXME - I'm sure this function would be about half
                // the size if anyone ever figures out how to improve
                // the DofMap interface... - RHS
                if (elem->p_refinement_flag() == Elem::JUST_REFINED)
                  {
                    assert (elem->p_level() > 0);
                    // P refinement of non-hierarchic bases will
                    // require a whole separate code path
                    assert (fe->is_hierarchic());
                    temp_fe_type = fe_type;
                    temp_fe_type.order =
                      static_cast<Order>(temp_fe_type.order - 1);
                    unsigned int old_index = 0, new_index = 0;
                    for (unsigned int n=0; n != n_nodes; ++n)
                      {
		        const unsigned int nc =
		          FEInterface::n_dofs_at_node (dim, temp_fe_type,
                                                       elem_type, n);
                        for (unsigned int i=0; i != nc; ++i)
                          {
                            Ue(new_index + i) =
                              old_vector(old_dof_indices[old_index++]);
                          }
                        new_index +=
		          FEInterface::n_dofs_at_node (dim, fe_type,
                                                       elem_type, n);
                      }
                  }
                else if (elem->p_refinement_flag() ==
                         Elem::JUST_COARSENED)
                  {
                    // P coarsening of non-hierarchic bases will
                    // require a whole separate code path
                    assert (fe->is_hierarchic());
                    temp_fe_type = fe_type;
                    temp_fe_type.order =
                      static_cast<Order>(temp_fe_type.order + 1);
                    unsigned int old_index = 0, new_index = 0;
                    for (unsigned int n=0; n != n_nodes; ++n)
                      {
		        const unsigned int nc =
		          FEInterface::n_dofs_at_node (dim, fe_type,
                                                       elem_type, n);
                        for (unsigned int i=0; i != nc; ++i)
                          {
                            Ue(new_index++) =
                              old_vector(old_dof_indices[old_index+i]);
                          }
                        old_index +=
		          FEInterface::n_dofs_at_node (dim, temp_fe_type,
                                                       elem_type, n);
                      }
                  }
                else
                  // If there's no p refinement, we can copy every DoF
		  for (unsigned int i=0; i<new_n_dofs; i++)
		    Ue(i) = old_vector(old_dof_indices[i]);
	      }
          } 
	  else { // fe type is Lagrange
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

            // We may have to clean up a parent's p_level
	    if (elem->refinement_flag() == Elem::JUST_REFINED)
              (const_cast<Elem *>(parent))->hack_p_level(old_parent_level);
          }  // end fe_type if()

	  for (unsigned int i = 0; i < new_n_dofs; i++) 
	    if (Ue(i) != 0.)
              new_vector.set(new_dof_indices[i], Ue(i));
        }  // end elem loop
    } // end variables loop

  new_vector.close();

  // If the old vector was parallel, we need to update it
  // and free the localized copies
  if (old_v.size() != old_v.local_size())
    {
      // We may have to set dof values that this processor doesn't
      // own in certain special cases, like LAGRANGE FIRST or
      // HERMITE THIRD elements on second-order meshes
      for (unsigned int i=0; i!=new_v.size(); i++)
        if (new_vector(i) != 0.0)
          new_v.set(i, new_vector(i));
      new_v.close();
    }

  dof_map.enforce_constraints_exactly(*this, &new_v);

#else

  // AMR is disabled: simply copy the vector
  new_vector = old_vector;
  
#endif // #ifdef ENABLE_AMR

  STOP_LOG("project_vector()", "System");
}



/**
 * This method projects an analytic function onto the solution via L2
 * projections and nodal interpolations on each element.
 */
void System::project_solution (Number fptr(const Point& p,
                                           const Parameters& parameters,
                                           const std::string& sys_name,
                                           const std::string& unknown_name),
                               Gradient gptr(const Point& p,
                                             const Parameters& parameters,
                                             const std::string& sys_name,
                                             const std::string& unknown_name),
                               Parameters& parameters) const
{
  this->project_vector(fptr, gptr, parameters, *solution);

  solution->localize(*current_local_solution);
}



/**
 * This method projects an analytic function via L2 projections and
 * nodal interpolations on each element.
 */
void System::project_vector (Number fptr(const Point& p,
                                         const Parameters& parameters,
                                         const std::string& sys_name,
                                         const std::string& unknown_name),
                             Gradient gptr(const Point& p,
                                           const Parameters& parameters,
                                           const std::string& sys_name,
                                           const std::string& unknown_name),
                             Parameters& parameters,
                             NumericVector<Number>& new_vector) const
{
  START_LOG ("project_vector()", "System");

  // We need data to project
  assert(fptr);

  /**
   * This method projects an analytic solution to the current
   * mesh.  The input function \p fptr gives the analytic solution,
   * while the \p new_vector (which should already be correctly sized)
   * gives the solution (to be computed) on the current mesh.
   */

  // The number of variables in this system
  const unsigned int n_variables = this->n_vars();

  // The dimensionality of the current mesh
  const unsigned int dim = _mesh.mesh_dimension();

  // The DofMap for this system
  const DofMap& dof_map = this->get_dof_map();

  // The element matrix and RHS for projections.
  // Note that Ke is always real-valued, whereas
  // Fe may be complex valued if complex number
  // support is enabled
  DenseMatrix<Real> Ke;
  DenseVector<Number> Fe;
  // The new element coefficients
  DenseVector<Number> Ue;


  // Loop over all the variables in the system
  for (unsigned int var=0; var<n_variables; var++)
    {
      // Get FE objects of the appropriate type
      const FEType& fe_type = dof_map.variable_type(var);     
      AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));      

      // Prepare variables for projection
      AutoPtr<QBase> qrule     (fe_type.default_quadrature_rule(dim));
      AutoPtr<QBase> qedgerule (fe_type.default_quadrature_rule(1));
      AutoPtr<QBase> qsiderule (fe_type.default_quadrature_rule(dim-1));

      // The values of the shape functions at the quadrature
      // points
      const std::vector<std::vector<Real> >& phi = fe->get_phi();

      // The gradients of the shape functions at the quadrature
      // points on the child element.
      const std::vector<std::vector<RealGradient> > *dphi = NULL;

      const FEContinuity cont = fe->get_continuity();

      if (cont == C_ONE)
        {
          // We'll need gradient data for a C1 projection
          assert(gptr);

          const std::vector<std::vector<RealGradient> >&
            ref_dphi = fe->get_dphi();
          dphi = &ref_dphi;
        }

      // The Jacobian * quadrature weight at the quadrature points
      const std::vector<Real>& JxW =
	fe->get_JxW();
     
      // The XYZ locations of the quadrature points
      const std::vector<Point>& xyz_values =
	fe->get_xyz();

      // The global DOF indices
      std::vector<unsigned int> dof_indices;
      // Side/edge DOF indices
      std::vector<unsigned int> side_dofs;
   
      // Iterators for the active elements on local processor
      MeshBase::element_iterator       elem_it =
		      _mesh.active_local_elements_begin();
      const MeshBase::element_iterator elem_end = 
		      _mesh.active_local_elements_end();
   
      for ( ; elem_it != elem_end; ++elem_it)
	{
	  const Elem* elem = *elem_it;

	  // Update the DOF indices for this element based on
          // the current mesh
	  dof_map.dof_indices (elem, dof_indices, var);

	  // The number of DOFs on the element
	  const unsigned int n_dofs = dof_indices.size();

          // Fixed vs. free DoFs on edge/face projections
          std::vector<char> dof_is_fixed(n_dofs, false); // bools
          std::vector<int> free_dof(n_dofs, 0);

	  // The element type
	  const ElemType elem_type = elem->type();

	  // The number of nodes on the new element
	  const unsigned int n_nodes = elem->n_nodes();

          // Zero the interpolated values
          Ue.resize (n_dofs); Ue.zero();

          // In general, we need a series of
          // projections to ensure a unique and continuous
          // solution.  We start by interpolating nodes, then
          // hold those fixed and project edges, then
          // hold those fixed and project faces, then
          // hold those fixed and project interiors

          // Interpolate node values first
          unsigned int current_dof = 0;
          for (unsigned int n=0; n!= n_nodes; ++n)
            {
              // FIXME: this should go through the DofMap,
              // not duplicate dof_indices code badly!
	      const unsigned int nc =
		FEInterface::n_dofs_at_node (dim, fe_type, elem_type,
                                             n);
              if (!elem->is_vertex(n))
                {
                  current_dof += nc;
                  continue;
                }
              if (cont == DISCONTINUOUS)
                {
                  assert(nc == 0);
                }
              // Assume that C_ZERO elements have a single nodal
              // value shape function
              else if (cont == C_ZERO)
                {
                  assert(nc == 1);
		  Ue(current_dof) = fptr(elem->point(n),
                                         parameters,
                                         this->name(),
                                         this->variable_name(var));
                  dof_is_fixed[current_dof] = true;
                  current_dof++;
                }
              // The hermite element vertex shape functions are weird
              else if (fe_type.family == HERMITE)
                {
                  Ue(current_dof) = fptr(elem->point(n),
                                         parameters,
                                         this->name(),
                                         this->variable_name(var));
                  dof_is_fixed[current_dof] = true;
                  current_dof++;
                  Gradient g = gptr(elem->point(n),
                                    parameters,
                                    this->name(),
                                    this->variable_name(var));
                  // x derivative
                  Ue(current_dof) = g(0);
                  dof_is_fixed[current_dof] = true;
                  current_dof++;
                  if (dim > 1)
                    {
                      // We'll finite difference mixed derivatives
                      Point nxminus = elem->point(n),
                            nxplus = elem->point(n);
                      nxminus(0) -= TOLERANCE;
                      nxplus(0) += TOLERANCE;
                      Gradient gxminus = gptr(nxminus,
                                              parameters,
                                              this->name(),
                                              this->variable_name(var));
                      Gradient gxplus = gptr(nxminus,
                                             parameters,
                                             this->name(),
                                             this->variable_name(var));
                      // y derivative
                      Ue(current_dof) = g(1);
                      dof_is_fixed[current_dof] = true;
                      current_dof++;
                      // xy derivative
                      Ue(current_dof) = (gxplus(1) - gxminus(1))
                                        / 2. / TOLERANCE;
                      dof_is_fixed[current_dof] = true;
                      current_dof++;

                      if (dim > 2)
                        {
                          // z derivative
                          Ue(current_dof) = g(2);
                          dof_is_fixed[current_dof] = true;
                          current_dof++;
                          // xz derivative
                          Ue(current_dof) = (gxplus(2) - gxminus(2))
                                            / 2. / TOLERANCE;
                          dof_is_fixed[current_dof] = true;
                          current_dof++;
                          // We need new points for yz
                          Point nyminus = elem->point(n),
                                nyplus = elem->point(n);
                          nyminus(1) -= TOLERANCE;
                          nyplus(1) += TOLERANCE;
                          Gradient gyminus = gptr(nyminus,
                                                  parameters,
                                                  this->name(),
                                                  this->variable_name(var));
                          Gradient gyplus = gptr(nyminus,
                                                 parameters,
                                                 this->name(),
                                                 this->variable_name(var));
                          // xz derivative
                          Ue(current_dof) = (gyplus(2) - gyminus(2))
                                            / 2. / TOLERANCE;
                          dof_is_fixed[current_dof] = true;
                          current_dof++;
                          // Getting a 2nd order xyz is more tedious
                          Point nxmym = elem->point(n),
                                nxmyp = elem->point(n),
                                nxpym = elem->point(n),
                                nxpyp = elem->point(n);
                          nxmym(0) -= TOLERANCE;
                          nxmym(1) -= TOLERANCE;
                          nxmyp(0) -= TOLERANCE;
                          nxmyp(1) += TOLERANCE;
                          nxpym(0) += TOLERANCE;
                          nxpym(1) -= TOLERANCE;
                          nxpyp(0) += TOLERANCE;
                          nxpyp(1) += TOLERANCE;
                          Gradient gxmym = gptr(nxmym,
                                                parameters,
                                                this->name(),
                                                this->variable_name(var));
                          Gradient gxmyp = gptr(nxmyp,
                                                parameters,
                                                this->name(),
                                                this->variable_name(var));
                          Gradient gxpym = gptr(nxpym,
                                                parameters,
                                                this->name(),
                                                this->variable_name(var));
                          Gradient gxpyp = gptr(nxpyp,
                                                parameters,
                                                this->name(),
                                                this->variable_name(var));
                          Number gxzplus = (gxpyp(2) - gxmyp(2))
                                         / 2. / TOLERANCE;
                          Number gxzminus = (gxpym(2) - gxmym(2))
                                          / 2. / TOLERANCE;
                          // xyz derivative
                          Ue(current_dof) = (gxzplus - gxzminus)
                                            / 2. / TOLERANCE;
                          dof_is_fixed[current_dof] = true;
                          current_dof++;
                        }
                    }
                }
              // Assume that other C_ONE elements have a single nodal
              // value shape function and nodal gradient component
              // shape functions
              else if (cont == C_ONE)
                {
                  assert(nc == 1 + dim);
		  Ue(current_dof) = fptr(elem->point(n),
                                         parameters,
                                         this->name(),
                                         this->variable_name(var));
                  dof_is_fixed[current_dof] = true;
                  current_dof++;
                  Gradient g = gptr(elem->point(n),
                                    parameters,
                                    this->name(),
                                    this->variable_name(var));
                  for (unsigned int i=0; i!= dim; ++i)
                    {
		      Ue(current_dof) = g(i);
                      dof_is_fixed[current_dof] = true;
                      current_dof++;
                    }
                }
              else
                error();
            }

          // In 3D, project any edge values next
          if (dim > 2 && cont != DISCONTINUOUS)
            for (unsigned int e=0; e != elem->n_edges(); ++e)
              {
		FEInterface::dofs_on_edge(elem, dim, fe_type, e,
                                          side_dofs);

                // Some edge dofs are on nodes and already
                // fixed, others are free to calculate
                unsigned int free_dofs = 0;
                for (unsigned int i=0; i != side_dofs.size(); ++i)
                  if (!dof_is_fixed[side_dofs[i]])
                    free_dof[free_dofs++] = i;
	        Ke.resize (free_dofs, free_dofs); Ke.zero();
	        Fe.resize (free_dofs); Fe.zero();
                // The new edge coefficients
                DenseVector<Number> Uedge(free_dofs);

                // Initialize FE data on the edge
                fe->attach_quadrature_rule (qedgerule.get());
	        fe->edge_reinit (elem, e);
	        const unsigned int n_qp = qedgerule->n_points();

	        // Loop over the quadrature points
	        for (unsigned int qp=0; qp<n_qp; qp++)
	          {
	            // solution at the quadrature point	      
	            Number fineval = fptr(xyz_values[qp],
                                          parameters,
                                          this->name(),
                                          this->variable_name(var));
	            // solution grad at the quadrature point	      
	            Gradient finegrad;
                    if (cont == C_ONE)
                      finegrad = gptr(xyz_values[qp], parameters,
                                      this->name(),
				      this->variable_name(var));

                    // Form edge projection matrix
                    for (unsigned int sidei=0, freei=0; 
                         sidei != side_dofs.size(); ++sidei)
                      {
                        unsigned int i = side_dofs[sidei];
                        // fixed DoFs aren't test functions
                        if (dof_is_fixed[i])
                          continue;
			for (unsigned int sidej=0, freej=0;
                             sidej != side_dofs.size(); ++sidej)
                          {
                            unsigned int j = side_dofs[sidej];
                            if (dof_is_fixed[j])
			      Fe(freei) -= phi[i][qp] * phi[j][qp] *
                                           JxW[qp] * Ue(j);
                            else
                              Ke(freei,freej) += phi[i][qp] *
						 phi[j][qp] * JxW[qp];
                            if (cont == C_ONE)
                              {
                                if (dof_is_fixed[j])
                                  Fe(freei) -= ((*dphi)[i][qp] *
					        (*dphi)[j][qp]) *
                                                JxW[qp] * Ue(j);
                                else
                                  Ke(freei,freej) += ((*dphi)[i][qp] *
                                                      (*dphi)[j][qp])
                                                      * JxW[qp];
                              }
                            if (!dof_is_fixed[j])
                              freej++;
                          }
                        Fe(freei) += phi[i][qp] * fineval * JxW[qp];
                        if (cont == C_ONE)
			  Fe(freei) += (finegrad * (*dphi)[i][qp]) *
                                       JxW[qp];
                        freei++;
                      }
	          }

                Ke.cholesky_solve(Fe, Uedge);

                // Transfer new edge solutions to element
		for (unsigned int i=0; i != free_dofs; ++i)
                  {
                    Number &ui = Ue(side_dofs[free_dof[i]]);
                    assert(std::abs(ui) < TOLERANCE ||
                           std::abs(ui - Uedge(i)) < TOLERANCE);
                    ui = Uedge(i);
                    dof_is_fixed[side_dofs[free_dof[i]]] = true;
                  }
              }
                 
	  // Project any side values (edges in 2D, faces in 3D)
          if (dim > 1 && cont != DISCONTINUOUS)
            for (unsigned int s=0; s != elem->n_sides(); ++s)
              {
		FEInterface::dofs_on_side(elem, dim, fe_type, s,
                                          side_dofs);

		// Some side dofs are on nodes/edges and already
                // fixed, others are free to calculate
                unsigned int free_dofs = 0;
                for (unsigned int i=0; i != side_dofs.size(); ++i)
                  if (!dof_is_fixed[side_dofs[i]])
                    free_dof[free_dofs++] = i;
	        Ke.resize (free_dofs, free_dofs); Ke.zero();
	        Fe.resize (free_dofs); Fe.zero();
                // The new side coefficients
                DenseVector<Number> Uside(free_dofs);

                // Initialize FE data on the side
                fe->attach_quadrature_rule (qsiderule.get());
	        fe->reinit (elem, s);
	        const unsigned int n_qp = qsiderule->n_points();

	        // Loop over the quadrature points
	        for (unsigned int qp=0; qp<n_qp; qp++)
	          {
	            // solution at the quadrature point	      
	            Number fineval = fptr(xyz_values[qp],
                                          parameters,
                                          this->name(),
                                          this->variable_name(var));
	            // solution grad at the quadrature point	      
	            Gradient finegrad;
                    if (cont == C_ONE)
                      finegrad = gptr(xyz_values[qp], parameters,
                                      this->name(),
				      this->variable_name(var));

                    // Form side projection matrix
                    for (unsigned int sidei=0, freei=0;
                         sidei != side_dofs.size(); ++sidei)
                      {
                        unsigned int i = side_dofs[sidei];
                        // fixed DoFs aren't test functions
                        if (dof_is_fixed[i])
                          continue;
			for (unsigned int sidej=0, freej=0;
                             sidej != side_dofs.size(); ++sidej)
                          {
                            unsigned int j = side_dofs[sidej];
                            if (dof_is_fixed[j])
			      Fe(freei) -= phi[i][qp] * phi[j][qp] *
                                           JxW[qp] * Ue(j);
                            else
			      Ke(freei,freej) += phi[i][qp] *
						 phi[j][qp] * JxW[qp];
                            if (cont == C_ONE)
                              {
                                if (dof_is_fixed[j])
                                  Fe(freei) -= ((*dphi)[i][qp] *
					        (*dphi)[j][qp]) *
                                               JxW[qp] * Ue(j);
                                else
                                  Ke(freei,freej) += ((*dphi)[i][qp] *
						      (*dphi)[j][qp])
                                                     * JxW[qp];
                              }
                            if (!dof_is_fixed[j])
                              freej++;
                          }
                        Fe(freei) += (fineval * phi[i][qp]) * JxW[qp];
                        if (cont == C_ONE)
			  Fe(freei) += (finegrad * (*dphi)[i][qp]) *
                                       JxW[qp];
                        freei++;
                      }
	          }

                Ke.cholesky_solve(Fe, Uside);

                // Transfer new side solutions to element
		for (unsigned int i=0; i != free_dofs; ++i)
                  {
                    Number &ui = Ue(side_dofs[free_dof[i]]);
                    assert(std::abs(ui) < TOLERANCE ||
                           std::abs(ui - Uside(i)) < TOLERANCE);
                    ui = Uside(i);
                    dof_is_fixed[side_dofs[free_dof[i]]] = true;
                  }
              }

	  // Project the interior values, finally

	  // Some interior dofs are on nodes/edges/sides and
          // already fixed, others are free to calculate
          unsigned int free_dofs = 0;
          for (unsigned int i=0; i != n_dofs; ++i)
            if (!dof_is_fixed[i])
              free_dof[free_dofs++] = i;
	  Ke.resize (free_dofs, free_dofs); Ke.zero();
	  Fe.resize (free_dofs); Fe.zero();
          // The new interior coefficients
          DenseVector<Number> Uint(free_dofs);

          // Initialize FE data
          fe->attach_quadrature_rule (qrule.get());
	  fe->reinit (elem);
	  const unsigned int n_qp = qrule->n_points();

	  // Loop over the quadrature points
	  for (unsigned int qp=0; qp<n_qp; qp++)
	    {
	      // solution at the quadrature point	      
	      Number fineval = fptr(xyz_values[qp],
                                    parameters,
                                    this->name(),
                                    this->variable_name(var));
	      // solution grad at the quadrature point	      
	      Gradient finegrad;
              if (cont == C_ONE)
                finegrad = gptr(xyz_values[qp], parameters,
                                this->name(),
				this->variable_name(var));

              // Form interior projection matrix
              for (unsigned int i=0, freei=0; i != n_dofs; ++i)
                {
                  // fixed DoFs aren't test functions
                  if (dof_is_fixed[i])
                    continue;
		  for (unsigned int j=0, freej=0; j != n_dofs; ++j)
                    {
                      if (dof_is_fixed[j])
			Fe(freei) -= phi[i][qp] * phi[j][qp] * JxW[qp]
                                     * Ue(j);
                      else
			Ke(freei,freej) += phi[i][qp] * phi[j][qp] *
                                           JxW[qp];
                      if (cont == C_ONE)
                        {
                          if (dof_is_fixed[j])
			    Fe(freei) -= ((*dphi)[i][qp] *
					 (*dphi)[j][qp]) * JxW[qp] *
                                         Ue(j);
                          else
			    Ke(freei,freej) += ((*dphi)[i][qp] *
						(*dphi)[j][qp]) *
                                               JxW[qp];
                        }
                      if (!dof_is_fixed[j])
                        freej++;
                    }
		  Fe(freei) += phi[i][qp] * fineval * JxW[qp];
                  if (cont == C_ONE)
                    Fe(freei) += (finegrad * (*dphi)[i][qp]) * JxW[qp];
                  freei++;
                }
	    }
          Ke.cholesky_solve(Fe, Uint);

          // Transfer new interior solutions to element
	  for (unsigned int i=0; i != free_dofs; ++i)
            {
              Number &ui = Ue(free_dof[i]);
              assert(std::abs(ui) < TOLERANCE ||
                     std::abs(ui - Uint(i)) < TOLERANCE);
              ui = Uint(i);
              dof_is_fixed[free_dof[i]] = true;
            }

          // Make sure every DoF got reached!
	  for (unsigned int i=0; i != n_dofs; ++i)
            assert(dof_is_fixed[i]);

	  const unsigned int
	    first = new_vector.first_local_index(),
	    last  = new_vector.last_local_index();
	  
          for (unsigned int i = 0; i < n_dofs; i++) 
	    if (Ue(i) != 0.)
	      if ((dof_indices[i] >= first) &&
		  (dof_indices[i] <  last)) new_vector.set(dof_indices[i], Ue(i));
        }  // end elem loop
    } // end variables loop

  new_vector.close();

  dof_map.enforce_constraints_exactly(*this, &new_vector);

  STOP_LOG("project_vector()", "System");
}
