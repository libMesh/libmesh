// $Id: fe_clough.C,v 1.6 2005-02-28 19:05:58 roystgnr Exp $

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



// Local includes
#include "dense_matrix.h"
#include "dense_vector.h"
#include "dof_map.h"
#include "elem.h"
#include "fe.h"
#include "fe_interface.h"
#include "quadrature_clough.h"




// ------------------------------------------------------------
// Hierarchic-specific implementations
template <unsigned int Dim, FEFamily T>
void FE<Dim,T>::nodal_soln(const Elem* elem,
			   const Order order,
			   const std::vector<Number>& elem_soln,
			   std::vector<Number>&       nodal_soln)
{
  const unsigned int n_nodes = elem->n_nodes();
  
  const ElemType type = elem->type();

  nodal_soln.resize(n_nodes);


  
  switch (order)
    {
      // Piecewise cubic shape functions only
    case THIRD:
      {

	const unsigned int n_sf =
	  FE<Dim,T>::n_shape_functions(type, order);
	
	for (unsigned int n=0; n<n_nodes; n++)
	  {
	    const Point mapped_point = FE<Dim,T>::inverse_map(elem,
							      elem->point(n));

	    assert (elem_soln.size() == n_sf);

	    // Zero before summation
	    nodal_soln[n] = 0;

	    // u_i = Sum (alpha_i phi_i)
	    for (unsigned int i=0; i<n_sf; i++)
	      nodal_soln[n] += elem_soln[i]*FE<Dim,T>::shape(elem,
							     order,
							     i,
							     mapped_point);	    
	  }

	return;
      }
      
    default:
      {
	error();
      }
    }
}



template <unsigned int Dim, FEFamily T>
unsigned int FE<Dim,T>::n_dofs(const ElemType t, const Order o)
{
  switch (o)
    {
      // Piecewise cubic Clough-Tocher element
    case THIRD:
      {
	switch (t)
	  {
	  case TRI6:
	    return 12;
	    
	  default:
	    {
#ifdef DEBUG
	      std::cerr << "ERROR: Bad ElemType = " << t
			<< " for " << o << "th order approximation!" 
			<< std::endl;
#endif
	      error();	    
	    }
	  }
      }
      
    default:
      {
	error();
      }
    }
  
  error();  
  return 0;
}



template <unsigned int Dim, FEFamily T>
unsigned int FE<Dim,T>::n_dofs_at_node(const ElemType t,
				       const Order o,
				       const unsigned int n)
{
  switch (o)
    {
      // The third-order hierarchic shape functions
    case THIRD:
      {
	switch (t)
	  {
	    // The 2D Clough-Tocher defined on a 6-noded triangle
	  case TRI6:
	    {
	      switch (n)
		{
		case 0:
		case 1:
		case 2:
		  return 3;

		case 3:
		case 4:
		case 5:
		  return 1;

		default:
		  error();
		}
	    }

	  default:
	    {
#ifdef DEBUG
	      std::cerr << "ERROR: Bad ElemType = " << t
			<< " for " << o << "th order approximation!" 
			<< std::endl;
#endif
	      error();	    
	    }
	    
	  }
      }
    default:
      {
	error();
      }
    }
  
  error();
  
  return 0;
}



template <unsigned int Dim, FEFamily T>
unsigned int FE<Dim,T>::n_dofs_per_elem(const ElemType t,
					const Order o)
{
  switch (o)
    {
      // The third-order Clough-Tocher shape functions
    case THIRD:
      {
	switch (t)
	  {
	    // The 2D hierarchic defined on a 6-noded triangle
	  case TRI6:
	    return 0;

	  default:
	    {
#ifdef DEBUG
	      std::cerr << "ERROR: Bad ElemType = " << t
			<< " for " << o << "th order approximation!" 
			<< std::endl;
#endif
	      error();	    
	    }
	    
	  }
      }
      // Otherwise no DOFS per element
    default:
      error();	    
      return 0;
    }
}



template <unsigned int Dim, FEFamily T>
void FE<Dim,T>::compute_constraints (std::map<unsigned int,
				     std::map<unsigned int, float> > &
				     constraints,
				     DofMap &dof_map,
				     const unsigned int variable_number,
				     const Elem* elem)
{
#ifdef ENABLE_AMR
  // Only constrain elements in 2,3D.
  if (Dim == 1)
    return;

  assert (elem != NULL);

  const FEType& fe_type = dof_map.variable_type(variable_number);

  AutoPtr<FEBase> my_fe (FEBase::build(Dim, fe_type));
  AutoPtr<FEBase> parent_fe (FEBase::build(Dim, fe_type));

  QClough my_qface(Dim-1, fe_type.default_quadrature_order());
  my_fe->attach_quadrature_rule (&my_qface);
  std::vector<Point> parent_qface;

  const std::vector<Real>& JxW = my_fe->get_JxW();
  const std::vector<Point>& q_point = my_fe->get_xyz();
  const std::vector<std::vector<Real> >& phi = my_fe->get_phi();
  const std::vector<std::vector<Real> >& parent_phi =
		  parent_fe->get_phi();
  const std::vector<Point>& face_normals = my_fe->get_normals();
  const std::vector<std::vector<RealGradient> >& dphi =
		  my_fe->get_dphi();
  const std::vector<std::vector<RealGradient> >& parent_dphi =
		  parent_fe->get_dphi();
  std::vector<unsigned int> child_dof_indices, parent_dof_indices;
  std::vector<unsigned int> my_side_dofs, parent_side_dofs;

  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;
  std::vector<DenseVector<Number> > Ue;

  // Look at the element faces.  Check to see if we need to
  // build constraints.
  for (unsigned int s=0; s<elem->n_sides(); s++)
    if (elem->neighbor(s) != NULL)
      // constrain dofs shared between
      // this element and ones coarser
      // than this element.
      if (elem->neighbor(s)->level() < elem->level()) 
        {
          // Get pointers to the elements of interest and its parent.
          const Elem* parent = elem->parent();
	  unsigned int s_parent = s;

          // This can't happen...  Only level-0 elements have NULL
          // parents, and no level-0 elements can be at a higher
          // level than their neighbors!
          assert (parent != NULL);

	  my_fe->reinit(elem, s);

	  dof_map.dof_indices (elem, child_dof_indices,
			       variable_number);
	  dof_map.dof_indices (parent, parent_dof_indices,
			       variable_number);

	  const unsigned int n_qp = my_qface.n_points();

	  FEInterface::inverse_map (Dim, fe_type, parent, q_point,
				    parent_qface);

	  parent_fe->reinit(parent, &parent_qface);

	  // We're only concerned with DOFs whose values and/or first
	  // derivatives are supported on side nodes
	  dofs_on_side(elem, fe_type.order, s, my_side_dofs);
	  dofs_on_side(parent, fe_type.order, s_parent, parent_side_dofs);
	  const unsigned int n_side_dofs = my_side_dofs.size();
	  assert(n_side_dofs == parent_side_dofs.size());

	  Ke.resize (n_side_dofs, n_side_dofs);
	  Ue.resize(n_side_dofs);

	  // Form the projection matrix, (inner product of fine basis
	  // functions against fine test functions)
	  for (unsigned int is = 0; is != n_side_dofs; ++is)
	    {
	      const unsigned int i = my_side_dofs[is];
	      for (unsigned int js = 0; js != n_side_dofs; ++js)
	        {
	          const unsigned int j = my_side_dofs[js];
		  for (unsigned int qp = 0; qp != n_qp; ++qp)
		    Ke(is,js) += JxW[qp] * (phi[i][qp] * phi[j][qp] +
					    (dphi[i][qp] *
					     face_normals[qp]) *
					    (dphi[j][qp] *
					     face_normals[qp]));
		}
	    }

	  // Form the right hand sides, (inner product of coarse basis
	  // functions against fine test functions)
	  for (unsigned int is = 0; is != n_side_dofs; ++is)
	    {
	      const unsigned int i = parent_side_dofs[is];
	      Fe.resize (n_side_dofs);
	      for (unsigned int js = 0; js != n_side_dofs; ++js)
		{
	          const unsigned int j = my_side_dofs[js];
	          for (unsigned int qp = 0; qp != n_qp; ++qp)
		    Fe(js) += JxW[qp] * (parent_phi[i][qp] *
					 phi[j][qp] +
					 (parent_dphi[i][qp] *
					  face_normals[qp]) *
					 (dphi[j][qp] *
					  face_normals[qp]));
		}
	      Ke.cholesky_solve(Fe, Ue[is]);
	    }
	  for (unsigned int is = 0; is != n_side_dofs; ++is)
	    {
	      const unsigned int i = parent_side_dofs[is];
	      const unsigned int their_dof_g = parent_dof_indices[i];
	      for (unsigned int js = 0; js != n_side_dofs; ++js)
	        {
	          const unsigned int j = my_side_dofs[js];
	          const unsigned int my_dof_g = child_dof_indices[j];
		  const Real their_dof_value = Ue[is](js);
		  if (their_dof_g == my_dof_g)
		    {
		      assert(std::abs(their_dof_value-1.) < 1.e-5);
		      for (unsigned int k = 0; k != n_side_dofs; ++k)
		        assert(k == is || std::abs(Ue[k](is)) < 1.e-5);
		      continue;
		    }
		  if (std::abs(their_dof_value) < 1.e-5)
		    continue;

		  DofMap::DofConstraintRow& constraint_row =
				  constraints[my_dof_g];

		  constraint_row.insert(std::make_pair(their_dof_g,
						       static_cast<float>(their_dof_value)));
	        }
	    }
	}
#endif
}



template <unsigned int Dim, FEFamily T>
bool FE<Dim,T>::shapes_need_reinit() const
{
  return true;
}


//--------------------------------------------------------------
// Explicit instantiations
template class FE<1,CLOUGH>;  // FIXME: 1D Not yet functional!
template class FE<2,CLOUGH>;
template class FE<3,CLOUGH>;  // FIXME: 2D Not yet functional!
