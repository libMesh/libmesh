// $Id: fe_boundary.C,v 1.15 2003-05-15 23:34:35 benkirk Exp $

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
#include <math.h>


// Local includes
#include "mesh_common.h"
#include "fe.h"
#include "quadrature.h"
#include "elem.h"
#include "fe_macro.h"



//-------------------------------------------------------
// Full specialization for 1D when this is a useless method
template <>
void FE<1,HIERARCHIC>::reinit(const Elem*,
			      const unsigned int)
{
  std::cerr << "ERROR: This method only makes sense for 2D, 3D elements!"
	    << std::endl;
  error();
}



//-------------------------------------------------------
// Full specialization for 1D when this is a useless method
template <>
void FE<1,LAGRANGE>::reinit(const Elem*,
			    const unsigned int)
{
  std::cerr << "ERROR: This method only makes sense for 2D, 3D elements!"
	    << std::endl;
  error();
}



//-------------------------------------------------------
// Full specialization for 1D when this is a useless method
template <>
void FE<1,MONOMIAL>::reinit(const Elem*,
			    const unsigned int)
{
  std::cerr << "ERROR: This method only makes sense for 2D, 3D elements!"
	    << std::endl;
  error();
}



//-------------------------------------------------------
// Method for 2D, 3D
template <unsigned int Dim, FEFamily T>
void FE<Dim,T>::reinit(const Elem* elem,
		       const unsigned int s)
{
  assert (elem  != NULL);
  assert (qrule != NULL);
  // We don't do this for 1D elements!
  assert (Dim != 1);
  
  const AutoPtr<Elem> side(elem->build_side(s));
      
  qrule->init(side->type());
  
  
  // Compute the Jacobian*Weight on the face for integration
  {
    elem_type = side->type();
    init_shape_functions (qrule->get_points(),  elem, s);
    compute_map          (qrule->get_weights(), elem, s);
  }

  // make a copy of the Jacobian for integration
  const std::vector<Real>  JxW_int(JxW);


  // Find where the integration points are located on the
  // full element.
  std::vector<Point> qp;
  
  inverse_map (elem, xyz, qp);
  

  // compute the shape function and derivative values
  {
    elem_type = elem->type();
    init_shape_functions    (qp, elem);
    compute_map             (qrule->get_weights(), elem);
    compute_shape_functions ();
  }  
  
  // copy back old data
  JxW = JxW_int;
}



template <unsigned int Dim, FEFamily T>
void FE<Dim,T>::init_shape_functions(const std::vector<Point>& qp,
				     const Elem* elem,
				     const unsigned int s)
{
  assert (elem  != NULL);
  
  const AutoPtr<Elem> side(elem->build_side(s));
  
  // The element type and order to use in
  // the map
  const Order    mapping_order     (side->default_order()); 
  const ElemType mapping_elem_type (side->type());

  // The number of quadrature points.
  const unsigned int n_qp = qp.size();
	
  const unsigned int n_mapping_shape_functions =
    n_shape_functions (mapping_elem_type,
		       mapping_order);
  
  // resize the vectors to hold current data
  // Psi are the shape functions used for the FE mapping
  {
    psi_map.resize        (n_mapping_shape_functions);
    dpsidxi_map.resize    (n_mapping_shape_functions);
    if (Dim == 3)
      dpsideta_map.resize (n_mapping_shape_functions);
    
    for (unsigned int i=0; i<n_mapping_shape_functions; i++)
      {
	psi_map[i].resize        (n_qp);
	dpsidxi_map[i].resize    (n_qp);
	if (Dim == 3)
	  dpsideta_map[i].resize (n_qp);
      }
  }
  
  
  // Compute the value of shape function i at quadrature point p
  // (Lagrange shape functions are used for the mapping)
  for (unsigned int i=0; i<n_mapping_shape_functions; i++)
    for (unsigned int p=0; p<n_qp; p++)
      {
	psi_map[i][p]        = FE<Dim-1,LAGRANGE>::shape       (mapping_elem_type, mapping_order, i,    qp[p]);
	dpsidxi_map[i][p]    = FE<Dim-1,LAGRANGE>::shape_deriv (mapping_elem_type, mapping_order, i, 0, qp[p]);
	if (Dim == 3)
	  dpsideta_map[i][p] = FE<Dim-1,LAGRANGE>::shape_deriv (mapping_elem_type, mapping_order, i, 1, qp[p]);
      }

  // compute the normal vectors at each quadrature point p
  tangents.resize(n_qp);
  normals.resize(n_qp);
  for (unsigned int p=0; p<n_qp; p++)
    {
      tangents[p].resize(2);
      
      if (Dim == 2)
	{
	  // compute d(xyz)/dxi at the
	  // quadrature point.
	  Point dxyzdxi_map;

	  for (unsigned int i=0; i<n_mapping_shape_functions; i++)
	    dxyzdxi_map.add_scaled(side->point(i), dpsidxi_map[i][p]);
	  
	  const Point n(dxyzdxi_map(1), -dxyzdxi_map(0), 0.);
	  
	  normals[p]     = n.unit();
	  tangents[p][0] = dxyzdxi_map.unit();
	  tangents[p][1] = dxyzdxi_map.cross(n).unit();
	}
	
      else if (Dim == 3)
	{
	  // compute d(xyz)/dxi and d(xyz)/deta at the
	  // quadrature point.
	  Point dxyzdxi_map;
	  Point dxyzdeta_map;

	  for (unsigned int i=0; i<n_mapping_shape_functions; i++)
	    {
	      dxyzdxi_map.add_scaled (side->point(i), dpsidxi_map[i][p]);
	      dxyzdeta_map.add_scaled(side->point(i), dpsideta_map[i][p]);
	    }
	  
	  const Point n  = dxyzdxi_map.cross(dxyzdeta_map);
	  normals[p]     = n.unit();
	  tangents[p][0] = dxyzdxi_map.unit();
	  tangents[p][1] = n.cross(dxyzdxi_map).unit();
	}
    }
}

  

void FEBase::compute_map(const std::vector<Real>& qw,
			 const Elem* elem,
			 const unsigned int s)
{
  assert (elem  != NULL);

  // Get an AutoPtr to the side.  Don't worry about deleting it.
  const AutoPtr<Elem> side(elem->build_side(s));       
  const unsigned int n_qp = qw.size();

  switch (dim)
    {
      
    case 2:
      {
	//-------------------------------------------------------------
	// Resize the vectors to hold data at the quadrature points
	{  
	  xyz.resize(n_qp);
	  dxyzdxi_map.resize(n_qp);
	  
	  JxW.resize(n_qp);
	}
	
	// Clear the entities that will be summed
	for (unsigned int p=0; p<n_qp; p++)
	  {
	    xyz[p].zero();
	    dxyzdxi_map[p].zero();
	  }
	
	// compute x, dxdxi at the quadrature points    
	for (unsigned int i=0; i<psi_map.size(); i++) // sum over the nodes
	  {
	    const Point& side_point = side->point(i);
	    
	    for (unsigned int p=0; p<n_qp; p++) // for each quadrature point...
	      {	  
		xyz[p].add_scaled        (side_point, psi_map[i][p]);
		dxyzdxi_map[p].add_scaled(side_point, dpsidxi_map[i][p]);
	      }
	  }
	  
	
	// compute the jacobian at the quadrature points
	for (unsigned int p=0; p<n_qp; p++)
	  {
	    const Real jac = sqrt(dxdxi_map(p)*dxdxi_map(p) +
				  dydxi_map(p)*dydxi_map(p));
	    
	    assert (jac > 0.);
	    
	    JxW[p] = jac*qw[p];
	  }
	
	// done computing the map	
	return;
      }


    case 3:
      {
	//----------------------------------------------------------------
	// Resize the vectors to hold data at the quadrature points
	{  
	  xyz.resize(n_qp);
	  dxyzdxi_map.resize(n_qp);
	  dxyzdeta_map.resize(n_qp);
	  
	  JxW.resize(n_qp);
	}
    
	// Clear the entities that will be summed
	for (unsigned int p=0; p<n_qp; p++)
	  {
	    xyz[p].zero();
	    dxyzdxi_map[p].zero();
	    dxyzdeta_map[p].zero();
	  }
	
	// compute x, dxdxi at the quadrature points    
	for (unsigned int i=0; i<psi_map.size(); i++) // sum over the nodes
	  {
	    const Point& side_point = side->point(i);
	    
	    for (unsigned int p=0; p<n_qp; p++) // for each quadrature point...
	      {
		xyz[p].add_scaled         (side_point, psi_map[i][p]);
		dxyzdxi_map[p].add_scaled (side_point, dpsidxi_map[i][p]);
		dxyzdeta_map[p].add_scaled(side_point, dpsideta_map[i][p]);
	      }
	  }
    
    
	
	// compute the jacobian at the quadrature points, see
	// http://sp81.msi.umn.edu:999/fluent/fidap/help/theory/thtoc.htm
	for (unsigned int p=0; p<n_qp; p++)
	  {
	    const Real g11 = (dxdxi_map(p)*dxdxi_map(p) +
			      dydxi_map(p)*dydxi_map(p) +
			      dzdxi_map(p)*dzdxi_map(p));
	    
	    const Real g12 = (dxdxi_map(p)*dxdeta_map(p) +
			      dydxi_map(p)*dydeta_map(p) +
			      dzdxi_map(p)*dzdeta_map(p));
	    
	    const Real g21 = g12;
	    
	    const Real g22 = (dxdeta_map(p)*dxdeta_map(p) +
			      dydeta_map(p)*dydeta_map(p) +
			      dzdeta_map(p)*dzdeta_map(p));
	    
	    
	    const Real jac = sqrt(g11*g22 - g12*g21);
	    
	    assert (jac > 0.);

	    JxW[p] = jac*qw[p];
	  }
	
	// done computing the map
	return;  
      }


    default:
      error();
      
    }
}




//--------------------------------------------------------------
// Explicit instantiations (doesn't make sense in 1D!) using fe_macro.h's macro

INSTANTIATE_FE(2);

INSTANTIATE_FE(3);

