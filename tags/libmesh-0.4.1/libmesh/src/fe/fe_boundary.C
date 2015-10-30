// $Id: fe_boundary.C,v 1.22 2003-11-05 22:26:45 benkirk Exp $

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
#include <math.h>


// Local includes
#include "libmesh_common.h"
#include "fe.h"
#include "quadrature.h"
#include "elem.h"
#include "fe_macro.h"
#include "libmesh_logging.h"




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

  // Build the side of interest 
  const AutoPtr<Elem> side(elem->build_side(s));

  // initialize quadrature rule
  qrule->init(side->type());

  // The last side we computed the shape functions for
  static unsigned int last_s = libMesh::invalid_uint;

  // We might not need to reinitialize the shape functions
  if ((this->get_type() != elem->type()) ||
      (s != last_s) ||
      this->shapes_need_reinit())
    {
      // Set the element type
      elem_type = elem->type();

      // Set the last_s
      last_s = s;
      
      // Initialize the face shape functions
      this->init_face_shape_functions (qrule->get_points(),  side.get());
    }
  
  // Compute the Jacobian*Weight on the face for integration
  this->compute_face_map (qrule->get_weights(), side.get());

  // make a copy of the Jacobian for integration
  const std::vector<Real> JxW_int(JxW);

  // Find where the integration points are located on the
  // full element.
  std::vector<Point> qp; this->inverse_map (elem, xyz, qp);
  
  // compute the shape function and derivative values
  // at the points qp
  this->reinit  (elem, &qp);
      
  // copy back old data
  JxW = JxW_int;
}



template <unsigned int Dim, FEFamily T>
void FE<Dim,T>::init_face_shape_functions(const std::vector<Point>& qp,
					  const Elem* side)
{
  assert (side  != NULL);
  
  /**
   * Start logging the shape function initialization
   */
  START_LOG("init_face_shape_functions()", "FE");

  // The element type and order to use in
  // the map
  const Order    mapping_order     (side->default_order()); 
  const ElemType mapping_elem_type (side->type());

  // The number of quadrature points.
  const unsigned int n_qp = qp.size();
	
  const unsigned int n_mapping_shape_functions =
    FE<Dim,LAGRANGE>::n_shape_functions (mapping_elem_type,
					 mapping_order);
  
  // resize the vectors to hold current data
  // Psi are the shape functions used for the FE mapping
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
      
  
      // Compute the value of shape function i at quadrature point p
      // (Lagrange shape functions are used for the mapping)
      for (unsigned int p=0; p<n_qp; p++)
	{
	  psi_map[i][p]        = FE<Dim-1,LAGRANGE>::shape       (mapping_elem_type, mapping_order, i,    qp[p]);
	  dpsidxi_map[i][p]    = FE<Dim-1,LAGRANGE>::shape_deriv (mapping_elem_type, mapping_order, i, 0, qp[p]);
	  if (Dim == 3)
	    dpsideta_map[i][p] = FE<Dim-1,LAGRANGE>::shape_deriv (mapping_elem_type, mapping_order, i, 1, qp[p]);
	}
    }

  
  /**
   * Stop logging the shape function initialization
   */
  STOP_LOG("init_face_shape_functions()", "FE");
}

  

void FEBase::compute_face_map(const std::vector<Real>& qw,
			      const Elem* side)
{
  assert (side  != NULL);

  START_LOG("compute_face_map()", "FE");

  // The number of quadrature points.
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
	  tangents.resize(n_qp);
	  normals.resize(n_qp);
	  
	  JxW.resize(n_qp);
	}
	
	// Clear the entities that will be summed
	// Compute the tangent & normal at the quadrature point
	for (unsigned int p=0; p<n_qp; p++)
	  {
	    tangents[p].resize(DIM-1); // 1 Tangent in 2D, 2 in 3D
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

	// Compute the tangent & normal at the quadrature point
	for (unsigned int p=0; p<n_qp; p++)
	  {
	    const Point n(dxyzdxi_map[p](1), -dxyzdxi_map[p](0), 0.);
	    
	    normals[p]     = n.unit();
	    tangents[p][0] = dxyzdxi_map[p].unit();
#if DIM == 3  // Only good in 3D space
	    tangents[p][1] = dxyzdxi_map[p].cross(n).unit();
#endif
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
	break;
      }


      
    case 3:
      {
	//----------------------------------------------------------------
	// Resize the vectors to hold data at the quadrature points
	{  
	  xyz.resize(n_qp);
	  dxyzdxi_map.resize(n_qp);
	  dxyzdeta_map.resize(n_qp);
	  tangents.resize(n_qp);
	  normals.resize(n_qp);

	  JxW.resize(n_qp);
	}
    
	// Clear the entities that will be summed
	for (unsigned int p=0; p<n_qp; p++)
	  {
	    tangents[p].resize(DIM-1); // 1 Tangent in 2D, 2 in 3D
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

	// Compute the tangents & normal at the quadrature point
	for (unsigned int p=0; p<n_qp; p++)
	  {	    
	    const Point n  = dxyzdxi_map[p].cross(dxyzdeta_map[p]);
	    normals[p]     = n.unit();
	    tangents[p][0] = dxyzdxi_map[p].unit();
	    tangents[p][1] = n.cross(dxyzdxi_map[p]).unit();
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
	break;
      }


    default:
      error();
      
    }
  STOP_LOG("compute_face_map()", "FE");
}




//--------------------------------------------------------------
// Explicit instantiations (doesn't make sense in 1D!) using fe_macro.h's macro

INSTANTIATE_FE(2);

INSTANTIATE_FE(3);

