// $Id: fe_boundary.C,v 1.3 2003-01-20 17:06:19 jwpeterson Exp $

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



//-------------------------------------------------------
// Full specialization for 1D when this is a useless method
template <>
void FE<1,HIERARCHIC>::reinit(QBase*,
			      const Elem*,
			      const unsigned int)
{
  std::cerr << "ERROR: This method only makes sense for 2D, 3D elements!"
	    << std::endl;
  error();
};



//-------------------------------------------------------
// Full specialization for 1D when this is a useless method
template <>
void FE<1,LAGRANGE>::reinit(QBase*,
			    const Elem*,
			    const unsigned int)
{
  std::cerr << "ERROR: This method only makes sense for 2D, 3D elements!"
	    << std::endl;
  error();
};



//-------------------------------------------------------
// Full specialization for 1D when this is a useless method
template <>
void FE<1,MONOMIAL>::reinit(QBase*,
			    const Elem*,
			    const unsigned int)
{
  std::cerr << "ERROR: This method only makes sense for 2D, 3D elements!"
	    << std::endl;
  error();
};



//-------------------------------------------------------
// Method for 2D, 3D
template <unsigned int Dim, FEFamily T>
void FE<Dim,T>::reinit(QBase* qside,
		       const Elem* elem,
		       const unsigned int s)
{
  assert (elem != NULL);
  assert (qside != NULL);
  // We don't do this for 1D elements!
  assert (Dim != 1);
  
  const AutoPtr<Elem> side(elem->build_side(s));
  
  assert (qside   != NULL);
  assert (qrule   != NULL);
  
  
  qrule->init(elem->type(), s);
  qside->init(side->type());

  assert (qrule->n_points() == qside->n_points());
  
  
  // Compute the Jacobian on the face for integration
  {
    elem_type = side->type();
    init_shape_functions (qside, elem, s);
    compute_map          (qside, elem, s);
  };
  
  // make a copy of the Jacobian for integration
  const std::vector<real>  JxW_int(JxW);
  const std::vector<Point> xyz_int(xyz);

  // compute the shape function and derivative values
  {
    elem_type = elem->type();
    init_shape_functions    (qrule, elem);
    compute_map             (qrule, elem);
    compute_shape_functions (qrule);
  };  
  
  // copy back old data
  {
    JxW = JxW_int;
    xyz = xyz_int;
  };
  
  return;
};



template <unsigned int Dim, FEFamily T>
void FE<Dim,T>::init_shape_functions(const QBase* qrule,
				     const Elem* elem,
				     const unsigned int s)
{
  assert (qrule != NULL);
  assert (elem  != NULL);
  
  const AutoPtr<Elem> side(elem->build_side(s));
  
  // The element type and order to use in
  // the map
  const Order    mapping_order     (side->default_order()); 
  const ElemType mapping_elem_type (side->type());
  
  const unsigned int      n_qp = qrule->n_points();
  const std::vector<Point>& qp = qrule->get_points();
	
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
      };
  };
  
  
  // Compute the value of shape function i at quadrature point p
  // (Lagrange shape functions are used for the mapping)
  for (unsigned int i=0; i<n_mapping_shape_functions; i++)
    for (unsigned int p=0; p<n_qp; p++)
      {
	psi_map[i][p]        = FE<Dim-1,LAGRANGE>::shape       (mapping_elem_type, mapping_order, i,    qp[p]);
	dpsidxi_map[i][p]    = FE<Dim-1,LAGRANGE>::shape_deriv (mapping_elem_type, mapping_order, i, 0, qp[p]);
	if (Dim == 3)
	  dpsideta_map[i][p] = FE<Dim-1,LAGRANGE>::shape_deriv (mapping_elem_type, mapping_order, i, 1, qp[p]);
      };

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
	    dxyzdxi_map  += side->point(i)*dpsidxi_map [i][p];
	  
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
	      dxyzdxi_map  += side->point(i)*dpsidxi_map [i][p];
	      dxyzdeta_map += side->point(i)*dpsideta_map[i][p];
	    };
	  
	  const Point n  = dxyzdxi_map.cross(dxyzdeta_map);
	  normals[p]     = n.unit();
	  tangents[p][0] = dxyzdxi_map.unit();
	  tangents[p][1] = n.cross(dxyzdxi_map).unit();
	};
    };
  
  return;
};

  

void FEBase::compute_map(const QBase* qrule,
			 const Elem* elem,
			 const unsigned int s)
{
  switch (dim)
    {
      
    case 2:
      {
	assert (qrule != NULL);
	assert (elem  != NULL);

	const AutoPtr<Elem> side(elem->build_side(s));
	
	const unsigned int      n_qp = qrule->n_points();
	const std::vector<real> & qw = qrule->get_weights();
	
	
	{
	  //-------------------------------------------------------------
	  // Resize the vectors to hold data at the quadrature points
	  {  
	    xyz.resize(n_qp);
	    dxyzdxi_map.resize(n_qp);
	    
	    jac.resize(n_qp);
	    JxW.resize(n_qp);
	  }
	  
	  // Clear the entities that will be summed
	  for (unsigned int p=0; p<n_qp; p++)
	    {
	      xyz[p].clear();
	      dxyzdxi_map[p].clear();
	    };
	  
	  // compute x, dxdxi at the quadrature points    
	  for (unsigned int i=0; i<psi_map.size(); i++) // sum over the nodes
	    for (unsigned int p=0; p<n_qp; p++) // for each quadrature point...
	      {	  
		xyz[p]         += side->point(i)*psi_map[i][p];
		dxyzdxi_map[p] += side->point(i)*dpsidxi_map[i][p];
	      };
	  
	  
	  // compute the jacobian at the quadrature points
	  for (unsigned int p=0; p<n_qp; p++)
	    {
	      
	      jac[p] = sqrt(dxdxi_map(p)*dxdxi_map(p) +
			    dydxi_map(p)*dydxi_map(p));
	      
	      assert (jac[p] > 0.);
	      
	      JxW[p] = jac[p]*qw[p];
	    };
	};
	// done computing the map
	
	return;
      }


    case 3:
      {
	assert (qrule != NULL);
	assert (elem  != NULL);

	const AutoPtr<Elem> side(elem->build_side(s));
  
	const unsigned int      n_qp = qrule->n_points();
	const std::vector<real> & qw = qrule->get_weights();


	{
	  //----------------------------------------------------------------
	  // Resize the vectors to hold data at the quadrature points
	  {  
	    xyz.resize(n_qp);
	    dxyzdxi_map.resize(n_qp);
	    dxyzdeta_map.resize(n_qp);
      
	    jac.resize(n_qp);
	    JxW.resize(n_qp);
	  }
    
	  // Clear the entities that will be summed
	  for (unsigned int p=0; p<n_qp; p++)
	    {
	      xyz[p].clear();
	      dxyzdxi_map[p].clear();
	      dxyzdeta_map[p].clear();
	    };
    
	  // compute x, dxdxi at the quadrature points    
	  for (unsigned int i=0; i<psi_map.size(); i++) // sum over the nodes
	    for (unsigned int p=0; p<n_qp; p++) // for each quadrature point...
	      {
		xyz[p]          += side->point(i)*psi_map[i][p];
		dxyzdxi_map[p]  += side->point(i)*dpsidxi_map[i][p];
		dxyzdeta_map[p] += side->point(i)*dpsideta_map[i][p];
	      }
    
    

	  // compute the jacobian at the quadrature points, see
	  // http://sp81.msi.umn.edu:999/fluent/fidap/help/theory/thtoc.htm
	  for (unsigned int p=0; p<n_qp; p++)
	    {

	      const real g11 = (dxdxi_map(p)*dxdxi_map(p) +
				dydxi_map(p)*dydxi_map(p) +
				dzdxi_map(p)*dzdxi_map(p));
	
	      const real g12 = (dxdxi_map(p)*dxdeta_map(p) +
				dydxi_map(p)*dydeta_map(p) +
				dzdxi_map(p)*dzdeta_map(p));
	
	      const real g21 = g12;
	
	      const real g22 = (dxdeta_map(p)*dxdeta_map(p) +
				dydeta_map(p)*dydeta_map(p) +
				dzdeta_map(p)*dzdeta_map(p));
	
	
	      jac[p] = sqrt(g11*g22 - g12*g21);
	
	      assert (jac[p] > 0.);

	      JxW[p] = jac[p]*qw[p];
	    };
	};
	// done computing the map

	return;  
      }


    default:
      error();

    };

  error();

  return;
};




//--------------------------------------------------------------
// Explicit instantiations (doesn't make sense in 1D!)
template class FE<2,HIERARCHIC>;
template class FE<3,HIERARCHIC>;

template class FE<2,LAGRANGE>;
template class FE<3,LAGRANGE>;

template class FE<2,MONOMIAL>;
template class FE<3,MONOMIAL>;
