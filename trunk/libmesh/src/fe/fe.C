// $Id: fe.C,v 1.3 2003-01-20 17:06:18 jwpeterson Exp $

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



// Local includes
#include "fe.h"
#include "quadrature.h"
#include "elem.h"




// ------------------------------------------------------------
// FEBase class members
AutoPtr<FEBase> FEBase::build (const unsigned int dim,
			       const FEType& fet)
{
  // The stupid AutoPtr<FEBase> ap(); return ap;
  // construct is required to satisfy IBM's xlC
  
  switch (dim)
    {
      // 1D
    case 1:
      {
	switch (fet.family)
	  {
	  case LAGRANGE:
	    {
	      AutoPtr<FEBase> ap(new FE<1,LAGRANGE>(fet));
	      return ap;
	    };
		   
	  case HIERARCHIC:
	    {
	      AutoPtr<FEBase> ap(new FE<1,HIERARCHIC>(fet));
	      return ap;
	    };
	    
	  case MONOMIAL:
	    {
	      AutoPtr<FEBase> ap(new FE<1,MONOMIAL>(fet));
	      return ap;
	    };
	    
	  default:
	    error();
	  };
      };

      
      // 2D
    case 2:
      {
	switch (fet.family)
	  {
	  case LAGRANGE:
	    {
	      AutoPtr<FEBase> ap(new FE<2,LAGRANGE>(fet));
	      return ap;
	    };
	    
	  case HIERARCHIC:
	    {
	      AutoPtr<FEBase> ap(new FE<2,HIERARCHIC>(fet));
	      return ap;
	    };
	    
	  case MONOMIAL:
	    {
	      AutoPtr<FEBase> ap(new FE<2,MONOMIAL>(fet));
	      return ap;
	    };
	    
	  default:
	    error();
	  };
      };

      
      // 3D
    case 3:
      {
	switch (fet.family)
	  {
	  case LAGRANGE:
	    {
	      AutoPtr<FEBase> ap(new FE<3,LAGRANGE>(fet));
	      return ap;
	    };
	    
	  case HIERARCHIC:
	    {
	      AutoPtr<FEBase> ap(new FE<3,HIERARCHIC>(fet));
	      return ap;
	    };
	    
	  case MONOMIAL:
	    {
	      AutoPtr<FEBase> ap(new FE<3,MONOMIAL>(fet));
	      return ap;
	    };
	    
	  default:
	    error();
	  };
      };

    default:
      error();
    };

  error();
  AutoPtr<FEBase> ap(NULL);
  return ap;
};



template <unsigned int Dim, FEFamily T>
void FE<Dim,T>::reinit(const Elem* elem)
{
  assert (qrule   != NULL);
  assert (elem    != NULL);

  qrule->init(elem->type());
    
  // update the type in accordance to the current cell
  // and reinit if the cell type has changed or (as in
  // the case of the hierarchics) the shape functions
  // depend on the particular element
  if ((get_type() != elem->type()) ||
      (T == HIERARCHIC))
    {
      elem_type = elem->type();
      init_shape_functions (qrule, elem);
    }
  
  // Compute the map for this element.  In the future we can specify
  // different types of maps
  compute_map (qrule, elem);

  // Compute the shape functions and the derivatives at all of the
  // quadrature points.  This part is dimension-independet
  compute_shape_functions (qrule);
};




template <unsigned int Dim, FEFamily T>
void FE<Dim,T>::init_shape_functions(const QBase* qrule,
				     const Elem* elem)
{
  assert (qrule != NULL);
  assert (elem  != NULL);

  // The element type and order to use in
  // the map
  const Order    mapping_order     (elem->default_order());
  const ElemType mapping_elem_type (elem->type());
    
  // Number of shape functions in the finite element approximation
  // space.
  const unsigned int n_approx_shape_functions = n_shape_functions();

  // Number of shape functions used to construt the map
  // (Lagrange shape functions are used for mapping)
  const unsigned int n_mapping_shape_functions =
    FE<Dim,LAGRANGE>::n_shape_functions (mapping_elem_type,
					 mapping_order);

  // The number and location of the quadrature points.
  const unsigned int        n_qp = qrule->n_points();
  const std::vector<Point>&   qp = qrule->get_points();

  
  // resize the vectors to hold current data
  // Phi are the shape functions used for the FE approximation
  // Phi_map are the shape functions used for the FE mapping
  {
    phi.resize     (n_approx_shape_functions);
    dphi.resize    (n_approx_shape_functions);
    dphidx.resize  (n_approx_shape_functions);
    dphidy.resize  (n_approx_shape_functions);
    dphidz.resize  (n_approx_shape_functions);
    dphidxi.resize (n_approx_shape_functions);
    
    if (Dim > 1)
      dphideta.resize      (n_approx_shape_functions);
    
    if (Dim == 3)
      dphidzeta.resize     (n_approx_shape_functions);
      

    
    phi_map.resize         (n_mapping_shape_functions);
    dphidxi_map.resize     (n_mapping_shape_functions);
    
    if (Dim > 1)
      dphideta_map.resize  (n_mapping_shape_functions);
    
    if (Dim == 3)
      dphidzeta_map.resize (n_mapping_shape_functions);
    
    for (unsigned int i=0; i<n_approx_shape_functions; i++)
      {
	phi[i].resize         (n_qp);
	dphi[i].resize        (n_qp);
	dphidx[i].resize      (n_qp);
	dphidy[i].resize      (n_qp);
	dphidz[i].resize      (n_qp);
	dphidxi[i].resize     (n_qp);
	   
	if (Dim > 1)
	  dphideta[i].resize  (n_qp);
	    
	   
	if (Dim == 3)	     
	  dphidzeta[i].resize (n_qp);
	     
      };
       
    for (unsigned int i=0; i<n_mapping_shape_functions; i++)
      {
	phi_map[i].resize         (n_qp);
	dphidxi_map[i].resize     (n_qp);
	   
	if (Dim > 1)
	  dphideta_map[i].resize  (n_qp);
	   
	if (Dim == 3)
	  dphidzeta_map[i].resize (n_qp);
      };
  };


  
  
  switch (Dim)
    {

      //------------------------------------------------------------
      // 1D
    case 1:
      {
	// Compute the value of the approximation shape function i at quadrature point p
	for (unsigned int i=0; i<n_approx_shape_functions; i++)
	  for (unsigned int p=0; p<n_qp; p++)
	    {
	      phi[i][p]      = FE<Dim,T>::shape       (elem, get_order(), i,    qp[p]);
	      dphidxi[i][p]  = FE<Dim,T>::shape_deriv (elem, get_order(), i, 0, qp[p]);
	    };
	
	// Compute the value of the mapping shape function i at quadrature point p
	// (Lagrange shape functions are used for mapping)
	for (unsigned int i=0; i<n_mapping_shape_functions; i++)
	  for (unsigned int p=0; p<n_qp; p++)
	    {
	      phi_map[i][p]      = FE<Dim,LAGRANGE>::shape       (mapping_elem_type, mapping_order, i,    qp[p]);
	      dphidxi_map[i][p]  = FE<Dim,LAGRANGE>::shape_deriv (mapping_elem_type, mapping_order, i, 0, qp[p]);
	    };
		
	return;
      };


      
      //------------------------------------------------------------
      // 2D
    case 2:
      {
	// Compute the value of the approximation shape function i at quadrature point p
	for (unsigned int i=0; i<n_approx_shape_functions; i++)
	  for (unsigned int p=0; p<n_qp; p++)
	    {
	      phi[i][p]      = FE<Dim,T>::shape       (elem, get_order(), i,    qp[p]);
	      dphidxi[i][p]  = FE<Dim,T>::shape_deriv (elem, get_order(), i, 0, qp[p]);
	      dphideta[i][p] = FE<Dim,T>::shape_deriv (elem, get_order(), i, 1, qp[p]);
	    };
	
	// Compute the value of the mapping shape function i at quadrature point p
	// (Lagrange shape functions are used for mapping)
	for (unsigned int i=0; i<n_mapping_shape_functions; i++)
	  for (unsigned int p=0; p<n_qp; p++)
	    {
	      phi_map[i][p]      = FE<Dim,LAGRANGE>::shape       (mapping_elem_type, mapping_order, i,    qp[p]);
	      dphidxi_map[i][p]  = FE<Dim,LAGRANGE>::shape_deriv (mapping_elem_type, mapping_order, i, 0, qp[p]);
	      dphideta_map[i][p] = FE<Dim,LAGRANGE>::shape_deriv (mapping_elem_type, mapping_order, i, 1, qp[p]);
	    };
			
       	return;
      };


      
      //------------------------------------------------------------
      // 3D
    case 3:
      {
	// Compute the value of the approximation shape function i at quadrature point p
	for (unsigned int i=0; i<n_approx_shape_functions; i++)
	  for (unsigned int p=0; p<n_qp; p++)
	    {
	      phi[i][p]       = FE<Dim,T>::shape       (elem, get_order(), i,    qp[p]);
	      dphidxi[i][p]   = FE<Dim,T>::shape_deriv (elem, get_order(), i, 0, qp[p]);
	      dphideta[i][p]  = FE<Dim,T>::shape_deriv (elem, get_order(), i, 1, qp[p]);
	      dphidzeta[i][p] = FE<Dim,T>::shape_deriv (elem, get_order(), i, 2, qp[p]);
	    };
	
	// Compute the value of the mapping shape function i at quadrature point p
	// (Lagrange shape functions are used for mapping)
	for (unsigned int i=0; i<n_mapping_shape_functions; i++)
	  for (unsigned int p=0; p<n_qp; p++)
	    {
	      phi_map[i][p]       = FE<Dim,LAGRANGE>::shape       (mapping_elem_type, mapping_order, i,    qp[p]);
	      dphidxi_map[i][p]   = FE<Dim,LAGRANGE>::shape_deriv (mapping_elem_type, mapping_order, i, 0, qp[p]);
	      dphideta_map[i][p]  = FE<Dim,LAGRANGE>::shape_deriv (mapping_elem_type, mapping_order, i, 1, qp[p]);
	      dphidzeta_map[i][p] = FE<Dim,LAGRANGE>::shape_deriv (mapping_elem_type, mapping_order, i, 2, qp[p]);
	    };
			
	return;
      };


    default:
      error();
    };

  error();
  return;
};




void FEBase::compute_shape_functions(const QBase* qrule)
{
  assert (qrule != NULL);
  
  const unsigned int n_qp = qrule->n_points();


  //-------------------------------------------------------------------------
  // Compute the shape function values (and derivatives)
  // at the Quadrature points.  Note that the actual values
  // have already been computed via init_shape_functions

  // Compute the value of the derivative shape function i at quadrature point p
  switch (dim)
    {
      
    case 1:
      {
	for (unsigned int i=0; i<phi.size(); i++)
	  for (unsigned int p=0; p<n_qp; p++)
	    {
	      // dphi/dx    = (dphi/dxi)*(dxi/dx)
	      dphi[i][p](0) =
		dphidx[i][p] = dphidxi[i][p]*dxidx_map[p];
	      
	      dphi[i][p](1) = dphidy[i][p] = 0.;
	      dphi[i][p](2) = dphidz[i][p] = 0.;
	    };
	  
	break;
      };

    case 2:
      {
	for (unsigned int i=0; i<phi.size(); i++)
	  for (unsigned int p=0; p<n_qp; p++)
	    {
	      // dphi/dx    = (dphi/dxi)*(dxi/dx) + (dphi/deta)*(deta/dx)
	      dphi[i][p](0) =
		dphidx[i][p] = (dphidxi[i][p]*dxidx_map[p] +
				dphideta[i][p]*detadx_map[p]);
	      
	      // dphi/dy    = (dphi/dxi)*(dxi/dy) + (dphi/deta)*(deta/dy)
	      dphi[i][p](1) =
		dphidy[i][p] = (dphidxi[i][p]*dxidy_map[p] +
				dphideta[i][p]*detady_map[p]);

	      dphi[i][p](2) = dphidz[i][p] = 0.;
	    };

	break;
      };
    
    case 3:
      {
	for (unsigned int i=0; i<phi.size(); i++)
	  for (unsigned int p=0; p<n_qp; p++)
	    {
	      // dphi/dx    = (dphi/dxi)*(dxi/dx) + (dphi/deta)*(deta/dx) + (dphi/dzeta)*(dzeta/dx);
	      dphi[i][p](0) =
		dphidx[i][p] = (dphidxi[i][p]*dxidx_map[p] +
				dphideta[i][p]*detadx_map[p] +
				dphidzeta[i][p]*dzetadx_map[p]);
		
	      // dphi/dy    = (dphi/dxi)*(dxi/dy) + (dphi/deta)*(deta/dy) + (dphi/dzeta)*(dzeta/dy);
	      dphi[i][p](1) =
		dphidy[i][p] = (dphidxi[i][p]*dxidy_map[p] +
				dphideta[i][p]*detady_map[p] +
				dphidzeta[i][p]*dzetady_map[p]);
		
	      // dphi/dz    = (dphi/dxi)*(dxi/dz) + (dphi/deta)*(deta/dz) + (dphi/dzeta)*(dzeta/dz);
	      dphi[i][p](2) =
		dphidz[i][p] = (dphidxi[i][p]*dxidz_map[p] +
				dphideta[i][p]*detadz_map[p] +
				dphidzeta[i][p]*dzetadz_map[p]);	      
	    };

	break;
      };

    default:
      {
	error();
      };
    };
};



bool FEBase::on_reference_element(const Point& p, const ElemType t, const real eps)
{
  assert (eps >= 0.);
  
  const real xi   = p(0);
  const real eta  = p(1);
  const real zeta = p(2);
  
  switch (t)
    {

    case EDGE2:
    case EDGE3:
    case EDGE4:
      {
	// The reference 1D element is [-1,1].
	if ((xi >= -1.-eps) &&
	    (xi <=  1.+eps))
	  return true;

	break;
      };

      
    case TRI3:
    case TRI6:
      {
	// The reference triangle is isocoles
	// and is bound by xi=0, eta=0, and xi+eta=1.
	if ((xi  >= 0.-eps) &&
	    (eta >= 0.-eps) &&
	    ((xi + eta) <= 1.+eps))
	  return true;

	break;
      };

      
    case QUAD4:
    case QUAD8:
    case QUAD9:
      {
	// The reference quadrilateral element is [-1,1]^2.
	if ((xi  >= -1.-eps) &&
	    (xi  <=  1.+eps) &&
	    (eta >= -1.-eps) &&
	    (eta <=  1.+eps))
	  return true;
		
	break;
      };


    case TET4:
    case TET10:
      {
	// The reference tetrahedral is isocoles
	// and is bound by xi=0, eta=0, zeta=0,
	// and xi+eta+zeta=1.
	if ((xi   >= 0.-eps) &&
	    (eta  >= 0.-eps) &&
	    (zeta >= 0.-eps) &&
	    ((xi + eta + zeta) <= 1.+eps))
	  return true;
		
	break;
      };

      
    case HEX8:
    case HEX20:
    case HEX27:
      {
	/*
	  if ((xi   >= -1.) &&
	  (xi   <=  1.) &&
	  (eta  >= -1.) &&
	  (eta  <=  1.) &&
	  (zeta >= -1.) &&
	  (zeta <=  1.))
	  return true;
	*/
	
	// The reference hexahedral element is [-1,1]^3.
	if ((xi   >= -1.-eps) &&
	    (xi   <=  1.+eps) &&
	    (eta  >= -1.-eps) &&
	    (eta  <=  1.+eps) &&
	    (zeta >= -1.-eps) &&
	    (zeta <=  1.+eps))
	  {
	    //	    std::cout << "Strange Point:\n";
	    //	    p.print();
	    return true;
	  }

	break;
      };

    case PRISM6:
    case PRISM18:
      {
	// Figure this one out...
	if ((xi   >= 0.-eps) &&
	    (eta  >= 0.-eps) &&
	    (zeta >= 0.-eps) &&
	    (zeta <= 1.+eps) &&
	    ((xi + eta) <= 1.+eps))
	  return true;

	break;
      };


    case PYRAMID5:
      {
	std::cerr << "BEN: Implement this you lazy bastard!"
		  << std::endl;
	error();

	break;
      };
      
    default:
      std::cerr << "ERROR: Unknown element type " << t << std::endl;
      error();
    };

  // If we get here then the point is _not_ in the
  // reference element.   Better return false.
  
  return false;
};




//--------------------------------------------------------------
// Explicit instantiations
template class FE<1,HIERARCHIC>;
template class FE<2,HIERARCHIC>;
template class FE<3,HIERARCHIC>;

template class FE<1,LAGRANGE>;
template class FE<2,LAGRANGE>;
template class FE<3,LAGRANGE>;

template class FE<1,MONOMIAL>;
template class FE<2,MONOMIAL>;
template class FE<3,MONOMIAL>;


