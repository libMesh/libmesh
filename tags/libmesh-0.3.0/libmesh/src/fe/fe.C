// $Id: fe.C,v 1.11 2003-02-13 22:56:09 benkirk Exp $

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
#include "inf_fe.h"
#include "quadrature.h"
#include "elem.h"


// ------------------------------------------------------------
// FE class members
template <unsigned int Dim, FEFamily T>
unsigned int FE<Dim,T>::n_quadrature_points () const
{ 
  assert (qrule != NULL);  
  return qrule->n_points(); 
}




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
}




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
	     
      }
       
    for (unsigned int i=0; i<n_mapping_shape_functions; i++)
      {
	phi_map[i].resize         (n_qp);
	dphidxi_map[i].resize     (n_qp);
	   
	if (Dim > 1)
	  dphideta_map[i].resize  (n_qp);
	   
	if (Dim == 3)
	  dphidzeta_map[i].resize (n_qp);
      }
  }


      
#ifdef ENABLE_INFINITE_ELEMENTS
  //------------------------------------------------------------
  // Initialize the data fields, which should only be used for infinite 
  // elements, to some sensible values, so that using a FE with the
  // variational formulation of an InfFE, correct element matrices are
  // returned

 {
    weight.resize  (n_qp);
    dweight.resize (n_qp);
    dphase.resize  (n_qp);
    
    for (unsigned int p=0; p<n_qp; p++)
      {
        weight[p] = 1.;
	dweight[p].zero();
	dphase[p].zero();
      }

 }
#endif // ifdef ENABLE_INFINITE_ELEMENTS


  
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
	    }
	
	// Compute the value of the mapping shape function i at quadrature point p
	// (Lagrange shape functions are used for mapping)
	for (unsigned int i=0; i<n_mapping_shape_functions; i++)
	  for (unsigned int p=0; p<n_qp; p++)
	    {
	      phi_map[i][p]      = FE<Dim,LAGRANGE>::shape       (mapping_elem_type, mapping_order, i,    qp[p]);
	      dphidxi_map[i][p]  = FE<Dim,LAGRANGE>::shape_deriv (mapping_elem_type, mapping_order, i, 0, qp[p]);
	    }
		
	return;
      }


      
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
	    }
	
	// Compute the value of the mapping shape function i at quadrature point p
	// (Lagrange shape functions are used for mapping)
	for (unsigned int i=0; i<n_mapping_shape_functions; i++)
	  for (unsigned int p=0; p<n_qp; p++)
	    {
	      phi_map[i][p]      = FE<Dim,LAGRANGE>::shape       (mapping_elem_type, mapping_order, i,    qp[p]);
	      dphidxi_map[i][p]  = FE<Dim,LAGRANGE>::shape_deriv (mapping_elem_type, mapping_order, i, 0, qp[p]);
	      dphideta_map[i][p] = FE<Dim,LAGRANGE>::shape_deriv (mapping_elem_type, mapping_order, i, 1, qp[p]);
	    }
			
       	return;
      }


      
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
	    }
	
	// Compute the value of the mapping shape function i at quadrature point p
	// (Lagrange shape functions are used for mapping)
	for (unsigned int i=0; i<n_mapping_shape_functions; i++)
	  for (unsigned int p=0; p<n_qp; p++)
	    {
	      phi_map[i][p]       = FE<Dim,LAGRANGE>::shape       (mapping_elem_type, mapping_order, i,    qp[p]);
	      dphidxi_map[i][p]   = FE<Dim,LAGRANGE>::shape_deriv (mapping_elem_type, mapping_order, i, 0, qp[p]);
	      dphideta_map[i][p]  = FE<Dim,LAGRANGE>::shape_deriv (mapping_elem_type, mapping_order, i, 1, qp[p]);
	      dphidzeta_map[i][p] = FE<Dim,LAGRANGE>::shape_deriv (mapping_elem_type, mapping_order, i, 2, qp[p]);
	    }
			
	return;
      }


    default:
      error();
    }
}




    
#ifdef ENABLE_INFINITE_ELEMENTS

template <unsigned int Dim, FEFamily T>
void FE<Dim,T>::init_base_shape_functions(const QBase* q, 
					  const Elem* e)
{ 
  elem_type = e->type(); 
  init_shape_functions(q, e); 
}

#endif


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


