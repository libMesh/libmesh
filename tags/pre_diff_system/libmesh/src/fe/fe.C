// $Id: fe.C,v 1.48 2006-05-27 10:55:11 roystgnr Exp $

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
#include "fe.h"
#include "fe_interface.h"
#include "quadrature_gauss.h"
#include "elem.h"
#include "libmesh_logging.h"
#include "fe_macro.h"



// ------------------------------------------------------------
// FE class members
template <unsigned int Dim, FEFamily T>
unsigned int FE<Dim,T>::n_shape_functions () const
{
  return FE<Dim,T>::n_dofs (elem_type,
           static_cast<Order>(fe_type.order + _p_level));
}


template <unsigned int Dim, FEFamily T>
void FE<Dim,T>::attach_quadrature_rule (QBase* q)
{
  assert (q != NULL); 
  qrule = q;
  // make sure we don't cache results from a previous quadrature rule
  elem_type = INVALID_ELEM;
  return;
}


template <unsigned int Dim, FEFamily T>
unsigned int FE<Dim,T>::n_quadrature_points () const
{ 
  assert (qrule != NULL);  
  return qrule->n_points(); 
}


template <unsigned int Dim, FEFamily T>
void FE<Dim,T>::dofs_on_side(const Elem* const elem,
			     const Order o,
			     unsigned int s,
			     std::vector<unsigned int>& di)
{
  assert(elem != NULL);
  assert(s < elem->n_sides());

  di.clear();
  unsigned int nodenum = 0;
  const unsigned int n_nodes = elem->n_nodes();
  for (unsigned int n = 0; n != n_nodes; ++n)
    {
      const unsigned int n_dofs = n_dofs_at_node(elem->type(),
						 static_cast<Order>(o + elem->p_level()), n);
      if (elem->is_node_on_side(n, s))
	for (unsigned int i = 0; i != n_dofs; ++i)
	  di.push_back(nodenum++);
      else
	nodenum += n_dofs;
    }
}



template <unsigned int Dim, FEFamily T>
void FE<Dim,T>::dofs_on_edge(const Elem* const elem,
			     const Order o,
			     unsigned int e,
			     std::vector<unsigned int>& di)
{
  assert(elem != NULL);
  assert(e < elem->n_edges());

  di.clear();
  unsigned int nodenum = 0;
  const unsigned int n_nodes = elem->n_nodes();
  for (unsigned int n = 0; n != n_nodes; ++n)
    {
      const unsigned int n_dofs = n_dofs_at_node(elem->type(),
						 static_cast<Order>(o + elem->p_level()), n);
      if (elem->is_node_on_edge(n, e))
	for (unsigned int i = 0; i != n_dofs; ++i)
	  di.push_back(nodenum++);
      else
	nodenum += n_dofs;
    }
}



template <unsigned int Dim, FEFamily T>
void FE<Dim,T>::reinit(const Elem* elem,
		       const std::vector<Point>* const pts)
{
  assert (elem != NULL);

  // Initialize the shape functions at the user-specified
  // points
  if (pts != NULL)
    {
      // Set the type and p level for this element
      elem_type = elem->type();
      _p_level = elem->p_level();

      // Initialize the shape functions
      this->init_shape_functions (*pts, elem);

      // The shape functions do not correspond to the qrule
      shapes_on_quadrature = false;
    }
  
  // If there are no user specified points, we use the
  // quadrature rule
  
  // update the type in accordance to the current cell
  // and reinit if the cell type has changed or (as in
  // the case of the hierarchics) the shape functions need
  // reinit, since they depend on the particular element
  else 
    {
      assert (qrule   != NULL);
      qrule->init(elem->type(), elem->p_level());

      if (elem_type != elem->type() ||
          _p_level != elem->p_level() ||
          !shapes_on_quadrature ||
	  this->shapes_need_reinit())
        {
          // Set the type and p level for this element
          elem_type = elem->type();
          _p_level = elem->p_level();
          // Initialize the shape functions
          this->init_shape_functions (qrule->get_points(), elem);
        }
      else
        {
          // Set the type and p level for this element
          elem_type = elem->type();
          _p_level = elem->p_level();
        }
      
      // The shape functions correspond to the qrule
      shapes_on_quadrature = true;
    }
  
  // Compute the map for this element.  In the future we can specify
  // different types of maps
  if (pts != NULL)
    {
      std::vector<Real> dummy_weights (pts->size(), 1.);
      
      this->compute_map (dummy_weights, elem);
    }
  else
    {
      this->compute_map (qrule->get_weights(), elem);
    }

  // Compute the shape functions and the derivatives at all of the
  // quadrature points.
  this->compute_shape_functions (elem);
}



template <unsigned int Dim, FEFamily T>
void FE<Dim,T>::init_shape_functions(const std::vector<Point>& qp,
				     const Elem* elem)
{
  assert (elem  != NULL);
  calculations_started = true;
  
  // Start logging the shape function initialization
  START_LOG("init_shape_functions()", "FE");

  
  // The number of quadrature points.
  const unsigned int n_qp = qp.size();

  // The element type and order to use in
  // the map
  const Order    mapping_order     (elem->default_order());
  const ElemType mapping_elem_type (elem->type());
    
  // Number of shape functions in the finite element approximation
  // space.
  const unsigned int n_approx_shape_functions =
    this->n_shape_functions(this->get_type(),
			    this->get_order());

  // Number of shape functions used to construt the map
  // (Lagrange shape functions are used for mapping)
  const unsigned int n_mapping_shape_functions =
    FE<Dim,LAGRANGE>::n_shape_functions (mapping_elem_type,
					 mapping_order);
  
  // resize the vectors to hold current data
  // Phi are the shape functions used for the FE approximation
  // Phi_map are the shape functions used for the FE mapping
  {
    if (calculate_phi)
      phi.resize     (n_approx_shape_functions);
    if (calculate_dphi)
      {
        dphi.resize    (n_approx_shape_functions);
        dphidx.resize  (n_approx_shape_functions);
        dphidy.resize  (n_approx_shape_functions);
        dphidz.resize  (n_approx_shape_functions);
        dphidxi.resize (n_approx_shape_functions);

        if (Dim > 1)
          dphideta.resize      (n_approx_shape_functions);
    
        if (Dim == 3)
          dphidzeta.resize     (n_approx_shape_functions);
      }
#ifdef ENABLE_SECOND_DERIVATIVES
    if (calculate_d2phi)
      {
        d2phi.resize     (n_approx_shape_functions);
        d2phidx2.resize  (n_approx_shape_functions);
        d2phidxdy.resize (n_approx_shape_functions);
        d2phidxdz.resize (n_approx_shape_functions);
        d2phidy2.resize  (n_approx_shape_functions);
        d2phidydz.resize (n_approx_shape_functions);
        d2phidz2.resize  (n_approx_shape_functions);
        d2phidxi2.resize (n_approx_shape_functions);
        if (Dim > 1)
          {
            d2phidxideta.resize (n_approx_shape_functions);
            d2phideta2.resize   (n_approx_shape_functions);
          }
        if (Dim > 2)
          {
            d2phidxidzeta.resize  (n_approx_shape_functions);
            d2phidetadzeta.resize (n_approx_shape_functions);
            d2phidzeta2.resize    (n_approx_shape_functions);
          }
      }
#endif // ifdef ENABLE_SECOND_DERIVATIVES
    
    phi_map.resize         (n_mapping_shape_functions);
    dphidxi_map.resize     (n_mapping_shape_functions);
#ifdef ENABLE_SECOND_DERIVATIVES
    d2phidxi2_map.resize   (n_mapping_shape_functions);
#endif // ifdef ENABLE_SECOND_DERIVATIVES
    
    if (Dim > 1)
      {
        dphideta_map.resize  (n_mapping_shape_functions);
#ifdef ENABLE_SECOND_DERIVATIVES
        d2phidxideta_map.resize   (n_mapping_shape_functions);
        d2phideta2_map.resize     (n_mapping_shape_functions);
#endif // ifdef ENABLE_SECOND_DERIVATIVES
      }
    
    if (Dim == 3)
      {
        dphidzeta_map.resize (n_mapping_shape_functions);
#ifdef ENABLE_SECOND_DERIVATIVES
        d2phidxidzeta_map.resize  (n_mapping_shape_functions);
        d2phidetadzeta_map.resize (n_mapping_shape_functions);
        d2phidzeta2_map.resize    (n_mapping_shape_functions);
#endif // ifdef ENABLE_SECOND_DERIVATIVES
      }
    
    for (unsigned int i=0; i<n_approx_shape_functions; i++)
      {
	if (calculate_phi)
	  phi[i].resize         (n_qp);
	if (calculate_dphi)
	  {
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
#ifdef ENABLE_SECOND_DERIVATIVES
	if (calculate_d2phi)
	  {
            d2phi[i].resize     (n_qp);
            d2phidx2[i].resize  (n_qp);
            d2phidxdy[i].resize (n_qp);
            d2phidxdz[i].resize (n_qp);
            d2phidy2[i].resize  (n_qp);
            d2phidydz[i].resize (n_qp);
            d2phidz2[i].resize  (n_qp);
            d2phidxi2[i].resize (n_qp);
            if (Dim > 1)
              {
                d2phidxideta[i].resize (n_qp);
                d2phideta2[i].resize   (n_qp);
              }
            if (Dim > 2)
              {
                d2phidxidzeta[i].resize  (n_qp);
                d2phidetadzeta[i].resize (n_qp);
                d2phidzeta2[i].resize    (n_qp);
              }
	  }
#endif // ifdef ENABLE_SECOND_DERIVATIVES
      }
       
    for (unsigned int i=0; i<n_mapping_shape_functions; i++)
      {
	phi_map[i].resize         (n_qp);
	dphidxi_map[i].resize     (n_qp);
#ifdef ENABLE_SECOND_DERIVATIVES
	d2phidxi2_map[i].resize     (n_qp);
	if (Dim > 1)
          {
            d2phidxideta_map[i].resize (n_qp);
            d2phideta2_map[i].resize (n_qp);
          }
        if (Dim > 2)
          {
            d2phidxidzeta_map[i].resize  (n_qp);
            d2phidetadzeta_map[i].resize (n_qp);
            d2phidzeta2_map[i].resize    (n_qp);
          }
#endif // ifdef ENABLE_SECOND_DERIVATIVES
	   
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

  // Optimize for the affine elements case:
  bool has_affine_map = elem->has_affine_map();
  
  switch (Dim)
    {

      //------------------------------------------------------------
      // 1D
    case 1:
      {
	// Compute the value of the approximation shape function i at quadrature point p
	if (calculate_phi)
	  for (unsigned int i=0; i<n_approx_shape_functions; i++)
	    for (unsigned int p=0; p<n_qp; p++)
	      phi[i][p]      = FE<Dim,T>::shape       (elem, fe_type.order, i,    qp[p]);
	if (calculate_dphi)
	  for (unsigned int i=0; i<n_approx_shape_functions; i++)
	    for (unsigned int p=0; p<n_qp; p++)
	      dphidxi[i][p]  = FE<Dim,T>::shape_deriv (elem, fe_type.order, i, 0, qp[p]);
#ifdef ENABLE_SECOND_DERIVATIVES
	if (calculate_d2phi)
	  for (unsigned int i=0; i<n_approx_shape_functions; i++)
	    for (unsigned int p=0; p<n_qp; p++)
              d2phidxi2[i][p] = FE<Dim,T>::shape_second_deriv (elem, fe_type.order, i, 0, qp[p]);
#endif // ifdef ENABLE_SECOND_DERIVATIVES
	
	// Compute the value of the mapping shape function i at quadrature point p
	// (Lagrange shape functions are used for mapping)
        if (has_affine_map)
          {
	    for (unsigned int i=0; i<n_mapping_shape_functions; i++)
              {
	        phi_map[i][0]      = FE<Dim,LAGRANGE>::shape       (mapping_elem_type, mapping_order, i,    qp[0]);
	        dphidxi_map[i][0]  = FE<Dim,LAGRANGE>::shape_deriv (mapping_elem_type, mapping_order, i, 0, qp[0]);
#ifdef ENABLE_SECOND_DERIVATIVES
                d2phidxi2_map[i][0] = FE<Dim,LAGRANGE>::shape_second_deriv (mapping_elem_type, mapping_order, i, 0, qp[0]);
#endif // ifdef ENABLE_SECOND_DERIVATIVES
	        for (unsigned int p=1; p<n_qp; p++)
                  {
	            phi_map[i][p]      = FE<Dim,LAGRANGE>::shape (mapping_elem_type, mapping_order, i,    qp[p]);
	            dphidxi_map[i][p]  = dphidxi_map[i][0];
#ifdef ENABLE_SECOND_DERIVATIVES
                    d2phidxi2_map[i][p] = d2phidxi2_map[i][0];
#endif // ifdef ENABLE_SECOND_DERIVATIVES
                  }
              }
          }
        else
	  for (unsigned int i=0; i<n_mapping_shape_functions; i++)
	    for (unsigned int p=0; p<n_qp; p++)
	      {
	        phi_map[i][p]      = FE<Dim,LAGRANGE>::shape       (mapping_elem_type, mapping_order, i,    qp[p]);
	        dphidxi_map[i][p]  = FE<Dim,LAGRANGE>::shape_deriv (mapping_elem_type, mapping_order, i, 0, qp[p]);
#ifdef ENABLE_SECOND_DERIVATIVES
                d2phidxi2_map[i][p] = FE<Dim,LAGRANGE>::shape_second_deriv (mapping_elem_type, mapping_order, i, 0, qp[p]);
#endif // ifdef ENABLE_SECOND_DERIVATIVES
	      }
		
	break;
      }


      
      //------------------------------------------------------------
      // 2D
    case 2:
      {
	// Compute the value of the approximation shape function i at quadrature point p
	if (calculate_phi)
	  for (unsigned int i=0; i<n_approx_shape_functions; i++)
	    for (unsigned int p=0; p<n_qp; p++)
	      phi[i][p]      = FE<Dim,T>::shape       (elem, fe_type.order, i,    qp[p]);
	if (calculate_dphi)
	  for (unsigned int i=0; i<n_approx_shape_functions; i++)
	    for (unsigned int p=0; p<n_qp; p++)
	      {
	        dphidxi[i][p]  = FE<Dim,T>::shape_deriv (elem, fe_type.order, i, 0, qp[p]);
	        dphideta[i][p] = FE<Dim,T>::shape_deriv (elem, fe_type.order, i, 1, qp[p]);
	      }
#ifdef ENABLE_SECOND_DERIVATIVES
	if (calculate_d2phi)
	  for (unsigned int i=0; i<n_approx_shape_functions; i++)
	    for (unsigned int p=0; p<n_qp; p++)
	      {
                d2phidxi2[i][p] = FE<Dim,T>::shape_second_deriv (elem, fe_type.order, i, 0, qp[p]);
                d2phidxideta[i][p] = FE<Dim,T>::shape_second_deriv (elem, fe_type.order, i, 1, qp[p]);
                d2phideta2[i][p] = FE<Dim,T>::shape_second_deriv (elem, fe_type.order, i, 2, qp[p]);
	      }
#endif // ifdef ENABLE_SECOND_DERIVATIVES
	
	// Compute the value of the mapping shape function i at quadrature point p
	// (Lagrange shape functions are used for mapping)
        if (has_affine_map)
          {
	    for (unsigned int i=0; i<n_mapping_shape_functions; i++)
              {
	        phi_map[i][0]      = FE<Dim,LAGRANGE>::shape       (mapping_elem_type, mapping_order, i,    qp[0]);
	        dphidxi_map[i][0]  = FE<Dim,LAGRANGE>::shape_deriv (mapping_elem_type, mapping_order, i, 0, qp[0]);
	        dphideta_map[i][0] = FE<Dim,LAGRANGE>::shape_deriv (mapping_elem_type, mapping_order, i, 1, qp[0]);
#ifdef ENABLE_SECOND_DERIVATIVES
                d2phidxi2_map[i][0] = FE<Dim,LAGRANGE>::shape_second_deriv (mapping_elem_type, mapping_order, i, 0, qp[0]);
                d2phidxideta_map[i][0] = FE<Dim,LAGRANGE>::shape_second_deriv (mapping_elem_type, mapping_order, i, 1, qp[0]);
                d2phideta2_map[i][0] = FE<Dim,LAGRANGE>::shape_second_deriv (mapping_elem_type, mapping_order, i, 2, qp[0]);
#endif // ifdef ENABLE_SECOND_DERIVATIVES
	        for (unsigned int p=1; p<n_qp; p++)
                  {
	            phi_map[i][p]      = FE<Dim,LAGRANGE>::shape (mapping_elem_type, mapping_order, i,    qp[p]);
	            dphidxi_map[i][p]  = dphidxi_map[i][0];
	            dphideta_map[i][p] = dphideta_map[i][0];
#ifdef ENABLE_SECOND_DERIVATIVES
                    d2phidxi2_map[i][p] = d2phidxi2_map[i][0];
                    d2phidxideta_map[i][p] = d2phidxideta_map[i][0];
                    d2phideta2_map[i][p] = d2phideta2_map[i][0];
#endif // ifdef ENABLE_SECOND_DERIVATIVES
                  }
              }
          }
        else
	  for (unsigned int i=0; i<n_mapping_shape_functions; i++)
	    for (unsigned int p=0; p<n_qp; p++)
	      {
	        phi_map[i][p]      = FE<Dim,LAGRANGE>::shape       (mapping_elem_type, mapping_order, i,    qp[p]);
	        dphidxi_map[i][p]  = FE<Dim,LAGRANGE>::shape_deriv (mapping_elem_type, mapping_order, i, 0, qp[p]);
	        dphideta_map[i][p] = FE<Dim,LAGRANGE>::shape_deriv (mapping_elem_type, mapping_order, i, 1, qp[p]);
#ifdef ENABLE_SECOND_DERIVATIVES
                d2phidxi2_map[i][p] = FE<Dim,LAGRANGE>::shape_second_deriv (mapping_elem_type, mapping_order, i, 0, qp[p]);
                d2phidxideta_map[i][p] = FE<Dim,LAGRANGE>::shape_second_deriv (mapping_elem_type, mapping_order, i, 1, qp[p]);
                d2phideta2_map[i][p] = FE<Dim,LAGRANGE>::shape_second_deriv (mapping_elem_type, mapping_order, i, 2, qp[p]);
#endif // ifdef ENABLE_SECOND_DERIVATIVES
	      }
			
       	break;
      }


      
      //------------------------------------------------------------
      // 3D
    case 3:
      {
	// Compute the value of the approximation shape function i at quadrature point p
	if (calculate_phi)
	  for (unsigned int i=0; i<n_approx_shape_functions; i++)
	    for (unsigned int p=0; p<n_qp; p++)
	      phi[i][p]       = FE<Dim,T>::shape       (elem, fe_type.order, i,    qp[p]);
	if (calculate_dphi)
	  for (unsigned int i=0; i<n_approx_shape_functions; i++)
	    for (unsigned int p=0; p<n_qp; p++)
	      {
	        dphidxi[i][p]   = FE<Dim,T>::shape_deriv (elem, fe_type.order, i, 0, qp[p]);
	        dphideta[i][p]  = FE<Dim,T>::shape_deriv (elem, fe_type.order, i, 1, qp[p]);
	        dphidzeta[i][p] = FE<Dim,T>::shape_deriv (elem, fe_type.order, i, 2, qp[p]);
	      }
#ifdef ENABLE_SECOND_DERIVATIVES
	if (calculate_d2phi)
	  for (unsigned int i=0; i<n_approx_shape_functions; i++)
	    for (unsigned int p=0; p<n_qp; p++)
	      {
                d2phidxi2[i][p] = FE<Dim,T>::shape_second_deriv (elem, fe_type.order, i, 0, qp[p]);
                d2phidxideta[i][p] = FE<Dim,T>::shape_second_deriv (elem, fe_type.order, i, 1, qp[p]);
                d2phideta2[i][p] = FE<Dim,T>::shape_second_deriv (elem, fe_type.order, i, 2, qp[p]);
                d2phidxidzeta[i][p] = FE<Dim,T>::shape_second_deriv (elem, fe_type.order, i, 3, qp[p]);
                d2phidetadzeta[i][p] = FE<Dim,T>::shape_second_deriv (elem, fe_type.order, i, 4, qp[p]);
                d2phidzeta2[i][p] = FE<Dim,T>::shape_second_deriv (elem, fe_type.order, i, 5, qp[p]);
	      }
#endif // ifdef ENABLE_SECOND_DERIVATIVES
	
	// Compute the value of the mapping shape function i at quadrature point p
	// (Lagrange shape functions are used for mapping)
        if (has_affine_map)
          {
	    for (unsigned int i=0; i<n_mapping_shape_functions; i++)
              {
	        phi_map[i][0]      = FE<Dim,LAGRANGE>::shape       (mapping_elem_type, mapping_order, i,    qp[0]);
	        dphidxi_map[i][0]  = FE<Dim,LAGRANGE>::shape_deriv (mapping_elem_type, mapping_order, i, 0, qp[0]);
	        dphideta_map[i][0] = FE<Dim,LAGRANGE>::shape_deriv (mapping_elem_type, mapping_order, i, 1, qp[0]);
	        dphidzeta_map[i][0] = FE<Dim,LAGRANGE>::shape_deriv (mapping_elem_type, mapping_order, i, 2, qp[0]);
#ifdef ENABLE_SECOND_DERIVATIVES
                d2phidxi2_map[i][0] = FE<Dim,LAGRANGE>::shape_second_deriv (mapping_elem_type, mapping_order, i, 0, qp[0]);
                d2phidxideta_map[i][0] = FE<Dim,LAGRANGE>::shape_second_deriv (mapping_elem_type, mapping_order, i, 1, qp[0]);
                d2phideta2_map[i][0] = FE<Dim,LAGRANGE>::shape_second_deriv (mapping_elem_type, mapping_order, i, 2, qp[0]);
                d2phidxideta_map[i][0] = FE<Dim,LAGRANGE>::shape_second_deriv (mapping_elem_type, mapping_order, i, 3, qp[0]);
                d2phidetadzeta_map[i][0] = FE<Dim,LAGRANGE>::shape_second_deriv (mapping_elem_type, mapping_order, i, 4, qp[0]);
                d2phidzeta2_map[i][0] = FE<Dim,LAGRANGE>::shape_second_deriv (mapping_elem_type, mapping_order, i, 5, qp[0]);
#endif // ifdef ENABLE_SECOND_DERIVATIVES
	        for (unsigned int p=1; p<n_qp; p++)
                  {
	            phi_map[i][p]      = FE<Dim,LAGRANGE>::shape (mapping_elem_type, mapping_order, i,    qp[p]);
	            dphidxi_map[i][p]  = dphidxi_map[i][0];
	            dphideta_map[i][p] = dphideta_map[i][0];
	            dphidzeta_map[i][p] = dphidzeta_map[i][0];
#ifdef ENABLE_SECOND_DERIVATIVES
                    d2phidxi2_map[i][p] = d2phidxi2_map[i][0];
                    d2phidxideta_map[i][p] = d2phidxideta_map[i][0];
                    d2phideta2_map[i][p] = d2phideta2_map[i][0];
                    d2phidxideta_map[i][p] = d2phidxideta_map[i][0];
                    d2phidetadzeta_map[i][p] = d2phidetadzeta_map[i][0];
                    d2phidzeta2_map[i][p] = d2phidzeta2_map[i][0];
#endif // ifdef ENABLE_SECOND_DERIVATIVES
                  }
              }
          }
        else
	  for (unsigned int i=0; i<n_mapping_shape_functions; i++)
	    for (unsigned int p=0; p<n_qp; p++)
	      {
	        phi_map[i][p]       = FE<Dim,LAGRANGE>::shape       (mapping_elem_type, mapping_order, i,    qp[p]);
	        dphidxi_map[i][p]   = FE<Dim,LAGRANGE>::shape_deriv (mapping_elem_type, mapping_order, i, 0, qp[p]);
	        dphideta_map[i][p]  = FE<Dim,LAGRANGE>::shape_deriv (mapping_elem_type, mapping_order, i, 1, qp[p]);
	        dphidzeta_map[i][p] = FE<Dim,LAGRANGE>::shape_deriv (mapping_elem_type, mapping_order, i, 2, qp[p]);
#ifdef ENABLE_SECOND_DERIVATIVES
                d2phidxi2_map[i][p] = FE<Dim,LAGRANGE>::shape_second_deriv (mapping_elem_type, mapping_order, i, 0, qp[p]);
                d2phidxideta_map[i][p] = FE<Dim,LAGRANGE>::shape_second_deriv (mapping_elem_type, mapping_order, i, 1, qp[p]);
                d2phideta2_map[i][p] = FE<Dim,LAGRANGE>::shape_second_deriv (mapping_elem_type, mapping_order, i, 2, qp[p]);
                d2phidxideta_map[i][p] = FE<Dim,LAGRANGE>::shape_second_deriv (mapping_elem_type, mapping_order, i, 3, qp[p]);
                d2phidetadzeta_map[i][p] = FE<Dim,LAGRANGE>::shape_second_deriv (mapping_elem_type, mapping_order, i, 4, qp[p]);
                d2phidzeta2_map[i][p] = FE<Dim,LAGRANGE>::shape_second_deriv (mapping_elem_type, mapping_order, i, 5, qp[p]);
#endif // ifdef ENABLE_SECOND_DERIVATIVES
	      }
			
	break;
      }


    default:
      error();
    }
  
  // Stop logging the shape function initialization
  STOP_LOG("init_shape_functions()", "FE");
}




template <unsigned int Dim, FEFamily T>
void FE<Dim,T>::compute_proj_constraints (DofConstraints &constraints,
				          DofMap &dof_map,
				          const unsigned int variable_number,
				          const Elem* elem)
{
#ifdef ENABLE_AMR
  // Only constrain elements in 2,3D.
  if (Dim == 1)
    return;

  // Only constrain active elements with this method
  if (!elem->active())
    return;

  assert (elem != NULL);

  const FEType& base_fe_type = dof_map.variable_type(variable_number);

  // Construct FE objects for this element and its neighbors.
  AutoPtr<FEBase> my_fe (FEBase::build(Dim, base_fe_type));
  const FEContinuity cont = my_fe->get_continuity();
  if (cont == DISCONTINUOUS)
    return;
  assert (cont == C_ZERO || cont == C_ONE);

  AutoPtr<FEBase> neigh_fe (FEBase::build(Dim, base_fe_type));

  QGauss my_qface(Dim-1, base_fe_type.default_quadrature_order());
  my_fe->attach_quadrature_rule (&my_qface);
  std::vector<Point> neigh_qface;

  const std::vector<Real>& JxW = my_fe->get_JxW();
  const std::vector<Point>& q_point = my_fe->get_xyz();
  const std::vector<std::vector<Real> >& phi = my_fe->get_phi();
  const std::vector<std::vector<Real> >& neigh_phi =
		  neigh_fe->get_phi();
  const std::vector<Point> *face_normals = NULL;
  const std::vector<std::vector<RealGradient> > *dphi = NULL;
  const std::vector<std::vector<RealGradient> > *neigh_dphi = NULL;
  std::vector<unsigned int> my_dof_indices, neigh_dof_indices;
  std::vector<unsigned int> my_side_dofs, neigh_side_dofs;

  if (cont != C_ZERO)
    {
      const std::vector<Point>& ref_face_normals =
        my_fe->get_normals();
      face_normals = &ref_face_normals;
      const std::vector<std::vector<RealGradient> >& ref_dphi =
	my_fe->get_dphi();
      dphi = &ref_dphi;
      const std::vector<std::vector<RealGradient> >& ref_neigh_dphi =
	neigh_fe->get_dphi();
      neigh_dphi = &ref_neigh_dphi;
    }

  DenseMatrix<Real> Ke;
  DenseVector<Real> Fe;
  std::vector<DenseVector<Real> > Ue;

  // Look at the element faces.  Check to see if we need to
  // build constraints.
  for (unsigned int s=0; s<elem->n_sides(); s++)
    if (elem->neighbor(s) != NULL)
      {
        // Get pointers to the element's neighbor.
        const Elem* neigh = elem->neighbor(s);

        // h refinement constraints:
        // constrain dofs shared between
        // this element and ones coarser
        // than this element.
        if (neigh->level() < elem->level()) 
          {
            const Elem *ancestor = elem;
            while (neigh->level() < ancestor->level())
              ancestor = ancestor->parent();
	    unsigned int s_neigh = neigh->which_neighbor_am_i(ancestor);

            // Find the minimum p level; we build the h constraint
            // matrix with this and then constrain away all higher p
            // DoFs.
            assert(neigh->active());
            const unsigned int min_p_level =
              std::min(elem->p_level(), neigh->p_level());

            // we may need to make the FE objects reinit with the
            // minimum shared p_level
            // FIXME - I hate using const_cast<> and avoiding
            // accessor functions; there's got to be a
            // better way to do this!
            const unsigned int old_elem_level = elem->p_level();
            if (old_elem_level != min_p_level)
              (const_cast<Elem *>(elem))->hack_p_level(min_p_level);
            const unsigned int old_neigh_level = neigh->p_level();
            if (old_neigh_level != min_p_level)
              (const_cast<Elem *>(neigh))->hack_p_level(min_p_level);

	    my_fe->reinit(elem, s);

	    dof_map.dof_indices (elem, my_dof_indices,
			         variable_number);
	    dof_map.dof_indices (neigh, neigh_dof_indices,
			         variable_number);

	    const unsigned int n_qp = my_qface.n_points();

	    FEInterface::inverse_map (Dim, base_fe_type, neigh,
                                      q_point, neigh_qface);

	    neigh_fe->reinit(neigh, &neigh_qface);

	    // We're only concerned with DOFs whose values (and/or first
	    // derivatives for C1 elements) are supported on side nodes
	    dofs_on_side(elem, base_fe_type.order, s, my_side_dofs);
	    dofs_on_side(neigh, base_fe_type.order, s_neigh, neigh_side_dofs);

            // We're done with functions that examine Elem::p_level(),
            // so let's unhack those levels
            if (elem->p_level() != old_elem_level)
              (const_cast<Elem *>(elem))->hack_p_level(old_elem_level);
            if (neigh->p_level() != old_neigh_level)
              (const_cast<Elem *>(neigh))->hack_p_level(old_neigh_level);

	    const unsigned int n_side_dofs = my_side_dofs.size();
	    assert(n_side_dofs == neigh_side_dofs.size());

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
                      {
		        Ke(is,js) += JxW[qp] * (phi[i][qp] * phi[j][qp]);
                        if (cont != C_ZERO)
		          Ke(is,js) += JxW[qp] * (((*dphi)[i][qp] *
					         (*face_normals)[qp]) *
					        ((*dphi)[j][qp] *
					         (*face_normals)[qp]));
                      }
		  }
	      }

	    // Form the right hand sides, (inner product of coarse basis
	    // functions against fine test functions)
	    for (unsigned int is = 0; is != n_side_dofs; ++is)
	      {
	        const unsigned int i = neigh_side_dofs[is];
	        Fe.resize (n_side_dofs);
	        for (unsigned int js = 0; js != n_side_dofs; ++js)
		  {
	            const unsigned int j = my_side_dofs[js];
	            for (unsigned int qp = 0; qp != n_qp; ++qp)
                      {
		        Fe(js) += JxW[qp] * (neigh_phi[i][qp] *
					     phi[j][qp]);
                        if (cont != C_ZERO)
		          Fe(js) += JxW[qp] * (((*neigh_dphi)[i][qp] *
					        (*face_normals)[qp]) *
					       ((*dphi)[j][qp] *
					        (*face_normals)[qp]));
                      }
		  }
	        Ke.cholesky_solve(Fe, Ue[is]);
	      }
	    for (unsigned int is = 0; is != n_side_dofs; ++is)
	      {
	        const unsigned int i = neigh_side_dofs[is];
	        const unsigned int their_dof_g = neigh_dof_indices[i];
                assert(their_dof_g != DofObject::invalid_id);
	        for (unsigned int js = 0; js != n_side_dofs; ++js)
	          {
	            const unsigned int j = my_side_dofs[js];
	            const unsigned int my_dof_g = my_dof_indices[j];
                    assert(my_dof_g != DofObject::invalid_id);
		    const Real their_dof_value = Ue[is](js);
		    if (their_dof_g == my_dof_g)
		      {
		        assert(std::abs(their_dof_value-1.) < 1.e-5);
		        for (unsigned int k = 0; k != n_side_dofs; ++k)
		          assert(k == is || std::abs(Ue[k](js)) < 1.e-5);
		        continue;
		      }
		    if (std::abs(their_dof_value) < 1.e-5)
		      continue;

		    DofConstraintRow& constraint_row =
                      constraints[my_dof_g];

		    constraint_row.insert(std::make_pair(their_dof_g,
						         their_dof_value));
	          }
	      }
	  }
        // p refinement constraints:
        // constrain dofs shared between
        // active elements and neighbors with
        // lower polynomial degrees
        const unsigned int min_p_level =
          neigh->min_p_level_by_neighbor(elem, elem->p_level());
        if (min_p_level < elem->p_level())
          {
            // Adaptive p refinement of non-hierarchic bases will
            // require more coding
            assert(my_fe->is_hierarchic());
            dof_map.constrain_p_dofs(variable_number, elem,
                                     s, min_p_level);
          }
      }
#endif
}


    
#ifdef ENABLE_INFINITE_ELEMENTS

template <unsigned int Dim, FEFamily T>
void FE<Dim,T>::init_base_shape_functions(const std::vector<Point>& qp,
					  const Elem* e)
{ 
  // I don't understand infinite elements well enough to risk
  // calculating too little.  :-(  RHS
  calculate_phi = calculate_dphi = calculate_d2phi = true;

  elem_type = e->type(); 
  init_shape_functions(qp, e); 
}

#endif
    


//--------------------------------------------------------------
// Explicit instantiations using macro from fe_macro.h

INSTANTIATE_FE(1);

INSTANTIATE_FE(2);

INSTANTIATE_FE(3);

