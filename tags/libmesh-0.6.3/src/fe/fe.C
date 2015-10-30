// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson
  
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
#include "elem.h"
#include "libmesh_logging.h"
#include "fe_macro.h"
#include "quadrature.h"



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
  libmesh_assert (q != NULL); 
  qrule = q;
  // make sure we don't cache results from a previous quadrature rule
  elem_type = INVALID_ELEM;
  return;
}


template <unsigned int Dim, FEFamily T>
unsigned int FE<Dim,T>::n_quadrature_points () const
{ 
  libmesh_assert (qrule != NULL);  
  return qrule->n_points(); 
}


template <unsigned int Dim, FEFamily T>
void FE<Dim,T>::dofs_on_side(const Elem* const elem,
			     const Order o,
			     unsigned int s,
			     std::vector<unsigned int>& di)
{
  libmesh_assert(elem != NULL);
  libmesh_assert(s < elem->n_sides());

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
  libmesh_assert(elem != NULL);
  libmesh_assert(e < elem->n_edges());

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
  libmesh_assert (elem != NULL);

  // Try to avoid calling init_shape_functions
  // even when shapes_need_reinit
  bool cached_nodes_still_fit = false;

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
  // reinit, since they depend on the particular element shape
  else 
    {
      libmesh_assert (qrule   != NULL);
      qrule->init(elem->type(), elem->p_level());

      if (elem_type != elem->type() ||
          _p_level != elem->p_level() ||
          !shapes_on_quadrature)
        {
          // Set the type and p level for this element
          elem_type = elem->type();
          _p_level = elem->p_level();
          // Initialize the shape functions
          this->init_shape_functions (qrule->get_points(), elem);

          if (this->shapes_need_reinit())
            {
              cached_nodes.resize(elem->n_nodes());
              for (unsigned int n = 0; n != elem->n_nodes(); ++n)
                {
                  cached_nodes[n] = elem->point(n);
                }
            }
        }
      else
        {
          libmesh_assert(elem->n_nodes() > 1);

          cached_nodes_still_fit = true;
          if (cached_nodes.size() != elem->n_nodes())
            cached_nodes_still_fit = false;
          else
            for (unsigned int n = 1; n < elem->n_nodes(); ++n)
              {
                if ((elem->point(n) - elem->point(0)) !=
                    (cached_nodes[n] - cached_nodes[0]))
                  {
                    cached_nodes_still_fit = false;
                    break;
                  }
              }

          if (this->shapes_need_reinit() && !cached_nodes_still_fit)
            {
              this->init_shape_functions (qrule->get_points(), elem);
              cached_nodes.resize(elem->n_nodes());
              for (unsigned int n = 0; n != elem->n_nodes(); ++n)
                cached_nodes[n] = elem->point(n);
            }
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
  if (!cached_nodes_still_fit)
    this->compute_shape_functions (elem);
}



template <unsigned int Dim, FEFamily T>
void FE<Dim,T>::init_shape_functions(const std::vector<Point>& qp,
				     const Elem* elem)
{
  libmesh_assert (elem  != NULL);
  calculations_started = true;

  // If the user forgot to request anything, we'll be safe and
  // calculate everything:
#ifdef ENABLE_SECOND_DERIVATIVES
  if (!calculate_phi && !calculate_dphi && !calculate_d2phi)
    calculate_phi = calculate_dphi = calculate_d2phi = true;
#else
  if (!calculate_phi && !calculate_dphi)
    calculate_phi = calculate_dphi = true;
#endif // ENABLE_SECOND_DERIVATIVES
  
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
      libmesh_error();
    }
  
  // Stop logging the shape function initialization
  STOP_LOG("init_shape_functions()", "FE");
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

#endif // ENABLE_INFINITE_ELEMENTS
    


//--------------------------------------------------------------
// Explicit instantiations using macro from fe_macro.h

INSTANTIATE_FE(1);

INSTANTIATE_FE(2);

INSTANTIATE_FE(3);

