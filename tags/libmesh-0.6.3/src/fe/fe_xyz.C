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
#include "libmesh_logging.h"
#include "fe.h"
#include "fe_macro.h"
#include "elem.h"




// ------------------------------------------------------------
// XYZ-specific implementations
template <unsigned int Dim>
void FEXYZ<Dim>::init_shape_functions(const std::vector<Point>& qp,
				      const Elem* elem)
{
  libmesh_assert (elem  != NULL);
  this->calculations_started = true;

  // If the user forgot to request anything, we'll be safe and
  // calculate everything:
#ifdef ENABLE_SECOND_DERIVATIVES
  if (!this->calculate_phi && !this->calculate_dphi && !this->calculate_d2phi)
    this->calculate_phi = this->calculate_dphi = this->calculate_d2phi = true;
#else
  if (!this->calculate_phi && !this->calculate_dphi)
    this->calculate_phi = this->calculate_dphi = true;
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
    // (note: GCC 3.4.0 requires the use of this-> here)
    if (this->calculate_phi)
      this->phi.resize     (n_approx_shape_functions);
    if (this->calculate_dphi)
      {
        this->dphi.resize    (n_approx_shape_functions);
        this->dphidx.resize  (n_approx_shape_functions);
        this->dphidy.resize  (n_approx_shape_functions);
        this->dphidz.resize  (n_approx_shape_functions);
      }
          
      if (Dim > 1)
        this->dphideta_map.resize  (n_mapping_shape_functions);
    
      if (Dim == 3)
        this->dphidzeta_map.resize (n_mapping_shape_functions);
#ifdef ENABLE_SECOND_DERIVATIVES
    if (this->calculate_d2phi)
      {
        this->d2phi.resize     (n_approx_shape_functions);
        this->d2phidx2.resize  (n_approx_shape_functions);
        this->d2phidxdy.resize (n_approx_shape_functions);
        this->d2phidxdz.resize (n_approx_shape_functions);
        this->d2phidy2.resize  (n_approx_shape_functions);
        this->d2phidydz.resize (n_approx_shape_functions);
        this->d2phidz2.resize  (n_approx_shape_functions);
        this->d2phidxi2.resize (n_approx_shape_functions);
      }
      if (Dim > 1)
        {
          this->d2phidxideta_map.resize (n_mapping_shape_functions);
          this->d2phideta2_map.resize   (n_mapping_shape_functions);
        }
      if (Dim > 2)
        {
          this->d2phidxidzeta_map.resize  (n_mapping_shape_functions);
          this->d2phidetadzeta_map.resize (n_mapping_shape_functions);
          this->d2phidzeta2_map.resize    (n_mapping_shape_functions);
        }
#endif // ifdef ENABLE_SECOND_DERIVATIVES
    
    this->phi_map.resize         (n_mapping_shape_functions);
    this->dphidxi_map.resize     (n_mapping_shape_functions);
#ifdef ENABLE_SECOND_DERIVATIVES
    this->d2phidxi2_map.resize   (n_mapping_shape_functions);
#endif // ifdef ENABLE_SECOND_DERIVATIVES
    
    for (unsigned int i=0; i<n_approx_shape_functions; i++)
      {
        if (this->calculate_phi)
	  this->phi[i].resize           (n_qp);
        if (this->calculate_dphi)
          {
	    this->dphi[i].resize        (n_qp);
	    this->dphidx[i].resize      (n_qp);
	    this->dphidy[i].resize      (n_qp);
	    this->dphidz[i].resize      (n_qp);
          }
#ifdef ENABLE_SECOND_DERIVATIVES
        if (this->calculate_d2phi)
          {
	    this->d2phi[i].resize       (n_qp);
	    this->d2phidx2[i].resize    (n_qp);
	    this->d2phidxdy[i].resize   (n_qp);
	    this->d2phidy2[i].resize    (n_qp);
	    this->d2phidydz[i].resize   (n_qp);
	    this->d2phidz2[i].resize    (n_qp);
          }
#endif // ifdef ENABLE_SECOND_DERIVATIVES
      }
       
    for (unsigned int i=0; i<n_mapping_shape_functions; i++)
      {
	this->phi_map[i].resize         (n_qp);
	this->dphidxi_map[i].resize     (n_qp);
#ifdef ENABLE_SECOND_DERIVATIVES
	this->d2phidxi2_map[i].resize   (n_qp);
        if (Dim > 1)
          {
	    this->d2phidxideta_map[i].resize   (n_qp);
	    this->d2phideta2_map[i].resize     (n_qp);
          }
	if (Dim > 2)
          {
	    this->d2phidxidzeta_map[i].resize  (n_qp);
	    this->d2phidetadzeta_map[i].resize (n_qp);
	    this->d2phidzeta2_map[i].resize    (n_qp);
          }
#endif // ifdef ENABLE_SECOND_DERIVATIVES
	   
	if (Dim > 1)
	  this->dphideta_map[i].resize  (n_qp);
	   
	if (Dim == 3)
	  this->dphidzeta_map[i].resize (n_qp);
      }
  }


      
#ifdef ENABLE_INFINITE_ELEMENTS
  //------------------------------------------------------------
  // Initialize the data fields, which should only be used for infinite 
  // elements, to some sensible values, so that using a FE with the
  // variational formulation of an InfFE, correct element matrices are
  // returned

 {
    this->weight.resize  (n_qp);
    this->dweight.resize (n_qp);
    this->dphase.resize  (n_qp);
    
    for (unsigned int p=0; p<n_qp; p++)
      {
        this->weight[p] = 1.;
	this->dweight[p].zero();
	this->dphase[p].zero();
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
        // Compute the value of the mapping shape function i at quadrature point p
        // (Lagrange shape functions are used for mapping)
        if (has_affine_map)
          {
            for (unsigned int i=0; i<n_mapping_shape_functions; i++)
              {
                this->phi_map[i][0]      = FE<Dim,LAGRANGE>::shape       (mapping_elem_type, mapping_order, i,    qp[0]);
                this->dphidxi_map[i][0]  = FE<Dim,LAGRANGE>::shape_deriv (mapping_elem_type, mapping_order, i, 0, qp[0]);
#ifdef ENABLE_SECOND_DERIVATIVES
                this->d2phidxi2_map[i][0] = FE<Dim,LAGRANGE>::shape_second_deriv (mapping_elem_type, mapping_order, i, 0, qp[0]);
#endif // ifdef ENABLE_SECOND_DERIVATIVES
                for (unsigned int p=1; p<n_qp; p++)
                  {
                    this->phi_map[i][p]      = FE<Dim,LAGRANGE>::shape (mapping_elem_type, mapping_order, i,    qp[p]);
                    this->dphidxi_map[i][p]  = this->dphidxi_map[i][0];
#ifdef ENABLE_SECOND_DERIVATIVES
                    this->d2phidxi2_map[i][p] = this->d2phidxi2_map[i][0];
#endif // ifdef ENABLE_SECOND_DERIVATIVES
                  }
              }
          }
        else
          for (unsigned int i=0; i<n_mapping_shape_functions; i++)
            for (unsigned int p=0; p<n_qp; p++)
              {
                this->phi_map[i][p]      = FE<Dim,LAGRANGE>::shape       (mapping_elem_type, mapping_order, i,    qp[p]);
                this->dphidxi_map[i][p]  = FE<Dim,LAGRANGE>::shape_deriv (mapping_elem_type, mapping_order, i, 0, qp[p]);
#ifdef ENABLE_SECOND_DERIVATIVES
                this->d2phidxi2_map[i][p] = FE<Dim,LAGRANGE>::shape_second_deriv (mapping_elem_type, mapping_order, i, 0, qp[p]);
#endif // ifdef ENABLE_SECOND_DERIVATIVES
              }

        break;

      }


      
      //------------------------------------------------------------
      // 2D
    case 2:
      {
 	// Compute the value of the mapping shape function i at quadrature point p
	// (Lagrange shape functions are used for mapping)
        if (has_affine_map)
          {
	    for (unsigned int i=0; i<n_mapping_shape_functions; i++)
              {
	        this->phi_map[i][0]      = FE<Dim,LAGRANGE>::shape       (mapping_elem_type, mapping_order, i,    qp[0]);
	        this->dphidxi_map[i][0]  = FE<Dim,LAGRANGE>::shape_deriv (mapping_elem_type, mapping_order, i, 0, qp[0]);
	        this->dphideta_map[i][0] = FE<Dim,LAGRANGE>::shape_deriv (mapping_elem_type, mapping_order, i, 1, qp[0]);
#ifdef ENABLE_SECOND_DERIVATIVES
                this->d2phidxi2_map[i][0] = FE<Dim,LAGRANGE>::shape_second_deriv (mapping_elem_type, mapping_order, i, 0, qp[0]);
                this->d2phidxideta_map[i][0] = FE<Dim,LAGRANGE>::shape_second_deriv (mapping_elem_type, mapping_order, i, 1, qp[0]);
                this->d2phideta2_map[i][0] = FE<Dim,LAGRANGE>::shape_second_deriv (mapping_elem_type, mapping_order, i, 2, qp[0]);
#endif // ifdef ENABLE_SECOND_DERIVATIVES
	        for (unsigned int p=1; p<n_qp; p++)
                  {
	            this->phi_map[i][p]      = FE<Dim,LAGRANGE>::shape (mapping_elem_type, mapping_order, i,    qp[p]);
	            this->dphidxi_map[i][p]  = this->dphidxi_map[i][0];
	            this->dphideta_map[i][p] = this->dphideta_map[i][0];
#ifdef ENABLE_SECOND_DERIVATIVES
                    this->d2phidxi2_map[i][p] = this->d2phidxi2_map[i][0];
                    this->d2phidxideta_map[i][p] = this->d2phidxideta_map[i][0];
                    this->d2phideta2_map[i][p] = this->d2phideta2_map[i][0];
#endif // ifdef ENABLE_SECOND_DERIVATIVES
                  }
              }
          }
        else
	  for (unsigned int i=0; i<n_mapping_shape_functions; i++)
	    for (unsigned int p=0; p<n_qp; p++)
	      {
	        this->phi_map[i][p]      = FE<Dim,LAGRANGE>::shape       (mapping_elem_type, mapping_order, i,    qp[p]);
	        this->dphidxi_map[i][p]  = FE<Dim,LAGRANGE>::shape_deriv (mapping_elem_type, mapping_order, i, 0, qp[p]);
	        this->dphideta_map[i][p] = FE<Dim,LAGRANGE>::shape_deriv (mapping_elem_type, mapping_order, i, 1, qp[p]);
#ifdef ENABLE_SECOND_DERIVATIVES
                this->d2phidxi2_map[i][p] = FE<Dim,LAGRANGE>::shape_second_deriv (mapping_elem_type, mapping_order, i, 0, qp[p]);
                this->d2phidxideta_map[i][p] = FE<Dim,LAGRANGE>::shape_second_deriv (mapping_elem_type, mapping_order, i, 1, qp[p]);
                this->d2phideta2_map[i][p] = FE<Dim,LAGRANGE>::shape_second_deriv (mapping_elem_type, mapping_order, i, 2, qp[p]);
#endif // ifdef ENABLE_SECOND_DERIVATIVES
	      }
			
       	break;
      }


      
      //------------------------------------------------------------
      // 3D
    case 3:
      {
	// Compute the value of the mapping shape function i at quadrature point p
	// (Lagrange shape functions are used for mapping)
        if (has_affine_map)
          {
	    for (unsigned int i=0; i<n_mapping_shape_functions; i++)
              {
	        this->phi_map[i][0]      = FE<Dim,LAGRANGE>::shape       (mapping_elem_type, mapping_order, i,    qp[0]);
	        this->dphidxi_map[i][0]  = FE<Dim,LAGRANGE>::shape_deriv (mapping_elem_type, mapping_order, i, 0, qp[0]);
	        this->dphideta_map[i][0] = FE<Dim,LAGRANGE>::shape_deriv (mapping_elem_type, mapping_order, i, 1, qp[0]);
	        this->dphidzeta_map[i][0] = FE<Dim,LAGRANGE>::shape_deriv (mapping_elem_type, mapping_order, i, 2, qp[0]);
#ifdef ENABLE_SECOND_DERIVATIVES
                this->d2phidxi2_map[i][0] = FE<Dim,LAGRANGE>::shape_second_deriv (mapping_elem_type, mapping_order, i, 0, qp[0]);
                this->d2phidxideta_map[i][0] = FE<Dim,LAGRANGE>::shape_second_deriv (mapping_elem_type, mapping_order, i, 1, qp[0]);
                this->d2phideta2_map[i][0] = FE<Dim,LAGRANGE>::shape_second_deriv (mapping_elem_type, mapping_order, i, 2, qp[0]);
                this->d2phidxideta_map[i][0] = FE<Dim,LAGRANGE>::shape_second_deriv (mapping_elem_type, mapping_order, i, 3, qp[0]);
                this->d2phidetadzeta_map[i][0] = FE<Dim,LAGRANGE>::shape_second_deriv (mapping_elem_type, mapping_order, i, 4, qp[0]);
                this->d2phidzeta2_map[i][0] = FE<Dim,LAGRANGE>::shape_second_deriv (mapping_elem_type, mapping_order, i, 5, qp[0]);
#endif // ifdef ENABLE_SECOND_DERIVATIVES
	        for (unsigned int p=1; p<n_qp; p++)
                  {
	            this->phi_map[i][p]      = FE<Dim,LAGRANGE>::shape (mapping_elem_type, mapping_order, i,    qp[p]);
	            this->dphidxi_map[i][p]  = this->dphidxi_map[i][0];
	            this->dphideta_map[i][p] = this->dphideta_map[i][0];
	            this->dphidzeta_map[i][p] = this->dphidzeta_map[i][0];
#ifdef ENABLE_SECOND_DERIVATIVES
                    this->d2phidxi2_map[i][p] = this->d2phidxi2_map[i][0];
                    this->d2phidxideta_map[i][p] = this->d2phidxideta_map[i][0];
                    this->d2phideta2_map[i][p] = this->d2phideta2_map[i][0];
                    this->d2phidxideta_map[i][p] = this->d2phidxideta_map[i][0];
                    this->d2phidetadzeta_map[i][p] = this->d2phidetadzeta_map[i][0];
                    this->d2phidzeta2_map[i][p] = this->d2phidzeta2_map[i][0];
#endif // ifdef ENABLE_SECOND_DERIVATIVES
                  }
              }
          }
        else
	  for (unsigned int i=0; i<n_mapping_shape_functions; i++)
	    for (unsigned int p=0; p<n_qp; p++)
	      {
	        this->phi_map[i][p]       = FE<Dim,LAGRANGE>::shape       (mapping_elem_type, mapping_order, i,    qp[p]);
	        this->dphidxi_map[i][p]   = FE<Dim,LAGRANGE>::shape_deriv (mapping_elem_type, mapping_order, i, 0, qp[p]);
	        this->dphideta_map[i][p]  = FE<Dim,LAGRANGE>::shape_deriv (mapping_elem_type, mapping_order, i, 1, qp[p]);
	        this->dphidzeta_map[i][p] = FE<Dim,LAGRANGE>::shape_deriv (mapping_elem_type, mapping_order, i, 2, qp[p]);
#ifdef ENABLE_SECOND_DERIVATIVES
                this->d2phidxi2_map[i][p] = FE<Dim,LAGRANGE>::shape_second_deriv (mapping_elem_type, mapping_order, i, 0, qp[p]);
                this->d2phidxideta_map[i][p] = FE<Dim,LAGRANGE>::shape_second_deriv (mapping_elem_type, mapping_order, i, 1, qp[p]);
                this->d2phideta2_map[i][p] = FE<Dim,LAGRANGE>::shape_second_deriv (mapping_elem_type, mapping_order, i, 2, qp[p]);
                this->d2phidxideta_map[i][p] = FE<Dim,LAGRANGE>::shape_second_deriv (mapping_elem_type, mapping_order, i, 3, qp[p]);
                this->d2phidetadzeta_map[i][p] = FE<Dim,LAGRANGE>::shape_second_deriv (mapping_elem_type, mapping_order, i, 4, qp[p]);
                this->d2phidzeta2_map[i][p] = FE<Dim,LAGRANGE>::shape_second_deriv (mapping_elem_type, mapping_order, i, 5, qp[p]);
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




template <unsigned int Dim>
void FEXYZ<Dim>::compute_shape_functions (const Elem* elem)
{
  libmesh_assert (elem != NULL);
  
  //-------------------------------------------------------------------------
  // Compute the shape function values (and derivatives)
  // at the Quadrature points.  Note that the actual values
  // have already been computed via init_shape_functions

  // Start logging the shape function computation
  START_LOG("compute_shape_functions()", "FE");

  const std::vector<Point>& xyz_qp = this->get_xyz();
  
  // Compute the value of the derivative shape function i at quadrature point p
  switch (this->dim)
    {
      
    case 1:
      {
        if (this->calculate_phi)
	  for (unsigned int i=0; i<this->phi.size(); i++)
	    for (unsigned int p=0; p<this->phi[i].size(); p++)
	      this->phi[i][p] = FE<Dim,XYZ>::shape (elem, this->fe_type.order, i, xyz_qp[p]);
        if (this->calculate_dphi)
	  for (unsigned int i=0; i<this->dphi.size(); i++)
	    for (unsigned int p=0; p<this->dphi[i].size(); p++)
	      {
	        this->dphi[i][p](0) =
		  this->dphidx[i][p] = FE<Dim,XYZ>::shape_deriv (elem, this->fe_type.order, i, 0, xyz_qp[p]);
	      
	        this->dphi[i][p](1) = this->dphidy[i][p] = 0.;
	        this->dphi[i][p](2) = this->dphidz[i][p] = 0.;
	      }
#ifdef ENABLE_SECOND_DERIVATIVES
        if (this->calculate_d2phi)
	  for (unsigned int i=0; i<this->d2phi.size(); i++)
	    for (unsigned int p=0; p<this->d2phi[i].size(); p++)
	      {
	        this->d2phi[i][p](0,0) =
		  this->d2phidx2[i][p] = FE<Dim,XYZ>::shape_second_deriv (elem, this->fe_type.order, i, 0, xyz_qp[p]);
	      
#if DIM>1
	        this->d2phi[i][p](0,1) = this->d2phidxdy[i][p] =
	        this->d2phi[i][p](1,0) = 0.;
	        this->d2phi[i][p](1,1) = this->d2phidy2[i][p] = 0.;
#if DIM>2
	        this->d2phi[i][p](0,2) = this->d2phidxdz[i][p] =
	        this->d2phi[i][p](2,0) = 0.;
	        this->d2phi[i][p](1,2) = this->d2phidydz[i][p] =
	        this->d2phi[i][p](2,1) = 0.;
	        this->d2phi[i][p](2,2) = this->d2phidz2[i][p] = 0.;
#endif
#endif
	      }
#endif // ifdef ENABLE_SECOND_DERIVATIVES

	// All done
	break;
      }

    case 2:
      {
        if (this->calculate_phi)
	  for (unsigned int i=0; i<this->phi.size(); i++)
	    for (unsigned int p=0; p<this->phi[i].size(); p++)
	      this->phi[i][p] = FE<Dim,XYZ>::shape (elem, this->fe_type.order, i, xyz_qp[p]);
        if (this->calculate_dphi)
	  for (unsigned int i=0; i<this->dphi.size(); i++)
	    for (unsigned int p=0; p<this->dphi[i].size(); p++)
	      {
	        this->dphi[i][p](0) =
		  this->dphidx[i][p] = FE<Dim,XYZ>::shape_deriv (elem, this->fe_type.order, i, 0, xyz_qp[p]);
	      
	        this->dphi[i][p](1) =
		  this->dphidy[i][p] = FE<Dim,XYZ>::shape_deriv (elem, this->fe_type.order, i, 1, xyz_qp[p]);
	      
#if DIM == 3  
	        this->dphi[i][p](2) = // can only assign to the Z component if DIM==3
#endif
		this->dphidz[i][p] = 0.;
	      }
#ifdef ENABLE_SECOND_DERIVATIVES
        if (this->calculate_d2phi)
	  for (unsigned int i=0; i<this->d2phi.size(); i++)
	    for (unsigned int p=0; p<this->d2phi[i].size(); p++)
	      {
	        this->d2phi[i][p](0,0) =
		  this->d2phidx2[i][p] = FE<Dim,XYZ>::shape_second_deriv (elem, this->fe_type.order, i, 0, xyz_qp[p]);
	      
	        this->d2phi[i][p](0,1) = this->d2phidxdy[i][p] =
	        this->d2phi[i][p](1,0) = FE<Dim,XYZ>::shape_second_deriv (elem, this->fe_type.order, i, 1, xyz_qp[p]);
	        this->d2phi[i][p](1,1) = 
                  this->d2phidy2[i][p] = FE<Dim,XYZ>::shape_second_deriv (elem, this->fe_type.order, i, 2, xyz_qp[p]);
#if DIM>2
	        this->d2phi[i][p](0,2) = this->d2phidxdz[i][p] =
	        this->d2phi[i][p](2,0) = 0.;
	        this->d2phi[i][p](1,2) = this->d2phidydz[i][p] =
	        this->d2phi[i][p](2,1) = 0.;
	        this->d2phi[i][p](2,2) = this->d2phidz2[i][p] = 0.;
#endif
	      }
#endif // ifdef ENABLE_SECOND_DERIVATIVES

	// All done
	break;
      }
    
    case 3:
      {
        if (this->calculate_dphi)
	  for (unsigned int i=0; i<this->phi.size(); i++)
	    for (unsigned int p=0; p<this->phi[i].size(); p++)
	      this->phi[i][p] = FE<Dim,XYZ>::shape (elem, this->fe_type.order, i, xyz_qp[p]);
	       
        if (this->calculate_dphi)
	  for (unsigned int i=0; i<this->dphi.size(); i++)
	    for (unsigned int p=0; p<this->dphi[i].size(); p++)
	      {
	        this->dphi[i][p](0) =
		  this->dphidx[i][p] = FE<Dim,XYZ>::shape_deriv (elem, this->fe_type.order, i, 0, xyz_qp[p]);
		
	        this->dphi[i][p](1) =
		  this->dphidy[i][p] = FE<Dim,XYZ>::shape_deriv (elem, this->fe_type.order, i, 1, xyz_qp[p]);
		
	        this->dphi[i][p](2) =
		  this->dphidz[i][p] = FE<Dim,XYZ>::shape_deriv (elem, this->fe_type.order, i, 2, xyz_qp[p]);	      
	      }
#ifdef ENABLE_SECOND_DERIVATIVES
        if (this->calculate_d2phi)
	  for (unsigned int i=0; i<this->d2phi.size(); i++)
	    for (unsigned int p=0; p<this->d2phi[i].size(); p++)
	      {
	        this->d2phi[i][p](0,0) =
		  this->d2phidx2[i][p] = FE<Dim,XYZ>::shape_second_deriv (elem, this->fe_type.order, i, 0, xyz_qp[p]);
	      
	        this->d2phi[i][p](0,1) = this->d2phidxdy[i][p] =
	        this->d2phi[i][p](1,0) = FE<Dim,XYZ>::shape_second_deriv (elem, this->fe_type.order, i, 1, xyz_qp[p]);
	        this->d2phi[i][p](1,1) = 
                  this->d2phidy2[i][p] = FE<Dim,XYZ>::shape_second_deriv (elem, this->fe_type.order, i, 2, xyz_qp[p]);
	        this->d2phi[i][p](0,2) = this->d2phidxdz[i][p] =
	        this->d2phi[i][p](2,0) = FE<Dim,XYZ>::shape_second_deriv (elem, this->fe_type.order, i, 3, xyz_qp[p]);
	        this->d2phi[i][p](1,2) = this->d2phidydz[i][p] =
	        this->d2phi[i][p](2,1) = FE<Dim,XYZ>::shape_second_deriv (elem, this->fe_type.order, i, 4, xyz_qp[p]);
	        this->d2phi[i][p](2,2) = this->d2phidz2[i][p] = FE<Dim,XYZ>::shape_second_deriv (elem, this->fe_type.order, i, 5, xyz_qp[p]);
	      }
#endif // ifdef ENABLE_SECOND_DERIVATIVES

	// All done
	break;
      }

    default:
      {
	libmesh_error();
      }
    }
  
  // Stop logging the shape function computation
  STOP_LOG("compute_shape_functions()", "FE");
}



template <unsigned int Dim, FEFamily T>
void FE<Dim,T>::nodal_soln(const Elem* elem,
			   const Order order,
			   const std::vector<Number>& elem_soln,
			   std::vector<Number>&       nodal_soln)
{
  const unsigned int n_nodes = elem->n_nodes();
  
  const ElemType type = elem->type();

  nodal_soln.resize(n_nodes);

  const Order totalorder = static_cast<Order>(order + elem->p_level());
  
  switch (totalorder)
    {
      // Constant shape functions
    case CONSTANT:
      {
	libmesh_assert (elem_soln.size() == 1);
	
	const Number val = elem_soln[0];
	
	for (unsigned int n=0; n<n_nodes; n++)
	  nodal_soln[n] = val;
	
	return;
      }


      // For other bases do interpolation at the nodes
      // explicitly.
    default:
      {

	const unsigned int n_sf =
	  FE<Dim,T>::n_shape_functions(type, totalorder);
	
	for (unsigned int n=0; n<n_nodes; n++)
	  {
	    libmesh_assert (elem_soln.size() == n_sf);

	    // Zero before summation
	    nodal_soln[n] = 0;

	    // u_i = Sum (alpha_i phi_i)
	    for (unsigned int i=0; i<n_sf; i++)
	      nodal_soln[n] += elem_soln[i]*FE<Dim,T>::shape(elem,
							     order,
							     i,
							     elem->point(n));
	  }

	return;
      }
    }
}



template <unsigned int Dim, FEFamily T>
unsigned int FE<Dim,T>::n_dofs(const ElemType t, const Order o)
{
  switch (o)
    {

      // constant shape functions
      // no matter what shape there is only one DOF.
    case CONSTANT:
      return 1;


      // Discontinuous linear shape functions
      // expressed in the XYZ monomials.
    case FIRST:
      {
	switch (t)
	  {
	  case EDGE2:
	  case EDGE3:
	  case EDGE4:
	    return 2;

	  case TRI3:
	  case TRI6:
	  case QUAD4:
	  case QUAD8:
	  case QUAD9:
	    return 3;

	  case TET4:
	  case TET10:
	  case HEX8:
	  case HEX20:
	  case HEX27:
	  case PRISM6:
	  case PRISM15:
	  case PRISM18:
	  case PYRAMID5:
	    return 4;
	    
	  default:
	    {
#ifdef DEBUG
	      std::cerr << "ERROR: Bad ElemType = " << t
			<< " for " << o << "th order approximation!" 
			<< std::endl;
#endif
	      libmesh_error();	    
	    }
	  }
      }


      // Discontinuous quadratic shape functions
      // expressed in the XYZ monomials.
    case SECOND:
      {
	switch (t)
	  {
	  case EDGE2:
	  case EDGE3:
	  case EDGE4:
	    return 3;

	  case TRI3:
	  case TRI6:
	  case QUAD4:
	  case QUAD8:
	  case QUAD9:
	    return 6;

	  case TET4:
	  case TET10:
	  case HEX8:
	  case HEX20:
	  case HEX27:
	  case PRISM6:
	  case PRISM15:
	  case PRISM18:
	  case PYRAMID5:
	    return 10;
	    
	  default:
	    {
#ifdef DEBUG
	      std::cerr << "ERROR: Bad ElemType = " << t
			<< " for " << o << "th order approximation!" 
			<< std::endl;
#endif
	      libmesh_error();	    
	    }
	  }
      }


      // Discontinuous cubic shape functions
      // expressed in the XYZ monomials.
    case THIRD:
      {
	switch (t)
	  {
	  case EDGE2:
	  case EDGE3:
	  case EDGE4:
	    return 4;

	  case TRI3:
	  case TRI6:
	  case QUAD4:
	  case QUAD8:
	  case QUAD9:
	    return 10;

	  case TET4:
	  case TET10:
	  case HEX8:
	  case HEX20:
	  case HEX27:
	  case PRISM6:
	  case PRISM15:
	  case PRISM18:
	  case PYRAMID5:
	    return 20;
	    
	  default:
	    {
#ifdef DEBUG
	      std::cerr << "ERROR: Bad ElemType = " << t
			<< " for " << o << "th order approximation!" 
			<< std::endl;
#endif
	      libmesh_error();	    
	    }
	  }
      }


      // Discontinuous quartic shape functions
      // expressed in the XYZ monomials.
    case FOURTH:
      {
	switch (t)
	  {
	  case EDGE2:
	  case EDGE3:
	    return 5;

	  case TRI3:
	  case TRI6:
	  case QUAD4:
	  case QUAD8:
	  case QUAD9:
	    return 15;

	  case TET4:
	  case TET10:
	  case HEX8:
	  case HEX20:
	  case HEX27:
	  case PRISM6:
	  case PRISM15:
	  case PRISM18:
	  case PYRAMID5:
	    return 35;
	    
	  default:
	    {
#ifdef DEBUG
	      std::cerr << "ERROR: Bad ElemType = " << t
			<< " for " << o << "th order approximation!" 
			<< std::endl;
#endif
	      libmesh_error();	    
	    }
	  }
      }

      
    default:
      {
        const unsigned int order = static_cast<unsigned int>(o);
	switch (t)
	  {
	  case EDGE2:
	  case EDGE3:
	    return (order+1);

	  case TRI3:
	  case TRI6:
	  case QUAD4:
	  case QUAD8:
	  case QUAD9:
	    return (order+1)*(order+2)/2;

	  case TET4:
	  case TET10:
	  case HEX8:
	  case HEX20:
	  case HEX27:
	  case PRISM6:
	  case PRISM15:
	  case PRISM18:
	  case PYRAMID5:
	    return (order+1)*(order+2)*(order+3)/6;
	    
	  default:
	    {
#ifdef DEBUG
	      std::cerr << "ERROR: Bad ElemType = " << t
			<< " for " << o << "th order approximation!" 
			<< std::endl;
#endif
	      libmesh_error();	    
	    }
	  }
      }
    }
  
  libmesh_error();
  
  return 0;
}



template <unsigned int Dim, FEFamily T>
unsigned int FE<Dim,T>::n_dofs_at_node(const ElemType,
				       const Order,
				       const unsigned int)
{
  // Monomials elements have no dofs at nodes
  // (just on the element)
  return 0;
}



template <unsigned int Dim, FEFamily T>
unsigned int FE<Dim,T>::n_dofs_per_elem(const ElemType t,
					const Order o)
{
  switch (o)
    {
      // constant shape functions always have 1 DOF per element
    case CONSTANT:
      return 1;

      
      // Discontinuous linear shape functions
      // expressed in the XYZ monomials.
    case FIRST:
      {
	switch (t)
	  {
	    // 1D linears have 2 DOFs per element
	  case EDGE2:
	  case EDGE3:
	  case EDGE4:
	    return 2;

	    // 2D linears have 3 DOFs per element
	  case TRI3:
	  case TRI6:
	  case QUAD4:
	  case QUAD8:
	  case QUAD9:
	    return 3;

	    // 3D linears have 4 DOFs per element
 	  case TET4:
	  case TET10:
	  case HEX8:
	  case HEX20:
	  case HEX27:
	  case PRISM6:
	  case PRISM15:
	  case PRISM18:
	  case PYRAMID5:
	    return 4;
	    
	  default:
	    {
#ifdef DEBUG
	      std::cerr << "ERROR: Bad ElemType = " << t
			<< " for " << o << "th order approximation!" 
			<< std::endl;
#endif
	      libmesh_error();	    
	    }
	  }
      }


      // Discontinuous quadratic shape functions
      // expressed in the XYZ monomials.
    case SECOND:
      {
	switch (t)
	  {
	    // 1D quadratics have 3 DOFs per element
	  case EDGE2:
	  case EDGE3:
	  case EDGE4:
	    return 3;

	    // 2D quadratics have 6 DOFs per element
	  case TRI3:
	  case TRI6:
	  case QUAD4:
	  case QUAD8:
	  case QUAD9:
	    return 6;

	    // 3D quadratics have 10 DOFs per element
	  case TET4:
	  case TET10:
	  case HEX8:
	  case HEX20:
	  case HEX27:
	  case PRISM6:
	  case PRISM15:
	  case PRISM18:
	  case PYRAMID5:
	    return 10;
	    
	  default:
	    {
#ifdef DEBUG
	      std::cerr << "ERROR: Bad ElemType = " << t
			<< " for " << o << "th order approximation!" 
			<< std::endl;
#endif
	      libmesh_error();	    
	    }
	  }
      }


      // Discontinuous cubic shape functions
      // expressed in the XYZ monomials.
    case THIRD:
      {
	switch (t)
	  {
	  case EDGE2:
	  case EDGE3:
	  case EDGE4:
	    return 4;

	  case TRI3:
	  case TRI6:
	  case QUAD4:
	  case QUAD8:
	  case QUAD9:
	    return 10;

	  case TET4:
	  case TET10:
	  case HEX8:
	  case HEX20:
	  case HEX27:
	  case PRISM6:
	  case PRISM15:
	  case PRISM18:
	  case PYRAMID5:
	    return 20;
	    
	  default:
	    {
#ifdef DEBUG
	      std::cerr << "ERROR: Bad ElemType = " << t
			<< " for " << o << "th order approximation!" 
			<< std::endl;
#endif
	      libmesh_error();	    
	    }
	  }
      }


      // Discontinuous quartic shape functions
      // expressed in the XYZ monomials.
    case FOURTH:
      {
	switch (t)
	  {
	  case EDGE2:
	  case EDGE3:
	  case EDGE4:
	    return 5;

	  case TRI3:
	  case TRI6:
	  case QUAD4:
	  case QUAD8:
	  case QUAD9:
	    return 15;

	  case TET4:
	  case TET10:
	  case HEX8:
	  case HEX20:
	  case HEX27:
	  case PRISM6:
	  case PRISM15:
	  case PRISM18:
	  case PYRAMID5:
	    return 35;
	    
	  default:
	    {
#ifdef DEBUG
	      std::cerr << "ERROR: Bad ElemType = " << t
			<< " for " << o << "th order approximation!" 
			<< std::endl;
#endif
	      libmesh_error();	    
	    }
	  }
      }
      
    default:
      {
        const unsigned int order = static_cast<unsigned int>(o);
	switch (t)
	  {
	  case EDGE2:
	  case EDGE3:
	    return (order+1);

	  case TRI3:
	  case TRI6:
	  case QUAD4:
	  case QUAD8:
	  case QUAD9:
	    return (order+1)*(order+2)/2;

	  case TET4:
	  case TET10:
	  case HEX8:
	  case HEX20:
	  case HEX27:
	  case PRISM6:
	  case PRISM15:
	  case PRISM18:
	  case PYRAMID5:
	    return (order+1)*(order+2)*(order+3)/6;
	    
	  default:
	    {
#ifdef DEBUG
	      std::cerr << "ERROR: Bad ElemType = " << t
			<< " for " << o << "th order approximation!" 
			<< std::endl;
#endif
	      libmesh_error();	    
	    }
	  }
      }
      return 0;
    }
}



template <unsigned int Dim, FEFamily T>
FEContinuity FE<Dim,T>::get_continuity() const
{
  return DISCONTINUOUS;
}



template <unsigned int Dim, FEFamily T>
bool FE<Dim,T>::is_hierarchic() const
{
  return true;
}



#ifdef ENABLE_AMR
template <unsigned int Dim, FEFamily T>
void FE<Dim,T>::compute_constraints (DofConstraints &,
				     DofMap &,
				     const unsigned int,
				     const Elem*)
{
  // Monomials are discontinuous...  No constraints.
  return;
}
#endif // #ifdef ENABLE_AMR



template <unsigned int Dim, FEFamily T>
bool FE<Dim,T>::shapes_need_reinit() const
{
  return false;
}



//--------------------------------------------------------------
// Explicit instantiation of member functions
INSTANTIATE_MBRF(1,XYZ);
INSTANTIATE_MBRF(2,XYZ);
INSTANTIATE_MBRF(3,XYZ);
template void  FEXYZ<1>::init_shape_functions(const std::vector<Point>&,
					    const Elem*);
template void  FEXYZ<2>::init_shape_functions(const std::vector<Point>&,
					    const Elem*);
template void  FEXYZ<3>::init_shape_functions(const std::vector<Point>&,
					    const Elem*);
template void  FEXYZ<1>::compute_shape_functions(const Elem*);
template void  FEXYZ<2>::compute_shape_functions(const Elem*);
template void  FEXYZ<3>::compute_shape_functions(const Elem*);

#ifdef ENABLE_AMR
template void FE<2,XYZ>::compute_constraints(DofConstraints&, DofMap&, 
					     const unsigned int, const Elem*);
template void FE<3,XYZ>::compute_constraints(DofConstraints&, DofMap&, 
					     const unsigned int, const Elem*);
#endif // #ifdef ENABLE_AMR
