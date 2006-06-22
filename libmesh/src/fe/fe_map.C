// $Id: fe_map.C,v 1.40 2006-06-22 14:46:08 benkirk Exp $

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
#include <cmath> // for std::sqrt, std::abs


// Local includes
#include "fe.h"
#include "elem.h"
#include "libmesh_logging.h"
#include "fe_macro.h"




void FEBase::compute_affine_map(const std::vector<Real>& qw,
			        const Elem* elem)
{
   // Start logging the map computation.
  START_LOG("compute_affine_map()", "FE");  

  assert (elem  != NULL);

  const unsigned int        n_qp = qw.size();

  switch (this->dim)
    {
      //--------------------------------------------------------------------
      // 1D
    case 1:
      {
	//------------------------------------------------------------------
	// Compute the values at the quadrature points,
	// the constant Jacobian
	
	// Resize the vectors to hold data at the quadrature points
	{  
	  xyz.resize(n_qp);
	  dxyzdxi_map.resize(n_qp);
	  dxidx_map.resize(n_qp);
#ifdef ENABLE_SECOND_DERIVATIVES
	  d2xyzdxi2_map.resize(n_qp);
#endif
	  
	  JxW.resize(n_qp);
	}
	
	// Clear the entities that will be summed
	for (unsigned int p=0; p<n_qp; p++)
	  xyz[p].zero();
	dxyzdxi_map[0].zero();
#ifdef ENABLE_SECOND_DERIVATIVES
	d2xyzdxi2_map[0].zero();
#endif
	
	
	// compute x at the quadrature points, dxdxi once
	for (unsigned int i=0; i<phi_map.size(); i++) // sum over the nodes
	  {
	    // Reference to the point, helps eliminate
	    // exessive temporaries in the inner loop
	    const Point& elem_point = elem->point(i);
	    
	    for (unsigned int p=0; p<n_qp; p++) // for each quadrature point
	      xyz[p].add_scaled        (elem_point, phi_map[i][p]    );
	      
	    dxyzdxi_map[0].add_scaled(elem_point, dphidxi_map[i][0]);
#ifdef ENABLE_SECOND_DERIVATIVES
	    d2xyzdxi2_map[0].add_scaled(elem_point, d2phidxi2_map[i][0]);
#endif
	  }

	for (unsigned int p=1; p<n_qp; p++) // for each extra quadrature point
          {
	    dxyzdxi_map[p] = dxyzdxi_map[0];
#ifdef ENABLE_SECOND_DERIVATIVES
	    d2xyzdxi2_map[p] = d2xyzdxi2_map[0];
#endif
          }

	/*
        // Test the inverse map
	for (unsigned int p=0; p<n_qp; p++)
	{
	const Point p_inv = inverse_map (elem, xyz[p]);
	    
	std::cout << "qp[p]   = ";
	qrule->qp(p).print();
	std::cout << "inv_map = ";
	p_inv.print();
	}
	*/

	// compute the jacobian once

	// Symbolically, the matrix determinant is
	//
	// jac = | dx/dxi | = dx/dxi
	//         
	
	// Compute the Jacobian.  This assumes the
	// 1D edge lives in 1D space.
	const Real jac = dxdxi_map(0);
	    
	if (jac <= 0.)
	  {
	    std::cerr << "ERROR: negative Jacobian: "
		      << jac
		      << std::endl;
	    error();
	  }
	    
	assert (dxdxi_map(0) != 0.);
	    
	dxidx_map[0] = 1./dxdxi_map(0);
	JxW[0] = jac*qw[0];
 
	for (unsigned int p=1; p<n_qp; p++)
          {
	    dxidx_map[p] = dxidx_map[0];
	    JxW[p] = jac*qw[p];
          }

	// done computing the map
	break;
      }

      
      //--------------------------------------------------------------------
      // 2D
    case 2:
      {
	//------------------------------------------------------------------
	// Compute the (x,y) values at the quadrature points,
	// the Jacobian at the quadrature points

	// Resize the vectors to hold data at the quadrature points
	{  
	  xyz.resize(n_qp);
	  dxyzdxi_map.resize(n_qp);
	  dxyzdeta_map.resize(n_qp);
	  dxidx_map.resize(n_qp);
	  dxidy_map.resize(n_qp);
	  dxidz_map.resize(n_qp);
	  detadx_map.resize(n_qp);
	  detady_map.resize(n_qp);
	  detadz_map.resize(n_qp);
#ifdef ENABLE_SECOND_DERIVATIVES
	  d2xyzdxi2_map.resize(n_qp);
	  d2xyzdxideta_map.resize(n_qp);
	  d2xyzdeta2_map.resize(n_qp);
#endif
	  
	  JxW.resize(n_qp);
	}
	
	// Clear the entities that will be summed
	for (unsigned int p=0; p<n_qp; p++)
	  xyz[p].zero();

	dxyzdxi_map[0].zero();
	dxyzdeta_map[0].zero();
#ifdef ENABLE_SECOND_D0RIVATIVES
	d2xyzdxi2_map[0].zero();
	d2xyzdxideta_map[0].zero();
	d2xyzdeta2_map[0].zero();
#endif
	
	
	// compute (x,y) at the quadrature points, derivatives once
	for (unsigned int i=0; i<phi_map.size(); i++) // sum over the nodes
	  {
	    // Reference to the point, helps eliminate
	    // exessive temporaries in the inner loop
	    const Point& elem_point = elem->point(i);
	    
	    for (unsigned int p=0; p<n_qp; p++) // for each quadrature point
	      xyz[p].add_scaled          (elem_point, phi_map[i][p]     );

	    dxyzdxi_map[0].add_scaled  (elem_point, dphidxi_map[i][0] );
	    dxyzdeta_map[0].add_scaled (elem_point, dphideta_map[i][0]);
#ifdef ENABLE_SECOND_DERIVATIVES
	    d2xyzdxi2_map[0].add_scaled    (elem_point,
					    d2phidxi2_map[i][0]);
	    d2xyzdxideta_map[0].add_scaled (elem_point,
					    d2phidxideta_map[i][0]);
	    d2xyzdeta2_map[0].add_scaled   (elem_point,
					    d2phideta2_map[i][0]);
#endif
	  }

	for (unsigned int p=1; p<n_qp; p++) // for each extra quadrature point
          {
	    dxyzdxi_map[p] = dxyzdxi_map[0];
	    dxyzdeta_map[p] = dxyzdeta_map[0];
#ifdef ENABLE_SECOND_DERIVATIVES
	    d2xyzdxi2_map[p] = d2xyzdxi2_map[0];
	    d2xyzdxideta_map[p] = d2xyzdxideta_map[0];
	    d2xyzdeta2_map[p] = d2xyzdeta2_map[0];
#endif
          }
	
	/*
        // Test the inverse map
	for (unsigned int p=0; p<n_qp; p++)
	{
	const Point p_inv = inverse_map (elem, xyz[p]);
	
	std::cout << "qp[p]   = ";
	qrule->qp(p).print();
	std::cout << "inv_map = ";
	p_inv.print();
	}
	*/
	
	// compute the jacobian once
	const Real dx_dxi = dxdxi_map(0), dx_deta = dxdeta_map(0),
	           dy_dxi = dydxi_map(0), dy_deta = dydeta_map(0),
	           dz_dxi = dzdxi_map(0), dz_deta = dzdeta_map(0);

#if DIM == 2
	// Compute the Jacobian.  This assumes the 2D face
	// lives in 2D space
	//
	// Symbolically, the matrix determinant is
	//
	//         | dx/dxi  dx/deta |
	// jac =   | dy/dxi  dy/deta |
	//         
	// jac = dx/dxi*dy/deta - dx/deta*dy/dxi 
	const Real jac = (dx_dxi*dy_deta - dx_deta*dy_dxi);
	    
	if (jac <= 0.)
	  {
	    std::cerr << "ERROR: negative Jacobian: "
		      << jac
		      << std::endl;
	    error();
	  }
	    
	JxW[0] = jac*qw[0];
	    
	// Compute the shape function derivatives wrt x,y at the
	// quadrature points
	const Real inv_jac = 1./jac;
	    
	dxidx_map[0]  =  dy_deta*inv_jac; //dxi/dx  =  (1/J)*dy/deta
	dxidy_map[0]  = -dx_deta*inv_jac; //dxi/dy  = -(1/J)*dx/deta
	detadx_map[0] = -dy_dxi* inv_jac; //deta/dx = -(1/J)*dy/dxi
	detady_map[0] =  dx_dxi* inv_jac; //deta/dy =  (1/J)*dx/dxi

	dxidz_map[0] = detadz_map[0] = 0.;
#else
	// Compute the Jacobian.  This assumes a 2D face in
	// 3D space.
	//
	// The transformation matrix T from local to global
	// coordinates is
	//
	//         | dx/dxi  dx/deta |
	//     T = | dy/dxi  dy/deta |
	//         | dz/dxi  dz/deta |
	// note det(T' T) = det(T')det(T) = det(T)det(T)
	// so det(T) = std::sqrt(det(T' T))
	//
	//----------------------------------------------
	// Notes:
	//
	//       dX = R dXi -> R'dX = R'R dXi
	// (R^-1)dX =   dXi    [(R'R)^-1 R']dX = dXi 
	//
	// so R^-1 = (R'R)^-1 R'
	//
	// and R^-1 R = (R'R)^-1 R'R = I.
	//
	const Real g11 = (dx_dxi*dx_dxi +
			  dy_dxi*dy_dxi +
			  dz_dxi*dz_dxi);
	    
	const Real g12 = (dx_dxi*dx_deta +
			  dy_dxi*dy_deta +
			  dz_dxi*dz_deta);
	    
	const Real g21 = g12;
	    
	const Real g22 = (dx_deta*dx_deta +
			  dy_deta*dy_deta +
			  dz_deta*dz_deta);

	const Real det = (g11*g22 - g12*g21);

	if (det <= 0.)
	  {
	    std::cerr << "ERROR: negative Jacobian! "
		      << std::endl;
	    error();
	  }
	      
	const Real inv_det = 1./det;
	const Real jac = std::sqrt(det);
	    
	JxW[0] = jac*qw[0];

	const Real g11inv =  g22*inv_det;
	const Real g12inv = -g12*inv_det;
	const Real g21inv = -g21*inv_det;
	const Real g22inv =  g11*inv_det;

	dxidx_map[0]  = g11inv*dx_dxi + g12inv*dx_deta;
	dxidy_map[0]  = g11inv*dy_dxi + g12inv*dy_deta;
	dxidz_map[0]  = g11inv*dz_dxi + g12inv*dz_deta;
	    
	detadx_map[0] = g21inv*dx_dxi + g22inv*dx_deta;
	detady_map[0] = g21inv*dy_dxi + g22inv*dy_deta;
	detadz_map[0] = g21inv*dz_dxi + g22inv*dz_deta;
	    	    	    
#endif
	for (unsigned int p=1; p<n_qp; p++) // for each extra quadrature point
          {
	    JxW[p] = jac*qw[p];

	    dxidx_map[p]  = dxidx_map[0];
	    dxidy_map[p]  = dxidy_map[0];
	    dxidz_map[p]  = dxidz_map[0];
	    
	    detadx_map[p] = detadx_map[0];
	    detady_map[p] = detady_map[0];
	    detadz_map[p] = detadz_map[0];
	  }
       
	// done computing the map
	break;
      }


      
      //--------------------------------------------------------------------
      // 3D
    case 3:
      {
	//------------------------------------------------------------------
	// Compute the (x,y,z) values at the quadrature points,
	// the Jacobian at the quadrature points

	// Resize the vectors to hold data at the quadrature points
	{  
	  xyz.resize           (n_qp);
	  dxyzdxi_map.resize   (n_qp);
	  dxyzdeta_map.resize  (n_qp);
	  dxyzdzeta_map.resize (n_qp);
	  dxidx_map.resize     (n_qp);
	  dxidy_map.resize     (n_qp);
	  dxidz_map.resize     (n_qp);
	  detadx_map.resize    (n_qp);
	  detady_map.resize    (n_qp);
	  detadz_map.resize    (n_qp);
	  dzetadx_map.resize   (n_qp);
	  dzetady_map.resize   (n_qp);
	  dzetadz_map.resize   (n_qp);
#ifdef ENABLE_SECOND_DERIVATIVES
	  d2xyzdxi2_map.resize(n_qp);
	  d2xyzdxideta_map.resize(n_qp);
	  d2xyzdxidzeta_map.resize(n_qp);
	  d2xyzdeta2_map.resize(n_qp);
	  d2xyzdetadzeta_map.resize(n_qp);
	  d2xyzdzeta2_map.resize(n_qp);
#endif
	  
	  JxW.resize (n_qp);
	}
    
	// Clear the entities that will be summed
	for (unsigned int p=0; p<n_qp; p++)
	  xyz[p].zero           ();
	dxyzdxi_map[0].zero   ();
	dxyzdeta_map[0].zero  ();
	dxyzdzeta_map[0].zero ();
#ifdef ENABLE_SECOND_DERIVATIVES
	d2xyzdxi2_map[0].zero();
	d2xyzdxideta_map[0].zero();
	d2xyzdxidzeta_map[0].zero();
	d2xyzdeta2_map[0].zero();
	d2xyzdetadzeta_map[0].zero();
	d2xyzdzeta2_map[0].zero();
#endif
	
	
	// compute (x,y,z) at the quadrature points,
        // dxdxi,   dydxi,   dzdxi,
	// dxdeta,  dydeta,  dzdeta,
	// dxdzeta, dydzeta, dzdzeta  all once
	for (unsigned int i=0; i<phi_map.size(); i++) // sum over the nodes
	  {
	    // Reference to the point, helps eliminate
	    // exessive temporaries in the inner loop
	    const Point& elem_point = elem->point(i);
	    
	    for (unsigned int p=0; p<n_qp; p++) // for each quadrature point
	      xyz[p].add_scaled           (elem_point, phi_map[i][p]      );
	    dxyzdxi_map[0].add_scaled   (elem_point, dphidxi_map[i][0]  );
	    dxyzdeta_map[0].add_scaled  (elem_point, dphideta_map[i][0] );
	    dxyzdzeta_map[0].add_scaled (elem_point, dphidzeta_map[i][0]);
#ifdef ENABLE_SECOND_DERIVATIVES
	    d2xyzdxi2_map[0].add_scaled      (elem_point,
					      d2phidxi2_map[i][0]);
	    d2xyzdxideta_map[0].add_scaled   (elem_point,
					      d2phidxideta_map[i][0]);
	    d2xyzdxidzeta_map[0].add_scaled  (elem_point,
					      d2phidxidzeta_map[i][0]);
	    d2xyzdeta2_map[0].add_scaled     (elem_point,
					      d2phideta2_map[i][0]);
	    d2xyzdetadzeta_map[0].add_scaled (elem_point,
					      d2phidetadzeta_map[i][0]);
	    d2xyzdzeta2_map[0].add_scaled    (elem_point,
					      d2phidzeta2_map[i][0]);
#endif
	    for (unsigned int p=0; p<n_qp; p++) // for each extra quadrature point
              {
	        dxyzdxi_map[p] = dxyzdxi_map[0];
	        dxyzdeta_map[p] = dxyzdeta_map[0];
	        dxyzdzeta_map[p] = dxyzdzeta_map[0];
#ifdef ENABLE_SECOND_DERIVATIVES
	        d2xyzdxi2_map[p] = d2xyzdxi2_map[0];
	        d2xyzdxideta_map[p] = d2xyzdxideta_map[0];
	        d2xyzdxidzeta_map[p] = d2xyzdxidzeta_map[0];
	        d2xyzdeta2_map[p] = d2xyzdeta2_map[0];
	        d2xyzdetadzeta_map[p] = d2xyzdetadzeta_map[0];
	        d2xyzdzeta2_map[p] = d2xyzdzeta2_map[0];
#endif
	      }
	  }
	
	/*
        // Test the inverse map
	for (unsigned int p=0; p<n_qp; p++)
	{
	const Point p_inv = inverse_map (elem, xyz[p]);
	    
	std::cout << "qp[p]   = ";
	qrule->qp(p).print();
	std::cout << "inv_map = ";
	p_inv.print();
	}
	*/
	  
	// compute the jacobian once
	const Real
	  dx_dxi   = dxdxi_map(0),   dy_dxi   = dydxi_map(0),   dz_dxi   = dzdxi_map(0),
	  dx_deta  = dxdeta_map(0),  dy_deta  = dydeta_map(0),  dz_deta  = dzdeta_map(0),
	  dx_dzeta = dxdzeta_map(0), dy_dzeta = dydzeta_map(0), dz_dzeta = dzdzeta_map(0);
	    
	// Symbolically, the matrix determinant is
	//
	//         | dx/dxi   dy/dxi   dz/dxi   |
	// jac =   | dx/deta  dy/deta  dz/deta  |
	//         | dx/dzeta dy/dzeta dz/dzeta |
	// 
	// jac = dx/dxi*(dy/deta*dz/dzeta - dz/deta*dy/dzeta) +
	//       dy/dxi*(dz/deta*dx/dzeta - dx/deta*dz/dzeta) +
	//       dz/dxi*(dx/deta*dy/dzeta - dy/deta*dx/dzeta)

	const Real jac = (dx_dxi*(dy_deta*dz_dzeta - dz_deta*dy_dzeta)  +
			  dy_dxi*(dz_deta*dx_dzeta - dx_deta*dz_dzeta)  +
			  dz_dxi*(dx_deta*dy_dzeta - dy_deta*dx_dzeta));
	    
	if (jac <= 0.)
	  {
	    std::cerr << "ERROR: negative Jacobian: "
		      << jac
		      << std::endl;
	    error();
	  }

	JxW[0] = jac*qw[0];
	    
	    // Compute the shape function derivatives wrt x,y at the
	    // quadrature points
	const Real inv_jac  = 1./jac;	    
	    
	dxidx_map[0]   = (dy_deta*dz_dzeta - dz_deta*dy_dzeta)*inv_jac;
	dxidy_map[0]   = (dz_deta*dx_dzeta - dx_deta*dz_dzeta)*inv_jac;
	dxidz_map[0]   = (dx_deta*dy_dzeta - dy_deta*dx_dzeta)*inv_jac;
	    
	detadx_map[0]  = (dz_dxi*dy_dzeta  - dy_dxi*dz_dzeta )*inv_jac;
	detady_map[0]  = (dx_dxi*dz_dzeta  - dz_dxi*dx_dzeta )*inv_jac;
	detadz_map[0]  = (dy_dxi*dx_dzeta  - dx_dxi*dy_dzeta )*inv_jac;
	    
	dzetadx_map[0] = (dy_dxi*dz_deta   - dz_dxi*dy_deta  )*inv_jac;
	dzetady_map[0] = (dz_dxi*dx_deta   - dx_dxi*dz_deta  )*inv_jac;
	dzetadz_map[0] = (dx_dxi*dy_deta   - dy_dxi*dx_deta  )*inv_jac;
	
	for (unsigned int p=0; p<n_qp; p++)
	  {
	    JxW[p] = jac*qw[p];
 
	    dxidx_map[p]   = dxidx_map[0];
	    dxidy_map[p]   = dxidy_map[0];
	    dxidz_map[p]   = dxidz_map[0];
   
	    detadx_map[p]  = detadx_map[0];
	    detady_map[p]  = detady_map[0];
	    detadz_map[p]  = detadz_map[0];

	    dzetadx_map[p] = dzetadx_map[0];
	    dzetady_map[p] = dzetady_map[0];
	    dzetadz_map[p] = dzetadz_map[0];
	  }
	// done computing the map
	break;
      }



    default:
      error();
    }
  
  STOP_LOG("compute_affine_map()", "FE");  
}



void FEBase::compute_map(const std::vector<Real>& qw,
			 const Elem* elem)
{
  if (elem->has_affine_map())
    {
      compute_affine_map(qw, elem);
      return;
    }
  
   // Start logging the map computation.
  START_LOG("compute_map()", "FE");

  assert (elem  != NULL);
  
  const unsigned int        n_qp = qw.size();

  switch (this->dim)
    {
      //--------------------------------------------------------------------
      // 1D
    case 1:
      {
	//------------------------------------------------------------------
	// Compute the values at the quadrature points,
	// the Jacobian at the quadrature points
	
	// Resize the vectors to hold data at the quadrature points
	{  
	  xyz.resize(n_qp);
	  dxyzdxi_map.resize(n_qp);
	  dxidx_map.resize(n_qp);
#ifdef ENABLE_SECOND_DERIVATIVES
	  d2xyzdxi2_map.resize(n_qp);
#endif
	  
	  JxW.resize(n_qp);
	}
	
	// Clear the entities that will be summed
	for (unsigned int p=0; p<n_qp; p++)
	  {
	    xyz[p].zero();
	    dxyzdxi_map[p].zero();
#ifdef ENABLE_SECOND_DERIVATIVES
	    d2xyzdxi2_map[p].zero();
#endif
	  }
	
	
	// compute x, dxdxi at the quadrature points    
	for (unsigned int i=0; i<phi_map.size(); i++) // sum over the nodes
	  {
	    // Reference to the point, helps eliminate
	    // exessive temporaries in the inner loop
	    const Point& elem_point = elem->point(i);
	    
	    for (unsigned int p=0; p<n_qp; p++) // for each quadrature point
	      {	  
		xyz[p].add_scaled        (elem_point, phi_map[i][p]    );
		dxyzdxi_map[p].add_scaled(elem_point, dphidxi_map[i][p]);
#ifdef ENABLE_SECOND_DERIVATIVES
		d2xyzdxi2_map[p].add_scaled(elem_point,
					    d2phidxi2_map[i][p]);
#endif
	      }
	  }

	/*
        // Test the inverse map
	for (unsigned int p=0; p<n_qp; p++)
	{
	const Point p_inv = inverse_map (elem, xyz[p]);
	    
	std::cout << "qp[p]   = ";
	qrule->qp(p).print();
	std::cout << "inv_map = ";
	p_inv.print();
	}
	*/

	// compute the jacobian at the quadrature points
	for (unsigned int p=0; p<n_qp; p++)
	  {
	    // Symbolically, the matrix determinant is
	    //
	    // jac = | dx/dxi | = dx/dxi
	    //         
	    
	    // Compute the Jacobian.  This assumes the
	    // 1D edge lives in 1D space.
	    const Real jac = dxdxi_map(p);
	    
	    if (jac <= 0.)
	      {
		std::cerr << "ERROR: negative Jacobian: "
			  << jac
			  << std::endl;
		error();
	      }
	    
	    assert (dxdxi_map(p) != 0.);
	    
	    dxidx_map[p] = 1./dxdxi_map(p);
	    
	    JxW[p] = jac*qw[p];
	  }

	// done computing the map
	break;
      }

      
      //--------------------------------------------------------------------
      // 2D
    case 2:
      {
	//------------------------------------------------------------------
	// Compute the (x,y) values at the quadrature points,
	// the Jacobian at the quadrature points

	// Resize the vectors to hold data at the quadrature points
	{  
	  xyz.resize(n_qp);
	  dxyzdxi_map.resize(n_qp);
	  dxyzdeta_map.resize(n_qp);
	  dxidx_map.resize(n_qp);
	  dxidy_map.resize(n_qp);
	  dxidz_map.resize(n_qp);
	  detadx_map.resize(n_qp);
	  detady_map.resize(n_qp);
	  detadz_map.resize(n_qp);
#ifdef ENABLE_SECOND_DERIVATIVES
	  d2xyzdxi2_map.resize(n_qp);
	  d2xyzdxideta_map.resize(n_qp);
	  d2xyzdeta2_map.resize(n_qp);
#endif
	  
	  JxW.resize(n_qp);
	}
	
	// Clear the entities that will be summed
	for (unsigned int p=0; p<n_qp; p++)
	  {
	    xyz[p].zero();
	    dxyzdxi_map[p].zero();
	    dxyzdeta_map[p].zero();
#ifdef ENABLE_SECOND_DERIVATIVES
	    d2xyzdxi2_map[p].zero();
	    d2xyzdxideta_map[p].zero();
	    d2xyzdeta2_map[p].zero();
#endif
	  }
	
	
	// compute (x,y), dxdxi, dydxi, dxdeta, dydeta at the quadrature points
	for (unsigned int i=0; i<phi_map.size(); i++) // sum over the nodes
	  {
	    // Reference to the point, helps eliminate
	    // exessive temporaries in the inner loop
	    const Point& elem_point = elem->point(i);
	    
	    for (unsigned int p=0; p<n_qp; p++) // for each quadrature point
	      {	  
		xyz[p].add_scaled          (elem_point, phi_map[i][p]     );
		dxyzdxi_map[p].add_scaled  (elem_point, dphidxi_map[i][p] );
		dxyzdeta_map[p].add_scaled (elem_point, dphideta_map[i][p]);
#ifdef ENABLE_SECOND_DERIVATIVES
		d2xyzdxi2_map[p].add_scaled    (elem_point,
						d2phidxi2_map[i][p]);
		d2xyzdxideta_map[p].add_scaled (elem_point,
						d2phidxideta_map[i][p]);
		d2xyzdeta2_map[p].add_scaled   (elem_point,
						d2phideta2_map[i][p]);
#endif
	      }
	  }
	
	/*
        // Test the inverse map
	for (unsigned int p=0; p<n_qp; p++)
	{
	const Point p_inv = inverse_map (elem, xyz[p]);
	
	std::cout << "qp[p]   = ";
	qrule->qp(p).print();
	std::cout << "inv_map = ";
	p_inv.print();
	}
	*/
	
	// compute the jacobian at the quadrature points
	for (unsigned int p=0; p<n_qp; p++)
	  {
	    const Real
	      dx_dxi = dxdxi_map(p), dx_deta = dxdeta_map(p),
	      dy_dxi = dydxi_map(p), dy_deta = dydeta_map(p),
	      dz_dxi = dzdxi_map(p), dz_deta = dzdeta_map(p);

#if DIM == 2
	    // Compute the Jacobian.  This assumes the 2D face
	    // lives in 2D space
	    //
	    // Symbolically, the matrix determinant is
	    //
	    //         | dx/dxi  dx/deta |
	    // jac =   | dy/dxi  dy/deta |
	    //         
	    // jac = dx/dxi*dy/deta - dx/deta*dy/dxi 
	    const Real jac = (dx_dxi*dy_deta - dx_deta*dy_dxi);
	    
	    if (jac <= 0.)
	      {
		std::cerr << "ERROR: negative Jacobian: "
			  << jac
			  << std::endl;
		error();
	      }
	    
	    JxW[p] = jac*qw[p];
	    
	    // Compute the shape function derivatives wrt x,y at the
	    // quadrature points
	    const Real
	      inv_jac = 1./jac;
	    
	    dxidx_map[p]  =  dy_deta*inv_jac; //dxi/dx  =  (1/J)*dy/deta
	    dxidy_map[p]  = -dx_deta*inv_jac; //dxi/dy  = -(1/J)*dx/deta
	    detadx_map[p] = -dy_dxi* inv_jac; //deta/dx = -(1/J)*dy/dxi
	    detady_map[p] =  dx_dxi* inv_jac; //deta/dy =  (1/J)*dx/dxi

	    dxidz_map[p] = detadz_map[p] = 0.;
#else
	    // Compute the Jacobian.  This assumes a 2D face in
	    // 3D space.
	    //
	    // The transformation matrix T from local to global
	    // coordinates is
	    //
	    //         | dx/dxi  dx/deta |
	    //     T = | dy/dxi  dy/deta |
	    //         | dz/dxi  dz/deta |
	    // note det(T' T) = det(T')det(T) = det(T)det(T)
	    // so det(T) = std::sqrt(det(T' T))
	    //
	    //----------------------------------------------
	    // Notes:
	    //
	    //       dX = R dXi -> R'dX = R'R dXi
	    // (R^-1)dX =   dXi    [(R'R)^-1 R']dX = dXi 
	    //
	    // so R^-1 = (R'R)^-1 R'
	    //
	    // and R^-1 R = (R'R)^-1 R'R = I.
	    //
	    const Real g11 = (dx_dxi*dx_dxi +
			      dy_dxi*dy_dxi +
			      dz_dxi*dz_dxi);
	    
	    const Real g12 = (dx_dxi*dx_deta +
			      dy_dxi*dy_deta +
			      dz_dxi*dz_deta);
	    
	    const Real g21 = g12;
	    
	    const Real g22 = (dx_deta*dx_deta +
			      dy_deta*dy_deta +
			      dz_deta*dz_deta);

	    const Real det = (g11*g22 - g12*g21);

	    if (det <= 0.)
	      {
		std::cerr << "ERROR: negative Jacobian! "
			  << std::endl;
		error();
	      }
	      
	    const Real inv_det = 1./det;
	    const Real jac = std::sqrt(det);
	    
	    JxW[p] = jac*qw[p];

	    const Real g11inv =  g22*inv_det;
	    const Real g12inv = -g12*inv_det;
	    const Real g21inv = -g21*inv_det;
	    const Real g22inv =  g11*inv_det;

	    dxidx_map[p]  = g11inv*dx_dxi + g12inv*dx_deta;
	    dxidy_map[p]  = g11inv*dy_dxi + g12inv*dy_deta;
	    dxidz_map[p]  = g11inv*dz_dxi + g12inv*dz_deta;
	    
	    detadx_map[p] = g21inv*dx_dxi + g22inv*dx_deta;
	    detady_map[p] = g21inv*dy_dxi + g22inv*dy_deta;
	    detadz_map[p] = g21inv*dz_dxi + g22inv*dz_deta;
	    	    	    
#endif
	  }
       
	// done computing the map
	break;
      }


      
      //--------------------------------------------------------------------
      // 3D
    case 3:
      {
	//------------------------------------------------------------------
	// Compute the (x,y,z) values at the quadrature points,
	// the Jacobian at the quadrature points

	// Resize the vectors to hold data at the quadrature points
	{  
	  xyz.resize           (n_qp);
	  dxyzdxi_map.resize   (n_qp);
	  dxyzdeta_map.resize  (n_qp);
	  dxyzdzeta_map.resize (n_qp);
	  dxidx_map.resize     (n_qp);
	  dxidy_map.resize     (n_qp);
	  dxidz_map.resize     (n_qp);
	  detadx_map.resize    (n_qp);
	  detady_map.resize    (n_qp);
	  detadz_map.resize    (n_qp);
	  dzetadx_map.resize   (n_qp);
	  dzetady_map.resize   (n_qp);
	  dzetadz_map.resize   (n_qp);
#ifdef ENABLE_SECOND_DERIVATIVES
	  d2xyzdxi2_map.resize(n_qp);
	  d2xyzdxideta_map.resize(n_qp);
	  d2xyzdxidzeta_map.resize(n_qp);
	  d2xyzdeta2_map.resize(n_qp);
	  d2xyzdetadzeta_map.resize(n_qp);
	  d2xyzdzeta2_map.resize(n_qp);
#endif
	  
	  JxW.resize (n_qp);
	}
    
	// Clear the entities that will be summed
	for (unsigned int p=0; p<n_qp; p++)
	  {
	    xyz[p].zero           ();
	    dxyzdxi_map[p].zero   ();
	    dxyzdeta_map[p].zero  ();
	    dxyzdzeta_map[p].zero ();
#ifdef ENABLE_SECOND_DERIVATIVES
	    d2xyzdxi2_map[p].zero();
	    d2xyzdxideta_map[p].zero();
	    d2xyzdxidzeta_map[p].zero();
	    d2xyzdeta2_map[p].zero();
	    d2xyzdetadzeta_map[p].zero();
	    d2xyzdzeta2_map[p].zero();
#endif
	  }
	
	
	// compute (x,y,z), dxdxi,   dydxi,   dzdxi,
	//                  dxdeta,  dydeta,  dzdeta,
	//                  dxdzeta, dydzeta, dzdzeta
	// at the quadrature points    
	for (unsigned int i=0; i<phi_map.size(); i++) // sum over the nodes
	  {
	    // Reference to the point, helps eliminate
	    // exessive temporaries in the inner loop
	    const Point& elem_point = elem->point(i);
	    
	    for (unsigned int p=0; p<n_qp; p++) // for each quadrature point
	      {	  
		xyz[p].add_scaled           (elem_point, phi_map[i][p]      );
		dxyzdxi_map[p].add_scaled   (elem_point, dphidxi_map[i][p]  );
		dxyzdeta_map[p].add_scaled  (elem_point, dphideta_map[i][p] );
		dxyzdzeta_map[p].add_scaled (elem_point, dphidzeta_map[i][p]);
#ifdef ENABLE_SECOND_DERIVATIVES
	        d2xyzdxi2_map[p].add_scaled      (elem_point,
					       d2phidxi2_map[i][p]);
	        d2xyzdxideta_map[p].add_scaled   (elem_point,
					       d2phidxideta_map[i][p]);
	        d2xyzdxidzeta_map[p].add_scaled  (elem_point,
					       d2phidxidzeta_map[i][p]);
	        d2xyzdeta2_map[p].add_scaled     (elem_point,
					       d2phideta2_map[i][p]);
	        d2xyzdetadzeta_map[p].add_scaled (elem_point,
					       d2phidetadzeta_map[i][p]);
	        d2xyzdzeta2_map[p].add_scaled    (elem_point,
					       d2phidzeta2_map[i][p]);
#endif
	      }
	  }
	
	/*
        // Test the inverse map
	for (unsigned int p=0; p<n_qp; p++)
	{
	const Point p_inv = inverse_map (elem, xyz[p]);
	    
	std::cout << "qp[p]   = ";
	qrule->qp(p).print();
	std::cout << "inv_map = ";
	p_inv.print();
	}
	*/
	  
	// compute the jacobian at the quadrature points
	for (unsigned int p=0; p<n_qp; p++)
	  {
	    const Real
	      dx_dxi   = dxdxi_map(p),   dy_dxi   = dydxi_map(p),   dz_dxi   = dzdxi_map(p),
	      dx_deta  = dxdeta_map(p),  dy_deta  = dydeta_map(p),  dz_deta  = dzdeta_map(p),
	      dx_dzeta = dxdzeta_map(p), dy_dzeta = dydzeta_map(p), dz_dzeta = dzdzeta_map(p);
	    
	    // Symbolically, the matrix determinant is
	    //
	    //         | dx/dxi   dy/dxi   dz/dxi   |
	    // jac =   | dx/deta  dy/deta  dz/deta  |
	    //         | dx/dzeta dy/dzeta dz/dzeta |
	    // 
	    // jac = dx/dxi*(dy/deta*dz/dzeta - dz/deta*dy/dzeta) +
	    //       dy/dxi*(dz/deta*dx/dzeta - dx/deta*dz/dzeta) +
	    //       dz/dxi*(dx/deta*dy/dzeta - dy/deta*dx/dzeta)
	    
	    const Real jac = (dx_dxi*(dy_deta*dz_dzeta - dz_deta*dy_dzeta)  +
			      dy_dxi*(dz_deta*dx_dzeta - dx_deta*dz_dzeta)  +
			      dz_dxi*(dx_deta*dy_dzeta - dy_deta*dx_dzeta));
	    
	    if (jac <= 0.)
	      {
		std::cerr << "ERROR: negative Jacobian: "
			  << jac
			  << std::endl;
		error();
	      }

	    JxW[p] = jac*qw[p];
	    
	    // Compute the shape function derivatives wrt x,y at the
	    // quadrature points
	    const Real
	      inv_jac  = 1./jac;	    
	    
	    dxidx_map[p]   = (dy_deta*dz_dzeta - dz_deta*dy_dzeta)*inv_jac;
	    dxidy_map[p]   = (dz_deta*dx_dzeta - dx_deta*dz_dzeta)*inv_jac;
	    dxidz_map[p]   = (dx_deta*dy_dzeta - dy_deta*dx_dzeta)*inv_jac;
	    
	    detadx_map[p]  = (dz_dxi*dy_dzeta  - dy_dxi*dz_dzeta )*inv_jac;
	    detady_map[p]  = (dx_dxi*dz_dzeta  - dz_dxi*dx_dzeta )*inv_jac;
	    detadz_map[p]  = (dy_dxi*dx_dzeta  - dx_dxi*dy_dzeta )*inv_jac;
	    
	    dzetadx_map[p] = (dy_dxi*dz_deta   - dz_dxi*dy_deta  )*inv_jac;
	    dzetady_map[p] = (dz_dxi*dx_deta   - dx_dxi*dz_deta  )*inv_jac;
	    dzetadz_map[p] = (dx_dxi*dy_deta   - dy_dxi*dx_deta  )*inv_jac;
	  }
	
	// done computing the map
	break;
      }



    default:
      error();
    }
  
  // Stop logging the map computation.
  STOP_LOG("compute_map()", "FE");  
}




template <unsigned int Dim, FEFamily T>
Point FE<Dim,T>::map (const Elem* elem,
		      const Point& reference_point)
{
  assert (elem != NULL);
    
  Point p;

  const ElemType type     = elem->type();
  const Order order       = elem->default_order();
  const unsigned int n_sf = FE<Dim,LAGRANGE>::n_shape_functions(type, order);

  // Lagrange basis functions are used for mapping
  for (unsigned int i=0; i<n_sf; i++)
    p.add_scaled (elem->point(i),
		  FE<Dim,LAGRANGE>::shape(type,
					  order,
					  i,
					  reference_point)
		  );

  return p;
}



template <unsigned int Dim, FEFamily T>
Point FE<Dim,T>::map_xi (const Elem* elem,
			 const Point& reference_point)
{
  assert (elem != NULL);
    
  Point p;

  const ElemType type     = elem->type();
  const Order order       = elem->default_order();
  const unsigned int n_sf = FE<Dim,LAGRANGE>::n_shape_functions(type, order);

  // Lagrange basis functions are used for mapping    
  for (unsigned int i=0; i<n_sf; i++)
    p.add_scaled (elem->point(i),
		  FE<Dim,LAGRANGE>::shape_deriv(type,
						order,
						i,
						0,
						reference_point)
		  );
    
  return p;
}



template <unsigned int Dim, FEFamily T>
Point FE<Dim,T>::map_eta (const Elem* elem,
			  const Point& reference_point)
{
  assert (elem != NULL);
    
  Point p;

  const ElemType type     = elem->type();
  const Order order       = elem->default_order();
  const unsigned int n_sf = FE<Dim,LAGRANGE>::n_shape_functions(type, order);
  
  // Lagrange basis functions are used for mapping
  for (unsigned int i=0; i<n_sf; i++)
    p.add_scaled (elem->point(i),
		  FE<Dim,LAGRANGE>::shape_deriv(type,
						order,
						i,
						1,
						reference_point)
		  );
    
  return p;
}



template <unsigned int Dim, FEFamily T>
Point FE<Dim,T>::map_zeta (const Elem* elem,
			   const Point& reference_point)
{
  assert (elem != NULL);
    
  Point p;

  const ElemType type     = elem->type();
  const Order order       = elem->default_order();
  const unsigned int n_sf = FE<Dim,LAGRANGE>::n_shape_functions(type, order);

  // Lagrange basis functions are used for mapping
  for (unsigned int i=0; i<n_sf; i++)
    p.add_scaled (elem->point(i),
		  FE<Dim,LAGRANGE>::shape_deriv(type,
						order,
						i,
						2,
						reference_point)
		  );
    
  return p;
}



template <unsigned int Dim, FEFamily T>
Point FE<Dim,T>::inverse_map (const Elem* elem,
			      const Point& physical_point,
			      const Real tolerance,
			      const bool secure)
{
  assert (elem != NULL);
  assert (tolerance >= 0.);

  
  // Start logging the map inversion.  
  START_LOG("inverse_map()", "FE");
  
  // How much did the point on the reference
  // element change by in this Newton step?
  Real error = 0.;
  
  //  The point on the reference element.  This is
  //  the "initial guess" for Newton's method.  The
  //  centroid seems like a good idea, but computing
  //  it is a little more intensive than, say taking
  //  the zero point.  
  // 
  //  Convergence should be insensitive of this choice
  //  for "good" elements.
  Point p; // the zero point.  No computation required

  //  The number of iterations in the map inversion process.
  unsigned int cnt = 0;




  //  Newton iteration loop.
  do
    {
      //  Where our current iterate \p p maps to.
      const Point physical_guess = FE<Dim,T>::map (elem, p);

      //  How far our current iterate is from the actual point.
      const Point delta = physical_point - physical_guess;

      //  Increment in current iterate \p p, will be computed.
      Point dp;


      //  The form of the map and how we invert it depends
      //  on the dimension that we are in.
      switch (Dim)
	{
	  
	  // ------------------------------------------------------------------
	  //  1D map inversion
	  // 
	  //  Here we find the point on a 1D reference element that maps to
	  //  the point \p physical_point in the domain.  This is a bit tricky
	  //  since we do not want to assume that the point \p physical_point
	  //  is also in a 1D domain.  In particular, this method might get
	  //  called on the edge of a 3D element, in which case
	  //  \p physical_point actually lives in 3D.
	case 1:
	  {
	    const Point dxi = FE<Dim,T>::map_xi (elem, p);
	    
	    //  Newton's method in this case looks like
	    // 
	    //  {X} - {X_n} = [J]*dp
	    // 
	    //  Where {X}, {X_n} are 3x1 vectors, [J] is a 3x1 matrix
	    //  d(x,y,z)/dxi, and we seek dp, a scalar.  Since the above
	    //  system is either overdetermined or rank-deficient, we will
	    //  solve the normal equations for this system
	    // 
	    //  [J]^T ({X} - {X_n}) = [J]^T [J] {dp}
	    // 
	    //  which involves the trivial inversion of the scalar
	    //  G = [J]^T [J]
	    const Real G = dxi*dxi;
	    
	    if (secure)
	      assert (G > 1.e-20);
	    
	    const Real Ginv = 1./G;
	    
	    const Real  dxidelta = dxi*delta;
	    
	    dp(0) = Ginv*dxidelta;

	    break;
	  }



	  // ------------------------------------------------------------------
	  //  2D map inversion
	  // 
	  //  Here we find the point on a 2D reference element that maps to
	  //  the point \p physical_point in the domain.  This is a bit tricky
	  //  since we do not want to assume that the point \p physical_point
	  //  is also in a 2D domain.  In particular, this method might get
	  //  called on the face of a 3D element, in which case
	  //  \p physical_point actually lives in 3D.
	case 2:
	  {
	    const Point dxi  = FE<Dim,T>::map_xi  (elem, p);
	    const Point deta = FE<Dim,T>::map_eta (elem, p);
	    
	    //  Newton's method in this case looks like
	    // 
	    //  {X} - {X_n} = [J]*{dp}
	    // 
	    //  Where {X}, {X_n} are 3x1 vectors, [J] is a 3x2 matrix
	    //  d(x,y,z)/d(xi,eta), and we seek {dp}, a 2x1 vector.  Since
	    //  the above system is either overdermined or rank-deficient,
	    //  we will solve the normal equations for this system
	    // 
	    //  [J]^T ({X} - {X_n}) = [J]^T [J] {dp}
	    // 
	    //  which involves the inversion of the 2x2 matrix
	    //  [G] = [J]^T [J]
	    const Real
	      G11 = dxi*dxi,  G12 = dxi*deta,
	      G21 = dxi*deta, G22 = deta*deta;
	    
	    
	    const Real det = (G11*G22 - G12*G21);
	    
	    if (secure)
	      assert (det > 1.e-20);

	    const Real inv_det = 1./det;
	    
	    const Real
	      Ginv11 =  G22*inv_det,
	      Ginv12 = -G12*inv_det,
	      
	      Ginv21 = -G21*inv_det,
	      Ginv22 =  G11*inv_det;
	    
	    
	    const Real  dxidelta  = dxi*delta;
	    const Real  detadelta = deta*delta;
	    
	    dp(0) = (Ginv11*dxidelta + Ginv12*detadelta);
	    dp(1) = (Ginv21*dxidelta + Ginv22*detadelta);

	    break;
	  }
	  

	  
	  // ------------------------------------------------------------------
	  //  3D map inversion
	  // 
	  //  Here we find the point in a 3D reference element that maps to
	  //  the point \p physical_point in a 3D domain. Nothing special
	  //  has to happen here, since (unless the map is singular because
	  //  you have a BAD element) the map will be invertable and we can
	  //  apply Newton's method directly.
	case 3:
	  {
       	    const Point dxi   = FE<Dim,T>::map_xi   (elem, p);
	    const Point deta  = FE<Dim,T>::map_eta  (elem, p);
	    const Point dzeta = FE<Dim,T>::map_zeta (elem, p);
	    
	    //  Newton's method in this case looks like
	    // 
	    //  {X} = {X_n} + [J]*{dp}
	    // 
	    //  Where {X}, {X_n} are 3x1 vectors, [J] is a 3x3 matrix
	    //  d(x,y,z)/d(xi,eta,zeta), and we seek {dp}, a 3x1 vector.
	    //  Since the above system is nonsingular for invertable maps
	    //  we will solve 
	    // 
	    //  {dp} = [J]^-1 ({X} - {X_n})
	    // 
	    //  which involves the inversion of the 3x3 matrix [J]
	    const Real
	      J11 = dxi(0), J12 = deta(0), J13 = dzeta(0),
	      J21 = dxi(1), J22 = deta(1), J23 = dzeta(1),
	      J31 = dxi(2), J32 = deta(2), J33 = dzeta(2);
	    
	    const Real det = (J11*(J22*J33 - J23*J32) +
			      J12*(J23*J31 - J21*J33) +
			      J13*(J21*J32 - J22*J31));
	    
	    if (secure)
	      assert (std::abs(det) > 1.e-20);

	    const Real inv_det = 1./det;
	    
	    const Real
	      Jinv11 =  (J22*J33 - J23*J32)*inv_det,
	      Jinv12 = -(J12*J33 - J13*J32)*inv_det,
	      Jinv13 =  (J12*J23 - J13*J22)*inv_det,
	      
	      Jinv21 = -(J21*J33 - J23*J31)*inv_det,
	      Jinv22 =  (J11*J33 - J13*J31)*inv_det,
	      Jinv23 = -(J11*J23 - J13*J21)*inv_det,
	      
	      Jinv31 =  (J21*J32 - J22*J31)*inv_det,
	      Jinv32 = -(J11*J32 - J12*J31)*inv_det,
	      Jinv33 =  (J11*J22 - J12*J21)*inv_det;
	    
	    dp(0) = (Jinv11*delta(0) +
		     Jinv12*delta(1) +
		     Jinv13*delta(2));
	    
	    dp(1) = (Jinv21*delta(0) +
		     Jinv22*delta(1) +
		     Jinv23*delta(2));
	    
	    dp(2) = (Jinv31*delta(0) +
		     Jinv32*delta(1) +
		     Jinv33*delta(2));

	    break;
	  }


	  //  Some other dimension?
	default:
	  error();
	} // end switch(Dim), dp now computed



      //  ||P_n+1 - P_n||
      error = dp.size();

      //  P_n+1 = P_n + dp
      p.add (dp);
      
      //  Increment the iteration count.
      cnt++;
      
      //  Watch for divergence of Newton's
      //  method.  Here's how it goes:
      //  (1) For good elements, we expect convergence in 10 iterations.
      //      - If called with (secure == true) and we have not yet converged
      //        print out a warning message.
      //      - If called with (secure == true) and we have not converged in
      //        20 iterations abort
      //  (2) This method may be called in cases when the target point is not
      //      inside the element and we have no business expecting convergence.
      //      For these cases if we have not converged in 10 iterations forget
      //      about it.
      if (cnt > 10)
	{
	  //  Do not bother about divergence when secure is false.
	  if (secure)
	    {
	      here();
	      std::cerr << "WARNING: Newton scheme has not converged in "
			<< cnt << " iterations:" << std::endl
			<< "   physical_point="
			<< physical_point
			<< "   physical_guess="
			<< physical_guess
			<< "   dp="
			<< dp
			<< "   p="
			<< p
			<< "   error=" << error
			<< std::endl;

	      if (cnt > 20)
		{
		  std::cerr << "ERROR: Newton scheme FAILED to converge in "
			    << cnt
			    << " iterations!"
			    << std::endl;

		  error();
		}
	    }
	  else
	    {
	      break;
	    }
	}
    }
  while (error > tolerance);



  //  If we are in debug mode do a sanity check.  Make sure
  //  the point \p p on the reference element actually does
  //  map to the point \p physical_point within a tolerance.
#ifdef DEBUG
  
  if (secure)
    {  
      const Point check = FE<Dim,T>::map (elem, p);
      const Point diff  = physical_point - check;
      
      if (diff.size() > tolerance)
	{
	  here();
	  std::cerr << "WARNING:  diff is "
		    << diff.size()
		    << std::endl
		    << " point="
		    << physical_point;
	  std::cerr << " local=" << check;
	  std::cerr << " lref= " << p;
	}
    }
  
#endif


  
  //  Stop logging the map inversion.
  STOP_LOG("inverse_map()", "FE");
  
  return p;
}



template <unsigned int Dim, FEFamily T>
void FE<Dim,T>::inverse_map (const Elem* elem,
			     const std::vector<Point>& physical_points,
			     std::vector<Point>&       reference_points)
{
  // The number of points to find the
  // inverse map of
  const unsigned int n_points = physical_points.size();

  // Resize the vector to hold the points
  // on the reference element
  reference_points.resize(n_points);

  // Find the coordinates on the reference
  // element of each point in physical space
  for (unsigned int p=0; p<n_points; p++)
    reference_points[p] =
      FE<Dim,T>::inverse_map (elem, physical_points[p]);
}



//--------------------------------------------------------------
// Explicit instantiations using the macro from fe_macro.h
INSTANTIATE_IMAP(1);
INSTANTIATE_IMAP(2);
INSTANTIATE_IMAP(3);

#if defined(__INTEL_COMPILER)
INSTANTIATE_MAP(1);
INSTANTIATE_MAP(2);
INSTANTIATE_MAP(3);
#endif
