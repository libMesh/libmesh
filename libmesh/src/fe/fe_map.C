// $Id: fe_map.C,v 1.6 2003-02-03 03:51:49 ddreyer Exp $

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
#include "fe.h"
#include "quadrature.h"
#include "elem.h"




void FEBase::compute_map(const QBase* qrule,
			 const Elem* elem)
{
  assert (qrule != NULL);
  assert (elem != NULL);
  
  const unsigned int        n_qp = qrule->n_points();
  const std::vector<Real> &   qw = qrule->get_weights();


  switch (dim)
    {


      //--------------------------------------------------------------------
      // 1D
    case 1:
      {
	//------------------------------------------------------------------
	// Compute the values at the quadrature points,
	// the Jacobian at the quadrature points, etc...
	
	// Resize the vectors to hold data at the quadrature points
	{  
	  xyz.resize(n_qp);
	  dxyzdxi_map.resize(n_qp);
	  dxidx_map.resize(n_qp);
	  
	  jac.resize(n_qp);
	  JxW.resize(n_qp);
	};
	
	// Clear the entities that will be summed
	for (unsigned int p=0; p<n_qp; p++)
	  {
	    xyz[p].clear();
	    dxyzdxi_map[p].clear();
	  };
	
	
	// compute x, dxdxi at the quadrature points    
	for (unsigned int i=0; i<phi_map.size(); i++) // sum over the nodes
	  for (unsigned int p=0; p<n_qp; p++) // for each quadrature point...
	    {	  
	      xyz[p]         += elem->point(i)*phi_map[i][p];
	      
	      dxyzdxi_map[p] += elem->point(i)*dphidxi_map[i][p];
	    };

	/*
        // Test the inverse map
	for (unsigned int p=0; p<n_qp; p++)
	{
	const Point p_inv = inverse_map (elem, xyz[p]);
	    
	std::cout << "qp[p]   = ";
	qrule->qp(p).print();
	std::cout << "inv_map = ";
	p_inv.print();
	};
	*/

	// compute the jacobian at the quadrature points
	for (unsigned int p=0; p<n_qp; p++)
	  {
	    jac[p] = dxdxi_map(p);
	    
	    if (jac[p] <= 0.)
	      {
		std::cerr << "ERROR: negative Jacobian: "
			  << jac[p]
			  << std::endl;
		error();
	      };
	    
	    assert (dxdxi_map(p) != 0.);
	    
	    dxidx_map[p] = 1./dxdxi_map(p);
	    
	    JxW[p] = jac[p]*qw[p];
	  };

	// done computing the map
	return;
      };


      
      //--------------------------------------------------------------------
      // 2D
    case 2:
      {
	//------------------------------------------------------------------
	// Compute the (x,y) values at the quadrature points,
	// the Jacobian at the quadrature points, etc..

	// Resize the vectors to hold data at the quadrature points
	{  
	  xyz.resize(n_qp);
	  dxyzdxi_map.resize(n_qp);
	  dxyzdeta_map.resize(n_qp);
	  dxidx_map.resize(n_qp);
	  dxidy_map.resize(n_qp);
	  detadx_map.resize(n_qp);
	  detady_map.resize(n_qp);
	  
	  jac.resize(n_qp);
	  JxW.resize(n_qp);
	};
	
	// Clear the entities that will be summed
	for (unsigned int p=0; p<n_qp; p++)
	  {
	    xyz[p].clear();
	    dxyzdxi_map[p].clear();
	    dxyzdeta_map[p].clear();
	  };
	
	
	// compute (x,y), dxdxi, dydxi, dxdeta, dydeta at the quadrature points    
	for (unsigned int i=0; i<phi_map.size(); i++) // sum over the nodes
	  for (unsigned int p=0; p<n_qp; p++) // for each quadrature point...
	    {	  
	      xyz[p]          += elem->point(i)*phi_map[i][p];
	      
	      dxyzdxi_map[p]  += elem->point(i)*dphidxi_map[i][p];
	      
	      dxyzdeta_map[p] += elem->point(i)*dphideta_map[i][p];
	    };
	
	/*
        // Test the inverse map
	for (unsigned int p=0; p<n_qp; p++)
	{
	const Point p_inv = inverse_map (elem, xyz[p]);
	    
	std::cout << "qp[p]   = ";
	qrule->qp(p).print();
	std::cout << "inv_map = ";
	p_inv.print();
	};
	*/
	
	// compute the jacobian at the quadrature points
	for (unsigned int p=0; p<n_qp; p++)
	  {
	    jac[p] = (dxdxi_map(p)*dydeta_map(p) - dxdeta_map(p)*dydxi_map(p));
	    
	    if (jac[p] <= 0.)
	      {
		std::cerr << "ERROR: negative Jacobian: "
			  << jac[p]
			  << std::endl;
		error();
	      };
	    
	    JxW[p] = jac[p]*qw[p];
	  };
	
	
	// Compute the shape function derivatives wrt x,y at the
	// quadrature points
	for (unsigned int p=0; p<n_qp; p++)
	  {
	    dxidx_map[p]  =  (dydeta_map(p))/jac[p];  // dxi/dx  =  (1/J)*dy/deta
	    dxidy_map[p]  = -(dxdeta_map(p))/jac[p];  // dxi/dy  = -(1/J)*dx/deta
	    detadx_map[p] = -(dydxi_map(p))/jac[p];   // deta/dx = -(1/J)*dy/dxi;
	    detady_map[p] =  (dxdxi_map(p))/jac[p];   // deta/dy =  (1/J)*dx/dxi;
	  };        

	// done computing the map
	return;
      };


      
      //--------------------------------------------------------------------
      // 3D
    case 3:
      {
	//------------------------------------------------------------------
	// Compute the (x,y,z) values at the quadrature points,
	// the Jacobian at the quadrature points, etc..

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
	  
	  jac.resize (n_qp);
	  JxW.resize (n_qp);
	};
    
	// Clear the entities that will be summed
	for (unsigned int p=0; p<n_qp; p++)
	  {
	    xyz[p].clear           ();
	    dxyzdxi_map[p].clear   ();
	    dxyzdeta_map[p].clear  ();
	    dxyzdzeta_map[p].clear ();
	  };
	
	
	// compute (x,y,z), dxdxi, dydxi, dxdeta, dydeta at the quadrature points    
	for (unsigned int i=0; i<phi_map.size(); i++) // sum over the nodes
	  for (unsigned int p=0; p<n_qp; p++) // for each quadrature point...
	    {	  
	      xyz[p]           += elem->point(i)*phi_map[i][p];
	      
	      dxyzdxi_map[p]   += elem->point(i)*dphidxi_map[i][p];
	      
	      dxyzdeta_map[p]  += elem->point(i)*dphideta_map[i][p];
	      
	      dxyzdzeta_map[p] += elem->point(i)*dphidzeta_map[i][p];
	    };
	
	/*
        // Test the inverse map
	for (unsigned int p=0; p<n_qp; p++)
	{
	const Point p_inv = inverse_map (elem, xyz[p]);
	    
	std::cout << "qp[p]   = ";
	qrule->qp(p).print();
	std::cout << "inv_map = ";
	p_inv.print();
	};
	*/
	  
	// compute the jacobian at the quadrature points
	for (unsigned int p=0; p<n_qp; p++)
	  {
	    // jac = dx/dxi*(dy/deta*dz/dzeta - dz/deta*dy/dzeta) +
	    //       dy/dxi*(dz/deta*dx/dzeta - dx/deta*dz/dzeta) +
	    //       dz/dxi*(dx/deta*dy/dzeta - dy/deta*dx/dzeta)
	    
	    jac[p] = (dxdxi_map(p)*(dydeta_map(p)*dzdzeta_map(p) - dzdeta_map(p)*dydzeta_map(p))  +
		      dydxi_map(p)*(dzdeta_map(p)*dxdzeta_map(p) - dxdeta_map(p)*dzdzeta_map(p))  +
		      dzdxi_map(p)*(dxdeta_map(p)*dydzeta_map(p) - dydeta_map(p)*dxdzeta_map(p)));
	    
	    if (jac[p] <= 0.)
	      {
		std::cerr << "ERROR: negative Jacobian: "
			  << jac[p]
			  << std::endl;
		error();
	      };

	    JxW[p] = jac[p]*qw[p];
	  };


	// Compute the shape function derivatives wrt x,y at the
	// quadrature points
	for (unsigned int p=0; p<n_qp; p++)
	  {
	    dxidx_map[p] = (dydeta_map(p)*dzdzeta_map(p) - dzdeta_map(p)*dydzeta_map(p))/jac[p];
	    dxidy_map[p] = (dzdeta_map(p)*dxdzeta_map(p) - dxdeta_map(p)*dzdzeta_map(p))/jac[p];
	    dxidz_map[p] = (dxdeta_map(p)*dydzeta_map(p) - dydeta_map(p)*dxdzeta_map(p))/jac[p];
	    
	    detadx_map[p] = (dzdxi_map(p)*dydzeta_map(p) - dydxi_map(p)*dzdzeta_map(p))/jac[p];
	    detady_map[p] = (dxdxi_map(p)*dzdzeta_map(p) - dzdxi_map(p)*dxdzeta_map(p))/jac[p];
	    detadz_map[p] = (dydxi_map(p)*dxdzeta_map(p) - dxdxi_map(p)*dydzeta_map(p))/jac[p];
	    
	    dzetadx_map[p] = (dydxi_map(p)*dzdeta_map(p) - dzdxi_map(p)*dydeta_map(p))/jac[p];
	    dzetady_map[p] = (dzdxi_map(p)*dxdeta_map(p) - dxdxi_map(p)*dzdeta_map(p))/jac[p];
	    dzetadz_map[p] = (dxdxi_map(p)*dydeta_map(p) - dydxi_map(p)*dxdeta_map(p))/jac[p];
	  };        

	// done computing the map
	return;
      };



    default:
      error();
    };

  error();

  return;      
};




template <unsigned int Dim, FEFamily T>
Point FE<Dim,T>::map (const Elem* elem,
		      const Point& reference_point)
{
  assert (elem != NULL);
    
  Point p;
    
  p.clear();

  const ElemType type     = elem->type();
  const Order order       = elem->default_order();
  const unsigned int n_sf = FE<Dim,LAGRANGE>::n_shape_functions(type, order);

  // Lagrange basis functions are used for mapping
  for (unsigned int i=0; i<n_sf; i++)
    p += elem->point(i)*FE<Dim,LAGRANGE>::shape(type,
						order,
						i,
						reference_point);

  return p;
};



template <unsigned int Dim, FEFamily T>
Point FE<Dim,T>::map_xi (const Elem* elem,
			 const Point& reference_point)
{
  assert (elem != NULL);
    
  Point p;
    
  p.clear();

  const ElemType type     = elem->type();
  const Order order       = elem->default_order();
  const unsigned int n_sf = FE<Dim,LAGRANGE>::n_shape_functions(type, order);

  // Lagrange basis functions are used for mapping    
  for (unsigned int i=0; i<n_sf; i++)
    p += elem->point(i)*FE<Dim,LAGRANGE>::shape_deriv(type,
						      order,
						      i,
						      0,
						      reference_point);
    
  return p;
};



template <unsigned int Dim, FEFamily T>
Point FE<Dim,T>::map_eta (const Elem* elem,
			  const Point& reference_point)
{
  assert (elem != NULL);
    
  Point p;
    
  p.clear();

  const ElemType type     = elem->type();
  const Order order       = elem->default_order();
  const unsigned int n_sf = FE<Dim,LAGRANGE>::n_shape_functions(type, order);
  
  // Lagrange basis functions are used for mapping
  for (unsigned int i=0; i<n_sf; i++)
    p += elem->point(i)*FE<Dim,LAGRANGE>::shape_deriv(type,
						  order,
						  i,
						  1,
						  reference_point);
    
  return p;
};



template <unsigned int Dim, FEFamily T>
Point FE<Dim,T>::map_zeta (const Elem* elem,
			   const Point& reference_point)
{
  assert (elem != NULL);
    
  Point p;
    
  p.clear();

  const ElemType type     = elem->type();
  const Order order       = elem->default_order();
  const unsigned int n_sf = FE<Dim,LAGRANGE>::n_shape_functions(type, order);

  // Lagrange basis functions are used for mapping
  for (unsigned int i=0; i<n_sf; i++)
    p += elem->point(i)*FE<Dim,LAGRANGE>::shape_deriv(type,
						      order,
						      i,
						      2,
						      reference_point);
    
  return p;
};



template <unsigned int Dim, FEFamily T>
Point FE<Dim,T>::inverse_map (const Elem* elem,
			      const Point& physical_point)
{
  assert (elem != NULL);

  switch (Dim)
    {
      
      //------------------------------------------------------------------
      // 1D map inversion
    case 1:
      {
	Real error = 0.;

	Point p;
	
	const Real
	  X = physical_point(0);

	unsigned int cnt = 0;
	
	do
	  {
	    // The actual update step
	    {
	      const Point physical_guess = FE<Dim,T>::map    (elem, p);
	      
	      const Point dxi            = FE<Dim,T>::map_xi (elem, p);
	      
	      const Real
		J = dxi(0);
	      
	      assert (J != 0.);
	      
	      const Real
		Jinv =  1./J;

	      
	      Point dp;
	      
	      dp(0) = Jinv*(-physical_guess(0) + X);
	      
	      error = dp.size();

	      p += dp;
	    };
	    
	    cnt++;
	    
	    if (cnt > 10)
	      {
		here();
		std::cerr << "WARNING: Newton scheme has not converged in "
			  << cnt << " iterations!"
			  << std::endl;
		
		if (cnt > 20)
		  {
		    std::cerr << "WARNING: Newton scheme FAILED to converge in "
			      << cnt << " iterations!"
			      << std::endl;
		    error();
		  };
	      };
	  }
	while (error > 1.e-6);

#ifdef DEBUG
	
	const Point check = FE<Dim,T>::map (elem, p);
	const Point diff  = physical_point - check;

	if (diff.size() > 1.e-6)
	  {
	    here();
	    std::cerr << "WARNING:  diff is "
		      << diff.size()
		      << std::endl;
	  };
       
#endif

	return p;
      };










      //------------------------------------------------------------------
      // 2D map inversion
    case 2:
      {
	Real error = 0.;

	Point p;
	
	const Real
	  X = physical_point(0), 
	  Y = physical_point(1);

	unsigned int cnt = 0;
	
	do
	  {
	    // The actual update step
	    {
	      const Point physical_guess = FE<Dim,T>::map (elem, p);
	      
	      const Point dxi   = FE<Dim,T>::map_xi   (elem, p);
	      const Point deta  = FE<Dim,T>::map_eta  (elem, p);
	      
	      const Real
		J11 = dxi(0), J12 = deta(0),
		J21 = dxi(1), J22 = deta(1);
	      
	      const Real det = (J11*J22 - J12*J21);
	      
	      assert (det > 0.);
	      assert (fabs(det) > 1.e-10);
	      
	      const Real inv_det = 1./det;
	      
	      const Real
		Jinv11 =  J22*inv_det,
		Jinv12 = -J12*inv_det,
		
		Jinv21 = -J21*inv_det,
		Jinv22 =  J11*inv_det;

	      
	      Point dp;
	      
	      dp(0) = (Jinv11*(-physical_guess(0) + X) +
		       Jinv12*(-physical_guess(1) + Y));
	      
	      dp(1) = (Jinv21*(-physical_guess(0) + X) +
		       Jinv22*(-physical_guess(1) + Y));
	      
	      error = dp.size();

	      p += dp;
	    };
	    
	    cnt++;

	    if (cnt > 10)
	      {
		here();
		std::cerr << "WARNING: Newton scheme has not converged in "
			  << cnt << " iterations!"
			  << std::endl;
		
		if (cnt > 20)
		  {
		    std::cerr << "WARNING: Newton scheme FAILED to converge in "
			      << cnt << " iterations!"
			      << std::endl;
		    error();
		  };
	      };
	  }
	while (error > 1.e-6);
	
#ifdef DEBUG
	
	const Point check = FE<Dim,T>::map (elem, p);
	const Point diff  = physical_point - check;

	if (diff.size() > 1.e-6)
	  {
	    here();
	    std::cerr << "WARNING:  diff is "
		      << diff.size()
		      << std::endl;
	  };
       
#endif

	return p;
      };










      //------------------------------------------------------------------
      // 3D map inversion
    case 3:
      {
	Real error = 0.;

	Point p;
	
	const Real
	  X = physical_point(0), 
	  Y = physical_point(1), 
	  Z = physical_point(2);

	unsigned int cnt = 0;
	
	do
	  {
	    // The actual update step
	    {
	      const Point physical_guess = FE<Dim,T>::map (elem, p);
	      
	      const Point dxi   = FE<Dim,T>::map_xi   (elem, p);
	      const Point deta  = FE<Dim,T>::map_eta  (elem, p);
	      const Point dzeta = FE<Dim,T>::map_zeta (elem, p);
	      
	      const Real
		J11 = dxi(0), J12 = deta(0), J13 = dzeta(0),
		J21 = dxi(1), J22 = deta(1), J23 = dzeta(1),
		J31 = dxi(2), J32 = deta(2), J33 = dzeta(2);
	      
	      const Real det = (J11*(J22*J33 - J23*J32) +
				J12*(J23*J31 - J21*J33) +
				J13*(J21*J32 - J22*J31));
	      
	      assert (det > 0.);
	      assert (fabs(det) > 1.e-10);
	      
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
	      
	      
	      Point dp;
	      
	      dp(0) = (Jinv11*(-physical_guess(0) + X) +
		       Jinv12*(-physical_guess(1) + Y) +
		       Jinv13*(-physical_guess(2) + Z));
	      
	      dp(1) = (Jinv21*(-physical_guess(0) + X) +
		       Jinv22*(-physical_guess(1) + Y) +
		       Jinv23*(-physical_guess(2) + Z));
	      
	      dp(2) = (Jinv31*(-physical_guess(0) + X) +
		       Jinv32*(-physical_guess(1) + Y) +
		       Jinv33*(-physical_guess(2) + Z));
	      
	      error = dp.size();

	      p += dp;
	    };
	    
	    cnt++;
	    
	    if (cnt > 10)
	      {
		here();
		std::cerr << "WARNING: Newton scheme has not converged in "
			  << cnt << " iterations!"
			  << std::endl;
		
		if (cnt > 20)
		  {
		    std::cerr << "ERROR: Newton scheme FAILED to converge in "
			      << cnt << " iterations!"
			      << std::endl;
		    return p;
		  };
	      };
	  }
	while (error > 1.e-6);

#ifdef DEBUG
	
	const Point check = FE<Dim,T>::map (elem, p);
	const Point diff  = physical_point - check;

	if (diff.size() > 1.e-6)
	  {
	    here();
	    std::cerr << "WARNING:  diff is "
		      << diff.size()
		      << std::endl;
	  };
       
#endif

	return p;
      };

      
    default:
      error();
      
    };

  error();

  Point p;
  
  return p;
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
