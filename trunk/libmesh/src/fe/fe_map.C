// $Id: fe_map.C,v 1.10 2003-02-07 04:00:42 jwpeterson Exp $

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
  assert (elem  != NULL);
  
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
	// the Jacobian at the quadrature points
	
	// Resize the vectors to hold data at the quadrature points
	{  
	  xyz.resize(n_qp);
	  dxyzdxi_map.resize(n_qp);
	  dxidx_map.resize(n_qp);
	  
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
	  {
	    // Reference to the point, helps eliminate
	    // exessive temporaries in the inner loop
	    const Point& elem_point = elem->point(i);
	    
	    for (unsigned int p=0; p<n_qp; p++) // for each quadrature point
	      {	  
		xyz[p].add_scaled        (elem_point, phi_map[i][p]    );
		dxyzdxi_map[p].add_scaled(elem_point, dphidxi_map[i][p]);
	      };
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
	      };
	    
	    assert (dxdxi_map(p) != 0.);
	    
	    dxidx_map[p] = 1./dxdxi_map(p);
	    
	    JxW[p] = jac*qw[p];
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
	// the Jacobian at the quadrature points

	// Resize the vectors to hold data at the quadrature points
	{  
	  xyz.resize(n_qp);
	  dxyzdxi_map.resize(n_qp);
	  dxyzdeta_map.resize(n_qp);
	  dxidx_map.resize(n_qp);
	  dxidy_map.resize(n_qp);
	  detadx_map.resize(n_qp);
	  detady_map.resize(n_qp);
	  
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
	  {
	    // Reference to the point, helps eliminate
	    // exessive temporaries in the inner loop
	    const Point& elem_point = elem->point(i);
	    
	    for (unsigned int p=0; p<n_qp; p++) // for each quadrature point
	      {	  
		xyz[p].add_scaled          (elem_point, phi_map[i][p]     );
		dxyzdxi_map[p].add_scaled  (elem_point, dphidxi_map[i][p] );
		dxyzdeta_map[p].add_scaled (elem_point, dphideta_map[i][p]);
	      };
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
	    const Real
	      dx_dxi  = dxdxi_map(p),  dy_dxi  = dydxi_map(p),
	      dx_deta = dxdeta_map(p), dy_deta = dydeta_map(p);
	    
	    // Symbolically, the matrix determinant is
	    //
	    //         | dx/dxi   dy/dxi  |
	    // jac =   | dx/deta  dy/deta |
	    //         
	    // jac = dx/dxi*dy/deta - dx/deta*dy/dxi 
	    
	    // Compute the Jacobian.  This assumes the 2D face
	    // lives in 2D space
	    const Real jac = (dx_dxi*dy_deta - dx_deta*dy_dxi);	    
	    
	    if (jac <= 0.)
	      {
		std::cerr << "ERROR: negative Jacobian: "
			  << jac
			  << std::endl;
		error();
	      };
	    
	    JxW[p] = jac*qw[p];
	    
	    // Compute the shape function derivatives wrt x,y at the
	    // quadrature points
	    const Real
	      inv_jac = 1./jac;
	    
	    dxidx_map[p]  =  dy_deta*inv_jac; //dxi/dx  =  (1/J)*dy/deta
	    dxidy_map[p]  = -dx_deta*inv_jac; //dxi/dy  = -(1/J)*dx/deta
	    detadx_map[p] = -dy_dxi* inv_jac; //deta/dx = -(1/J)*dy/dxi
	    detady_map[p] =  dx_dxi* inv_jac; //deta/dy =  (1/J)*dx/dxi
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
	      };
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
	      };

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
	  };
	
	// done computing the map
	return;
      };



    default:
      error();
    };
};




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
};



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
};



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
};



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
};



template <unsigned int Dim, FEFamily T>
Point FE<Dim,T>::inverse_map (const Elem* elem,
			      const Point& physical_point,
			      const Real tolerance)
{
  assert (elem != NULL);
  assert (tolerance >= 0.);

  /**
   * How much did the point on the reference
   * element change by in this Newton step?
   */
  Real error = 0.;

  /**
   * The point on the reference element.  This is
   * the "initial guess" for Newton's method.  The
   * centroid seems like a good idea, but computing
   * it is a little more intensive than, say taking
   * the zero point.  
   *
   * Convergence should be insensitive of this choice
   * for "good" elements.
   */
  //Point p = elem->centroid(); // A reasonable guess.  Requires computation
  Point p; // the zero point.  No computation required

  /**
   * The number of iterations in the map inversion process.
   */
  unsigned int cnt = 0;




  /**
   * Newton iteration loop.
   */
  do
    {
      /**
       * Where our current iterate \p p maps to.
       */
      const Point physical_guess = FE<Dim,T>::map (elem, p);

      /**
       * How far our current iterate is from the actual point.
       */
      const Point delta = physical_point - physical_guess;

      /**
       * Increment in current iterate \p p, will be computed.
       */
      Point dp;


      /**
       * The form of the map and how we invert it depends
       * on the dimension that we are in.
       */      
      switch (Dim)
	{
	  
	  /**	 
	   *------------------------------------------------------------------
	   * 1D map inversion
	   *
	   * Here we find the point on a 1D reference element that maps to
	   * the point \p physical_point in the domain.  This is a bit tricky
	   * since we do not want to assume that the point \p physical_point
	   * is also in a 1D domain.  In particular, this method might get
	   * called on the edge of a 3D element, in which case \p physical_point
	   * actually lives in 3D.
	   */
	case 1:
	  {
	    const Point dxi            = FE<Dim,T>::map_xi (elem, p);
	    
	    /**
	     * Newton's method in this case looks like
	     *
	     * {X} - {X_n} = [J]*dp
	     *
	     * Where {X}, {X_n} are 3x1 vectors, [J] is a 3x1 matrix
	     * d(x,y,z)/dxi, and we seek dp, a scalar.  Since the above
	     * system is either overdermined or rank-deficient, we will
	     * solve the normal equations for this system
	     *
	     * [J]^T ({X} - {X_n}) = [J]^T [J] {dp}
	     *
	     * which involves the trivial inversion of the scalar
	     * G = [J]^T [J]
	     */	    
	    const Real G = dxi*dxi;
	    
	    assert (G > 0.);
	    
	    const Real Ginv = 1./G;
	    
	    const Real  dxidelta = dxi*delta;
	    
	    dp(0) = Ginv*dxidelta;

	    break;
	  };



	  /**	 
	   *------------------------------------------------------------------
	   * 2D map inversion
	   *
	   * Here we find the point on a 2D reference element that maps to
	   * the point \p physical_point in the domain.  This is a bit tricky
	   * since we do not want to assume that the point \p physical_point
	   * is also in a 2D domain.  In particular, this method might get
	   * called on the face of a 3D element, in which case \p physical_point
	   * actually lives in 3D.
	   */
	case 2:
	  {
	    const Point dxi            = FE<Dim,T>::map_xi  (elem, p);
	    const Point deta           = FE<Dim,T>::map_eta (elem, p);
	    
	    /**
	     * Newton's method in this case looks like
	     *
	     * {X} - {X_n} = [J]*{dp}
	     *
	     * Where {X}, {X_n} are 3x1 vectors, [J] is a 3x2 matrix
	     * d(x,y,z)/d(xi,eta), and we seek {dp}, a 2x1 vector.  Since
	     * the above system is either overdermined or rank-deficient,
	     * we will solve the normal equations for this system
	     *
	     * [J]^T ({X} - {X_n}) = [J]^T [J] {dp}
	     *
	     * which involves the inversion of the 2x2 matrix
	     * [G] = [J]^T [J]
	     */
	    const Real
	      G11 = dxi*dxi,  G12 = dxi*deta,
	      G21 = dxi*deta, G22 = deta*deta;
	    
	    
	    const Real det = (G11*G22 - G12*G21);
	    
	    assert (det > 0.);
	    assert (fabs(det) > 1.e-10);
	    
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
	  };
	  

	  
	  /**	 
	   *------------------------------------------------------------------
	   * 3D map inversion
	   *
	   * Here we find the point in a 3D reference element that maps to
	   * the point \p physical_point in a 3D domain. Nothing special
	   * has to happen here, since (unless the map is singular because
	   * you have a BAD element) the map will be invertable and we can
	   * apply Newton's method directly.
	   */
	case 3:
	  {
       	    const Point dxi   = FE<Dim,T>::map_xi   (elem, p);
	    const Point deta  = FE<Dim,T>::map_eta  (elem, p);
	    const Point dzeta = FE<Dim,T>::map_zeta (elem, p);
	    
	    /**
	     * Newton's method in this case looks like
	     *
	     * {X} - {X_n} = [J]*{dp}
	     *
	     * Where {X}, {X_n} are 3x1 vectors, [J] is a 3x3 matrix
	     * d(x,y,z)/dxi, and we seek {dp}, a 3x1 vector. Since the above
	     * system is nonsingular for invertable maps we will solve 
	     *
	     * ({X} - {X_n}) = [J] {dp}
	     *
	     * which involves the inversion of the 3x3 matrix [J]
	     */	    
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
	  };


	  /**
	   * Some other dimension?
	   */
	default:
	  error();
	}; // end switch(Dim), dp now computed



      /**
       * ||P_n+1 - P_n||
       */
      error = dp.size();

      /**
       * P_n+1 = P_n + dp
       */
      p.add (dp);

      /**
       * Increment the iteration count.
       */
      cnt++;

      /**
       * Watch for divergence of Newton's
       * method.
       */
      if (cnt > 10)
	{
	  here();
	  std::cerr << "WARNING: Newton scheme has not converged in "
		    << cnt << " iterations!"
		    << std::endl;
	  
	  if (cnt > 20)
	    {
	      std::cerr << "ERROR: Newton scheme FAILED to converge in "
			<< cnt << " iterations!" << std::endl
			<< "p="; 
	      error();
	    };
	};
    }
  while (error > tolerance);



  /**
   * If we are in debug mode do a sanity check.  Make sure
   * the point \p p on the reference element actually does
   * map to the point \p physical_point within a tolerance.
   */ 
#ifdef DEBUG
	
  const Point check = FE<Dim,T>::map (elem, p);
  const Point diff  = physical_point - check;
  
  if (diff.size() > tolerance)
    {
      here();
      std::cerr << "WARNING:  diff is "
		<< diff.size()
		<< std::endl;
    };
  
#endif


  
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
