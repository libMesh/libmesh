// $Id: inf_fe.C,v 1.12 2003-02-20 12:54:57 spetersen Exp $

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
#include "mesh_config.h"
#ifdef ENABLE_INFINITE_ELEMENTS
#include "inf_fe.h"
#include "quadrature_gauss.h" /* this also includes "quadrature.h" */
#include "elem.h"
#include "fe.h"



// ------------------------------------------------------------
// InfFE class members



// Constructor
template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
InfFE<Dim,T_radial,T_map>::InfFE (const FEType& fet) :
  FEBase       (Dim, fet),

  _n_total_approx_sf (0),
  _n_total_qp        (0),

  base_qrule   (NULL),
  radial_qrule (NULL),
  base_elem    (NULL),
  base_fe      (NULL),
  /* 
   * initialize the current_fe_type to all the same
   * values as \p fet (since the FE families and coordinate
   * map type should @e not change), but use an invalid order
   * for the radial part (since this is the only order
   * that may change!).
   * the data structures like \p phi etc are not initialized 
   * through the constructor, but throught reinit()
   */
  current_fe_type ( FEType(fet.order, 
			   fet.family, 
			   INVALID_ORDER, 
			   fet.radial_family,      
			   fet.inf_map) )

{
  // Sanity check.  Make sure the family and
  // map specified in the template instantiation
  // matches the one in the FEType object
  assert (T_radial == fe_type.radial_family);
  assert (T_map    == fe_type.inf_map);

  // build the base_fe object, handle the AutoPtr
  AutoPtr<FEBase> ap_fb(FEBase::build(Dim-1, fet));
  base_fe = ap_fb.release();
}




// Desctructor
template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
InfFE<Dim,T_radial,T_map>::~InfFE ()
{
  // delete pointers, if necessary
  if (base_qrule != NULL)
    delete base_qrule;

  if (radial_qrule != NULL)
    delete radial_qrule;

  if (base_elem != NULL)
    delete base_elem;

  if (base_fe != NULL)
    delete base_fe;
}





template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
void InfFE<Dim,T_radial,T_map>:: attach_quadrature_rule (QBase* q)
{
  assert (q != NULL);

  Order base_int_order   = q->get_order();
//  Order radial_int_order = static_cast<Order>(
//      static_cast<unsigned int>(fe_type.radial_order) + 2 );
  Order radial_int_order = static_cast<Order>(
      2*(static_cast<unsigned int>(fe_type.radial_order) + 1) );
  // radial order rather conservative, may also work with other values...?
  // check this radial order again!!!!

  if (Dim != 1)
  {
    // build a Dim-1 quadrature rule of the type that we received
    AutoPtr<QBase> apq( QBase::build(q->type(), Dim-1, base_int_order) );
    base_qrule = apq.release();
    base_fe->attach_quadrature_rule(base_qrule);
  }

  // in radial direction, always use Gauss quadrature
  radial_qrule = new QGauss(1, radial_int_order);

  /* currently not used. But maybe helpful to store the QBase*
   * with which we initialized our own quadrature rules */
  qrule = q;
}





template <unsigned int Dim, FEFamily T_radial, InfMapType T_base>
void InfFE<Dim,T_radial,T_base>::update_base_elem (const Elem* inf_elem)
{
  if (base_elem != NULL)
    delete base_elem;
  base_elem = Base::build_elem(inf_elem);
}






template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
void InfFE<Dim,T_radial,T_map>::reinit(const Elem* inf_elem)
{
  // check InfFE quadrature rule and the new Elem* pointer
  assert (base_fe        != NULL);
  assert (base_fe->qrule != NULL);
  assert (base_fe->qrule == base_qrule);
  assert (radial_qrule   != NULL);
  assert (inf_elem       != NULL);


  // -----------------------------------------------------------------
  // update the type in accordance to the current cell
  // and reinit if the cell type has changed or (as in
  // the case of the hierarchics) the shape functions
  // depend on the particular element
  //
  // this is only necessary for the base! -- the radial
  // data only have to be reinitialized when the radial
  // order changes
  //
  // (note that FEBase::elem_type is initialized to INVALID_ELEM,
  //  so the following if statement is safe)
  if (  ( Dim != 1) &&
	(  (get_type() != inf_elem->type())  ||  
	   (fe_type.family == HIERARCHIC)  )  )
    {
      // store the new element type
      elem_type = inf_elem->type();


      // initialize the AutoPtr<Elem> for the base_fe
      update_base_elem(inf_elem);


      // initialize the base quadrature rule for the new element
      base_qrule->init( base_elem->type() );

      // initialize the shape functions in the base
      base_fe->init_base_shape_functions( base_fe->qrule, base_elem );

    }


  // -----------------------------------------------------------------
  // the radial part only needs to be re-initialized with
  // init_shape_functions() when the radial order changed
  if (current_fe_type.radial_order != fe_type.radial_order)
    {
      // update to the new radial order
      current_fe_type.radial_order = fe_type.radial_order;

      // the new element type has already been stored.
      // proceed with quadrature rule

      /*
       * initialize the radial quadrature rule with EDGE2
       * (the only ElemType that you should _not_ use is
       * INVALID_ELEM, since then, QBase::init() would think
       * it has already done the work!
       * But when the quadrature rule is initialized _once_,
       * even when the inf_elem->type() may change from
       * INFPRISM6 to INFHEX8, the radial data remains unchanged!
       */
      radial_qrule->init(EDGE2);

      // initialize the radial shape functions
      init_shape_functions (inf_elem);
    }


  // -----------------------------------------------------------------
  // Now that both the base and radial parts are properly initialized,
  // we only have to throw them together so that FEBase's::compute_map()
  // may directly be applied
  combine_base_radial();


  // Compute the map for this element.  In the future we can specify
  // different types of maps
  compute_map (qrule, inf_elem);


  // Compute the shape functions and the derivatives at all of the
  // quadrature points.  This part is dimension-independent
  compute_shape_functions ();


}






template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
void InfFE<Dim,T_radial,T_map>::init_shape_functions(const Elem* inf_elem)
{
  assert (qrule    != NULL);
  assert (inf_elem != NULL);

  // most of the data in here is separated as follows:
  //    base          everything in the base (handled by base_fe)
  //    radial        things in radial direction only
  //    total         the total of both, generally base x radial = total



  // -----------------------------------------------------------------
  // initialize most of the things related to mapping
  
  // The element type and order to use in the base map  
  const Order    base_mapping_order     ( base_elem->default_order() );
  const ElemType base_mapping_elem_type ( base_elem->type()          );

  // The order to use in the radial map (currently independent of the element type)
  const Order    radial_mapping_order   ( Radial::mapping_order() );    

  // Number of base shape functions used to construct the map
  // (Lagrange shape functions are used for mapping in the base)
  unsigned int n_base_mapping_shape_functions;

  // Note that the test used to be
  //if (Dim > 1)
  //  n_base_mapping_shape_functions = FE<1,LAGRANGE>::n_shape_functions (base_mapping_elem_type,
  //									base_mapping_order);
  //else
  //  n_base_mapping_shape_functions = 1;
  //
  // But that causes some compilers (icc 7.0, in particular) to complain.  Specifically,
  // it doesn't do the if-test at instantiation time, and there are undefined references
  // to FE<0,0> at link time.  Instead use the following (more redundant :-( ) test.  
  if (Dim == 1)
    n_base_mapping_shape_functions = 1;
  
  else if (Dim == 2)
    n_base_mapping_shape_functions = FE<1,LAGRANGE>::n_shape_functions (base_mapping_elem_type,
									base_mapping_order);
  else if (Dim == 3)
    n_base_mapping_shape_functions = FE<2,LAGRANGE>::n_shape_functions (base_mapping_elem_type,
									base_mapping_order);
  else
    error();
  

  // Note that Radial::n_mapping_shape_functions() is independent of the
  // element type
  const unsigned int n_radial_mapping_shape_functions =
	Radial::n_dofs(elem_type, radial_mapping_order);

  const unsigned int n_total_mapping_shape_functions = 
	n_radial_mapping_shape_functions * n_base_mapping_shape_functions;

  // Note that we don't use some variables concerning the quadrature
  // points in the base -- we don't need them, base_fe->init_shape_functions()
  // only used them
  



  // -----------------------------------------------------------------
  // initialize most of the things related to physical approximation
  
  // The order to use in radial approximation
  const Order    radial_approx_order ( fe_type.radial_order );

  unsigned int n_base_approx_shape_functions;
  if (Dim > 1)
    n_base_approx_shape_functions = base_fe->n_shape_functions();
  else
    n_base_approx_shape_functions = 1;

  const unsigned int n_radial_approx_shape_functions =
      Radial::n_dofs(elem_type, radial_approx_order);

  const unsigned int n_total_approx_shape_functions = 
      n_radial_approx_shape_functions * n_base_approx_shape_functions;

  // update class member field
  _n_total_approx_sf = n_total_approx_shape_functions;



  // The number and location of the radial quadrature points.
  // const unsigned int        n_radial_qp = qrule->n_points();
  // const std::vector<Point>&   radial_qp = qrule->get_points();

  const unsigned int        n_radial_qp = radial_qrule->n_points();
  const std::vector<Point>&   radial_qp = radial_qrule->get_points();

  const unsigned int        n_base_qp =  base_qrule->n_points();

  // The total number of quadrature points.
  const unsigned int        n_total_qp =  n_radial_qp * n_base_qp;

  
  // update class member field
  _n_total_qp = n_total_qp;


  // -----------------------------------------------------------------
  // resize the _base_ data fields

  // most of this has already been done in base_fe->init_shape_functions().
  // but e.g. we need the radial distance from the origin for _each_
  // base mapping node
  dist.resize(n_base_mapping_shape_functions);

 


  // -----------------------------------------------------------------
  // resize the _radial_ data fields
  
  // these fields are only a function of v (radial direction)
  mode.resize      (n_radial_approx_shape_functions);  // the radial polynomials (eval)
  dmodedv.resize   (n_radial_approx_shape_functions);

  som.resize       (n_radial_qp);  // the (1-v)/2 weight
  dsomdv.resize    (n_radial_qp);


  radial_map.resize    (n_radial_mapping_shape_functions);
  dradialdv_map.resize (n_radial_mapping_shape_functions);
  



  // -----------------------------------------------------------------
  // resize the _total_ data fields
  
  // the phase term varies with xi, eta and zeta(v): store it for _all_ qp
  //
  // when computing the phase, we need the base approximations
  // therefore, initialize the phase here, but evaluate it 
  // in combine_base_radial().
  //
  // the weight, though, is only needed at the radial quadrature points, n_radial_qp.
  // but for a uniform interface to the protected data fields
  // the weight data field (which are accessible from the outside) are expanded to n_total_qp.
  weight.resize      (n_total_qp);
  dweightdv.resize   (n_total_qp);
  dweight.resize     (n_total_qp);

  dphase.resize      (n_total_qp);
  dphasedxi.resize   (n_total_qp);
  dphasedeta.resize  (n_total_qp);
  dphasedzeta.resize (n_total_qp);
  



  // -----------------------------------------------------------------
  // InfFE's data fields phi, dphi, dphidx, phi_map etc hold the _total_
  // shape and mapping functions, respectively
  {
    phi.resize     (n_total_approx_shape_functions);
    dphi.resize    (n_total_approx_shape_functions);
    dphidx.resize  (n_total_approx_shape_functions);
    dphidy.resize  (n_total_approx_shape_functions);
    dphidz.resize  (n_total_approx_shape_functions);
    dphidxi.resize (n_total_approx_shape_functions);
    
    if (Dim > 1)
      dphideta.resize      (n_total_approx_shape_functions);
    
    if (Dim == 3)
      dphidzeta.resize     (n_total_approx_shape_functions);
      
    
    phi_map.resize         (n_total_mapping_shape_functions);
    dphidxi_map.resize     (n_total_mapping_shape_functions);
    
    if (Dim > 1)
      dphideta_map.resize  (n_total_mapping_shape_functions);
    
    if (Dim == 3)
      dphidzeta_map.resize (n_total_mapping_shape_functions);
  }



  // -----------------------------------------------------------------
  // collect all the for loops, where inner vectors are 
  // resized to the appropriate number of quadrature points
  {    
    for (unsigned int i=0; i<n_total_approx_shape_functions; i++)
      {
	phi[i].resize         (n_total_qp);
	dphi[i].resize        (n_total_qp);
	dphidx[i].resize      (n_total_qp);
	dphidy[i].resize      (n_total_qp);
	dphidz[i].resize      (n_total_qp);
	dphidxi[i].resize     (n_total_qp);
	   
	if (Dim > 1)
	  dphideta[i].resize  (n_total_qp);
	    
	   
	if (Dim == 3)	     
	  dphidzeta[i].resize (n_total_qp);
	     
      }
       
    for (unsigned int i=0; i<n_total_mapping_shape_functions; i++)
      {
	phi_map[i].resize         (n_total_qp);
	dphidxi_map[i].resize     (n_total_qp);
	   
	if (Dim > 1)
	  dphideta_map[i].resize  (n_total_qp);
	   
	if (Dim == 3)
	  dphidzeta_map[i].resize (n_total_qp);
      }
  }


  // these mapping shapes are required at the radial quadrature points
  for (unsigned int i=0; i<n_radial_mapping_shape_functions; i++)
    {
      radial_map[i].resize    (n_radial_qp);
      dradialdv_map[i].resize (n_radial_qp);
    }


  // these approximation shapes are required at the radial quadrature points
  for (unsigned int i=0; i<n_radial_approx_shape_functions; i++)
    {
      mode[i].resize    (n_radial_qp);
      dmodedv[i].resize (n_radial_qp);
    }


  // zero  the phase, since it is to be summed up
  for (unsigned int p=0; p<n_total_qp; p++)
    {
      dphasedxi[p]   = 0.;
      dphasedeta[p]  = 0.;
      dphasedzeta[p] = 0.;
    }




  // -----------------------------------------------------------------
  // -----------------------------------------------------------------
  // start computing the values -- note that we can only 
  // compute values which are _independent_ of the base shapes




  // -----------------------------------------------------------------
  // compute the radial distances
  for (unsigned int i=0; i<n_base_mapping_shape_functions; i++)
    {
      // this works, since the _base_ nodes are numbered in the 
      // same manner for the base element as for the infinite element


      dist[i] = Point( inf_elem->point(i) 
		       - inf_elem->point(i+n_base_mapping_shape_functions) ).size();


    }
  


  // -----------------------------------------------------------------
  // compute scalar values at radial quadrature points
  for (unsigned int p=0; p<n_radial_qp; p++)
    {
      som[p]       = Radial::decay       (radial_qp[p](0)); 
      dsomdv[p]    = Radial::decay_deriv (radial_qp[p](0)); 
    }



  // -----------------------------------------------------------------
  // compute scalar values at _all_ quadrature points  -- for uniform 
  // access from the outside to these fields
  for (unsigned int rp=0; rp<n_radial_qp; rp++)
    for (unsigned int bp=0; bp<n_base_qp; bp++)
      {
        weight   [ bp+rp*n_base_qp ] = Radial::D       (radial_qp[rp](0)); 
        dweightdv[ bp+rp*n_base_qp ] = Radial::D_deriv (radial_qp[rp](0));
      }
  
  

  // -----------------------------------------------------------------
  // evaluate the mode shapes in radial direction at radial quadrature points
  for (unsigned int i=0; i<n_radial_approx_shape_functions; i++)
    for (unsigned int p=0; p<n_radial_qp; p++)
      {
        mode[i][p]    = InfFE<Dim,T_radial,T_map>::eval       (radial_qp[p](0), radial_approx_order, i);
        dmodedv[i][p] = InfFE<Dim,T_radial,T_map>::eval_deriv (radial_qp[p](0), radial_approx_order, i);
      }



  // -----------------------------------------------------------------
  // evaluate the mapping functions in radial direction at radial quadrature points
  for (unsigned int i=0; i<n_radial_mapping_shape_functions; i++)
    for (unsigned int p=0; p<n_radial_qp; p++)
      {
	radial_map[i][p]    = InfFE<Dim,INFINITE_MAP,T_map>::eval       (radial_qp[p](0), radial_mapping_order, i);
	dradialdv_map[i][p] = InfFE<Dim,INFINITE_MAP,T_map>::eval_deriv (radial_qp[p](0), radial_mapping_order, i);
      }

 
}









template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
void InfFE<Dim,T_radial,T_map>::combine_base_radial()
{

  switch (Dim)
    {

      //------------------------------------------------------------
      // 1D
    case 1:
      {
	std::cout << "ERROR: Not implemented." << std::endl;
	error();
	return;
      }


      
      //------------------------------------------------------------
      // 2D
    case 2:
      {
	std::cout << "ERROR: Not implemented." << std::endl;
	error();
     	return;
      }


      
      //------------------------------------------------------------
      // 3D
    case 3:
      {
	// fast access to the approximation and mapping shapes of base_fe
	const std::vector<std::vector<Real> >& S  = base_fe->phi;
	const std::vector<std::vector<Real> >& Ss = base_fe->dphidxi;
	const std::vector<std::vector<Real> >& St = base_fe->dphideta;
	const std::vector<std::vector<Real> >& S_map  = base_fe->phi_map;
	const std::vector<std::vector<Real> >& Ss_map = base_fe->dphidxi_map;
	const std::vector<std::vector<Real> >& St_map = base_fe->dphideta_map;

	const unsigned int n_radial_qp = radial_qrule->n_points();
	const unsigned int n_base_qp   = base_qrule->  n_points();

	const unsigned int n_base_mapping_sf   = dist.size();
	const unsigned int n_radial_mapping_sf = radial_map.size();

	const unsigned int n_base_approx_sf   = base_fe->n_shape_functions();
	const unsigned int n_radial_approx_sf = Radial::n_dofs(elem_type, fe_type.radial_order);


	// compute the phase term derivatives
	{
  	  unsigned int tp=0;
	  for (unsigned int rp=0; rp<n_radial_qp; rp++)  // over radial qp's
	    for (unsigned int bp=0; bp<n_base_qp; bp++)  // over base qp's
	      {
	        // sum over all base shapes, to get the average distance
    	        for (unsigned int i=0; i<n_base_mapping_sf; i++)
	          {
	            dphasedxi[tp]   += Ss_map[i][bp] * dist[i] * radial_map   [1][rp];
		    dphasedeta[tp]  += St_map[i][bp] * dist[i] * radial_map   [1][rp];
		    dphasedzeta[tp] += S_map [i][bp] * dist[i] * dradialdv_map[1][rp];

#ifdef DEBUG
		    if (tp != (bp+rp*n_base_qp) )
		    {
		      std::cout << "ERROR: cannot count..." << std::endl
				<< " tp = " << tp 
				<< " bp = " << bp
				<< " rp = " << rp
				<< std::endl;	      
		      error();
		    }
#endif		    

		  }

		tp++;

	      } // loop radial and base qp's

	}


	// compute the overall approximation shape functions
	for (unsigned int rp=0; rp<n_radial_qp; rp++)  // over radial qp's
	  for (unsigned int bp=0; bp<n_base_qp; bp++)  // over base qp's
	    for (unsigned int ri=0; ri<n_radial_approx_sf; ri++)  // over radial approx shapes
	      for (unsigned int bi=0; bi<n_base_approx_sf; bi++)  // over base   approx shapes
	        {
		  // form the total shape function data fields
		  phi      [ bi+ri*n_base_approx_sf ][ bp+rp*n_base_qp ] = S [bi][bp] * mode[ri][rp] * som[ri];
		  dphidxi  [ bi+ri*n_base_approx_sf ][ bp+rp*n_base_qp ] = Ss[bi][bp] * mode[ri][rp] * som[ri];
		  dphideta [ bi+ri*n_base_approx_sf ][ bp+rp*n_base_qp ] = St[bi][bp] * mode[ri][rp] * som[ri];
		  dphidzeta[ bi+ri*n_base_approx_sf ][ bp+rp*n_base_qp ] = S [bi][bp] 
		      * (dmodedv[ri][rp] * som[ri] + mode[ri][rp] * dsomdv[ri]);
		}

	
	// compute the overall mapping functions
	for (unsigned int rp=0; rp<n_radial_qp; rp++)  // over radial qp's
	  for (unsigned int bp=0; bp<n_base_qp; bp++)  // over base qp's
	    for (unsigned int ri=0; ri<n_radial_mapping_sf; ri++)  // over radial mapping shapes
	      for (unsigned int bi=0; bi<n_base_mapping_sf; bi++)  // over base   mapping shapes
	        {
		  // form the total shape function data fields
		  phi_map      [ bi+ri*n_base_approx_sf ][ bp+rp*n_base_qp ] = S_map [bi][bp] * radial_map   [ri][rp];
		  dphidxi_map  [ bi+ri*n_base_approx_sf ][ bp+rp*n_base_qp ] = Ss_map[bi][bp] * radial_map   [ri][rp];
		  dphideta_map [ bi+ri*n_base_approx_sf ][ bp+rp*n_base_qp ] = St_map[bi][bp] * radial_map   [ri][rp];
		  dphidzeta_map[ bi+ri*n_base_approx_sf ][ bp+rp*n_base_qp ] = S_map [bi][bp] * dradialdv_map[ri][rp];
		}
			

	return;
      }


    default:
      error();
    }

  error();
  return;

}






template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
void InfFE<Dim,T_radial,T_map>::compute_shape_functions()
{
  assert (radial_qrule != NULL);

  const unsigned int n_total_qp  = _n_total_qp;


  //-------------------------------------------------------------------------
  // Compute the shape function values (and derivatives)
  // at the Quadrature points.  Note that the actual values
  // have already been computed via init_shape_functions

  // Compute the value of the derivative shape function i at quadrature point p
  switch (dim)
    {
      
    case 1:
      {
	std::cout << "ERROR: Not implemented." << std::endl;
	error();
	break;
      }

    case 2:
      {
	std::cout << "ERROR: Not implemented." << std::endl;
	error();
	break;
      }
    
    case 3:
      {
	// These are _all_ shape functions of this infinite element
	for (unsigned int i=0; i<phi.size(); i++)
	  for (unsigned int p=0; p<n_total_qp; p++)
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
	    }


	// This is the derivative of the phase term of this infinite element
	for (unsigned int p=0; p<n_total_qp; p++)
	  {
	    // the derivative of the phase term
	    dphase[p](0) = (dphasedxi[p]   * dxidx_map[p] +
			    dphasedeta[p]  * detadx_map[p] +
			    dphasedzeta[p] * dzetadx_map[p]);
		
	    dphase[p](1) = (dphasedxi[p]   * dxidy_map[p] +
			    dphasedeta[p]  * detady_map[p] +
			    dphasedzeta[p] * dzetady_map[p]);
		
	    dphase[p](2) = (dphasedxi[p]   * dxidz_map[p] +
			    dphasedeta[p]  * detadz_map[p] +
			    dphasedzeta[p] * dzetadz_map[p]);

	    // the derivative of the radial weight - varies only in radial direction,
	    // therefore dweightdxi = dweightdeta = 0.
	    dweight[p](0) = dweightdv[p] * dzetadx_map[p];
		
	    dweight[p](1) = dweightdv[p] * dzetady_map[p];
		
	    dweight[p](2) = dweightdv[p] * dzetadz_map[p];

	  }

	break;
      }



    default:
      {
	error();
      }
    }
}





//--------------------------------------------------------------
// Explicit instantiations
#include "inf_fe_instantiate_1D.h"
#include "inf_fe_instantiate_2D.h"
#include "inf_fe_instantiate_3D.h"



#endif //ifdef ENABLE_INFINITE_ELEMENTS

