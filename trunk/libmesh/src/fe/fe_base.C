// $Id: fe_base.C,v 1.15 2003-05-15 23:34:35 benkirk Exp $

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
#include "libmesh.h"
#include "quadrature.h"
#include "inf_fe.h"
#include "mesh_logging.h"



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
	    }
		   
	  case HIERARCHIC:
	    {
	      AutoPtr<FEBase> ap(new FE<1,HIERARCHIC>(fet));
	      return ap;
	    }
	    
	  case MONOMIAL:
	    {
	      AutoPtr<FEBase> ap(new FE<1,MONOMIAL>(fet));
	      return ap;
	    }
	    
	  default:
	    std::cout << "ERROR: Bad FEType.family= " << fet.family << std::endl;
	    error();
	  }
      }

      
      // 2D
    case 2:
      {
	switch (fet.family)
	  {
	  case LAGRANGE:
	    {
	      AutoPtr<FEBase> ap(new FE<2,LAGRANGE>(fet));
	      return ap;
	    }
	    
	  case HIERARCHIC:
	    {
	      AutoPtr<FEBase> ap(new FE<2,HIERARCHIC>(fet));
	      return ap;
	    }
	    
	  case MONOMIAL:
	    {
	      AutoPtr<FEBase> ap(new FE<2,MONOMIAL>(fet));
	      return ap;
	    }
	    
	  default:
	    std::cout << "ERROR: Bad FEType.family= " << fet.family << std::endl;
	    error();
	  }
      }

      
      // 3D
    case 3:
      {
	switch (fet.family)
	  {
	  case LAGRANGE:
	    {
	      AutoPtr<FEBase> ap(new FE<3,LAGRANGE>(fet));
	      return ap;
	    }
	    
	  case HIERARCHIC:
	    {
	      AutoPtr<FEBase> ap(new FE<3,HIERARCHIC>(fet));
	      return ap;
	    }
	    
	  case MONOMIAL:
	    {
	      AutoPtr<FEBase> ap(new FE<3,MONOMIAL>(fet));
	      return ap;
	    }
	    
	  default:
	    std::cout << "ERROR: Bad FEType.family= " << fet.family << std::endl;
	    error();
	  }
      }

    default:
      error();
    }

  error();
  AutoPtr<FEBase> ap(NULL);
  return ap;
}







#ifdef ENABLE_INFINITE_ELEMENTS


AutoPtr<FEBase> FEBase::build_InfFE (const unsigned int dim,
				     const FEType& fet)
{
  // The stupid AutoPtr<FEBase> ap(); return ap;
  // construct is required to satisfy IBM's xlC

  switch (dim)
    {

      // 1D
    case 1:
      {
	switch (fet.radial_family)
	  {
	  case INFINITE_MAP:
	    {
	      std::cerr << "ERROR: Don't build an infinite element " << std::endl
			<< " with FEFamily = " << fet.radial_family << std::endl;
	      error();
	    }

	  case JACOBI_20_00:
	    {
  	      switch (fet.inf_map)
	        {
		  case CARTESIAN:
		    {
		      AutoPtr<FEBase> ap(new InfFE<1,JACOBI_20_00,CARTESIAN>(fet));
		      return ap;
		    }
		  default:
		    std::cerr << "ERROR: Don't build an infinite element " << std::endl
			      << " with InfMapType = " << fet.inf_map << std::endl;
		    error();
		}
	    }

	  case JACOBI_30_00:
	    {
  	      switch (fet.inf_map)
	        {
		  case CARTESIAN:
		    {
		      AutoPtr<FEBase> ap(new InfFE<1,JACOBI_30_00,CARTESIAN>(fet));
		      return ap;
		    }
		  default:
		    std::cerr << "ERROR: Don't build an infinite element " << std::endl
			      << " with InfMapType = " << fet.inf_map << std::endl;
		    error();
		}
	    }

	  case LEGENDRE:
	    {
  	      switch (fet.inf_map)
	        {
		  case CARTESIAN:
		    {
		      AutoPtr<FEBase> ap(new InfFE<1,LEGENDRE,CARTESIAN>(fet));
		      return ap;
		    }
		  default:
		    std::cerr << "ERROR: Don't build an infinite element " << std::endl
			      << " with InfMapType = " << fet.inf_map << std::endl;
		    error();
		}
	    }

	  case LAGRANGE:
	    {
  	      switch (fet.inf_map)
	        {
		  case CARTESIAN:
		    {
		      AutoPtr<FEBase> ap(new InfFE<1,LAGRANGE,CARTESIAN>(fet));
		      return ap;
		    }
		  default:
		    std::cerr << "ERROR: Don't build an infinite element " << std::endl
			      << " with InfMapType = " << fet.inf_map << std::endl;
		    error();
		}
	    }


	    
	  default:
	    std::cerr << "ERROR: Bad FEType.radial_family= " << fet.radial_family << std::endl;
	    error();
	  }

      }

      


      // 2D
    case 2:
      {
	switch (fet.radial_family)
	  {
	  case INFINITE_MAP:
	    {
	      std::cerr << "ERROR: Don't build an infinite element " << std::endl
			<< " with FEFamily = " << fet.radial_family << std::endl;
	      error();
	    }

	  case JACOBI_20_00:
	    {
  	      switch (fet.inf_map)
	        {
		  case CARTESIAN:
		    {
		      AutoPtr<FEBase> ap(new InfFE<2,JACOBI_20_00,CARTESIAN>(fet));
		      return ap;
		    }
		  default:
		    std::cerr << "ERROR: Don't build an infinite element " << std::endl
			      << " with InfMapType = " << fet.inf_map << std::endl;
		    error();
		}
	    }

	  case JACOBI_30_00:
	    {
  	      switch (fet.inf_map)
	        {
		  case CARTESIAN:
		    {
		      AutoPtr<FEBase> ap(new InfFE<2,JACOBI_30_00,CARTESIAN>(fet));
		      return ap;
		    }
		  default:
		    std::cerr << "ERROR: Don't build an infinite element " << std::endl
			      << " with InfMapType = " << fet.inf_map << std::endl;
		    error();
		}
	    }

	  case LEGENDRE:
	    {
  	      switch (fet.inf_map)
	        {
		  case CARTESIAN:
		    {
		      AutoPtr<FEBase> ap(new InfFE<2,LEGENDRE,CARTESIAN>(fet));
		      return ap;
		    }
		  default:
		    std::cerr << "ERROR: Don't build an infinite element " << std::endl
			      << " with InfMapType = " << fet.inf_map << std::endl;
		    error();
		}
	    }

	  case LAGRANGE:
	    {
  	      switch (fet.inf_map)
	        {
		  case CARTESIAN:
		    {
		      AutoPtr<FEBase> ap(new InfFE<2,LAGRANGE,CARTESIAN>(fet));
		      return ap;
		    }
		  default:
		    std::cerr << "ERROR: Don't build an infinite element " << std::endl
			      << " with InfMapType = " << fet.inf_map << std::endl;
		    error();
		}
	    }


	    
	  default:
	    std::cerr << "ERROR: Bad FEType.radial_family= " << fet.radial_family << std::endl;
	    error();
	  }

      }

      


      // 3D
    case 3:
      {
	switch (fet.radial_family)
	  {
	  case INFINITE_MAP:
	    {
	      std::cerr << "ERROR: Don't build an infinite element " << std::endl
			<< " with FEFamily = " << fet.radial_family << std::endl;
	      error();
	    }

	  case JACOBI_20_00:
	    {
  	      switch (fet.inf_map)
	        {
		  case CARTESIAN:
		    {
		      AutoPtr<FEBase> ap(new InfFE<3,JACOBI_20_00,CARTESIAN>(fet));
		      return ap;
		    }
		  default:
		    std::cerr << "ERROR: Don't build an infinite element " << std::endl
			      << " with InfMapType = " << fet.inf_map << std::endl;
		    error();
		}
	    }

	  case JACOBI_30_00:
	    {
  	      switch (fet.inf_map)
	        {
		  case CARTESIAN:
		    {
		      AutoPtr<FEBase> ap(new InfFE<3,JACOBI_30_00,CARTESIAN>(fet));
		      return ap;
		    }
		  default:
		    std::cerr << "ERROR: Don't build an infinite element " << std::endl
			      << " with InfMapType = " << fet.inf_map << std::endl;
		    error();
		}
	    }

	  case LEGENDRE:
	    {
  	      switch (fet.inf_map)
	        {
		  case CARTESIAN:
		    {
		      AutoPtr<FEBase> ap(new InfFE<3,LEGENDRE,CARTESIAN>(fet));
		      return ap;
		    }
		  default:
		    std::cerr << "ERROR: Don't build an infinite element " << std::endl
			      << " with InfMapType = " << fet.inf_map << std::endl;
		    error();
		}
	    }

	  case LAGRANGE:
	    {
  	      switch (fet.inf_map)
	        {
		  case CARTESIAN:
		    {
		      AutoPtr<FEBase> ap(new InfFE<3,LAGRANGE,CARTESIAN>(fet));
		      return ap;
		    }
		  default:
		    std::cerr << "ERROR: Don't build an infinite element " << std::endl
			      << " with InfMapType = " << fet.inf_map << std::endl;
		    error();
		}
	    }


	    
	  default:
	    std::cerr << "ERROR: Bad FEType.radial_family= " << fet.radial_family << std::endl;
	    error();
	  }
      }

    default:
      error();
    }

  error();
  AutoPtr<FEBase> ap(NULL);
  return ap;
}




#endif // ifdef ENABLE_INFINITE_ELEMENTS











void FEBase::compute_shape_functions ()
{
  //-------------------------------------------------------------------------
  // Compute the shape function values (and derivatives)
  // at the Quadrature points.  Note that the actual values
  // have already been computed via init_shape_functions

  /**
   * Start logging the shape function computation
   */
  START_LOG("compute_shape_functions()", "FE");

  // Compute the value of the derivative shape function i at quadrature point p
  switch (dim)
    {
      
    case 1:
      {
	for (unsigned int i=0; i<phi.size(); i++)
	  for (unsigned int p=0; p<phi[i].size(); p++)
	    {
	      // dphi/dx    = (dphi/dxi)*(dxi/dx)
	      dphi[i][p](0) =
		dphidx[i][p] = dphidxi[i][p]*dxidx_map[p];
	      
	      dphi[i][p](1) = dphidy[i][p] = 0.;
	      dphi[i][p](2) = dphidz[i][p] = 0.;
	    }

	// All done
	break;
      }

    case 2:
      {
	for (unsigned int i=0; i<phi.size(); i++)
	  for (unsigned int p=0; p<phi[i].size(); p++)
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
	    }

	// All done
	break;
      }
    
    case 3:
      {
	for (unsigned int i=0; i<phi.size(); i++)
	  for (unsigned int p=0; p<phi[i].size(); p++)
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

	// All done
	break;
      }

    default:
      {
	error();
      }
    }
  
  /**
   * Stop logging the shape function computation
   */
  STOP_LOG("compute_shape_functions()", "FE");
}



bool FEBase::on_reference_element(const Point& p, const ElemType t, const Real eps)
{
  assert (eps >= 0.);
  
  const Real xi   = p(0);
  const Real eta  = p(1);
  const Real zeta = p(2);
  
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

	return false;
      }

      
    case TRI3:
    case TRI6:
      {
	// The reference triangle is isocoles
	// and is bound by xi=0, eta=0, and xi+eta=1.
	if ((xi  >= 0.-eps) &&
	    (eta >= 0.-eps) &&
	    ((xi + eta) <= 1.+eps))
	  return true;

	return false;
      }

      
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
		
	return false;
      }


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
		
	return false;
      }

      
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

	return false;
      }

    case PRISM6:
    case PRISM15:
    case PRISM18:
      {
	// Figure this one out...
	// inside the reference triange with zeta in [-1,1]
	if ((xi   >=  0.-eps) &&
	    (eta  >=  0.-eps) &&
	    (zeta >= -1.-eps) &&
	    (zeta <=  1.+eps) &&
	    ((xi + eta) <= 1.+eps))
	  return true;

	return false;
      }


    case PYRAMID5:
      {
	std::cerr << "BEN: Implement this you lazy bastard!"
		  << std::endl;
	error();

	return false;
      }

#ifdef ENABLE_INFINITE_ELEMENTS
    case INFHEX8:
      {      
	// The reference infhex8 is a [-1,1]^3.
	if ((xi   >= -1.-eps) &&
	    (xi   <=  1.+eps) &&
	    (eta  >= -1.-eps) &&
	    (eta  <=  1.+eps) &&
	    (zeta >= -1.-eps) &&
	    (zeta <=  1.+eps))
	  {
	    std::cout << "Found this point in infinite domain:\n";
	    p.print();
	    return true;
	  }

	return false;
      }
#endif

    default:
      std::cerr << "ERROR: Unknown element type " << t << std::endl;
      error();
    }

  // If we get here then the point is _not_ in the
  // reference element.   Better return false.
  
  return false;
}


