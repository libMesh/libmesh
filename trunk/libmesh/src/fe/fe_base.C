// $Id: fe_base.C,v 1.32 2005-07-22 16:31:52 roystgnr Exp $

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
#include "fe.h"
#include "inf_fe.h"
#include "libmesh_logging.h"



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
	  case CLOUGH:
	    {
	      AutoPtr<FEBase> ap(new FE<1,CLOUGH>(fet));
	      return ap;
	    }
	    
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
	    
#ifdef ENABLE_HIGHER_ORDER_SHAPES
	  case SZABAB:
	    {
	      AutoPtr<FEBase> ap(new FE<1,SZABAB>(fet));
	      return ap;
	    }

	  case BERNSTEIN:
	    {
	      AutoPtr<FEBase> ap(new FE<1,BERNSTEIN>(fet));
	      return ap;
	    }
#endif

	  case XYZ:
	    {
	      AutoPtr<FEBase> ap(new FEXYZ<1>(fet));
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
	  case CLOUGH:
	    {
	      AutoPtr<FEBase> ap(new FE<2,CLOUGH>(fet));
	      return ap;
	    }
	    
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
	    
#ifdef ENABLE_HIGHER_ORDER_SHAPES
	  case SZABAB:
	    {
	      AutoPtr<FEBase> ap(new FE<2,SZABAB>(fet));
	      return ap;
	    }

	  case BERNSTEIN:
	    {
	      AutoPtr<FEBase> ap(new FE<2,BERNSTEIN>(fet));
	      return ap;
	    }
#endif

	  case XYZ:
	    {
	      AutoPtr<FEBase> ap(new FEXYZ<2>(fet));
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
	  case CLOUGH:
	    {
	      std::cout << "ERROR: Clough-Tocher elements currently only support 1D and 2D" <<
		      std::endl;
	      error();
	    }
	    
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
	    
#ifdef ENABLE_HIGHER_ORDER_SHAPES
	  case SZABAB:
	    {
	      AutoPtr<FEBase> ap(new FE<3,SZABAB>(fet));
	      return ap;
	    }

	  case BERNSTEIN:
	    {
	      AutoPtr<FEBase> ap(new FE<3,BERNSTEIN>(fet));
	      return ap;
	    }
#endif

	  case XYZ:
	    {
	      AutoPtr<FEBase> ap(new FEXYZ<3>(fet));
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











void FEBase::compute_shape_functions (const Elem*)
{
  //-------------------------------------------------------------------------
  // Compute the shape function values (and derivatives)
  // at the Quadrature points.  Note that the actual values
  // have already been computed via init_shape_functions

  // Start logging the shape function computation
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
#ifdef ENABLE_SECOND_DERIVATIVES
	      d2phi[i][p](0,0) = d2phidx2[i][p] = 
		d2phidxi2[i][p]*dxidx_map[p]*dxidx_map[p];
#if DIM>1
	      d2phi[i][p](0,1) = d2phidxdy[i][p] = 
		d2phi[i][p](1,0) = 0.;
	      d2phi[i][p](1,1) = d2phidy2[i][p] = 0.;
#if DIM>2
	      d2phi[i][p](0,2) = d2phidxdz[i][p] =
		d2phi[i][p](2,0) = 0.;
	      d2phi[i][p](1,2) = d2phidydz[i][p] = 
		d2phi[i][p](2,1) = 0.;
	      d2phi[i][p](2,2) = d2phidz2[i][p] = 0.;
#endif
#endif
#endif
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
	      
	      // dphi/dz    = (dphi/dxi)*(dxi/dz) + (dphi/deta)*(deta/dz)
#if DIM == 3  
	      dphi[i][p](2) = // can only assign to the Z component if DIM==3
#endif
		dphidz[i][p] = (dphidxi[i][p]*dxidz_map[p] +
				dphideta[i][p]*detadz_map[p]);

#ifdef ENABLE_SECOND_DERIVATIVES
	      d2phi[i][p](0,0) = d2phidx2[i][p] = 
		d2phidxi2[i][p]*dxidx_map[p]*dxidx_map[p] +
		2*d2phidxideta[i][p]*dxidx_map[p]*detadx_map[p] +
		d2phideta2[i][p]*detadx_map[p]*detadx_map[p];
	      d2phi[i][p](0,1) = d2phidxdy[i][p] =
		d2phi[i][p](1,0) = 
		d2phidxi2[i][p]*dxidx_map[p]*dxidy_map[p] +
		d2phidxideta[i][p]*dxidx_map[p]*detady_map[p] +
		d2phideta2[i][p]*detadx_map[p]*detady_map[p] +
		d2phidxideta[i][p]*detadx_map[p]*dxidy_map[p];
	      d2phi[i][p](1,1) = d2phidy2[i][p] =
		d2phidxi2[i][p]*dxidy_map[p]*dxidy_map[p] +
		2*d2phidxideta[i][p]*dxidy_map[p]*detady_map[p] +
		d2phideta2[i][p]*detady_map[p]*detady_map[p];
#if DIM == 3  
	      d2phi[i][p](0,2) = d2phidxdz[i][p] = 
		d2phi[i][p](2,0) = 
		d2phidxi2[i][p]*dxidx_map[p]*dxidz_map[p] +
		d2phidxideta[i][p]*dxidx_map[p]*detadz_map[p] +
		d2phideta2[i][p]*detadx_map[p]*detadz_map[p] +
		d2phidxideta[i][p]*detadx_map[p]*dxidz_map[p];
	      d2phi[i][p](1,2) = d2phidydz[i][p] = 
		d2phi[i][p](2,1) =
		d2phidxi2[i][p]*dxidy_map[p]*dxidz_map[p] +
		d2phidxideta[i][p]*dxidy_map[p]*detadz_map[p] +
		d2phideta2[i][p]*detady_map[p]*detadz_map[p] +
		d2phidxideta[i][p]*detady_map[p]*dxidz_map[p];
	      d2phi[i][p](2,2) = d2phidz2[i][p] =
		d2phidxi2[i][p]*dxidz_map[p]*dxidz_map[p] +
		2*d2phidxideta[i][p]*dxidz_map[p]*detadz_map[p] +
		d2phideta2[i][p]*detadz_map[p]*detadz_map[p];
#endif
#endif
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

#ifdef ENABLE_SECOND_DERIVATIVES
	      d2phi[i][p](0,0) = d2phidx2[i][p] = 
		d2phidxi2[i][p]*dxidx_map[p]*dxidx_map[p] +
		2*d2phidxideta[i][p]*dxidx_map[p]*detadx_map[p] +
		2*d2phidxidzeta[i][p]*dxidx_map[p]*dzetadx_map[p] +
		2*d2phidetadzeta[i][p]*detadx_map[p]*dzetadx_map[p] +
		d2phideta2[i][p]*detadx_map[p]*detadx_map[p] +
		d2phidzeta2[i][p]*dzetadx_map[p]*dzetadx_map[p];
	      d2phi[i][p](0,1) = d2phidxdy[i][p] =
		d2phi[i][p](1,0) = 
		d2phidxi2[i][p]*dxidx_map[p]*dxidy_map[p] +
		d2phidxideta[i][p]*dxidx_map[p]*detady_map[p] +
		d2phidxidzeta[i][p]*dxidx_map[p]*dzetady_map[p] +
		d2phideta2[i][p]*detadx_map[p]*detady_map[p] +
		d2phidxideta[i][p]*detadx_map[p]*dxidy_map[p] +
		d2phidetadzeta[i][p]*detadx_map[p]*dzetady_map[p] +
		d2phidzeta2[i][p]*dzetadx_map[p]*dzetady_map[p] +
		d2phidxidzeta[i][p]*dzetadx_map[p]*dxidy_map[p] +
		d2phidetadzeta[i][p]*dzetadx_map[p]*detady_map[p];
	      d2phi[i][p](0,2) = d2phidxdz[i][p] = 
		d2phi[i][p](2,0) = 
		d2phidxi2[i][p]*dxidx_map[p]*dxidz_map[p] +
		d2phidxideta[i][p]*dxidx_map[p]*detadz_map[p] +
		d2phidxidzeta[i][p]*dxidx_map[p]*dzetadz_map[p] +
		d2phideta2[i][p]*detadx_map[p]*detadz_map[p] +
		d2phidxideta[i][p]*detadx_map[p]*dxidz_map[p] +
		d2phidetadzeta[i][p]*detadx_map[p]*dzetadz_map[p] +
		d2phidzeta2[i][p]*dzetadx_map[p]*dzetadz_map[p] +
		d2phidxidzeta[i][p]*dzetadx_map[p]*dxidz_map[p] +
		d2phidetadzeta[i][p]*dzetadx_map[p]*detadz_map[p];
	      d2phi[i][p](1,1) = d2phidy2[i][p] =
		d2phidxi2[i][p]*dxidy_map[p]*dxidy_map[p] +
		2*d2phidxideta[i][p]*dxidy_map[p]*detady_map[p] +
		2*d2phidxidzeta[i][p]*dxidy_map[p]*dzetady_map[p] +
		2*d2phidetadzeta[i][p]*detady_map[p]*dzetady_map[p] +
		d2phideta2[i][p]*detady_map[p]*detady_map[p] +
		d2phidzeta2[i][p]*dzetady_map[p]*dzetady_map[p];
	      d2phi[i][p](1,2) = d2phidydz[i][p] = 
		d2phi[i][p](2,1) =
		d2phidxi2[i][p]*dxidy_map[p]*dxidz_map[p] +
		d2phidxideta[i][p]*dxidy_map[p]*detadz_map[p] +
		d2phidxidzeta[i][p]*dxidy_map[p]*dzetadz_map[p] +
		d2phideta2[i][p]*detady_map[p]*detadz_map[p] +
		d2phidxideta[i][p]*detady_map[p]*dxidz_map[p] +
		d2phidetadzeta[i][p]*detady_map[p]*dzetadz_map[p] +
		d2phidzeta2[i][p]*dzetady_map[p]*dzetadz_map[p] +
		d2phidxidzeta[i][p]*dzetady_map[p]*dxidz_map[p] +
		d2phidetadzeta[i][p]*dzetady_map[p]*detadz_map[p];
	      d2phi[i][p](2,2) = d2phidz2[i][p] =
		d2phidxi2[i][p]*dxidz_map[p]*dxidz_map[p] +
		2*d2phidxideta[i][p]*dxidz_map[p]*detadz_map[p] +
		2*d2phidxidzeta[i][p]*dxidz_map[p]*dzetadz_map[p] +
		2*d2phidetadzeta[i][p]*detadz_map[p]*dzetadz_map[p] +
		d2phideta2[i][p]*detadz_map[p]*detadz_map[p] +
		d2phidzeta2[i][p]*dzetadz_map[p]*dzetadz_map[p];
#endif


	    }

	// All done
	break;
      }

    default:
      {
	error();
      }
    }
  
  // Stop logging the shape function computation
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
	    return true;
	  }
	return false;
      }

    case INFPRISM6:
      {      
	// inside the reference triange with zeta in [-1,1]
	if ((xi   >=  0.-eps) &&
	    (eta  >=  0.-eps) &&
	    (zeta >= -1.-eps) &&
	    (zeta <=  1.+eps) &&
	    ((xi + eta) <= 1.+eps))
	  {
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




void FEBase::print_JxW(std::ostream& os) const
{
  for (unsigned int i=0; i<JxW.size(); ++i)
    os << JxW[i] << std::endl;
}




void FEBase::print_phi(std::ostream& os) const
{
  for (unsigned int i=0; i<phi.size(); ++i)
    for (unsigned int j=0; j<phi[i].size(); ++j)
      os << " phi[" << i << "][" << j << "]=" << phi[i][j] << std::endl;
}




void FEBase::print_dphi(std::ostream& os) const
{
  for (unsigned int i=0; i<dphi.size(); ++i)
    for (unsigned int j=0; j<dphi[i].size(); ++j)
      os << " dphi[" << i << "][" << j << "]=" << dphi[i][j];
}



#ifdef ENABLE_SECOND_DERIVATIVES


void FEBase::print_d2phi(std::ostream& os) const
{
  for (unsigned int i=0; i<dphi.size(); ++i)
    for (unsigned int j=0; j<dphi[i].size(); ++j)
      os << " d2phi[" << i << "][" << j << "]=" << d2phi[i][j];
}

#endif




void FEBase::print_xyz(std::ostream& os) const
{
  for (unsigned int i=0; i<xyz.size(); ++i)
    os << xyz[i];
}




void FEBase::print_info(std::ostream& os) const
{
  os << "Shape functions at the Gauss pts." << std::endl;
  this->print_phi(os);
  
  os << "Shape function gradients at the Gauss pts." << std::endl;
  this->print_dphi(os);
  
  os << "XYZ locations of the Gauss pts." << std::endl;
  this->print_xyz(os);
  
  os << "Values of JxW at the Gauss pts." << std::endl;
  this->print_JxW(os);
}




std::ostream& operator << (std::ostream& os, const FEBase& fe)
{
  fe.print_info(os);
  return os;
}



