// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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
// For projection code:
#include "boundary_info.h"
#include "mesh_base.h"
#include "dense_matrix.h"
#include "dense_vector.h"
#include "dof_map.h"
#include "elem.h"
#include "fe_interface.h"
#include "numeric_vector.h"
#include "quadrature.h"
#include "quadrature_gauss.h"
#include "threads.h"

namespace libMesh
{



// ------------------------------------------------------------
// FEBase class members
AutoPtr<FEBase> FEBase::build (const unsigned int dim,
			       const FEType& fet)
{
  // The stupid AutoPtr<FEBase> ap(); return ap;
  // construct is required to satisfy IBM's xlC

  switch (dim)
    {
      // 0D
    case 0:
      {
	switch (fet.family)
	  {
	  case CLOUGH:
	    {
	      AutoPtr<FEBase> ap(new FE<0,CLOUGH>(fet));
	      return ap;
	    }
	    
	  case HERMITE:
	    {
	      AutoPtr<FEBase> ap(new FE<0,HERMITE>(fet));
	      return ap;
	    }
	    
	  case LAGRANGE:
	    {
	      AutoPtr<FEBase> ap(new FE<0,LAGRANGE>(fet));
	      return ap;
	    }
		   
	  case HIERARCHIC:
	    {
	      AutoPtr<FEBase> ap(new FE<0,HIERARCHIC>(fet));
	      return ap;
	    }
	    
	  case MONOMIAL:
	    {
	      AutoPtr<FEBase> ap(new FE<0,MONOMIAL>(fet));
	      return ap;
	    }
	    
#ifdef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES
	  case SZABAB:
	    {
	      AutoPtr<FEBase> ap(new FE<0,SZABAB>(fet));
	      return ap;
	    }

	  case BERNSTEIN:
	    {
	      AutoPtr<FEBase> ap(new FE<0,BERNSTEIN>(fet));
	      return ap;
	    }
#endif

	  case XYZ:
	    {
	      AutoPtr<FEBase> ap(new FEXYZ<0>(fet));
	      return ap;
	    }

          case SCALAR:
          {
	      AutoPtr<FEBase> ap(new FEScalar<0>(fet));
	      return ap;
          }

	  default:
	    libMesh::out << "ERROR: Bad FEType.family= " << fet.family << std::endl;
	    libmesh_error();
	  }
      }
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
	    
	  case HERMITE:
	    {
	      AutoPtr<FEBase> ap(new FE<1,HERMITE>(fet));
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
	    
#ifdef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES
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

          case SCALAR:
          {
	      AutoPtr<FEBase> ap(new FEScalar<1>(fet));
	      return ap;
          }

	  default:
	    libMesh::out << "ERROR: Bad FEType.family= " << fet.family << std::endl;
	    libmesh_error();
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
	    
	  case HERMITE:
	    {
	      AutoPtr<FEBase> ap(new FE<2,HERMITE>(fet));
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
	    
#ifdef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES
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

          case SCALAR:
          {
	      AutoPtr<FEBase> ap(new FEScalar<2>(fet));
	      return ap;
          }

	  default:
	    libMesh::out << "ERROR: Bad FEType.family= " << fet.family << std::endl;
	    libmesh_error();
	  }
      }

      
      // 3D
    case 3:
      {
	switch (fet.family)
	  {
	  case CLOUGH:
	    {
	      libMesh::out << "ERROR: Clough-Tocher elements currently only support 1D and 2D"
                            << std::endl;
	      libmesh_error();
	    }
	    
	  case HERMITE:
	    {
	      AutoPtr<FEBase> ap(new FE<3,HERMITE>(fet));
	      return ap;
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
	    
#ifdef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES
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

          case SCALAR:
          {
	      AutoPtr<FEBase> ap(new FEScalar<3>(fet));
	      return ap;
          }

	  default:
	    libMesh::out << "ERROR: Bad FEType.family= " << fet.family << std::endl;
	    libmesh_error();
	  }
      }

    default:
      libmesh_error();
    }

  libmesh_error();
  AutoPtr<FEBase> ap(NULL);
  return ap;
}







#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS


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
	      libMesh::err << "ERROR: Don't build an infinite element " << std::endl
			    << " with FEFamily = " << fet.radial_family << std::endl;
	      libmesh_error();
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
		    libMesh::err << "ERROR: Don't build an infinite element " << std::endl
			          << " with InfMapType = " << fet.inf_map << std::endl;
		    libmesh_error();
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
		    libMesh::err << "ERROR: Don't build an infinite element " << std::endl
			          << " with InfMapType = " << fet.inf_map << std::endl;
		    libmesh_error();
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
		    libMesh::err << "ERROR: Don't build an infinite element " << std::endl
			          << " with InfMapType = " << fet.inf_map << std::endl;
		    libmesh_error();
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
		    libMesh::err << "ERROR: Don't build an infinite element " << std::endl
			          << " with InfMapType = " << fet.inf_map << std::endl;
		    libmesh_error();
		}
	    }


	    
	  default:
	    libMesh::err << "ERROR: Bad FEType.radial_family= " << fet.radial_family << std::endl;
	    libmesh_error();
	  }

      }

      


      // 2D
    case 2:
      {
	switch (fet.radial_family)
	  {
	  case INFINITE_MAP:
	    {
	      libMesh::err << "ERROR: Don't build an infinite element " << std::endl
			    << " with FEFamily = " << fet.radial_family << std::endl;
	      libmesh_error();
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
		    libMesh::err << "ERROR: Don't build an infinite element " << std::endl
			          << " with InfMapType = " << fet.inf_map << std::endl;
		    libmesh_error();
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
		    libMesh::err << "ERROR: Don't build an infinite element " << std::endl
			          << " with InfMapType = " << fet.inf_map << std::endl;
		    libmesh_error();
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
		    libMesh::err << "ERROR: Don't build an infinite element " << std::endl
			          << " with InfMapType = " << fet.inf_map << std::endl;
		    libmesh_error();
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
		    libMesh::err << "ERROR: Don't build an infinite element " << std::endl
			          << " with InfMapType = " << fet.inf_map << std::endl;
		    libmesh_error();
		}
	    }


	    
	  default:
	    libMesh::err << "ERROR: Bad FEType.radial_family= " << fet.radial_family << std::endl;
	    libmesh_error();
	  }

      }

      


      // 3D
    case 3:
      {
	switch (fet.radial_family)
	  {
	  case INFINITE_MAP:
	    {
	      libMesh::err << "ERROR: Don't build an infinite element " << std::endl
			    << " with FEFamily = " << fet.radial_family << std::endl;
	      libmesh_error();
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
		    libMesh::err << "ERROR: Don't build an infinite element " << std::endl
			          << " with InfMapType = " << fet.inf_map << std::endl;
		    libmesh_error();
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
		    libMesh::err << "ERROR: Don't build an infinite element " << std::endl
			          << " with InfMapType = " << fet.inf_map << std::endl;
		    libmesh_error();
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
		    libMesh::err << "ERROR: Don't build an infinite element " << std::endl
			          << " with InfMapType = " << fet.inf_map << std::endl;
		    libmesh_error();
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
		    libMesh::err << "ERROR: Don't build an infinite element " << std::endl
			          << " with InfMapType = " << fet.inf_map << std::endl;
		    libmesh_error();
		}
	    }


	    
	  default:
	    libMesh::err << "ERROR: Bad FEType.radial_family= " << fet.radial_family << std::endl;
	    libmesh_error();
	  }
      }

    default:
      libmesh_error();
    }

  libmesh_error();
  AutoPtr<FEBase> ap(NULL);
  return ap;
}




#endif // ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS











void FEBase::compute_shape_functions (const Elem*)
{
  //-------------------------------------------------------------------------
  // Compute the shape function values (and derivatives)
  // at the Quadrature points.  Note that the actual values
  // have already been computed via init_shape_functions

  // Start logging the shape function computation
  START_LOG("compute_shape_functions()", "FE");

  calculations_started = true;

  // If the user forgot to request anything, we'll be safe and
  // calculate everything:
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  if (!calculate_phi && !calculate_dphi && !calculate_d2phi)
    calculate_phi = calculate_dphi = calculate_d2phi = true;
#else
  if (!calculate_phi && !calculate_dphi)
    calculate_phi = calculate_dphi = true;
#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES

  // Compute the value of the derivative shape function i at quadrature point p
  switch (dim)
    {

    case 0: // No derivatives in 0D
      break;

    case 1:
      {
	if (calculate_dphi)
	  for (unsigned int i=0; i<dphi.size(); i++)
	    for (unsigned int p=0; p<dphi[i].size(); p++)
	      {
	        // dphi/dx    = (dphi/dxi)*(dxi/dx)
	        dphi[i][p](0) =
		  dphidx[i][p] = dphidxi[i][p]*dxidx_map[p];
	      
#if LIBMESH_DIM>1
	        dphi[i][p](1) = dphidy[i][p] = 0.;
#endif
#if LIBMESH_DIM>2
	        dphi[i][p](2) = dphidz[i][p] = 0.;
#endif
	      }
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
	if (calculate_d2phi)
	  for (unsigned int i=0; i<d2phi.size(); i++)
	    for (unsigned int p=0; p<d2phi[i].size(); p++)
	      {
	        d2phi[i][p](0,0) = d2phidx2[i][p] = 
		  d2phidxi2[i][p]*dxidx_map[p]*dxidx_map[p];
#if LIBMESH_DIM>1
	        d2phi[i][p](0,1) = d2phidxdy[i][p] = 
		  d2phi[i][p](1,0) = 0.;
	        d2phi[i][p](1,1) = d2phidy2[i][p] = 0.;
#endif
#if LIBMESH_DIM>2
	        d2phi[i][p](0,2) = d2phidxdz[i][p] =
		  d2phi[i][p](2,0) = 0.;
	        d2phi[i][p](1,2) = d2phidydz[i][p] = 
		  d2phi[i][p](2,1) = 0.;
	        d2phi[i][p](2,2) = d2phidz2[i][p] = 0.;
#endif
	      }
#endif

	// All done
	break;
      }

    case 2:
      {
	if (calculate_dphi)
	  for (unsigned int i=0; i<dphi.size(); i++)
	    for (unsigned int p=0; p<dphi[i].size(); p++)
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
#if LIBMESH_DIM == 3  
	        dphi[i][p](2) = // can only assign to the Z component if LIBMESH_DIM==3
#endif
		dphidz[i][p] = (dphidxi[i][p]*dxidz_map[p] +
				dphideta[i][p]*detadz_map[p]);
	      }

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
	if (calculate_d2phi)
	  for (unsigned int i=0; i<d2phi.size(); i++)
	    for (unsigned int p=0; p<d2phi[i].size(); p++)
	      {
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
#if LIBMESH_DIM == 3  
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
	      }
#endif

	// All done
	break;
      }
    
    case 3:
      {
	if (calculate_dphi)
	  for (unsigned int i=0; i<dphi.size(); i++)
	    for (unsigned int p=0; p<dphi[i].size(); p++)
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

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
	if (calculate_d2phi)
	  for (unsigned int i=0; i<d2phi.size(); i++)
	    for (unsigned int p=0; p<d2phi[i].size(); p++)
	      {
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
	      }
#endif
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



bool FEBase::on_reference_element(const Point& p, const ElemType t, const Real eps)
{
  libmesh_assert (eps >= 0.);
  
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
	    //	    libMesh::out << "Strange Point:\n";
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
	libMesh::err << "BEN: Implement this you lazy bastard!"
		      << std::endl;
	libmesh_error();

	return false;
      }

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
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
      libMesh::err << "ERROR: Unknown element type " << t << std::endl;
      libmesh_error();
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



#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES


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



#ifdef LIBMESH_ENABLE_AMR

void FEBase::coarsened_dof_values(const NumericVector<Number> &old_vector,
                                  const DofMap &dof_map,
                                  const Elem *elem,
                                  DenseVector<Number> &Ue,
                                  const unsigned int var,
                                  const bool use_old_dof_indices)
{
  // Side/edge DOF indices
  std::vector<unsigned int> new_side_dofs, old_side_dofs;

  // FIXME: what about 2D shells in 3D space?
  unsigned int dim = elem->dim();

  // We use local FE objects for now
  // FIXME: we should use more, external objects instead for efficiency
  const FEType& base_fe_type = dof_map.variable_type(var);
  AutoPtr<FEBase> fe (FEBase::build(dim, base_fe_type));
  AutoPtr<FEBase> fe_coarse (FEBase::build(dim, base_fe_type));

  AutoPtr<QBase> qrule     (base_fe_type.default_quadrature_rule(dim));
  AutoPtr<QBase> qedgerule (base_fe_type.default_quadrature_rule(1));
  AutoPtr<QBase> qsiderule (base_fe_type.default_quadrature_rule(dim-1));
  std::vector<Point> coarse_qpoints;

  // The values of the shape functions at the quadrature
  // points
  const std::vector<std::vector<Real> >& phi_values =
    fe->get_phi();
  const std::vector<std::vector<Real> >& phi_coarse =
    fe_coarse->get_phi();

  // The gradients of the shape functions at the quadrature
  // points on the child element.
  const std::vector<std::vector<RealGradient> > *dphi_values =
    NULL;
  const std::vector<std::vector<RealGradient> > *dphi_coarse =
    NULL;

  const FEContinuity cont = fe->get_continuity();

  if (cont == C_ONE)
    {
      const std::vector<std::vector<RealGradient> >&
        ref_dphi_values = fe->get_dphi();
      dphi_values = &ref_dphi_values;
      const std::vector<std::vector<RealGradient> >&
        ref_dphi_coarse = fe_coarse->get_dphi();
      dphi_coarse = &ref_dphi_coarse;
    }

      // The Jacobian * quadrature weight at the quadrature points
      const std::vector<Real>& JxW =
        fe->get_JxW();

      // The XYZ locations of the quadrature points on the
      // child element
      const std::vector<Point>& xyz_values =
        fe->get_xyz();



  FEType fe_type = base_fe_type, temp_fe_type;
  const ElemType elem_type = elem->type();
  fe_type.order = static_cast<Order>(fe_type.order +
                                     elem->max_descendant_p_level());

  // Number of nodes on parent element
  const unsigned int n_nodes = elem->n_nodes();

  // Number of dofs on parent element
  const unsigned int new_n_dofs =
    FEInterface::n_dofs(dim, fe_type, elem_type);

  // Fixed vs. free DoFs on edge/face projections
  std::vector<char> dof_is_fixed(new_n_dofs, false); // bools
  std::vector<int> free_dof(new_n_dofs, 0);

  DenseMatrix<Real> Ke;
  DenseVector<Number> Fe;
  Ue.resize(new_n_dofs); Ue.zero();


  // When coarsening, in general, we need a series of
  // projections to ensure a unique and continuous
  // solution.  We start by interpolating nodes, then
  // hold those fixed and project edges, then
  // hold those fixed and project faces, then
  // hold those fixed and project interiors

  // Copy node values first
  {
  std::vector<unsigned int> node_dof_indices;
  if (use_old_dof_indices)
    dof_map.old_dof_indices (elem, node_dof_indices, var);
  else
    dof_map.dof_indices (elem, node_dof_indices, var);

  unsigned int current_dof = 0;
  for (unsigned int n=0; n!= n_nodes; ++n)
    {
      // FIXME: this should go through the DofMap,
      // not duplicate dof_indices code badly!
      const unsigned int my_nc =
        FEInterface::n_dofs_at_node (dim, fe_type,
                                     elem_type, n);
      if (!elem->is_vertex(n))
        {
          current_dof += my_nc;
          continue;
        }

      temp_fe_type = base_fe_type;
      // We're assuming here that child n shares vertex n,
      // which is wrong on non-simplices right now
      // ... but this code isn't necessary except on elements
      // where p refinement creates more vertex dofs; we have
      // no such elements yet.
/*
      if (elem->child(n)->p_level() < elem->p_level())
        {
          temp_fe_type.order = 
            static_cast<Order>(temp_fe_type.order +
                               elem->child(n)->p_level());
        }
*/
      const unsigned int nc =
        FEInterface::n_dofs_at_node (dim, temp_fe_type,
                                     elem_type, n);
      for (unsigned int i=0; i!= nc; ++i)
        {
          Ue(current_dof) =
            old_vector(node_dof_indices[current_dof]);
          dof_is_fixed[current_dof] = true;
          current_dof++;
        }
    }
  }

  // In 3D, project any edge values next
  if (dim > 2 && cont != DISCONTINUOUS)
    for (unsigned int e=0; e != elem->n_edges(); ++e)
      {
        FEInterface::dofs_on_edge(elem, dim, fe_type,
                                  e, new_side_dofs);

        // Some edge dofs are on nodes and already
        // fixed, others are free to calculate
        unsigned int free_dofs = 0;
        for (unsigned int i=0; i !=
             new_side_dofs.size(); ++i)
          if (!dof_is_fixed[new_side_dofs[i]])
            free_dof[free_dofs++] = i;
        Ke.resize (free_dofs, free_dofs); Ke.zero();
        Fe.resize (free_dofs); Fe.zero();
        // The new edge coefficients
        DenseVector<Number> Uedge(free_dofs);

        // Add projection terms from each child sharing
        // this edge
        for (unsigned int c=0; c != elem->n_children();
             ++c)
          {
            if (!elem->is_child_on_edge(c,e))
              continue;
            Elem *child = elem->child(c);

            std::vector<unsigned int> child_dof_indices;
            if (use_old_dof_indices)
              dof_map.old_dof_indices (child,
                child_dof_indices, var);
            else
              dof_map.dof_indices (child,
                child_dof_indices, var);
            const unsigned int child_n_dofs = child_dof_indices.size();

            temp_fe_type = base_fe_type;
            temp_fe_type.order = 
              static_cast<Order>(temp_fe_type.order +
                                 child->p_level());

            FEInterface::dofs_on_edge(child, dim,
              temp_fe_type, e, old_side_dofs);

            // Initialize both child and parent FE data
            // on the child's edge
            fe->attach_quadrature_rule (qedgerule.get());
            fe->edge_reinit (child, e);
            const unsigned int n_qp = qedgerule->n_points();

            FEInterface::inverse_map (dim, fe_type, elem,
                            xyz_values, coarse_qpoints);

            fe_coarse->reinit(elem, &coarse_qpoints);

            // Loop over the quadrature points
            for (unsigned int qp=0; qp<n_qp; qp++)
              {
                // solution value at the quadrature point
                Number fineval = libMesh::zero;
                // solution grad at the quadrature point
                Gradient finegrad;

                // Sum the solution values * the DOF
                // values at the quadrature point to
                // get the solution value and gradient.
                for (unsigned int i=0; i<child_n_dofs;
                     i++)
                  {
                    fineval +=
                      (old_vector(child_dof_indices[i])*
                      phi_values[i][qp]);
                    if (cont == C_ONE)
                      finegrad.add_scaled((*dphi_values)[i][qp],
                                          old_vector(child_dof_indices[i]));
                  }

                // Form edge projection matrix
                for (unsigned int sidei=0, freei=0; 
                     sidei != new_side_dofs.size();
                     ++sidei)
                  {
                    unsigned int i = new_side_dofs[sidei];
                    // fixed DoFs aren't test functions
                    if (dof_is_fixed[i])
                      continue;
                    for (unsigned int sidej=0, freej=0;
                         sidej != new_side_dofs.size();
                         ++sidej)
                      {
                        unsigned int j =
                          new_side_dofs[sidej];
                        if (dof_is_fixed[j])
                          Fe(freei) -=
                            phi_coarse[i][qp] *
                            phi_coarse[j][qp] * JxW[qp] *
                            Ue(j);
                        else
                          Ke(freei,freej) +=
                            phi_coarse[i][qp] *
                            phi_coarse[j][qp] * JxW[qp];
                        if (cont == C_ONE)
                          {
                            if (dof_is_fixed[j])
                              Fe(freei) -=
                                ((*dphi_coarse)[i][qp] *
                                 (*dphi_coarse)[j][qp]) *
                                JxW[qp] *
                                Ue(j);
                            else
                              Ke(freei,freej) +=
                                ((*dphi_coarse)[i][qp] *
                                 (*dphi_coarse)[j][qp])
                                * JxW[qp];
                          }
                        if (!dof_is_fixed[j])
                          freej++;
                      }
                    Fe(freei) += phi_coarse[i][qp] *
                                 fineval * JxW[qp];
                    if (cont == C_ONE)
                      Fe(freei) +=
                        (finegrad * (*dphi_coarse)[i][qp]) * JxW[qp];
                    freei++;
                  }
              }
          }
        Ke.cholesky_solve(Fe, Uedge);

        // Transfer new edge solutions to element
        for (unsigned int i=0; i != free_dofs; ++i)
          {
            Number &ui = Ue(new_side_dofs[free_dof[i]]);
            libmesh_assert(std::abs(ui) < TOLERANCE ||
                   std::abs(ui - Uedge(i)) < TOLERANCE);
            ui = Uedge(i);
            dof_is_fixed[new_side_dofs[free_dof[i]]] =
              true;
          }
      }
   
  // Project any side values (edges in 2D, faces in 3D)
  if (dim > 1 && cont != DISCONTINUOUS)
    for (unsigned int s=0; s != elem->n_sides(); ++s)
      {
        FEInterface::dofs_on_side(elem, dim, fe_type,
                                  s, new_side_dofs);

        // Some side dofs are on nodes/edges and already
        // fixed, others are free to calculate
        unsigned int free_dofs = 0;
        for (unsigned int i=0; i !=
             new_side_dofs.size(); ++i)
          if (!dof_is_fixed[new_side_dofs[i]])
            free_dof[free_dofs++] = i;
        Ke.resize (free_dofs, free_dofs); Ke.zero();
        Fe.resize (free_dofs); Fe.zero();
        // The new side coefficients
        DenseVector<Number> Uside(free_dofs);

        // Add projection terms from each child sharing
        // this side
        for (unsigned int c=0; c != elem->n_children();
             ++c)
          {
            if (!elem->is_child_on_side(c,s))
              continue;
            Elem *child = elem->child(c);

            std::vector<unsigned int> child_dof_indices;
            if (use_old_dof_indices)
              dof_map.old_dof_indices (child,
                child_dof_indices, var);
            else
              dof_map.dof_indices (child,
                child_dof_indices, var);
            const unsigned int child_n_dofs = child_dof_indices.size();

            temp_fe_type = base_fe_type;
            temp_fe_type.order = 
              static_cast<Order>(temp_fe_type.order +
                                 child->p_level());

            FEInterface::dofs_on_side(child, dim,
              temp_fe_type, s, old_side_dofs);

            // Initialize both child and parent FE data
            // on the child's side
            fe->attach_quadrature_rule (qsiderule.get());
            fe->reinit (child, s);
            const unsigned int n_qp = qsiderule->n_points();

            FEInterface::inverse_map (dim, fe_type, elem,
                            xyz_values, coarse_qpoints);

            fe_coarse->reinit(elem, &coarse_qpoints);

            // Loop over the quadrature points
            for (unsigned int qp=0; qp<n_qp; qp++)
              {
                // solution value at the quadrature point
                Number fineval = libMesh::zero;
                // solution grad at the quadrature point
                Gradient finegrad;

                // Sum the solution values * the DOF
                // values at the quadrature point to
                // get the solution value and gradient.
                for (unsigned int i=0; i<child_n_dofs;
                     i++)
                  {
                    fineval +=
                      (old_vector(child_dof_indices[i])*
                      phi_values[i][qp]);
                    if (cont == C_ONE)
                      finegrad.add_scaled((*dphi_values)[i][qp],
                                          old_vector(child_dof_indices[i]));
                  }

                // Form side projection matrix
                for (unsigned int sidei=0, freei=0;
                     sidei != new_side_dofs.size();
                     ++sidei)
                  {
                    unsigned int i = new_side_dofs[sidei];
                    // fixed DoFs aren't test functions
                    if (dof_is_fixed[i])
                      continue;
                    for (unsigned int sidej=0, freej=0;
                         sidej != new_side_dofs.size();
                         ++sidej)
                      {
                        unsigned int j =
                          new_side_dofs[sidej];
                        if (dof_is_fixed[j])
                          Fe(freei) -=
                            phi_coarse[i][qp] *
                            phi_coarse[j][qp] * JxW[qp] *
                            Ue(j);
                        else
                          Ke(freei,freej) +=
                            phi_coarse[i][qp] *
                            phi_coarse[j][qp] * JxW[qp];
                        if (cont == C_ONE)
                          {
                            if (dof_is_fixed[j])
                              Fe(freei) -=
                                ((*dphi_coarse)[i][qp] *
                                 (*dphi_coarse)[j][qp]) *
                                JxW[qp] *
                                Ue(j);
                            else
                              Ke(freei,freej) +=
                                ((*dphi_coarse)[i][qp] *
                                 (*dphi_coarse)[j][qp])
                                * JxW[qp];
                          }
                        if (!dof_is_fixed[j])
                          freej++;
                      }
                    Fe(freei) += (fineval * phi_coarse[i][qp]) * JxW[qp];
                    if (cont == C_ONE)
                      Fe(freei) +=
                        (finegrad * (*dphi_coarse)[i][qp]) * JxW[qp];
                    freei++;
                  }
              }
          }
        Ke.cholesky_solve(Fe, Uside);

        // Transfer new side solutions to element
        for (unsigned int i=0; i != free_dofs; ++i)
          {
            Number &ui = Ue(new_side_dofs[free_dof[i]]);
            libmesh_assert(std::abs(ui) < TOLERANCE ||
                   std::abs(ui - Uside(i)) < TOLERANCE);
            ui = Uside(i);
            dof_is_fixed[new_side_dofs[free_dof[i]]] =
              true;
          }
      }

  // Project the interior values, finally

  // Some interior dofs are on nodes/edges/sides and
  // already fixed, others are free to calculate
  unsigned int free_dofs = 0;
  for (unsigned int i=0; i != new_n_dofs; ++i)
    if (!dof_is_fixed[i])
      free_dof[free_dofs++] = i;
  Ke.resize (free_dofs, free_dofs); Ke.zero();
  Fe.resize (free_dofs); Fe.zero();
  // The new interior coefficients
  DenseVector<Number> Uint(free_dofs);

  // Add projection terms from each child
  for (unsigned int c=0; c != elem->n_children(); ++c)
    {
      Elem *child = elem->child(c);

      std::vector<unsigned int> child_dof_indices;
      if (use_old_dof_indices)
        dof_map.old_dof_indices (child,
          child_dof_indices, var);
      else
        dof_map.dof_indices (child,
          child_dof_indices, var);
      const unsigned int child_n_dofs = child_dof_indices.size();

      // Initialize both child and parent FE data
      // on the child's quadrature points
      fe->attach_quadrature_rule (qrule.get());
      fe->reinit (child);
      const unsigned int n_qp = qrule->n_points();

      FEInterface::inverse_map (dim, fe_type, elem,
        xyz_values, coarse_qpoints);

      fe_coarse->reinit(elem, &coarse_qpoints);

      // Loop over the quadrature points
      for (unsigned int qp=0; qp<n_qp; qp++)
        {
          // solution value at the quadrature point              
          Number fineval = libMesh::zero;
          // solution grad at the quadrature point              
          Gradient finegrad;

          // Sum the solution values * the DOF
          // values at the quadrature point to
          // get the solution value and gradient.
          for (unsigned int i=0; i<child_n_dofs; i++)
            {
              fineval +=
                (old_vector(child_dof_indices[i])*
                 phi_values[i][qp]);
              if (cont == C_ONE)
                finegrad.add_scaled((*dphi_values)[i][qp],
                                    old_vector(child_dof_indices[i]));
            }

          // Form interior projection matrix
          for (unsigned int i=0, freei=0;
               i != new_n_dofs; ++i)
            {
              // fixed DoFs aren't test functions
              if (dof_is_fixed[i])
                continue;
              for (unsigned int j=0, freej=0; j !=
                   new_n_dofs; ++j)
                {
                  if (dof_is_fixed[j])
                    Fe(freei) -=
                      phi_coarse[i][qp] *
                      phi_coarse[j][qp] * JxW[qp] *
                      Ue(j);
                  else
                    Ke(freei,freej) +=
                      phi_coarse[i][qp] *
                      phi_coarse[j][qp] * JxW[qp];
                  if (cont == C_ONE)
                    {
                      if (dof_is_fixed[j])
                        Fe(freei) -=
                          ((*dphi_coarse)[i][qp] *
                           (*dphi_coarse)[j][qp]) *
                          JxW[qp] * Ue(j);
                      else
                        Ke(freei,freej) +=
                          ((*dphi_coarse)[i][qp] *
                           (*dphi_coarse)[j][qp]) * JxW[qp];
                    }
                  if (!dof_is_fixed[j])
                    freej++;
                }
              Fe(freei) += phi_coarse[i][qp] * fineval *
                           JxW[qp];
              if (cont == C_ONE)
                Fe(freei) += (finegrad * (*dphi_coarse)[i][qp]) * JxW[qp];
              freei++;
            }
        }
    }
  Ke.cholesky_solve(Fe, Uint);

  // Transfer new interior solutions to element
  for (unsigned int i=0; i != free_dofs; ++i)
    {
      Number &ui = Ue(free_dof[i]);
      libmesh_assert(std::abs(ui) < TOLERANCE ||
             std::abs(ui - Uint(i)) < TOLERANCE);
      ui = Uint(i);
      dof_is_fixed[free_dof[i]] = true;
    }

  // Make sure every DoF got reached!
  for (unsigned int i=0; i != new_n_dofs; ++i)
    libmesh_assert(dof_is_fixed[i]);
}



void FEBase::compute_proj_constraints (DofConstraints &constraints,
				       DofMap &dof_map,
				       const unsigned int variable_number,
				       const Elem* elem)
{
  libmesh_assert (elem != NULL);

  const unsigned int Dim = elem->dim();

  // Only constrain elements in 2,3D.
  if (Dim == 1)
    return;

  // Only constrain active elements with this method
  if (!elem->active())
    return;

  const FEType& base_fe_type = dof_map.variable_type(variable_number);

  // Construct FE objects for this element and its neighbors.
  AutoPtr<FEBase> my_fe (FEBase::build(Dim, base_fe_type));
  const FEContinuity cont = my_fe->get_continuity();

  // We don't need to constrain discontinuous elements
  if (cont == DISCONTINUOUS)
    return;
  libmesh_assert (cont == C_ZERO || cont == C_ONE);

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
	    unsigned int s_neigh = neigh->which_neighbor_am_i(elem);
            libmesh_assert (s_neigh < neigh->n_neighbors());

            // Find the minimum p level; we build the h constraint
            // matrix with this and then constrain away all higher p
            // DoFs.
            libmesh_assert(neigh->active());
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
	    
	    // This function gets called element-by-element, so there
	    // will be a lot of memory allocation going on.  We can 
	    // at least minimize this for the case of the dof indices
	    // by efficiently preallocating the requisite storage.
	    // n_nodes is not necessarily n_dofs, but it is better
	    // than nothing!
	    my_dof_indices.reserve    (elem->n_nodes());
	    neigh_dof_indices.reserve (neigh->n_nodes());

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
	    FEInterface::dofs_on_side(elem,  Dim, base_fe_type, s,       my_side_dofs);
	    FEInterface::dofs_on_side(neigh, Dim, base_fe_type, s_neigh, neigh_side_dofs);

            // We're done with functions that examine Elem::p_level(),
            // so let's unhack those levels
            if (elem->p_level() != old_elem_level)
              (const_cast<Elem *>(elem))->hack_p_level(old_elem_level);
            if (neigh->p_level() != old_neigh_level)
              (const_cast<Elem *>(neigh))->hack_p_level(old_neigh_level);

	    const unsigned int n_side_dofs = my_side_dofs.size();
	    libmesh_assert(n_side_dofs == neigh_side_dofs.size());

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
                libmesh_assert(their_dof_g != DofObject::invalid_id);
	        for (unsigned int js = 0; js != n_side_dofs; ++js)
	          {
	            const unsigned int j = my_side_dofs[js];
	            const unsigned int my_dof_g = my_dof_indices[j];
                    libmesh_assert(my_dof_g != DofObject::invalid_id);
		    const Real their_dof_value = Ue[is](js);
		    if (their_dof_g == my_dof_g)
		      {
		        libmesh_assert(std::abs(their_dof_value-1.) < 1.e-5);
		        for (unsigned int k = 0; k != n_side_dofs; ++k)
		          libmesh_assert(k == is || std::abs(Ue[k](js)) < 1.e-5);
		        continue;
		      }
		    if (std::abs(their_dof_value) < 1.e-5)
		      continue;

		    // since we may be running this method concurretly 
		    // on multiple threads we need to acquire a lock 
		    // before modifying the shared constraint_row object.
		    {
		      Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);

		      DofConstraintRow& constraint_row =
			constraints[my_dof_g];
		      
		      constraint_row.insert(std::make_pair(their_dof_g,
							   their_dof_value));
		    }
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
            libmesh_assert(my_fe->is_hierarchic());
            dof_map.constrain_p_dofs(variable_number, elem,
                                     s, min_p_level);
          }
      }
}

#endif // #ifdef LIBMESH_ENABLE_AMR



#ifdef LIBMESH_ENABLE_PERIODIC

void FEBase::compute_periodic_constraints (DofConstraints &constraints,
				           DofMap &dof_map,
                                           PeriodicBoundaries &boundaries,
                                           const MeshBase &mesh,
				           const unsigned int variable_number,
				           const Elem* elem)
{
  // Only bother if we truly have periodic boundaries
  if (boundaries.empty())
    return;

  libmesh_assert (elem != NULL);
  
  // Only constrain active elements with this method
  if (!elem->active())
    return;

  const unsigned int Dim = elem->dim();
  
  const FEType& base_fe_type = dof_map.variable_type(variable_number);

  // Construct FE objects for this element and its pseudo-neighbors.
  AutoPtr<FEBase> my_fe (FEBase::build(Dim, base_fe_type));
  const FEContinuity cont = my_fe->get_continuity();

  // We don't need to constrain discontinuous elements
  if (cont == DISCONTINUOUS)
    return;
  libmesh_assert (cont == C_ZERO || cont == C_ONE);

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
    {
      if (elem->neighbor(s))
        continue;

      unsigned int boundary_id = mesh.boundary_info->boundary_id(elem, s);
      PeriodicBoundary *periodic = boundaries.boundary(boundary_id);
      if (periodic)
        {
          // Get pointers to the element's neighbor.
          const Elem* neigh = boundaries.neighbor(boundary_id, mesh, elem, s);

          // h refinement constraints:
          // constrain dofs shared between
          // this element and ones as coarse
          // as or coarser than this element.
          if (neigh->level() <= elem->level()) 
            {
	      unsigned int s_neigh = 
                mesh.boundary_info->side_with_boundary_id (neigh, periodic->pairedboundary);
              libmesh_assert(s_neigh != libMesh::invalid_uint);

#ifdef LIBMESH_ENABLE_AMR
              // Find the minimum p level; we build the h constraint
              // matrix with this and then constrain away all higher p
              // DoFs.
              libmesh_assert(neigh->active());
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
#endif // #ifdef LIBMESH_ENABLE_AMR

	      my_fe->reinit(elem, s);

	      dof_map.dof_indices (elem, my_dof_indices,
			           variable_number);
	      dof_map.dof_indices (neigh, neigh_dof_indices,
			           variable_number);

	      const unsigned int n_qp = my_qface.n_points();

              // Translate the quadrature points over to the
              // neighbor's boundary
              std::vector<Point> neigh_point = q_point;
              for (unsigned int i=0; i != neigh_point.size(); ++i)
                neigh_point[i] += periodic->translation_vector;

	      FEInterface::inverse_map (Dim, base_fe_type, neigh,
                                        neigh_point, neigh_qface);

	      neigh_fe->reinit(neigh, &neigh_qface);

	      // We're only concerned with DOFs whose values (and/or first
	      // derivatives for C1 elements) are supported on side nodes
	      FEInterface::dofs_on_side(elem, Dim, base_fe_type, s, my_side_dofs);
	      FEInterface::dofs_on_side(neigh, Dim, base_fe_type, s_neigh, neigh_side_dofs);

              // We're done with functions that examine Elem::p_level(),
              // so let's unhack those levels
#ifdef LIBMESH_ENABLE_AMR
              if (elem->p_level() != old_elem_level)
                (const_cast<Elem *>(elem))->hack_p_level(old_elem_level);
              if (neigh->p_level() != old_neigh_level)
                (const_cast<Elem *>(neigh))->hack_p_level(old_neigh_level);
#endif // #ifdef LIBMESH_ENABLE_AMR

	      const unsigned int n_side_dofs = my_side_dofs.size();
	      libmesh_assert(n_side_dofs == neigh_side_dofs.size());

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

              // Make sure we're not adding recursive constraints
              // due to the redundancy in the way we add periodic
              // boundary constraints
              std::vector<bool> recursive_constraint(n_side_dofs, false);

	      for (unsigned int is = 0; is != n_side_dofs; ++is)
	        {
	          const unsigned int i = neigh_side_dofs[is];
	          const unsigned int their_dof_g = neigh_dof_indices[i];
                  libmesh_assert(their_dof_g != DofObject::invalid_id);

                  if (!dof_map.is_constrained_dof(their_dof_g))
                    continue;

                  DofConstraintRow& their_constraint_row =
                    constraints[their_dof_g];

	          for (unsigned int js = 0; js != n_side_dofs; ++js)
	            {
	              const unsigned int j = my_side_dofs[js];
	              const unsigned int my_dof_g = my_dof_indices[j];
                      libmesh_assert(my_dof_g != DofObject::invalid_id);

                      if (their_constraint_row.count(my_dof_g))
                        recursive_constraint[js] = true;
	            }
                }
	      for (unsigned int js = 0; js != n_side_dofs; ++js)
	        {
                  if (recursive_constraint[js])
                    continue;

	          const unsigned int j = my_side_dofs[js];
	          const unsigned int my_dof_g = my_dof_indices[j];
                  libmesh_assert(my_dof_g != DofObject::invalid_id);

                  if (dof_map.is_constrained_dof(my_dof_g))
                    continue;

	          for (unsigned int is = 0; is != n_side_dofs; ++is)
	            {
	              const unsigned int i = neigh_side_dofs[is];
	              const unsigned int their_dof_g = neigh_dof_indices[i];
                      libmesh_assert(their_dof_g != DofObject::invalid_id);

		      const Real their_dof_value = Ue[is](js);
		      if (their_dof_g == my_dof_g)
		        {
		          libmesh_assert(std::abs(their_dof_value-1.) < 1.e-5);
		          for (unsigned int k = 0; k != n_side_dofs; ++k)
		            libmesh_assert(k == is || std::abs(Ue[k](js)) < 1.e-5);
		          continue;
		        }
		      if (std::abs(their_dof_value) < 1.e-5)
		        continue;

		      // since we may be running this method concurretly 
		      // on multiple threads we need to acquire a lock 
		      // before modifying the shared constraint_row object.
		      {
			Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);

			DofConstraintRow& constraint_row =
			  constraints[my_dof_g];

			constraint_row.insert(std::make_pair(their_dof_g,
							     their_dof_value));
		      }
		    }
	        }
	    }
          // p refinement constraints:
          // constrain dofs shared between
          // active elements and neighbors with
          // lower polynomial degrees
#ifdef LIBMESH_ENABLE_AMR
          const unsigned int min_p_level =
            neigh->min_p_level_by_neighbor(elem, elem->p_level());
          if (min_p_level < elem->p_level())
            {
              // Adaptive p refinement of non-hierarchic bases will
              // require more coding
              libmesh_assert(my_fe->is_hierarchic());
              dof_map.constrain_p_dofs(variable_number, elem,
                                       s, min_p_level);
            }
#endif // #ifdef LIBMESH_ENABLE_AMR
        }
    }
}

#endif // LIBMESH_ENABLE_PERIODIC


} // namespace libMesh
