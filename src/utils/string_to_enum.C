// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include <algorithm>
#include <map>

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/enum_elem_quality.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/enum_eigen_solver_type.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_inf_map_type.h"
#include "libmesh/enum_io_package.h"
#include "libmesh/enum_norm_type.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_parallel_type.h"
#include "libmesh/enum_point_locator_type.h"
#include "libmesh/enum_preconditioner_type.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/enum_solver_package.h"
#include "libmesh/enum_solver_type.h"
#include "libmesh/enum_subset_solve_mode.h"
#include "libmesh/enum_xdr_mode.h"
#include "libmesh/elem.h"

namespace libMesh
{


// ------------------------------------------------------------
// Anonymous namespace to hold local data & methods
namespace {


  // Reverse a map
  template <typename MapIter, class MapType>
  inline
  void build_reverse_map (MapIter it, MapIter end, MapType& reverse)
  {
    reverse.clear();

    for (; it != end; ++it)
      {
	// If the forward map is not invertible, we might already have
	// found a preimage of it->second.  Choose the "largest"
	// preimage according to operator<; for std::string this will
	// give us the longest, hopefully most specific name
	// corresponding to an enum.
	typename MapType::iterator preimage = reverse.find(it->second);
        if (preimage == reverse.end())
          reverse.insert (std::make_pair(it->second, it->first));
        else if (preimage->second < it->first)
          preimage->second = it->first;
      }
  }

#define INSTANTIATE_ENUM_MAPS(ENUM_NAME,VAR_NAME) \
  std::map<std::string, ENUM_NAME> VAR_NAME##_to_enum; \
 \
  std::map<ENUM_NAME, std::string> enum_to_##VAR_NAME; \
 \
  void init_##VAR_NAME##_to_enum (); \
 \
  /* Initialize the enum_to_elem_type on first call */ \
  void init_enum_to_##VAR_NAME () \
  { \
    /* Build reverse map */ \
    if (enum_to_##VAR_NAME .empty()) \
      { \
	/* Initialize elem_type_to_enum on first call */ \
	init_##VAR_NAME##_to_enum(); \
 \
	build_reverse_map (VAR_NAME##_to_enum.begin(), \
			   VAR_NAME##_to_enum.end(), \
			   enum_to_##VAR_NAME); \
      } \
  }

INSTANTIATE_ENUM_MAPS(ElemType, elem_type)

  //----------------------------------------------------

  // Initialize elem_type_to_enum on first call
  void init_elem_type_to_enum ()
  {
    if (elem_type_to_enum.empty())
      {
	elem_type_to_enum["EDGE"      ]=EDGE2;
	elem_type_to_enum["EDGE2"     ]=EDGE2;
	elem_type_to_enum["EDGE3"     ]=EDGE3;
	elem_type_to_enum["EDGE4"     ]=EDGE4;

	elem_type_to_enum["TRI"       ]=TRI3;
	elem_type_to_enum["TRI3"      ]=TRI3;
	elem_type_to_enum["TRI6"      ]=TRI6;

	elem_type_to_enum["QUAD"      ]=QUAD4;
	elem_type_to_enum["QUAD4"     ]=QUAD4;
	elem_type_to_enum["QUAD8"     ]=QUAD8;
	elem_type_to_enum["QUAD9"     ]=QUAD9;

	elem_type_to_enum["TET"       ]=TET4;
	elem_type_to_enum["TET4"      ]=TET4;
	elem_type_to_enum["TET10"     ]=TET10;

	elem_type_to_enum["HEX"       ]=HEX8;
	elem_type_to_enum["HEX8"      ]=HEX8;
	elem_type_to_enum["HEX20"     ]=HEX20;
	elem_type_to_enum["HEX27"     ]=HEX27;

	elem_type_to_enum["PRISM"     ]=PRISM6;
	elem_type_to_enum["PRISM6"    ]=PRISM6;
	elem_type_to_enum["PRISM15"   ]=PRISM15;
	elem_type_to_enum["PRISM18"   ]=PRISM18;

	elem_type_to_enum["PYRAMID"   ]=PYRAMID5;
	elem_type_to_enum["PYRAMID5"  ]=PYRAMID5;

	elem_type_to_enum["INFEDGE"   ]=INFEDGE2;
	elem_type_to_enum["INFEDGE2"  ]=INFEDGE2;

	elem_type_to_enum["INFQUAD"   ]=INFQUAD4;
	elem_type_to_enum["INFQUAD4"  ]=INFQUAD4;
	elem_type_to_enum["INFQUAD6"  ]=INFQUAD6;

	elem_type_to_enum["INFHEX"    ]=INFHEX8;
	elem_type_to_enum["INFHEX8"   ]=INFHEX8;
	elem_type_to_enum["INFHEX16"  ]=INFHEX16;
	elem_type_to_enum["INFHEX18"  ]=INFHEX18;

	elem_type_to_enum["INFPRISM"  ]=INFPRISM6;
	elem_type_to_enum["INFPRISM6" ]=INFPRISM6;
	elem_type_to_enum["INFPRISM12"]=INFPRISM12;

	elem_type_to_enum["NODE"      ]=NODEELEM;
	elem_type_to_enum["NODEELEM"  ]=NODEELEM;
      }
  }


INSTANTIATE_ENUM_MAPS(Order, order)

  // Initialize order_to_enum on first call
  void init_order_to_enum ()
  {
    if (order_to_enum.empty())
      {
	order_to_enum["CONSTANT"     ]=CONSTANT;
	order_to_enum["FIRST"        ]=FIRST;
	order_to_enum["SECOND"       ]=SECOND;
	order_to_enum["THIRD"        ]=THIRD;
	order_to_enum["FOURTH"       ]=FOURTH;
	order_to_enum["FIFTH"        ]=FIFTH;
	order_to_enum["SIXTH"        ]=SIXTH;
	order_to_enum["SEVENTH"      ]=SEVENTH;
	order_to_enum["EIGHTH"       ]=EIGHTH;
	order_to_enum["NINTH"        ]=NINTH;
	order_to_enum["TENTH"        ]=TENTH;

	order_to_enum["ELEVENTH"     ]=ELEVENTH;
	order_to_enum["TWELFTH"      ]=TWELFTH;
	order_to_enum["THIRTEENTH"   ]=THIRTEENTH;
	order_to_enum["FOURTEENTH"   ]=FOURTEENTH;
	order_to_enum["FIFTEENTH"    ]=FIFTEENTH;
	order_to_enum["SIXTEENTH"    ]=SIXTEENTH;
	order_to_enum["SEVENTEENTH"  ]=SEVENTEENTH;
	order_to_enum["EIGHTTEENTH"  ]=EIGHTTEENTH;
	order_to_enum["NINTEENTH"    ]=NINTEENTH;
	order_to_enum["TWENTIETH"    ]=TWENTIETH;

	order_to_enum["TWENTYFIRST"  ]=TWENTYFIRST;
	order_to_enum["TWENTYSECOND" ]=TWENTYSECOND;
	order_to_enum["TWENTYTHIRD"  ]=TWENTYTHIRD;
	order_to_enum["TWENTYFOURTH" ]=TWENTYFOURTH;
	order_to_enum["TWENTYFIFTH"  ]=TWENTYFIFTH;
	order_to_enum["TWENTYSIXTH"  ]=TWENTYSIXTH;
	order_to_enum["TWENTYSEVENTH"]=TWENTYSEVENTH;
	order_to_enum["TWENTYEIGHTH" ]=TWENTYEIGHTH;
	order_to_enum["TWENTYNINTH"  ]=TWENTYNINTH;
	order_to_enum["THIRTIETH"    ]=THIRTIETH;

	order_to_enum["THIRTYFIRST"  ]=THIRTYFIRST;
	order_to_enum["THIRTYSECOND" ]=THIRTYSECOND;
	order_to_enum["THIRTYTHIRD"  ]=THIRTYTHIRD;
	order_to_enum["THIRTYFOURTH" ]=THIRTYFOURTH;
	order_to_enum["THIRTYFIFTH"  ]=THIRTYFIFTH;
	order_to_enum["THIRTYSIXTH"  ]=THIRTYSIXTH;
	order_to_enum["THIRTYSEVENTH"]=THIRTYSEVENTH;
	order_to_enum["THIRTYEIGHTH" ]=THIRTYEIGHTH;
	order_to_enum["THIRTYNINTH"  ]=THIRTYNINTH;
	order_to_enum["FORTIETH"    ]=FORTIETH;

	order_to_enum["FORTYFIRST"  ]=FORTYFIRST;
	order_to_enum["FORTYSECOND" ]=FORTYSECOND;
	order_to_enum["FORTYTHIRD"  ]=FORTYTHIRD;
      }
  }



INSTANTIATE_ENUM_MAPS(FEFamily, fefamily)

  // Initialize fefamily_to_enum on first call
  void init_fefamily_to_enum ()
  {
    if (fefamily_to_enum.empty())
      {
	fefamily_to_enum["LAGRANGE"     ]=LAGRANGE;
	fefamily_to_enum["LAGRANGE_VEC" ]=LAGRANGE_VEC;
	fefamily_to_enum["L2_LAGRANGE"  ]=L2_LAGRANGE;
	fefamily_to_enum["HIERARCHIC"   ]=HIERARCHIC;
	fefamily_to_enum["L2_HIERARCHIC"]=L2_HIERARCHIC;
	fefamily_to_enum["MONOMIAL"     ]=MONOMIAL;
	fefamily_to_enum["SCALAR"       ]=SCALAR;
	fefamily_to_enum["XYZ"          ]=XYZ;
	fefamily_to_enum["BERNSTEIN"    ]=BERNSTEIN;
	fefamily_to_enum["SZABAB"       ]=SZABAB;
	fefamily_to_enum["INFINITE_MAP" ]=INFINITE_MAP;
	fefamily_to_enum["JACOBI_20_00" ]=JACOBI_20_00;
	fefamily_to_enum["JACOBI_30_00" ]=JACOBI_30_00;
	fefamily_to_enum["LEGENDRE"     ]=LEGENDRE;
	fefamily_to_enum["CLOUGH"       ]=CLOUGH;
	fefamily_to_enum["HERMITE"      ]=HERMITE;
	fefamily_to_enum["NEDELEC_ONE"  ]=NEDELEC_ONE;
      }

  }



INSTANTIATE_ENUM_MAPS(InfMapType, inf_map_type)

  // Initialize inf_map_type_to_enum on first call
  void init_inf_map_type_to_enum ()
  {
    if (inf_map_type_to_enum.empty())
      {
	inf_map_type_to_enum["CARTESIAN"  ]=CARTESIAN;
	inf_map_type_to_enum["SPHERICAL"  ]=SPHERICAL;
	inf_map_type_to_enum["ELLIPSOIDAL"]=ELLIPSOIDAL;
      }
  }


INSTANTIATE_ENUM_MAPS(QuadratureType, quadrature_type)

  // Initialize quadrature_type_to_enum on first call
  void init_quadrature_type_to_enum ()
  {
    if (quadrature_type_to_enum.empty())
      {
	quadrature_type_to_enum["QGAUSS"     ]=QGAUSS;
	quadrature_type_to_enum["QJACOBI_1_0"]=QJACOBI_1_0;
	quadrature_type_to_enum["QJACOBI_2_0"]=QJACOBI_2_0;
	quadrature_type_to_enum["QSIMPSON"   ]=QSIMPSON;
	quadrature_type_to_enum["QTRAP"      ]=QTRAP;
	quadrature_type_to_enum["QGRID"      ]=QGRID;
	quadrature_type_to_enum["QCLOUGH"    ]=QCLOUGH;
      }
  }


INSTANTIATE_ENUM_MAPS(PreconditionerType, preconditioner_type)

  // Initialize preconditioner_type_to_enum on first call
  void init_preconditioner_type_to_enum ()
  {
    if (preconditioner_type_to_enum.empty())
      {
	preconditioner_type_to_enum["IDENTITY_PRECOND"      ]=IDENTITY_PRECOND;
	preconditioner_type_to_enum["JACOBI_PRECOND"	    ]=JACOBI_PRECOND;
	preconditioner_type_to_enum["BLOCK_JACOBI_PRECOND"  ]=BLOCK_JACOBI_PRECOND;
	preconditioner_type_to_enum["SOR_PRECOND"           ]=SOR_PRECOND;
	preconditioner_type_to_enum["SSOR_PRECOND"          ]=SSOR_PRECOND;
	preconditioner_type_to_enum["EISENSTAT_PRECOND"	    ]=EISENSTAT_PRECOND;
	preconditioner_type_to_enum["ASM_PRECOND"	    ]=ASM_PRECOND;
	preconditioner_type_to_enum["CHOLESKY_PRECOND"	    ]=CHOLESKY_PRECOND;
	preconditioner_type_to_enum["ICC_PRECOND"	    ]=ICC_PRECOND;
	preconditioner_type_to_enum["ILU_PRECOND"           ]=ILU_PRECOND;
	preconditioner_type_to_enum["LU_PRECOND"            ]=LU_PRECOND;
	preconditioner_type_to_enum["USER_PRECOND"          ]=USER_PRECOND;
	preconditioner_type_to_enum["SHELL_PRECOND"         ]=SHELL_PRECOND;
	preconditioner_type_to_enum["AMG_PRECOND"           ]=AMG_PRECOND;
	preconditioner_type_to_enum["INVALID_PRECONDITIONER"]=INVALID_PRECONDITIONER;

        //shorter
      	preconditioner_type_to_enum["IDENTITY"    ]=IDENTITY_PRECOND;
	preconditioner_type_to_enum["JACOBI"	  ]=JACOBI_PRECOND;
	preconditioner_type_to_enum["BLOCK_JACOBI"]=BLOCK_JACOBI_PRECOND;
	preconditioner_type_to_enum["SOR"         ]=SOR_PRECOND;
	preconditioner_type_to_enum["SSOR"        ]=SSOR_PRECOND;
	preconditioner_type_to_enum["EISENSTAT"	  ]=EISENSTAT_PRECOND;
	preconditioner_type_to_enum["ASM"	  ]=ASM_PRECOND;
	preconditioner_type_to_enum["CHOLESKY"	  ]=CHOLESKY_PRECOND;
	preconditioner_type_to_enum["ICC"	  ]=ICC_PRECOND;
	preconditioner_type_to_enum["ILU"         ]=ILU_PRECOND;
	preconditioner_type_to_enum["LU"          ]=LU_PRECOND;
	preconditioner_type_to_enum["USER"        ]=USER_PRECOND;
	preconditioner_type_to_enum["SHELL"       ]=SHELL_PRECOND;
	preconditioner_type_to_enum["AMG"         ]=AMG_PRECOND;
	preconditioner_type_to_enum["INVALID"     ]=INVALID_PRECONDITIONER;
      }
  }


#ifdef LIBMESH_ENABLE_AMR

INSTANTIATE_ENUM_MAPS(Elem::RefinementState, refinementstate_type)

  // Initialize refinementstate_type_to_enum on first call
  void init_refinementstate_type_to_enum ()
  {
    if (refinementstate_type_to_enum.empty())
      {
	refinementstate_type_to_enum["COARSEN"                ]=Elem::COARSEN;
	refinementstate_type_to_enum["DO_NOTHING"             ]=Elem::DO_NOTHING;
	refinementstate_type_to_enum["REFINE"                 ]=Elem::REFINE;
	refinementstate_type_to_enum["JUST_REFINED"           ]=Elem::JUST_REFINED;
	refinementstate_type_to_enum["JUST_COARSENED"         ]=Elem::JUST_COARSENED;
	refinementstate_type_to_enum["INACTIVE"               ]=Elem::INACTIVE;
	refinementstate_type_to_enum["COARSEN_INACTIVE"       ]=Elem::COARSEN_INACTIVE;
	refinementstate_type_to_enum["INVALID_REFINEMENTSTATE"]=Elem::INVALID_REFINEMENTSTATE;
      }
  }
#endif // LIBMESH_ENABLE_AMR


INSTANTIATE_ENUM_MAPS(EigenSolverType, eigensolvertype)

  // Initialize eigensolvertype_to_enum on first call
  void init_eigensolvertype_to_enum ()
  {
    if (eigensolvertype_to_enum.empty())
      {
	eigensolvertype_to_enum["POWER"              ]=POWER;
	eigensolvertype_to_enum["LAPACK"             ]=LAPACK;
	eigensolvertype_to_enum["SUBSPACE"           ]=SUBSPACE;
	eigensolvertype_to_enum["ARNOLDI"            ]=ARNOLDI;
	eigensolvertype_to_enum["LANCZOS"            ]=LANCZOS;
	eigensolvertype_to_enum["KRYLOVSCHUR"        ]=KRYLOVSCHUR;
	eigensolvertype_to_enum["INVALID_EIGENSOLVER"]=INVALID_EIGENSOLVER;
      }
  }


INSTANTIATE_ENUM_MAPS(SolverType, solvertype)

  // Initialize solvertype_to_enum on first call
  void init_solvertype_to_enum ()
  {
    if (solvertype_to_enum.empty())
      {
	solvertype_to_enum["CG"            ]=CG;
	solvertype_to_enum["CGN"           ]=CGN;
	solvertype_to_enum["CGS"           ]=CGS;
	solvertype_to_enum["CR"            ]=CR;
	solvertype_to_enum["QMR"           ]=QMR;
	solvertype_to_enum["TCQMR"         ]=TCQMR;
	solvertype_to_enum["TFQMR"         ]=TFQMR;
	solvertype_to_enum["BICG"          ]=BICG;
	solvertype_to_enum["MINRES"        ]=MINRES;
	solvertype_to_enum["GMRES"         ]=GMRES;
	solvertype_to_enum["LSQR"          ]=LSQR;
	solvertype_to_enum["JACOBI"        ]=JACOBI;
	solvertype_to_enum["SOR_FORWARD"   ]=SOR_FORWARD;
	solvertype_to_enum["SOR_BACKWARD"  ]=SOR_BACKWARD;
	solvertype_to_enum["SSOR"          ]=SSOR;
	solvertype_to_enum["RICHARDSON"    ]=RICHARDSON;
	solvertype_to_enum["CHEBYSHEV"     ]=CHEBYSHEV;
	solvertype_to_enum["INVALID_SOLVER"]=INVALID_SOLVER;
      }
  }


INSTANTIATE_ENUM_MAPS(ElemQuality, elemquality)

  // Initialize elemquality_to_enum on first call
  void init_elemquality_to_enum ()
  {
    if (elemquality_to_enum.empty())
      {
	    elemquality_to_enum["ASPECT_RATIO"       ]=ASPECT_RATIO;
	    elemquality_to_enum["SKEW"               ]=SKEW;
	    elemquality_to_enum["SHEAR"              ]=SHEAR;
	    elemquality_to_enum["SHAPE"              ]=SHAPE;
	    elemquality_to_enum["MAX_ANGLE"          ]=MAX_ANGLE;
	    elemquality_to_enum["MIN_ANGLE"          ]=MIN_ANGLE;
	    elemquality_to_enum["CONDITION"          ]=CONDITION;
	    elemquality_to_enum["DISTORTION"         ]=DISTORTION;
	    elemquality_to_enum["TAPER"              ]=TAPER;
	    elemquality_to_enum["WARP"               ]=WARP;
	    elemquality_to_enum["STRETCH"            ]=STRETCH;
	    elemquality_to_enum["DIAGONAL"           ]=DIAGONAL;
	    elemquality_to_enum["ASPECT_RATIO_BETA"  ]=ASPECT_RATIO_BETA;
	    elemquality_to_enum["ASPECT_RATIO_GAMMA" ]=ASPECT_RATIO_GAMMA;
	    elemquality_to_enum["SIZE"               ]=SIZE;
	    elemquality_to_enum["JACOBIAN"           ]=JACOBIAN;
      }
  }


INSTANTIATE_ENUM_MAPS(IOPackage, iopackage)

  // Initialize iopackage_to_enum on first call
  void init_iopackage_to_enum ()
  {
    if (iopackage_to_enum.empty())
      {
	    iopackage_to_enum["TECPLOT" ]=TECPLOT;
	    iopackage_to_enum["GMV"     ]=GMV;
	    iopackage_to_enum["GMSH"    ]=GMSH;
	    iopackage_to_enum["VTK"     ]=VTK;
	    iopackage_to_enum["DIVA"    ]=DIVA;
	    iopackage_to_enum["TETGEN"  ]=TETGEN;
	    iopackage_to_enum["UCD"     ]=UCD;
	    iopackage_to_enum["LIBMESH" ]=LIBMESH;
      }
  }


INSTANTIATE_ENUM_MAPS(FEMNormType, norm_type)

  // Initialize norm_type_to_enum on first call
  void init_norm_type_to_enum ()
  {
    if (norm_type_to_enum.empty())
      {
	    norm_type_to_enum["L2" ]=L2;
	    norm_type_to_enum["H1" ]=H1;
	    norm_type_to_enum["H2" ]=H2;
	    norm_type_to_enum["HCURL" ]=HCURL;
	    norm_type_to_enum["HDIV" ]=HDIV;

	    norm_type_to_enum["L1" ]=L1;
	    norm_type_to_enum["L_INF" ]=L_INF;

	    norm_type_to_enum["H1_SEMINORM" ]=H1_SEMINORM;
	    norm_type_to_enum["H2_SEMINORM" ]=H2_SEMINORM;
	    norm_type_to_enum["HCURL_SEMINORM" ]=HCURL_SEMINORM;
	    norm_type_to_enum["HDIV_SEMINORM" ]=HDIV_SEMINORM;

	    norm_type_to_enum["W1_INF_SEMINORM" ]=W1_INF_SEMINORM;
	    norm_type_to_enum["W2_INF_SEMINORM" ]=W2_INF_SEMINORM;

	    norm_type_to_enum["DISCRETE_L1" ]=DISCRETE_L1;
	    norm_type_to_enum["DISCRETE_L2" ]=DISCRETE_L2;
	    norm_type_to_enum["DISCRETE_L_INF" ]=DISCRETE_L_INF;

	    norm_type_to_enum["H1_X_SEMINORM" ]=H1_X_SEMINORM;
	    norm_type_to_enum["H1_Y_SEMINORM" ]=H1_Y_SEMINORM;
	    norm_type_to_enum["H1_Z_SEMINORM" ]=H1_Z_SEMINORM;

	    norm_type_to_enum["INVALID_NORM" ]=INVALID_NORM;
      }
  }


INSTANTIATE_ENUM_MAPS(ParallelType, parallel_type)

  // Initialize parallel_type_to_enum on first call
  void init_parallel_type_to_enum ()
  {
    if (parallel_type_to_enum.empty())
      {
	parallel_type_to_enum["AUTOMATIC" ]=AUTOMATIC;
	parallel_type_to_enum["SERIAL"    ]=SERIAL;
	parallel_type_to_enum["PARALLEL"  ]=PARALLEL;
	parallel_type_to_enum["GHOSTED"   ]=GHOSTED;
	parallel_type_to_enum["INVALID_PARALLELIZATION" ]=INVALID_PARALLELIZATION;
      }
  }


INSTANTIATE_ENUM_MAPS(PointLocatorType, point_locator_type)

  // Initialize point_locator_type_to_enum on first call
  void init_point_locator_type_to_enum ()
  {
    if (point_locator_type_to_enum.empty())
      {
	point_locator_type_to_enum["TREE" ]=TREE;
	point_locator_type_to_enum["LIST" ]=LIST;
	point_locator_type_to_enum["INVALID_LOCATOR" ]=INVALID_LOCATOR;
      }
  }


INSTANTIATE_ENUM_MAPS(SolverPackage, solverpackage_type)

  // Initialize solverpackage_type_to_enum on first call
  void init_solverpackage_type_to_enum ()
  {
    if (solverpackage_type_to_enum.empty())
      {
	solverpackage_type_to_enum["PETSC_SOLVERS"    ]=PETSC_SOLVERS;
	solverpackage_type_to_enum["TRILINOS_SOLVERS" ]=TRILINOS_SOLVERS;
	solverpackage_type_to_enum["LASPACK_SOLVERS"  ]=LASPACK_SOLVERS;
	solverpackage_type_to_enum["SLEPC_SOLVERS"    ]=SLEPC_SOLVERS;
	solverpackage_type_to_enum["EIGEN_SOLVERS"    ]=EIGEN_SOLVERS;
	solverpackage_type_to_enum["INVALID_SOLVER_PACKAGE" ]=INVALID_SOLVER_PACKAGE;
      }
  }


INSTANTIATE_ENUM_MAPS(SubsetSolveMode, subset_solve_mode)

  // Initialize subset_solve_mode_to_enum on first call
  void init_subset_solve_mode_to_enum ()
  {
    if (subset_solve_mode_to_enum.empty())
      {
	subset_solve_mode_to_enum["SUBSET_ZERO" ]=SUBSET_ZERO;
	subset_solve_mode_to_enum["SUBSET_COPY_RHS" ]=SUBSET_COPY_RHS;
	subset_solve_mode_to_enum["SUBSET_DONT_TOUCH" ]=SUBSET_DONT_TOUCH;
      }
  }


INSTANTIATE_ENUM_MAPS(XdrMODE, xdr_mode)

  // Initialize xdr_mode_to_enum on first call
  void init_xdr_mode_to_enum ()
  {
    if (xdr_mode_to_enum.empty())
      {
	xdr_mode_to_enum["UNKNOWN" ]=UNKNOWN;
	xdr_mode_to_enum["ENCODE"  ]=ENCODE;
	xdr_mode_to_enum["DECODE"  ]=DECODE;
	xdr_mode_to_enum["WRITE"   ]=WRITE;
	xdr_mode_to_enum["READ"    ]=READ;
      }
  }





#undef INSTANTIATE_ENUM_MAPS

} // end anonymous namespace



// ------------------------------------------------------
// Utility::string_to_enum<> & Utility::enum_to_string<>
// full specializations
namespace Utility {

#define INSTANTIATE_STRING_TO_ENUM(ENUM_NAME,VAR_NAME) \
  template <> \
  ENUM_NAME string_to_enum<ENUM_NAME> (const std::string& s) \
  { \
    init_##VAR_NAME##_to_enum(); \
 \
    std::string upper(s); \
    std::transform(upper.begin(), upper.end(), upper.begin(), ::toupper); \
 \
    if (!VAR_NAME##_to_enum.count(upper)) \
      { \
        libmesh_error_msg("No " #ENUM_NAME " named " + s + " found."); \
      } \
 \
    return VAR_NAME##_to_enum[upper]; \
  } \
 \
  template <> \
  std::string enum_to_string<ENUM_NAME> (const ENUM_NAME e) \
  { \
    init_enum_to_##VAR_NAME (); \
 \
    if (!enum_to_##VAR_NAME .count(e)) \
      libmesh_error(); \
 \
    return enum_to_##VAR_NAME [e]; \
  }


INSTANTIATE_STRING_TO_ENUM(ElemType,elem_type)
INSTANTIATE_STRING_TO_ENUM(Order,order)
INSTANTIATE_STRING_TO_ENUM(FEFamily,fefamily)
INSTANTIATE_STRING_TO_ENUM(InfMapType,inf_map_type)
INSTANTIATE_STRING_TO_ENUM(QuadratureType,quadrature_type)
INSTANTIATE_STRING_TO_ENUM(PreconditionerType,preconditioner_type)

#ifdef LIBMESH_ENABLE_AMR
INSTANTIATE_STRING_TO_ENUM(Elem::RefinementState,refinementstate_type)
#endif // LIBMESH_ENABLE_AMR

INSTANTIATE_STRING_TO_ENUM(SolverType,solvertype)
INSTANTIATE_STRING_TO_ENUM(EigenSolverType,eigensolvertype)
INSTANTIATE_STRING_TO_ENUM(ElemQuality,elemquality)
INSTANTIATE_STRING_TO_ENUM(IOPackage,iopackage)
INSTANTIATE_STRING_TO_ENUM(FEMNormType, norm_type)
INSTANTIATE_STRING_TO_ENUM(ParallelType, parallel_type)
INSTANTIATE_STRING_TO_ENUM(PointLocatorType, point_locator_type)
INSTANTIATE_STRING_TO_ENUM(SolverPackage,solverpackage_type)
INSTANTIATE_STRING_TO_ENUM(SubsetSolveMode,subset_solve_mode)
INSTANTIATE_STRING_TO_ENUM(XdrMODE,xdr_mode)

#undef INSTANTIATE_STRING_TO_ENUM

} // namespace Utility

} // namespace libMesh
