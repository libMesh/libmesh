// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/enum_to_string.h"
#include "libmesh/enum_convergence_flags.h"
#include "libmesh/enum_elem_quality.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/enum_eigen_solver_type.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_inf_map_type.h"
#include "libmesh/enum_io_package.h"
#include "libmesh/enum_norm_type.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_parallel_type.h"
#include "libmesh/enum_partitioner_type.h"
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
template <class MapType>
inline
std::map<typename MapType::mapped_type, typename MapType::key_type>
build_reverse_map (const MapType & forward)
{
  std::map<typename MapType::mapped_type, typename MapType::key_type> reverse;

  for (auto & [key, val] : forward)
    {
      // If the forward map is not invertible, we might already have
      // found a preimage of val.  Choose the "largest"
      // preimage according to operator<; for std::string this will
      // give us the longest, hopefully most specific name
      // corresponding to an enum.
      if (auto preimage = reverse.find(val);
          preimage == reverse.end())
        reverse.emplace (val, key);
      else if (preimage->second < key)
        preimage->second = key;
    }

  return reverse;
}

std::map<std::string, ElemType> elem_type_to_enum {
   {"EDGE"           , EDGE2},
   {"EDGE2"          , EDGE2},
   {"EDGE2"          , EDGE2},
   {"EDGE3"          , EDGE3},
   {"EDGE4"          , EDGE4},

   {"TRI"            , TRI3},
   {"TRI3"           , TRI3},
   {"TRISHELL3"      , TRISHELL3},
   {"TRI3SUBDIVISION", TRI3SUBDIVISION},
   {"TRI6"           , TRI6},
   {"TRI7"           , TRI7},

   {"QUAD"           , QUAD4},
   {"QUAD4"          , QUAD4},
   {"QUADSHELL4"     , QUADSHELL4},
   {"QUAD8"          , QUAD8},
   {"QUADSHELL8"     , QUADSHELL8},
   {"QUAD9"          , QUAD9},
   {"QUADSHELL9"     , QUADSHELL9},

   {"TET"            , TET4},
   {"TET4"           , TET4},
   {"TET10"          , TET10},
   {"TET14"          , TET14},

   {"HEX"            , HEX8},
   {"HEX8"           , HEX8},
   {"HEX20"          , HEX20},
   {"HEX27"          , HEX27},

   {"PRISM"          , PRISM6},
   {"PRISM6"         , PRISM6},
   {"PRISM15"        , PRISM15},
   {"PRISM18"        , PRISM18},
   {"PRISM20"        , PRISM20},
   {"PRISM21"        , PRISM21},

   {"PYRAMID"        , PYRAMID5},
   {"PYRAMID5"       , PYRAMID5},
   {"PYRAMID13"      , PYRAMID13},
   {"PYRAMID14"      , PYRAMID14},
   {"PYRAMID18"      , PYRAMID18},

   {"INFEDGE"        , INFEDGE2},
   {"INFEDGE2"       , INFEDGE2},

   {"INFQUAD"        , INFQUAD4},
   {"INFQUAD4"       , INFQUAD4},
   {"INFQUAD6"       , INFQUAD6},

   {"INFHEX"         , INFHEX8},
   {"INFHEX8"        , INFHEX8},
   {"INFHEX16"       , INFHEX16},
   {"INFHEX18"       , INFHEX18},

   {"INFPRISM"       , INFPRISM6},
   {"INFPRISM6"      , INFPRISM6},
   {"INFPRISM12"     , INFPRISM12},

   {"NODE"           , NODEELEM},
   {"NODEELEM"       , NODEELEM},

   {"INVALID_ELEM"   , INVALID_ELEM}
  };

std::map<ElemType, std::string> enum_to_elem_type =
  build_reverse_map(elem_type_to_enum);


std::map<std::string, ElemMappingType> elem_mapping_type_to_enum {
   {"LAGRANGE_MAP"          , LAGRANGE_MAP},
   {"RATIONAL_BERNSTEIN_MAP", RATIONAL_BERNSTEIN_MAP},
   {"INVALID_MAP"           , INVALID_MAP}
  };

std::map<ElemMappingType, std::string> enum_to_elem_mapping_type =
  build_reverse_map(elem_mapping_type_to_enum);


std::map<std::string, Order> order_to_enum {
   {"CONSTANT"     , CONSTANT},
   {"FIRST"        , FIRST},
   {"SECOND"       , SECOND},
   {"THIRD"        , THIRD},
   {"FOURTH"       , FOURTH},
   {"FIFTH"        , FIFTH},
   {"SIXTH"        , SIXTH},
   {"SEVENTH"      , SEVENTH},
   {"EIGHTH"       , EIGHTH},
   {"NINTH"        , NINTH},
   {"TENTH"        , TENTH},

   {"ELEVENTH"     , ELEVENTH},
   {"TWELFTH"      , TWELFTH},
   {"THIRTEENTH"   , THIRTEENTH},
   {"FOURTEENTH"   , FOURTEENTH},
   {"FIFTEENTH"    , FIFTEENTH},
   {"SIXTEENTH"    , SIXTEENTH},
   {"SEVENTEENTH"  , SEVENTEENTH},
   {"EIGHTTEENTH"  , EIGHTTEENTH},
   {"NINETEENTH"   , NINETEENTH},
   {"TWENTIETH"    , TWENTIETH},

   {"TWENTYFIRST"  , TWENTYFIRST},
   {"TWENTYSECOND" , TWENTYSECOND},
   {"TWENTYTHIRD"  , TWENTYTHIRD},
   {"TWENTYFOURTH" , TWENTYFOURTH},
   {"TWENTYFIFTH"  , TWENTYFIFTH},
   {"TWENTYSIXTH"  , TWENTYSIXTH},
   {"TWENTYSEVENTH", TWENTYSEVENTH},
   {"TWENTYEIGHTH" , TWENTYEIGHTH},
   {"TWENTYNINTH"  , TWENTYNINTH},
   {"THIRTIETH"    , THIRTIETH},

   {"THIRTYFIRST"  , THIRTYFIRST},
   {"THIRTYSECOND" , THIRTYSECOND},
   {"THIRTYTHIRD"  , THIRTYTHIRD},
   {"THIRTYFOURTH" , THIRTYFOURTH},
   {"THIRTYFIFTH"  , THIRTYFIFTH},
   {"THIRTYSIXTH"  , THIRTYSIXTH},
   {"THIRTYSEVENTH", THIRTYSEVENTH},
   {"THIRTYEIGHTH" , THIRTYEIGHTH},
   {"THIRTYNINTH"  , THIRTYNINTH},
   {"FORTIETH"     , FORTIETH},

   {"FORTYFIRST"   , FORTYFIRST},
   {"FORTYSECOND"  , FORTYSECOND},
   {"FORTYTHIRD"   , FORTYTHIRD},
   {"MAXIMUM"      , MAXIMUM}
  };

std::map<Order, std::string> enum_to_order =
  build_reverse_map(order_to_enum);


std::map<std::string, FEFamily> fefamily_to_enum {
   {"LAGRANGE"          , LAGRANGE},
   {"LAGRANGE_VEC"      , LAGRANGE_VEC},
   {"L2_LAGRANGE"       , L2_LAGRANGE},
   {"L2_LAGRANGE_VEC"   , L2_LAGRANGE_VEC},
   {"HIERARCHIC"        , HIERARCHIC},
   {"HIERARCHIC_VEC"    , HIERARCHIC_VEC},
   {"L2_HIERARCHIC"     , L2_HIERARCHIC},
   {"L2_HIERARCHIC_VEC" , L2_HIERARCHIC_VEC},
   {"SIDE_HIERARCHIC"   , SIDE_HIERARCHIC},
   {"MONOMIAL"          , MONOMIAL},
   {"MONOMIAL_VEC"      , MONOMIAL_VEC},
   {"SCALAR"            , SCALAR},
   {"XYZ"               , XYZ},
   {"BERNSTEIN"         , BERNSTEIN},
   {"RATIONAL_BERNSTEIN", RATIONAL_BERNSTEIN},
   {"SZABAB"            , SZABAB},
   {"INFINITE_MAP"      , INFINITE_MAP},
   {"JACOBI_20_00"      , JACOBI_20_00},
   {"JACOBI_30_00"      , JACOBI_30_00},
   {"LEGENDRE"          , LEGENDRE},
   {"CLOUGH"            , CLOUGH},
   {"HERMITE"           , HERMITE},
   {"SUBDIVISION"       , SUBDIVISION},
   {"NEDELEC_ONE"       , NEDELEC_ONE},
   {"RAVIART_THOMAS"    , RAVIART_THOMAS},
   {"L2_RAVIART_THOMAS" , L2_RAVIART_THOMAS}
  };

std::map<FEFamily, std::string> enum_to_fefamily =
  build_reverse_map(fefamily_to_enum);


std::map<std::string, InfMapType> inf_map_type_to_enum {
   {"CARTESIAN"  , CARTESIAN},
   {"SPHERICAL"  , SPHERICAL},
   {"ELLIPSOIDAL", ELLIPSOIDAL}
  };

std::map<InfMapType, std::string> enum_to_inf_map_type =
  build_reverse_map(inf_map_type_to_enum);


std::map<std::string, QuadratureType> quadrature_type_to_enum {
   {"QCLOUGH"          , QCLOUGH},
   {"QCOMPOSITE"       , QCOMPOSITE},
   {"QCONICAL"         , QCONICAL},
   {"QGAUSS"           , QGAUSS},
   {"QGAUSS_LOBATTO"   , QGAUSS_LOBATTO},
   {"QGRID"            , QGRID},
   {"QGRUNDMANN_MOLLER", QGRUNDMANN_MOLLER},
   {"QJACOBI_1_0"      , QJACOBI_1_0},
   {"QJACOBI_2_0"      , QJACOBI_2_0},
   {"QMONOMIAL"        , QMONOMIAL},
   {"QNODAL"           , QNODAL},
   {"QSIMPSON"         , QSIMPSON},
   {"QTRAP"            , QTRAP}
  };

std::map<QuadratureType, std::string> enum_to_quadrature_type =
  build_reverse_map(quadrature_type_to_enum);


std::map<std::string, PartitionerType> partitioner_type_to_enum {
   {"CENTROID_PARTITIONER"        , CENTROID_PARTITIONER},
   {"LINEAR_PARTITIONER"          , LINEAR_PARTITIONER},
   {"SFC_PARTITIONER"             , SFC_PARTITIONER},
   {"HILBERT_SFC_PARTITIONER"     , HILBERT_SFC_PARTITIONER},
   {"MORTON_SFC_PARTITIONER"      , MORTON_SFC_PARTITIONER},
   {"METIS_PARTITIONER"           , METIS_PARTITIONER},
   {"PARMETIS_PARTITIONER"        , PARMETIS_PARTITIONER},
   {"SUBDOMAIN_PARTITIONER"       , SUBDOMAIN_PARTITIONER},
   {"MAPPED_SUBDOMAIN_PARTITIONER", MAPPED_SUBDOMAIN_PARTITIONER},

      //shorter
   {"CENTROID"                    , CENTROID_PARTITIONER},
   {"LINEAR"                      , LINEAR_PARTITIONER},
   {"SFC"                         , SFC_PARTITIONER},
   {"HILBERT_SFC"                 , HILBERT_SFC_PARTITIONER},
   {"MORTON_SFC"                  , MORTON_SFC_PARTITIONER},
   {"METIS"                       , METIS_PARTITIONER},
   {"PARMETIS"                    , PARMETIS_PARTITIONER},
   {"SUBDOMAIN"                   , SUBDOMAIN_PARTITIONER},
   {"MAPPED_SUBDOMAIN"            , MAPPED_SUBDOMAIN_PARTITIONER},
  };

std::map<PartitionerType, std::string> enum_to_partitioner_type =
  build_reverse_map(partitioner_type_to_enum);


std::map<std::string, PreconditionerType> preconditioner_type_to_enum {
   {"IDENTITY_PRECOND"      , IDENTITY_PRECOND},
   {"JACOBI_PRECOND"        , JACOBI_PRECOND},
   {"BLOCK_JACOBI_PRECOND"  , BLOCK_JACOBI_PRECOND},
   {"SOR_PRECOND"           , SOR_PRECOND},
   {"SSOR_PRECOND"          , SSOR_PRECOND},
   {"EISENSTAT_PRECOND"     , EISENSTAT_PRECOND},
   {"ASM_PRECOND"           , ASM_PRECOND},
   {"CHOLESKY_PRECOND"      , CHOLESKY_PRECOND},
   {"ICC_PRECOND"           , ICC_PRECOND},
   {"ILU_PRECOND"           , ILU_PRECOND},
   {"LU_PRECOND"            , LU_PRECOND},
   {"USER_PRECOND"          , USER_PRECOND},
   {"SHELL_PRECOND"         , SHELL_PRECOND},
   {"AMG_PRECOND"           , AMG_PRECOND},
   {"SVD_PRECOND"           , SVD_PRECOND},
   {"INVALID_PRECONDITIONER", INVALID_PRECONDITIONER},

      //shorter
   {"IDENTITY"    , IDENTITY_PRECOND},
   {"JACOBI"      , JACOBI_PRECOND},
   {"BLOCK_JACOBI", BLOCK_JACOBI_PRECOND},
   {"SOR"         , SOR_PRECOND},
   {"SSOR"        , SSOR_PRECOND},
   {"EISENSTAT"   , EISENSTAT_PRECOND},
   {"ASM"         , ASM_PRECOND},
   {"CHOLESKY"    , CHOLESKY_PRECOND},
   {"ICC"         , ICC_PRECOND},
   {"ILU"         , ILU_PRECOND},
   {"LU"          , LU_PRECOND},
   {"USER"        , USER_PRECOND},
   {"SHELL"       , SHELL_PRECOND},
   {"AMG"         , AMG_PRECOND},
   {"SVD"         , SVD_PRECOND},
   {"INVALID"     , INVALID_PRECONDITIONER},
  };

std::map<PreconditionerType, std::string> enum_to_preconditioner_type =
  build_reverse_map(preconditioner_type_to_enum);


#ifdef LIBMESH_ENABLE_AMR
std::map<std::string, Elem::RefinementState> refinementstate_type_to_enum {
   {"COARSEN"                , Elem::COARSEN},
   {"DO_NOTHING"             , Elem::DO_NOTHING},
   {"REFINE"                 , Elem::REFINE},
   {"JUST_REFINED"           , Elem::JUST_REFINED},
   {"JUST_COARSENED"         , Elem::JUST_COARSENED},
   {"INACTIVE"               , Elem::INACTIVE},
   {"COARSEN_INACTIVE"       , Elem::COARSEN_INACTIVE},
   {"INVALID_REFINEMENTSTATE", Elem::INVALID_REFINEMENTSTATE},
  };

std::map<Elem::RefinementState, std::string> enum_to_refinementstate_type =
  build_reverse_map(refinementstate_type_to_enum);
#endif // LIBMESH_ENABLE_AMR


std::map<std::string, EigenSolverType> eigensolvertype_to_enum {
   {"POWER"              , POWER},
   {"LAPACK"             , LAPACK},
   {"SUBSPACE"           , SUBSPACE},
   {"ARNOLDI"            , ARNOLDI},
   {"LANCZOS"            , LANCZOS},
   {"KRYLOVSCHUR"        , KRYLOVSCHUR},
   {"INVALID_EIGENSOLVER", INVALID_EIGENSOLVER},
  };

std::map<EigenSolverType, std::string> enum_to_eigensolvertype =
  build_reverse_map(eigensolvertype_to_enum);


std::map<std::string, SolverType> solvertype_to_enum {
   {"CG"            , CG},
   {"CGN"           , CGN},
   {"CGS"           , CGS},
   {"CR"            , CR},
   {"QMR"           , QMR},
   {"TCQMR"         , TCQMR},
   {"TFQMR"         , TFQMR},
   {"BICG"          , BICG},
   {"BICGSTAB"      , BICGSTAB},
   {"MINRES"        , MINRES},
   {"GMRES"         , GMRES},
   {"LSQR"          , LSQR},
   {"JACOBI"        , JACOBI},
   {"SOR_FORWARD"   , SOR_FORWARD},
   {"SOR_BACKWARD"  , SOR_BACKWARD},
   {"SSOR"          , SSOR},
   {"RICHARDSON"    , RICHARDSON},
   {"CHEBYSHEV"     , CHEBYSHEV},
   {"SPARSELU"      , SPARSELU},
   {"INVALID_SOLVER", INVALID_SOLVER},
  };

std::map<SolverType, std::string> enum_to_solvertype =
  build_reverse_map(solvertype_to_enum);


std::map<std::string, ElemQuality> elemquality_to_enum {
   {"ASPECT_RATIO"       , ASPECT_RATIO},
   {"SKEW"               , SKEW},
   {"SHEAR"              , SHEAR},
   {"SHAPE"              , SHAPE},
   {"MAX_ANGLE"          , MAX_ANGLE},
   {"MIN_ANGLE"          , MIN_ANGLE},
   {"MAX_DIHEDRAL_ANGLE" , MAX_DIHEDRAL_ANGLE},
   {"MIN_DIHEDRAL_ANGLE" , MIN_DIHEDRAL_ANGLE},
   {"CONDITION"          , CONDITION},
   {"DISTORTION"         , DISTORTION},
   {"TAPER"              , TAPER},
   {"WARP"               , WARP},
   {"STRETCH"            , STRETCH},
   {"DIAGONAL"           , DIAGONAL},
   {"ASPECT_RATIO_BETA"  , ASPECT_RATIO_BETA},
   {"ASPECT_RATIO_GAMMA" , ASPECT_RATIO_GAMMA},
   {"SIZE"               , SIZE},
   {"JACOBIAN"           , JACOBIAN},
   {"TWIST"              , TWIST},
  };

std::map<ElemQuality, std::string> enum_to_elemquality =
  build_reverse_map(elemquality_to_enum);


std::map<std::string, IOPackage> iopackage_to_enum {
   {"TECPLOT" , TECPLOT},
   {"GMV"     , GMV},
   {"GMSH"    , GMSH},
   {"VTK"     , VTK},
   {"DIVA"    , DIVA},
   {"TETGEN"  , TETGEN},
   {"UCD"     , UCD},
   {"LIBMESH" , LIBMESH},
  };

std::map<IOPackage, std::string> enum_to_iopackage =
  build_reverse_map(iopackage_to_enum);


std::map<std::string, FEMNormType> norm_type_to_enum {
   {"L2"             , L2},
   {"H1"             , H1},
   {"H2"             , H2},
   {"HCURL"          , HCURL},
   {"HDIV"           , HDIV},

   {"L1"             , L1},
   {"L_INF"          , L_INF},

   {"H1_SEMINORM"    , H1_SEMINORM},
   {"H2_SEMINORM"    , H2_SEMINORM},
   {"HCURL_SEMINORM" , HCURL_SEMINORM},
   {"HDIV_SEMINORM"  , HDIV_SEMINORM},

   {"W1_INF_SEMINORM", W1_INF_SEMINORM},
   {"W2_INF_SEMINORM", W2_INF_SEMINORM},

   {"DISCRETE_L1"    , DISCRETE_L1},
   {"DISCRETE_L2"    , DISCRETE_L2},
   {"DISCRETE_L_INF" , DISCRETE_L_INF},

   {"H1_X_SEMINORM"  , H1_X_SEMINORM},
   {"H1_Y_SEMINORM"  , H1_Y_SEMINORM},
   {"H1_Z_SEMINORM"  , H1_Z_SEMINORM},

   {"INVALID_NORM"   , INVALID_NORM},
  };

std::map<FEMNormType, std::string> enum_to_norm_type =
  build_reverse_map(norm_type_to_enum);


std::map<std::string, ParallelType> parallel_type_to_enum {
   {"AUTOMATIC"               , AUTOMATIC},
   {"SERIAL"                  , SERIAL},
   {"PARALLEL"                , PARALLEL},
   {"GHOSTED"                 , GHOSTED},
   {"INVALID_PARALLELIZATION" , INVALID_PARALLELIZATION},
  };

std::map<ParallelType, std::string> enum_to_parallel_type =
  build_reverse_map(parallel_type_to_enum);


std::map<std::string, PointLocatorType> point_locator_type_to_enum {
   {"TREE"            , TREE},
   {"INVALID_LOCATOR" , INVALID_LOCATOR},
  };

std::map<PointLocatorType, std::string> enum_to_point_locator_type =
  build_reverse_map(point_locator_type_to_enum);


std::map<std::string, SolverPackage> solverpackage_type_to_enum {
   {"PETSC_SOLVERS"          , PETSC_SOLVERS},
   {"TRILINOS_SOLVERS"       , TRILINOS_SOLVERS},
   {"LASPACK_SOLVERS"        , LASPACK_SOLVERS},
   {"SLEPC_SOLVERS"          , SLEPC_SOLVERS},
   {"EIGEN_SOLVERS"          , EIGEN_SOLVERS},
   {"NLOPT_SOLVERS"          , NLOPT_SOLVERS},
   {"INVALID_SOLVER_PACKAGE" , INVALID_SOLVER_PACKAGE},
  };

std::map<SolverPackage, std::string> enum_to_solverpackage_type =
  build_reverse_map(solverpackage_type_to_enum);


std::map<std::string, SubsetSolveMode> subset_solve_mode_to_enum {
   {"SUBSET_ZERO"       , SUBSET_ZERO},
   {"SUBSET_COPY_RHS"   , SUBSET_COPY_RHS},
   {"SUBSET_DONT_TOUCH" , SUBSET_DONT_TOUCH},
  };

std::map<SubsetSolveMode, std::string> enum_to_subset_solve_mode =
  build_reverse_map(subset_solve_mode_to_enum);


std::map<std::string, XdrMODE> xdr_mode_to_enum {
   {"UNKNOWN" , UNKNOWN},
   {"ENCODE"  , ENCODE},
   {"DECODE"  , DECODE},
   {"WRITE"   , WRITE},
   {"READ"    , READ},
  };

std::map<XdrMODE, std::string> enum_to_xdr_mode =
  build_reverse_map(xdr_mode_to_enum);


std::map<std::string, LinearConvergenceReason> linear_convergence_reason_to_enum {
   {"CONVERGED_RTOL_NORMAL",       CONVERGED_RTOL_NORMAL},
   {"CONVERGED_ATOL_NORMAL",       CONVERGED_ATOL_NORMAL},
   {"CONVERGED_RTOL",              CONVERGED_RTOL},
   {"CONVERGED_ATOL",              CONVERGED_ATOL},
   {"CONVERGED_ITS",               CONVERGED_ITS},
   {"CONVERGED_CG_NEG_CURVE",      CONVERGED_CG_NEG_CURVE},
   {"CONVERGED_CG_CONSTRAINED",    CONVERGED_CG_CONSTRAINED},
   {"CONVERGED_STEP_LENGTH",       CONVERGED_STEP_LENGTH},
   {"CONVERGED_HAPPY_BREAKDOWN",   CONVERGED_HAPPY_BREAKDOWN},
   {"DIVERGED_NULL",               DIVERGED_NULL},
   {"DIVERGED_ITS",                DIVERGED_ITS},
   {"DIVERGED_DTOL",               DIVERGED_DTOL},
   {"DIVERGED_BREAKDOWN",          DIVERGED_BREAKDOWN},
   {"DIVERGED_BREAKDOWN_BICG",     DIVERGED_BREAKDOWN_BICG},
   {"DIVERGED_NONSYMMETRIC",       DIVERGED_NONSYMMETRIC},
   {"DIVERGED_INDEFINITE_PC",      DIVERGED_INDEFINITE_PC},
   {"DIVERGED_NAN",                DIVERGED_NAN},
   {"DIVERGED_INDEFINITE_MAT",     DIVERGED_INDEFINITE_MAT},
   {"DIVERGED_PCSETUP_FAILED",     DIVERGED_PCSETUP_FAILED},
   {"CONVERGED_ITERATING",         CONVERGED_ITERATING},
   {"UNKNOWN_FLAG",                UNKNOWN_FLAG},
  };

std::map<LinearConvergenceReason, std::string> enum_to_linear_convergence_reason =
  build_reverse_map(linear_convergence_reason_to_enum);

} // end anonymous namespace



// ------------------------------------------------------
// Utility::string_to_enum<> & Utility::enum_to_string<>
// full specializations
namespace Utility {

#define INSTANTIATE_STRING_TO_ENUM(ENUM_NAME,VAR_NAME)                  \
  template <>                                                           \
  ENUM_NAME string_to_enum<ENUM_NAME> (std::string_view s)              \
  {                                                                     \
    std::string upper(s);                                               \
    std::transform(upper.begin(), upper.end(), upper.begin(), ::toupper); \
                                                                        \
    if (!VAR_NAME##_to_enum.count(upper))                               \
      {                                                                 \
        libmesh_error_msg("No " #ENUM_NAME " named " << s << " found.");  \
      }                                                                 \
                                                                        \
    return VAR_NAME##_to_enum[upper];                                   \
  }                                                                     \
                                                                        \
  template <>                                                           \
  ENUM_NAME string_to_enum<ENUM_NAME> (const std::string & s)           \
  {                                                                     \
    return string_to_enum<ENUM_NAME>(std::string_view(s));              \
  }                                                                     \
                                                                        \
  template <>                                                           \
  ENUM_NAME string_to_enum<ENUM_NAME> (const char * s)                  \
  {                                                                     \
    return string_to_enum<ENUM_NAME>(std::string_view(s));              \
  }                                                                     \
                                                                        \
  template <>                                                           \
  std::string enum_to_string<ENUM_NAME> (const ENUM_NAME e)             \
  {                                                                     \
    if (!enum_to_##VAR_NAME .count(e))                                  \
      libmesh_error_msg("No " #ENUM_NAME " with enumeration " << e << " found."); \
                                                                        \
    return enum_to_##VAR_NAME [e];                                      \
  }



INSTANTIATE_STRING_TO_ENUM(ElemType,elem_type)
INSTANTIATE_STRING_TO_ENUM(ElemMappingType,elem_mapping_type)
INSTANTIATE_STRING_TO_ENUM(Order,order)
INSTANTIATE_STRING_TO_ENUM(FEFamily,fefamily)
INSTANTIATE_STRING_TO_ENUM(InfMapType,inf_map_type)
INSTANTIATE_STRING_TO_ENUM(QuadratureType,quadrature_type)
INSTANTIATE_STRING_TO_ENUM(PartitionerType,partitioner_type)
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
INSTANTIATE_STRING_TO_ENUM(LinearConvergenceReason, linear_convergence_reason)

#undef INSTANTIATE_STRING_TO_ENUM

} // namespace Utility

} // namespace libMesh
