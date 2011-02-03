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



// C++ includes
#include <algorithm>
#include <map>

// Local includes
#include "libmesh_common.h"
#include "string_to_enum.h"
#include "enum_elem_type.h"
#include "enum_order.h"
#include "enum_fe_family.h"
#include "enum_inf_map_type.h"
#include "enum_quadrature_type.h"
#include "enum_preconditioner_type.h"

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


  //----------------------------------------------------  
  std::map<std::string, ElemType> elem_type_to_enum;
  
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
      }
  }


  
  std::map<ElemType, std::string> enum_to_elem_type;

  // Initialize the enum_to_elem_type on first call
  void init_enum_to_elem_type ()
  {
    // Build reverse map
    if (enum_to_elem_type.empty())
      {
	// Initialize elem_type_to_enum on first call
	init_elem_type_to_enum();

	build_reverse_map (elem_type_to_enum.begin(),
			   elem_type_to_enum.end(),
			   enum_to_elem_type);
      }
  }



  
  //---------------------------------------------
  std::map<std::string, Order> order_to_enum;
  
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



  std::map<Order, std::string> enum_to_order;
  
  // Initialize the enum_to_order on first call
  void init_enum_to_order ()
  {
    // Build reverse map
    if (enum_to_order.empty())
      {
	// Initialize order_to_enum on first call
	init_order_to_enum();

	build_reverse_map (order_to_enum.begin(),
			   order_to_enum.end(),
			   enum_to_order);
      }
  }



  //---------------------------------------------------  
  std::map<std::string, FEFamily> fefamily_to_enum;

  // Initialize fefamily_to_enum on first call
  void init_fefamily_to_enum ()
  {
    if (fefamily_to_enum.empty())
      {
	fefamily_to_enum["LAGRANGE"    ]=LAGRANGE;
	fefamily_to_enum["HIERARCHIC"  ]=HIERARCHIC;
	fefamily_to_enum["MONOMIAL"    ]=MONOMIAL;
	fefamily_to_enum["SCALAR"      ]=SCALAR;
	fefamily_to_enum["XYZ"         ]=XYZ;
	fefamily_to_enum["BERNSTEIN"   ]=BERNSTEIN;
	fefamily_to_enum["SZABAB"      ]=SZABAB;
	fefamily_to_enum["INFINITE_MAP"]=INFINITE_MAP;
	fefamily_to_enum["JACOBI_20_00"]=JACOBI_20_00;
	fefamily_to_enum["JACOBI_30_00"]=JACOBI_30_00;
	fefamily_to_enum["LEGENDRE"    ]=LEGENDRE;
	fefamily_to_enum["CLOUGH"      ]=CLOUGH;
	fefamily_to_enum["HERMITE"     ]=HERMITE;
      }
    
  }


  std::map<FEFamily, std::string> enum_to_fefamily;
  
  // Initialize the enum_to_fefamily on first call
  void init_enum_to_fefamily ()
  {
    // Build reverse map
    if (enum_to_fefamily.empty())
      {
	// Initialize fefamily_to_enum on first call
	init_fefamily_to_enum();

	build_reverse_map (fefamily_to_enum.begin(),
			   fefamily_to_enum.end(),
			   enum_to_fefamily);
      }
  }



  //---------------------------------------------------  
  std::map<std::string, InfMapType> inf_map_type_to_enum;

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


  std::map<InfMapType, std::string> enum_to_inf_map_type;
  
  // Initialize the enum_to_inf_map_type on first call
  void init_enum_to_inf_map_type ()
  {
    // Build reverse map
    if (enum_to_inf_map_type.empty())
      {
	// Initialize inf_map_type_to_enum on first call
	init_inf_map_type_to_enum();

	build_reverse_map (inf_map_type_to_enum.begin(),
			   inf_map_type_to_enum.end(),
			   enum_to_inf_map_type);
      }
  }



  //---------------------------------------------------  
  std::map<std::string, QuadratureType> quadrature_type_to_enum;

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


  std::map<QuadratureType, std::string> enum_to_quadrature_type;
  
  // Initialize the enum_to_quadrature_type on first call
  void init_enum_to_quadrature_type ()
  {
    // Build reverse map
    if (enum_to_quadrature_type.empty())
      {
	// Initialize inf_map_type_to_enum on first call
	init_quadrature_type_to_enum();

	build_reverse_map (quadrature_type_to_enum.begin(),
			   quadrature_type_to_enum.end(),
			   enum_to_quadrature_type);
      }
  }


  //---------------------------------------------------  
  std::map<std::string, PreconditionerType> preconditioner_type_to_enum;

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


  std::map<PreconditionerType, std::string> enum_to_preconditioner_type;
  
  // Initialize the enum_to_preconditioner_type on first call
  void init_enum_to_preconditioner_type ()
  {
    // Build reverse map
    if (enum_to_preconditioner_type.empty())
      {
	// Initialize inf_map_type_to_enum on first call
	init_preconditioner_type_to_enum();

	build_reverse_map (preconditioner_type_to_enum.begin(),
			   preconditioner_type_to_enum.end(),
			   enum_to_preconditioner_type);
      }
  }
} // end anonymous namespace



// ------------------------------------------------------
// Utility::string_to_enum<> & Utility::enum_to_string<>
// full specializations
namespace Utility {

  //------------------------------------------------------
  // ElemType specialization
  template <>
  ElemType string_to_enum<ElemType> (const std::string& s)
  {
    init_elem_type_to_enum();

    std::string upper(s);
    std::transform(upper.begin(), upper.end(), upper.begin(), ::toupper);
    
    if (!elem_type_to_enum.count(upper))
      libmesh_error();
    
    return elem_type_to_enum[upper];
  }


  
  template <>
  std::string enum_to_string<ElemType> (const ElemType e)
  {
    init_enum_to_elem_type();

    if (!enum_to_elem_type.count(e))
      libmesh_error();

    return enum_to_elem_type[e];
  }


  
  //------------------------------------------------
  // Order specialization
  template <>
  Order string_to_enum<Order> (const std::string& s)
  {
    init_order_to_enum();
    
    std::string upper(s);
    std::transform(upper.begin(), upper.end(), upper.begin(), ::toupper);
    
    if (!order_to_enum.count(upper))
      libmesh_error();
    
    return order_to_enum[upper];
  }


  
  template <>
  std::string enum_to_string<Order> (const Order o)
  {
    init_enum_to_order();

    if (!enum_to_order.count(o))
      libmesh_error();

    return enum_to_order[o];
  }



  //------------------------------------------------------
  // FEFamily specialization
  template <>
  FEFamily string_to_enum<FEFamily> (const std::string& s)
  {
    init_fefamily_to_enum();
    
    std::string upper(s);
    std::transform(upper.begin(), upper.end(), upper.begin(), ::toupper);
    
    if (!fefamily_to_enum.count(upper))
      {
	libMesh::err << "ERROR: could not convert '" << upper << "' to enum." << std::endl;
	libmesh_error();
      }
    
    return fefamily_to_enum[upper];
  }


  
  template <>
  std::string enum_to_string<FEFamily> (const FEFamily f)
  {
    init_enum_to_fefamily();

    if (!enum_to_fefamily.count(f))
      libmesh_error();

    return enum_to_fefamily[f];
  }



  //------------------------------------------------------
  // InfMapType specialization
  template <>
  InfMapType string_to_enum<InfMapType> (const std::string& s)
  {
    init_inf_map_type_to_enum();
    
    std::string upper(s);
    std::transform(upper.begin(), upper.end(), upper.begin(), ::toupper);
    
    if (!inf_map_type_to_enum.count(upper))
      libmesh_error();
    
    return inf_map_type_to_enum[upper];
  }


  
  template <>
  std::string enum_to_string<InfMapType> (const InfMapType i)
  {
    init_enum_to_inf_map_type();

    if (!enum_to_inf_map_type.count(i))
      libmesh_error();

    return enum_to_inf_map_type[i];
  }



  //------------------------------------------------------
  // QuadratureType specialization
  template <>
  QuadratureType string_to_enum<QuadratureType> (const std::string& s)
  {
    init_quadrature_type_to_enum();
    
    std::string upper(s);
    std::transform(upper.begin(), upper.end(), upper.begin(), ::toupper);
    
    if (!quadrature_type_to_enum.count(upper))
      libmesh_error();
    
    return quadrature_type_to_enum[upper];
  }


  
  template <>
  std::string enum_to_string<QuadratureType> (const QuadratureType i)
  {
    init_enum_to_quadrature_type();

    if (!enum_to_quadrature_type.count(i))
      libmesh_error();

    return enum_to_quadrature_type[i];
  }


  //------------------------------------------------------
  // PreconditionerType specialization
  template <>
  PreconditionerType string_to_enum<PreconditionerType> (const std::string& s)
  {
    init_preconditioner_type_to_enum();
    
    std::string upper(s);
    std::transform(upper.begin(), upper.end(), upper.begin(), ::toupper);
    
    if (!preconditioner_type_to_enum.count(upper))
      libmesh_error();
    
    return preconditioner_type_to_enum[upper];
  }


  
  template <>
  std::string enum_to_string<PreconditionerType> (const PreconditionerType i)
  {
    init_enum_to_preconditioner_type();

    if (!enum_to_preconditioner_type.count(i))
      libmesh_error();

    return enum_to_preconditioner_type[i];
  }
  
}

} // namespace libMesh



