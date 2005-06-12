// $Id: string_to_enum.C,v 1.1 2005-06-12 14:20:17 benkirk Exp $

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
#include <map>

// Local includes
#include "libmesh_common.h"
#include "string_to_enum.h"
#include "enum_elem_type.h"
#include "enum_order.h"
#include "enum_fe_family.h"



// ------------------------------------------------------------
// Utility::string_to_enum<> full specializations
namespace Utility {

  //------------------------------------------------------
  // ElemType specialization
  template <>
  ElemType string_to_enum<ElemType> (const std::string& s)
  {
    static std::map<std::string, ElemType> map;
    
    // Initialize map on first call
    if (map.empty())
      {
	map["EDGE2"     ]=EDGE2;
	map["EDGE3"     ]=EDGE3;
	map["EDGE4"     ]=EDGE4;
	
	map["TRI3"      ]=TRI3;
	map["TRI6"      ]=TRI6;
	
	map["QUAD4"     ]=QUAD4;
	map["QUAD8"     ]=QUAD8;
	map["QUAD9"     ]=QUAD9;
	
	map["TET4"      ]=TET4;
	map["TET10"     ]=TET10;
	
	map["HEX8"      ]=HEX8;
	map["HEX20"     ]=HEX20;
	map["HEX27"     ]=HEX27;
	
	map["PRISM6"    ]=PRISM6;
	map["PRISM15"   ]=PRISM15;
	map["PRISM18"   ]=PRISM18;
	
	map["PYRAMID5"  ]=PYRAMID5;
	
#ifdef ENABLE_INFINITE_ELEMENTS
	map["INFEDGE2"  ]=INFEDGE2;
	
	map["INFQUAD4"  ]=INFQUAD4;
	map["INFQUAD6"  ]=INFQUAD6;
	
	map["INFHEX8"   ]=INFHEX8;
	map["INFHEX16"  ]=INFHEX16;
	map["INFHEX18"  ]=INFHEX18;
	
	map["INFPRISM6" ]=INFPRISM6;
	map["INFPRISM12"]=INFPRISM12;
#endif
      }
    
    if (!map.count(s))
      error();
    
    return map[s];
  }


  
  //------------------------------------------------
  // Order specialization
  template <>
  Order string_to_enum<Order> (const std::string& s)
  {
    static std::map<std::string, Order> map;
    
    // Initialize map on first call
    if (map.empty())
      {
	map["CONSTANT"     ]=CONSTANT;
	map["FIRST"        ]=FIRST;
	map["SECOND"       ]=SECOND;
	map["THIRD"        ]=THIRD;
	map["FOURTH"       ]=FOURTH;
	map["FIFTH"        ]=FIFTH;
	map["SIXTH"        ]=SIXTH;
	map["SEVENTH"      ]=SEVENTH;
	map["EIGHTH"       ]=EIGHTH;
	map["NINTH"        ]=NINTH;
	map["TENTH"        ]=TENTH;
			    			  
	map["ELEVENTH"     ]=ELEVENTH;
	map["TWELFTH"      ]=TWELFTH;
	map["THIRTEENTH"   ]=THIRTEENTH;
	map["FOURTEENTH"   ]=FOURTEENTH;
	map["FIFTEENTH"    ]=FIFTEENTH;
	map["SIXTEENTH"    ]=SIXTEENTH;
	map["SEVENTEENTH"  ]=SEVENTEENTH;
	map["EIGHTTEENTH"  ]=EIGHTTEENTH;
	map["NINTEENTH"    ]=NINTEENTH;
	map["TWENTIETH"    ]=TWENTIETH;
			    			  
	map["TWENTYFIRST"  ]=TWENTYFIRST;
	map["TWENTYSECOND" ]=TWENTYSECOND;
	map["TWENTYTHIRD"  ]=TWENTYTHIRD;
	map["TWENTYFOURTH" ]=TWENTYFOURTH;
	map["TWENTYFIFTH"  ]=TWENTYFIFTH;
	map["TWENTYSIXTH"  ]=TWENTYSIXTH;
	map["TWENTYSEVENTH"]=TWENTYSEVENTH;
	map["TWENTYEIGHTH" ]=TWENTYEIGHTH;
	map["TWENTYNINTH"  ]=TWENTYNINTH;
	map["THIRTIETH"    ]=THIRTIETH;
			    			  
	map["THIRTYFIRST"  ]=THIRTYFIRST;
	map["THIRTYSECOND" ]=THIRTYSECOND;
	map["THIRTYTHIRD"  ]=THIRTYTHIRD;
	map["THIRTYFOURTH" ]=THIRTYFOURTH;
	map["THIRTYFIFTH"  ]=THIRTYFIFTH;
	map["THIRTYSIXTH"  ]=THIRTYSIXTH;
	map["THIRTYSEVENTH"]=THIRTYSEVENTH;
	map["THIRTYEIGHTH" ]=THIRTYEIGHTH;
	map["THIRTYNINTH"  ]=THIRTYNINTH;
	map["FOURTIETH"    ]=FOURTIETH;
			    			  
	map["FOURTYFIRST"  ]=FOURTYFIRST;
	map["FOURTYSECOND" ]=FOURTYSECOND;
	map["FOURTYTHIRD"  ]=FOURTYTHIRD;
      }
    
    if (!map.count(s))
      error();
    
    return map[s];
  }



  //------------------------------------------------------
  // FEFamily specialization
  template <>
  FEFamily string_to_enum<FEFamily> (const std::string& s)
  {
    static std::map<std::string, FEFamily> map;
    
    // Initialize map on first call
    if (map.empty())
      {
	map["LAGRANGE"  ]=LAGRANGE;
	map["HIERARCHIC"]=HIERARCHIC;
	map["MONOMIAL"  ]=MONOMIAL;
	map["XYZ"       ]=XYZ;
#ifdef ENABLE_HIGHER_ORDER_SHAPES
	map["BERNSTEIN" ]=BERNSTEIN;
	map["SZABAB"    ]=SZABAB;
#endif
#ifdef ENABLE_INFINITE_ELEMENTS
	//TODO:[BSK] which ones?
#endif
	map["CLOUGH"    ]=CLOUGH;
      }
    
    if (!map.count(s))
      error();
    
    return map[s];
  }
}



