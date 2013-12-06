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
#include <iostream>

// Local includes
#include "libmesh/elem_type.h"

namespace libMesh
{

// ------------------------------------------------------------
// Element type definitions


std::string ElementTypes::basic_name (const ElemType t)
{
  std::string its_name;
  switch (t)
    {
    case EDGE2:
    case EDGE3:
    case EDGE4:
      {
	its_name = "Edge";
	break;
      }

    case TRI3:
    case TRI6:
      {
	its_name = "Triangle";
	break;
      }

    case QUAD4:
    case QUAD8:
    case QUAD9:
      {
	its_name = "Quadrilateral";
	break;
      }

    case TET4:
    case TET10:
      {
	its_name = "Tetrahedron";
	break;
      }

    case HEX8:
    case HEX20:
    case HEX27:
      {
	its_name = "Hexahedron";
	break;
      }

    case PRISM6:
    case PRISM18:
      {
	its_name = "Prism";
	break;
      }

    case PYRAMID5:
    case PYRAMID14:
      {
	its_name = "Pyramid";
	break;
      }



#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

    // infinite elements
    case INFEDGE2:
      {
	its_name = "Infinite Edge";
	break;
      }

    case INFQUAD4:
    case INFQUAD6:
      {
	its_name = "Infinite Quadrilateral";
	break;
      }

    case INFHEX8:
    case INFHEX16:
    case INFHEX18:
      {
 	its_name = "Infinite Hexahedron";
	break;
      }

    case INFPRISM6:
    case INFPRISM12:
      {
	its_name = "Infinite Prism";
	break;
      }

#endif


    default:
      {
	libMesh::out << "Undefined element type!." << std::endl;
	libmesh_error();
      }
    }
  return its_name;
}


std::string ElementTypes::name(const ElemType t)
{
  std::string its_name;
  switch (t)
    {
    case EDGE2:
      {
	its_name = "Edge 2";
	break;
      }

    case EDGE3:
      {
	its_name = "Edge 3";
	break;
      }

    case EDGE4:
      {
	its_name = "Edge 4";
	break;
      }

    case TRI3:
      {
	its_name = "Tri 3";
	break;
      }

    case TRI6:
      {
	its_name = "Tri 6";
	break;
      }

    case QUAD4:
      {
	its_name = "Quad 4";
	break;
      }

    case QUAD8:
      {
	its_name = "Quad 8";
	break;
      }

    case QUAD9:
      {
	its_name = "Quad 9";
	break;
      }

    case TET4:
      {
	its_name = "Tet 4";
	break;
      }

    case TET10:
      {
	its_name = "Tet 10";
	break;
      }

    case HEX8:
      {
	its_name = "Hex 8";
	break;
      }

    case HEX20:
      {
	its_name = "Hex 20";
	break;
      }

    case HEX27:
      {
	its_name = "Hex 27";
	break;
      }

    case PRISM6:
      {
	its_name = "Prism 6";
	break;
      }

    case PRISM18:
      {
	its_name = "Prism 8";
	break;
      }

    case PYRAMID5:
      {
	its_name = "Pyramid 5";
	break;
      }

    case PYRAMID14:
      {
	its_name = "Pyramid 14";
	break;
      }



#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

    case INFEDGE2:
      {
	its_name = "Infinite Edge 2";
	break;
      }

    case INFQUAD4:
      {
	its_name = "Infinite Quad 4";
	break;
      }

    case INFQUAD6:
      {
	its_name = "Infinite Quad 6";
	break;
      }

    case INFHEX8:
      {
 	its_name = "Infinite Hex 8";
	break;
      }

    case INFHEX16:
      {
 	its_name = "Infinite Hex 16";
	break;
      }

    case INFHEX18:
      {
 	its_name = "Infinite Hex 18";
	break;
      }

    case INFPRISM6:
      {
	its_name = "Infinite Prism 6";
	break;
      }

    case INFPRISM12:
      {
	its_name = "Infinite Prism 12";
	break;
      }

#endif



    default:
      {
	libMesh::err << "Undefined element type!." << std::endl;
	libmesh_error();
      }
    }
  return its_name;
}

} // namespace libMesh
