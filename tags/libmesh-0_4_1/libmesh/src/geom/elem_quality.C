// $Id: elem_quality.C,v 1.7 2003-09-25 21:46:56 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002-2003  Benjamin S. Kirk, John W. Peterson
  
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
#include <sstream> 

// Local includes
#include "libmesh_common.h"
#include "enum_elem_quality.h"
#include "elem_quality.h"

// ------------------------------------------------------------
// Quality function definitions

/**
 * This function returns a string containing some name
 * for q.  Useful for asking the enum what its name is.
 * I added this since you may want a simple way to attach
 * a name or description to the ElemQuality enums.
 * It can be removed if it is found to be useless.
 */
std::string Quality::name (const ElemQuality q) 
{
  std::string its_name;
  
  switch (q)
    {
      
    case ASPECT_RATIO:
      its_name = "Aspect Ratio";
      break;
      
    case SKEW:
      its_name = "Skew";
      break;
      
    case SHEAR:
      its_name = "Shear";
      break;
      
    case SHAPE:
      its_name = "Shape";
      break;
      
    case MAX_ANGLE:
      its_name = "Maximum Angle";
      break;

    case MIN_ANGLE:
      its_name = "Minimum Angle";
      break;
      
    case CONDITION:
      its_name = "Condition Number";
      break;
      
    case DISTORTION:
      its_name = "Distortion";
      break;
      
    case TAPER:
      its_name = "Taper";
      break;
      
    case WARP:
      its_name = "Warp";
      break;
      
    case STRETCH:
      its_name = "Stretch";
      break;
      
    case DIAGONAL:
      its_name = "Diagonal";
      break;

    case ASPECT_RATIO_BETA:
      its_name = "AR Beta";
      break;
      
    case ASPECT_RATIO_GAMMA:
      its_name = "AR Gamma";
      break;

    case SIZE:
      its_name = "Size";
      break;
      
    case JACOBIAN:
      its_name = "Jacobian";
      break;

    default:
      its_name = "Unknown";
      break;
    }
  
  return its_name;
}





/**
 * This function returns a string containing a short 
 * description of q.  Useful for asking the enum what
 * it computes.
 */
std::string Quality::describe (const ElemQuality q) 
{
  
  std::ostringstream desc; 
  
  switch (q)
    {
      
    case ASPECT_RATIO:
      desc << "Max edge length ratio" << std::endl
	   << "at element center." << std::endl
	   << std::endl
	   << "Suggested ranges:" << std::endl
	   << "Hexes: (1 -> 4)" << std::endl
	   << "Quads: (1 -> 4)";
      break;
      
    case SKEW:
      desc << "Maximum |cos A|, where A" << std::endl
	   << "is the angle between edges" << std::endl
	   << "at element center." << std::endl
	   << std::endl
	   << "Suggested ranges:" << std::endl
	   << "Hexes: (0 -> 0.5)" << std::endl
	   << "Quads: (0 -> 0.5)";
      break;
      
    case SHEAR:
      desc << "DIM / K(Js)" << std::endl
	   << std::endl
	   << "DIM   = element dimension." << std::endl
	   << "K(Js) = Condition number of " << std::endl
	   << "        Jacobian skew matrix." << std::endl
	   << std::endl
	   << "Suggested ranges:" << std::endl
	   << "Hexes(DIM=3): (0.3 -> 1)" << std::endl
	   << "Quads(DIM=2): (0.3 -> 1)";
      break;
      
    case SHAPE:
      desc << "DIM / K(Jw)" << std::endl
	   << std::endl
	   << "DIM   = element dimension." << std::endl
	   << "K(Jw) = Condition number of " << std::endl
	   << "        weighted Jacobian" << std::endl
	   << "        matrix." << std::endl
	   << std::endl
	   << "Suggested ranges:" << std::endl
	   << "Hexes(DIM=3): (0.3 -> 1)" << std::endl
	   << "Tets(DIM=3): (0.2 -> 1)" << std::endl
	   << "Quads(DIM=2): (0.3 -> 1).";
      break;
      
    case MAX_ANGLE:
      desc << "Largest included angle." << std::endl
	   << std::endl
	   << "Suggested ranges:" << std::endl
	   << "Quads: (90 -> 135)" << std::endl
	   << "Triangles: (60 -> 90)";
      break;

    case MIN_ANGLE:
      desc << "Smallest included angle." << std::endl
	   << std::endl
	   << "Suggested ranges:" << std::endl
	   << "Quads: (45 -> 90)" << std::endl
	   << "Triangles: (30 -> 60)";
      break;
      
    case CONDITION:
      desc << "Condition number of the" << std::endl
	   << "Jacobian matrix." << std::endl
	   << std::endl
	   << "Suggested ranges:" << std::endl
	   << "Quads: (1 -> 4)" << std::endl 
	   << "Hexes: (1 -> 8)" << std::endl
	   << "Tris: (1 -> 1.3)" << std::endl
	   << "Tets: (1 -> 3)";
      break;
      
    case DISTORTION:
      desc << "min |J| * A / <A>" << std::endl
	   << std::endl
	   << "|J| = norm of Jacobian matrix" << std::endl
	   << " A  = actual area" << std::endl
	   << "<A> = reference area" << std::endl
	   << std::endl
	   << "Suggested ranges:" << std::endl
	   << "Quads: (0.6 -> 1), <A>=4" << std::endl 
	   << "Hexes: (0.6 -> 1), <A>=8" << std::endl
	   << "Tris: (0.6 -> 1), <A>=1/2" << std::endl
	   << "Tets: (0.6 -> 1), <A>=1/6";
      break;
      
    case TAPER:
      desc << "Maximum ratio of lengths" << std::endl
	   << "derived from opposited edges." << std::endl
	   << std::endl
	   << "Suggested ranges:" << std::endl
	   << "Quads: (0.7 -> 1)" << std::endl
	   << "Hexes: (0.4 -> 1)";
      break;
      
    case WARP:
      desc << "cos D" << std::endl
	   << std::endl
	   << "D = minimum dihedral angle" << std::endl
	   << "    formed by diagonals." << std::endl
	   << std::endl
	   << "Suggested ranges:" << std::endl
	   << "Quads: (0.9 -> 1)";
      break;
      
    case STRETCH:
      desc << "Sqrt(3) * L_min / L_max" << std::endl
	   << std::endl
	   << "L_min = minimum edge length." << std::endl
	   << "L_max = maximum edge length." << std::endl
	   << std::endl
	   << "Suggested ranges:" << std::endl
	   << "Quads: (0.25 -> 1)" << std::endl
	   << "Hexes: (0.25 -> 1)";
      break;
      
    case DIAGONAL:
      desc << "D_min / D_max" << std::endl
	   << std::endl
	   << "D_min = minimum diagonal." << std::endl
	   << "D_max = maximum diagonal." << std::endl
	   << std::endl
	   << "Suggested ranges:" << std::endl
	   << "Hexes: (0.65 -> 1)";
      break;

    case ASPECT_RATIO_BETA:
      desc << "CR / (3 * IR)" << std::endl
	   << std::endl
	   << "CR = circumsphere radius" << std::endl
	   << "IR = inscribed sphere radius" << std::endl
	   << std::endl
	   << "Suggested ranges:" << std::endl
	   << "Tets: (1 -> 3)";
      break;
      
    case ASPECT_RATIO_GAMMA:
      desc << "S^(3/2) / 8.479670 * V" << std::endl
	   << std::endl
	   << "S = sum(si*si/6)" << std::endl
	   << "si = edge length" << std::endl
	   << "V = volume" << std::endl
	   << std::endl
	   << "Suggested ranges:" << std::endl
	   << "Tets: (1 -> 3)";
      break;

    case SIZE:
      desc << "min (|J|, |1/J|)" << std::endl
	   << std::endl
	   << "|J| = norm of Jacobian matrix." << std::endl
	   << std::endl
	   << "Suggested ranges:" << std::endl
	   << "Quads: (0.3 -> 1)" << std::endl 
	   << "Hexes: (0.5 -> 1)" << std::endl
	   << "Tris: (0.25 -> 1)" << std::endl
	   << "Tets: (0.2 -> 1)";
      break;
      
    case JACOBIAN:
      desc << "Minimum Jacobian divided by" << std::endl
	   << "the lengths of the DIM" << std::endl
	   << "largest edge vectors." << std::endl
	   << std::endl
	   << "DIM = element dimension." << std::endl
	   << std::endl
	   << "Suggested ranges:" << std::endl
	   << "Quads: (0.5 -> 1)" << std::endl 
	   << "Hexes: (0.5 -> 1)" << std::endl
	   << "Tris: (0.5 -> 1.155)" << std::endl
	   << "Tets: (0.5 -> 1.414)";
      break;

    default:
      desc << "Unknown";
      break;
    }
  
  return desc.str();
}


/**
 * Returns all valid quality metrics for
 * element type t.
 */
std::vector<ElemQuality> Quality::valid(const ElemType t)
{
  std::vector<ElemQuality> v;
  
  switch (t)
    {
    case EDGE2:
    case EDGE3:
    case EDGE4:
      {
	// None yet
	break;
      }

    case TRI3:
    case TRI6:
      {
	v.resize(7);
	v[0] = MAX_ANGLE;
	v[1] = MIN_ANGLE;
	v[2] = CONDITION;
	v[3] = JACOBIAN;
	v[4] = SIZE;
	v[5] = SHAPE;
	v[6] = DISTORTION;
	break;
      }

    case QUAD4:
    case QUAD8:
    case QUAD9:
      {
	v.resize(13);
	v[0]  = ASPECT_RATIO;
	v[1]  = SKEW;
	v[2]  = TAPER;
	v[3]  = WARP;
	v[4]  = STRETCH;
	v[5]  = MIN_ANGLE;
	v[6]  = MAX_ANGLE;
	v[7]  = CONDITION;
	v[8]  = JACOBIAN;
	v[9]  = SHEAR;
	v[10] = SHAPE;
	v[11] = SIZE;
	v[12] = DISTORTION;
	break;
      }

    case TET4:
    case TET10:
      {
	v.resize(7);
	v[0]  = ASPECT_RATIO_BETA;
	v[1]  = ASPECT_RATIO_GAMMA;
	v[2]  = CONDITION;
	v[3]  = JACOBIAN;
	v[4]  = SHAPE;
	v[5]  = SIZE;
	v[6]  = DISTORTION;
	break;
      }

    case HEX8:
    case HEX20:
    case HEX27:
      {
	v.resize(11);
	v[0]  = ASPECT_RATIO;
	v[1]  = SKEW;
	v[2]  = SHEAR;
	v[3] = SHAPE;
	v[4]  = CONDITION;
	v[5]  = JACOBIAN;
	v[6]  = DISTORTION;
	v[7]  = TAPER;
	v[8]  = STRETCH;
	v[9]  = DIAGONAL;
	v[10]  = SIZE;
	break;
      }

    case PRISM6:
    case PRISM18:
      {
	// None yet
	break;
      }

    case PYRAMID5:
      {
	// None yet
	break;
      }



#ifdef ENABLE_INFINITE_ELEMENTS

    case INFEDGE2:
      {
	// None yet
	break;
      }

    case INFQUAD4:
    case INFQUAD6:
      {
        // None yet
	break;
      }

    case INFHEX8:
    case INFHEX16:
    case INFHEX18:
      {
        // None yet
	break;
      }

    case INFPRISM6:
    case INFPRISM12:
      {
        // None yet
	break;
      }

#endif

      
    default:
      {
	std::cout << "Undefined element type!." << std::endl;
	error();
      }
    }
  
  return v;
}
