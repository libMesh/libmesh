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
#include <iostream>
#include <sstream>

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/elem_quality.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/enum_elem_quality.h"


namespace libMesh
{

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

    case MAX_DIHEDRAL_ANGLE:
      its_name = "Maximum Dihedral Angle";
      break;

    case MIN_DIHEDRAL_ANGLE:
      its_name = "Minimum Dihedral Angle";
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

    case SCALED_JACOBIAN:
      its_name = "Scaled Jacobian";
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

    case EDGE_LENGTH_RATIO:
    case ASPECT_RATIO:
      desc << "Max edge length ratio\n"
           << "at element center.\n"
           << '\n'
           << "Suggested ranges:\n"
           << "Hexes: (1 -> 4)\n"
           << "Quads: (1 -> 4)";
      break;

    case SKEW:
      desc << "Maximum |cos A|, where A\n"
           << "is the angle between edges\n"
           << "at element center.\n"
           << '\n'
           << "Suggested ranges:\n"
           << "Hexes: (0 -> 0.5)\n"
           << "Quads: (0 -> 0.5)";
      break;

    case SHEAR:
      desc << "LIBMESH_DIM / K(Js)\n"
           << '\n'
           << "LIBMESH_DIM   = element dimension.\n"
           << "K(Js) = Condition number of \n"
           << "        Jacobian skew matrix.\n"
           << '\n'
           << "Suggested ranges:\n"
           << "Hexes(LIBMESH_DIM=3): (0.3 -> 1)\n"
           << "Quads(LIBMESH_DIM=2): (0.3 -> 1)";
      break;

    case SHAPE:
      desc << "LIBMESH_DIM / K(Jw)\n"
           << '\n'
           << "LIBMESH_DIM   = element dimension.\n"
           << "K(Jw) = Condition number of \n"
           << "        weighted Jacobian\n"
           << "        matrix.\n"
           << '\n'
           << "Suggested ranges:\n"
           << "Hexes(LIBMESH_DIM=3): (0.3 -> 1)\n"
           << "Tets(LIBMESH_DIM=3): (0.2 -> 1)\n"
           << "Quads(LIBMESH_DIM=2): (0.3 -> 1).";
      break;

    case MAX_ANGLE:
      desc << "Largest angle between all adjacent pairs of edges (in 2D, sides).\n"
           << '\n'
           << "Suggested ranges:\n"
           << "Quads: (90 -> 135)\n"
           << "Triangles: (60 -> 90)";
      break;

    case MIN_ANGLE:
      desc << "Smallest angle between all adjacent pairs of edges (in 2D, sides).\n"
           << '\n'
           << "Suggested ranges:\n"
           << "Quads: (45 -> 90)\n"
           << "Triangles: (30 -> 60)";
      break;

    case MAX_DIHEDRAL_ANGLE:
      desc << "Largest angle between all adjacent pairs of sides (in 2D, equivalent to MAX_ANGLE).\n"
           << '\n'
           << "Suggested ranges:\n"
           << "Quads: (90 -> 135)\n"
           << "Triangles: (60 -> 90)";
      break;

    case MIN_DIHEDRAL_ANGLE:
      desc << "Smallest angle between all adjacent pairs of sides (in 2D, equivalent to MIN_ANGLE).\n"
           << '\n'
           << "Suggested ranges:\n"
           << "Quads: (45 -> 90)\n"
           << "Triangles: (30 -> 60)";
      break;

    case CONDITION:
      desc << "Condition number of the\n"
           << "Jacobian matrix.\n"
           << '\n'
           << "Suggested ranges:\n"
           << "Quads: (1 -> 4)\n"
           << "Hexes: (1 -> 8)\n"
           << "Tris: (1 -> 1.3)\n"
           << "Tets: (1 -> 3)";
      break;

    case DISTORTION:
      desc << "min |J| * A / <A>\n"
           << '\n'
           << "|J| = norm of Jacobian matrix\n"
           << " A  = actual area\n"
           << "<A> = reference area\n"
           << '\n'
           << "Suggested ranges:\n"
           << "Quads: (0.6 -> 1), <A>=4\n"
           << "Hexes: (0.6 -> 1), <A>=8\n"
           << "Tris: (0.6 -> 1), <A>=1/2\n"
           << "Tets: (0.6 -> 1), <A>=1/6";
      break;

    case TAPER:
      desc << "Maximum ratio of lengths\n"
           << "derived from opposite edges.\n"
           << '\n'
           << "Suggested ranges:\n"
           << "Quads: (0.7 -> 1)\n"
           << "Hexes: (0.4 -> 1)";
      break;

    case WARP:
      desc << "cos D\n"
           << '\n'
           << "D = minimum dihedral angle\n"
           << "    formed by diagonals.\n"
           << '\n'
           << "Suggested ranges:\n"
           << "Quads: (0.9 -> 1)";
      break;

    case STRETCH:
      desc << "Sqrt(3) * L_min / L_max\n"
           << '\n'
           << "L_min = minimum edge length.\n"
           << "L_max = maximum edge length.\n"
           << '\n'
           << "Suggested ranges:\n"
           << "Quads: (0.25 -> 1)\n"
           << "Hexes: (0.25 -> 1)";
      break;

    case DIAGONAL:
      desc << "D_min / D_max\n"
           << '\n'
           << "D_min = minimum diagonal.\n"
           << "D_max = maximum diagonal.\n"
           << '\n'
           << "Suggested ranges:\n"
           << "Hexes: (0.65 -> 1)";
      break;

    case ASPECT_RATIO_BETA:
      desc << "CR / (3 * IR)\n"
           << '\n'
           << "CR = circumsphere radius\n"
           << "IR = inscribed sphere radius\n"
           << '\n'
           << "Suggested ranges:\n"
           << "Tets: (1 -> 3)";
      break;

    case ASPECT_RATIO_GAMMA:
      desc << "S^(3/2) / 8.479670 * V\n"
           << '\n'
           << "S = sum(si*si/6)\n"
           << "si = edge length\n"
           << "V = volume\n"
           << '\n'
           << "Suggested ranges:\n"
           << "Tets: (1 -> 3)";
      break;

    case SIZE:
      desc << "min (|J|, |1/J|)\n"
           << '\n'
           << "|J| = norm of Jacobian matrix.\n"
           << '\n'
           << "Suggested ranges:\n"
           << "Quads: (0.3 -> 1)\n"
           << "Hexes: (0.5 -> 1)\n"
           << "Tris: (0.25 -> 1)\n"
           << "Tets: (0.2 -> 1)";
      break;

    case JACOBIAN:
    case SCALED_JACOBIAN:
      desc << "Minimum nodal Jacobian.\n"
           << "The nodal Jacobians are computed by taking the cross product (2D) or scalar product (3D) of the adjacent edges that meet at that node.\n"
           << "In the SCALED_JACOBIAN case, we also then divide by the lengths of each of the associated edges.\n"
           << "For Pyramid elements where four edges meet at the apex node, special handling is required.\n"
           << '\n'
           << "Suggested acceptable ranges (from Cubit documentation) for SCALED_JACOBIAN metric:\n"
           << "Quads/Hexes: (0.5 -> 1)\n"
           << "Tris/Tets: (0.2 -> 1.0)";
      break;

    default:
      desc << "Unknown";
      break;
    }

  return desc.str();
}


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
    case TRISHELL3:
    case TRI6:
    case TRI7:
      {
        v = {
          CONDITION,
          DISTORTION,
          EDGE_LENGTH_RATIO,
          JACOBIAN,
          SCALED_JACOBIAN,
          MAX_ANGLE,
          MIN_ANGLE,
          MAX_DIHEDRAL_ANGLE,
          MIN_DIHEDRAL_ANGLE,
          SHAPE,
          SIZE
        };

        break;
      }

    case QUAD4:
    case QUADSHELL4:
    case QUAD8:
    case QUADSHELL8:
    case QUAD9:
    case QUADSHELL9:
      {
        v = {
          ASPECT_RATIO,
          CONDITION,
          DISTORTION,
          EDGE_LENGTH_RATIO,
          JACOBIAN,
          SCALED_JACOBIAN,
          MAX_ANGLE,
          MIN_ANGLE,
          MAX_DIHEDRAL_ANGLE,
          MIN_DIHEDRAL_ANGLE,
          SHAPE,
          SHEAR,
          SIZE,
          SKEW,
          STRETCH,
          TAPER,
          WARP
        };

        break;
      }

    case TET4:
    case TET10:
    case TET14:
      {
        v = {
          ASPECT_RATIO_BETA,
          ASPECT_RATIO_GAMMA,
          CONDITION,
          DISTORTION,
          JACOBIAN,
          SCALED_JACOBIAN,
          MAX_ANGLE,
          MIN_ANGLE,
          MAX_DIHEDRAL_ANGLE,
          MIN_DIHEDRAL_ANGLE,
          SHAPE,
          SIZE
        };

        break;
      }

    case HEX8:
    case HEX20:
    case HEX27:
      {
        v = {
          ASPECT_RATIO,
          CONDITION,
          DIAGONAL,
          DISTORTION,
          JACOBIAN,
          SCALED_JACOBIAN,
          MAX_ANGLE,
          MIN_ANGLE,
          MAX_DIHEDRAL_ANGLE,
          MIN_DIHEDRAL_ANGLE,
          SHAPE,
          SHEAR,
          SIZE,
          SKEW,
          STRETCH,
          TAPER
        };

        break;
      }

    case PRISM6:
    case PRISM18:
    case PRISM20:
    case PRISM21:
      {
        v = {
          EDGE_LENGTH_RATIO,
          MAX_ANGLE,
          MIN_ANGLE,
          MAX_DIHEDRAL_ANGLE,
          MIN_DIHEDRAL_ANGLE,
        };

        break;
      }

    case PYRAMID5:
    case PYRAMID13:
    case PYRAMID14:
    case PYRAMID18:
      {
        v = {
          EDGE_LENGTH_RATIO,
          MAX_ANGLE,
          MIN_ANGLE,
          MAX_DIHEDRAL_ANGLE,
          MIN_DIHEDRAL_ANGLE,
        };

        break;
      }



#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

    case INFEDGE2:
      {
        // None yet
        break;
      }

    case INFQUAD4:
    case INFQUAD6:
    case INFHEX8:
    case INFHEX16:
    case INFHEX18:
    case INFPRISM6:
    case INFPRISM12:
      {
        v = {
          MAX_ANGLE,
          MIN_ANGLE,
          MAX_DIHEDRAL_ANGLE,
          MIN_DIHEDRAL_ANGLE,
        };

        break;
      }

#endif


    default:
      libmesh_error_msg("Undefined element type!");
    }

  return v;
}

} // namespace libMesh
