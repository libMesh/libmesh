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
#include <iostream>

// Local includes
#include "libmesh_common.h"
#include "quadrature_rules.h"

namespace libMesh
{



std::string QuadratureRules::name (const QuadratureType q) 
{
  std::string its_name;
  
  switch (q)
    {
      
    case QGAUSS:
      its_name = "Gauss-Legendre Quadrature";
      break;

    case QJACOBI_1_0:
      its_name = "Jacobi(1,0)-Gauss Quadrature";
      break;

    case QJACOBI_2_0:
      its_name = "Jacobi(2,0)-Gauss Quadrature";
      break;

    case QSIMPSON:
      its_name = "Simpson Rule";
      break;

    case QTRAP:
      its_name = "Trapezoidal Rule";
      break;


    default:
      { 
	libMesh::err << "ERROR: Bad qt=" << q << std::endl;
	libmesh_error();
      }
    }
  
  return its_name;
}

} // namespace libMesh





