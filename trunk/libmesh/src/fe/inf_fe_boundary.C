// $Id: inf_fe_boundary.C,v 1.9 2005-02-22 22:17:38 jwpeterson Exp $

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


// Local includes
#include "libmesh_config.h"
#ifdef ENABLE_INFINITE_ELEMENTS
#include "inf_fe.h"




//-------------------------------------------------------
// Method for 2D, 3D -- see inf_fe_1D.C for a 1D version of this
template <unsigned int Dim, FEFamily T_radial, InfMapType T_base>
void InfFE<Dim,T_radial,T_base>::reinit(const Elem*,
					const unsigned int)
{
  // We don't do this for 1D elements!
  //assert (Dim != 1);

  std::cerr << "ERROR: Boundary conditions for infinite elements "
	    << "not implemented!" << std::endl;
  error();
}




template <unsigned int Dim, FEFamily T_radial, InfMapType T_base>
void InfFE<Dim,T_radial,T_base>::init_face_shape_functions(const std::vector<Point>&,
							   const Elem*)
{
  // We don't do this for 1D elements!
  //assert (Dim != 1);
  
  std::cerr << "ERROR: Boundary conditions for infinite elements "
	    << "not implemented!" << std::endl;
  error();
}




//--------------------------------------------------------------
// Explicit instantiations - doesn't make sense in 1D, but as 
// long as we only return errors, we are fine... ;-) 
#include "inf_fe_instantiate_1D.h"
#include "inf_fe_instantiate_2D.h"
#include "inf_fe_instantiate_3D.h"



#endif //ifdef ENABLE_INFINITE_ELEMENTS

