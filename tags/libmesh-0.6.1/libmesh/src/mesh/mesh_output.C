// $Id: mesh_output.C,v 1.4 2007-10-21 20:48:50 benkirk Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson
  
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
#include "mesh_output.h"
#include "equation_systems.h"
#include "unstructured_mesh.h"

template <class MT>
void MeshOutput<MT>::
_build_variable_names_and_solution_vector (const EquationSystems& es,
					   std::vector<Number>& soln,
					   std::vector<std::string>& names)
{
  es.build_variable_names  (names);
  es.build_solution_vector (soln);
}



// Instantiate for our Mesh types.  If this becomes too cumbersome later,
// move any functions in this file to the header file instead.
template class MeshOutput<MeshBase>;
template class MeshOutput<UnstructuredMesh>;
