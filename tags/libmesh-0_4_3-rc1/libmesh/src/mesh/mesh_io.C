// $Id: mesh_io.C,v 1.2 2004-04-07 21:42:38 benkirk Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2004  Benjamin S. Kirk, John W. Peterson
  
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
#include "mesh_base.h"
#include "mesh.h"
#include "mesh_io.h"
#include "equation_systems.h"



// ------------------------------------------------------------
// MeshIO members
template <class MT>
void MeshIO<MT>::write_equation_systems (const std::string& fname,
					 const EquationSystems& es)
{
  // Build the nodal solution values & get the variable
  // names from the EquationSystems object
  std::vector<Number>      soln;
  std::vector<std::string> names;

  es.build_variable_names  (names);
  es.build_solution_vector (soln);

  this->write_nodal_data (fname, soln, names);  
}



template <class MT>
void MeshIO<MT>::skip_comment_lines (std::istream &in,
				     const char comment_start)
{    
  char c, line[256];
  
  while (in.get(c), c==comment_start) 
    in.getline (line, 255);
  
  // put back first character of
  // first non-comment line
  in.putback (c);
}


//--------------------------------------------------------------
// template instantiations
template class MeshIO<MeshBase>;
template class MeshIO<Mesh>;
