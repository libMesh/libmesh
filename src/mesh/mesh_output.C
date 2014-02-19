// The libMesh Finite Element Library.
// Copyright (C) 2002-2014 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/equation_systems.h"
#include "libmesh/mesh_output.h"
#include "libmesh/parallel.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/unstructured_mesh.h"

namespace libMesh
{

template <class MT>
void MeshOutput<MT>::
_build_variable_names_and_solution_vector (const EquationSystems& es,
                                           std::vector<Number>& soln,
                                           std::vector<std::string>& names,
                                           const std::set<std::string>* system_names)
{
  if(!_is_parallel_format)
  {
    // We need a serial mesh for MeshOutput for now
    const_cast<EquationSystems&>(es).allgather();
  }

  es.build_variable_names  (names, NULL, system_names);
  es.build_solution_vector (soln, system_names);

  // For now, if we're doing a parallel format we're going to broadcast the vector from processor 0
  // to all of the processors to mimic what build_solution_vector used to do.
  // this is TERRIBLE and WASTEFUL but it's only temporary until we redesign the output of build_solution_vector
  // and the inputs to the I/O... both of which should actually be NumericVectors....
  if(_is_parallel_format)
  {
    size_t size = soln.size();
    _obj->comm().broadcast(size);

    if(_obj->comm().rank())
      soln.resize(size);

    _obj->comm().broadcast(soln);
  }
}



// Instantiate for our Mesh types.  If this becomes too cumbersome later,
// move any functions in this file to the header file instead.
template class MeshOutput<MeshBase>;
template class MeshOutput<UnstructuredMesh>;
template class MeshOutput<ParallelMesh>;

} // namespace libMesh
