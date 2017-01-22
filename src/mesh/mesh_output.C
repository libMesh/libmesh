// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/distributed_mesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/mesh_output.h"
#include "libmesh/parallel.h"
#include "libmesh/unstructured_mesh.h"
#include "libmesh/numeric_vector.h"

namespace libMesh
{

template <class MT>
void MeshOutput<MT>::write_equation_systems (const std::string & fname,
                                             const EquationSystems & es,
                                             const std::set<std::string> * system_names)
{
  LOG_SCOPE("write_equation_systems()", "MeshOutput");

  // We may need to gather and/or renumber a DistributedMesh to output
  // it, making that const qualifier in our constructor a dirty lie
  MT & my_mesh = const_cast<MT &>(*_obj);

  // If we're asked to write data that's associated with a different
  // mesh, output files full of garbage are the result.
  libmesh_assert_equal_to(&es.get_mesh(), _obj);

  // A non-renumbered mesh may not have a contiguous numbering, and
  // that needs to be fixed before we can build a solution vector.
  if (my_mesh.max_elem_id() != my_mesh.n_elem() ||
      my_mesh.max_node_id() != my_mesh.n_nodes())
    {
      // If we were allowed to renumber then we should have already
      // been properly renumbered...
      libmesh_assert(!my_mesh.allow_renumbering());

      libmesh_do_once(libMesh::out <<
                      "Warning:  This MeshOutput subclass only supports meshes which are contiguously renumbered!"
                      << std::endl;);

      my_mesh.allow_renumbering(true);

      my_mesh.renumber_nodes_and_elements();

      // Not sure what good going back to false will do here, the
      // renumbering horses have already left the barn...
      my_mesh.allow_renumbering(false);
    }

  MeshSerializer serialize(const_cast<MT &>(*_obj), !_is_parallel_format, _serial_only_needed_on_proc_0);

  // Build the list of variable names that will be written.
  std::vector<std::string> names;
  es.build_variable_names  (names, libmesh_nullptr, system_names);

  if (!_is_parallel_format)
    {
      // Build the nodal solution values & get the variable
      // names from the EquationSystems object
      std::vector<Number> soln;
      es.build_solution_vector (soln, system_names);

      this->write_nodal_data (fname, soln, names);
    }
  else // _is_parallel_format
    {
      UniquePtr<NumericVector<Number> > parallel_soln =
        es.build_parallel_solution_vector(system_names);

      this->write_nodal_data (fname, *parallel_soln, names);
    }
}



template <class MT>
void MeshOutput<MT>::write_nodal_data (const std::string & fname,
                                       const NumericVector<Number> & parallel_soln,
                                       const std::vector<std::string> & names)
{
  // This is the fallback implementation for parallel I/O formats that
  // do not yet implement proper writing in parallel, and instead rely
  // on the full solution vector being available on all processors.
  std::vector<Number> soln;
  parallel_soln.localize(soln);
  this->write_nodal_data(fname, soln, names);
}



// Instantiate for our Mesh types.  If this becomes too cumbersome later,
// move any functions in this file to the header file instead.
template class MeshOutput<MeshBase>;
template class MeshOutput<UnstructuredMesh>;
template class MeshOutput<DistributedMesh>;

} // namespace libMesh
