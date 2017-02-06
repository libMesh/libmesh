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



#include <sstream>
#include <fstream>
#include <string>
#include <iomanip>

#include "libmesh/dof_map.h"
#include "libmesh/ensight_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe_interface.h"
#include "libmesh/libmesh.h"
#include "libmesh/system.h"
#include "libmesh/elem.h"

namespace libMesh
{

// Initialize the static data members by calling the static build functions.
std::map<ElemType, std::string> EnsightIO::_element_map = EnsightIO::build_element_map();

// Static function used to build the _element_map.
std::map<ElemType, std::string> EnsightIO::build_element_map()
{
  std::map<ElemType, std::string> ret;
  ret[EDGE2] = "bar2";
  ret[EDGE3] = "bar3";
  ret[QUAD4] = "quad4";
  ret[QUAD8] = "quad8";
  // ret[QUAD9] = "quad9"; // not supported
  ret[TRI3] = "tria3";
  ret[TRI6] = "tria6";
  ret[TET4] = "tetra4";
  ret[TET10] = "tetra10";
  ret[HEX8] = "hexa8";
  ret[HEX20] = "hexa20";
  // ret[HEX27] = "HEX27"; // not supported
  ret[PYRAMID5] = "pyramid5";
  return ret;
}


EnsightIO::EnsightIO (const std::string & filename,
                      const EquationSystems & eq) :
  MeshOutput<MeshBase> (eq.get_mesh()),
  _equation_systems(eq)
{
  if (_equation_systems.n_processors() == 1)
    _ensight_file_name = filename;
  else
    {
      std::ostringstream tmp_file;
      tmp_file << filename << "_rank" << _equation_systems.processor_id();
      _ensight_file_name = tmp_file.str();
    }
}



void EnsightIO::add_vector (const std::string & system_name,
                            const std::string & vec_description,
                            const std::string & u,
                            const std::string & v)
{
  libmesh_assert (_equation_systems.has_system(system_name));
  libmesh_assert (_equation_systems.get_system(system_name).has_variable(u));
  libmesh_assert (_equation_systems.get_system(system_name).has_variable(v));

  Vectors vec;
  vec.description = vec_description;
  vec.components.push_back(u);
  vec.components.push_back(v);

  _system_vars_map[system_name].EnsightVectors.push_back(vec);
}



void EnsightIO::add_vector (const std::string & system_name,
                            const std::string & vec_name,
                            const std::string & u,
                            const std::string & v,
                            const std::string & w)
{
  libmesh_assert(_equation_systems.has_system(system_name));
  libmesh_assert(_equation_systems.get_system(system_name).has_variable(u));
  libmesh_assert(_equation_systems.get_system(system_name).has_variable(v));
  libmesh_assert(_equation_systems.get_system(system_name).has_variable(w));

  Vectors vec;
  vec.description = vec_name;
  vec.components.push_back(u);
  vec.components.push_back(v);
  vec.components.push_back(w);
  _system_vars_map[system_name].EnsightVectors.push_back(vec);
}



void EnsightIO::add_scalar(const std::string & system_name,
                           const std::string & scl_description,
                           const std::string & s)
{
  libmesh_assert(_equation_systems.has_system(system_name));
  libmesh_assert(_equation_systems.get_system(system_name).has_variable(s));

  Scalars scl;
  scl.description = scl_description;
  scl.scalar_name = s;

  _system_vars_map[system_name].EnsightScalars.push_back(scl);
}



// This method must be implemented as it is pure virtual in
// the MeshOutput base class.
void EnsightIO::write (const std::string & name)
{
  // We may need to gather a DistributedMesh to output it, making that
  // const qualifier in our constructor a dirty lie
  MeshSerializer serialize(const_cast<MeshBase &>(this->mesh()), !_is_parallel_format);

  _ensight_file_name = name;
  this->write();
}



void EnsightIO::write (Real time)
{
  this->write_ascii(time);
  this->write_case();
}



void EnsightIO::write_ascii (Real time)
{
  _time_steps.push_back(time);

  this->write_geometry_ascii();
  this->write_solution_ascii();
}



void EnsightIO::write_geometry_ascii()
{
  std::ostringstream file;
  file << _ensight_file_name
       << ".geo"
       << std::setw(3)
       << std::setprecision(0)
       << std::setfill('0')
       << std::right
       << _time_steps.size()-1;

  // Open a stream to write the mesh
  std::ofstream mesh_stream(file.str().c_str());

  mesh_stream << "EnSight Gold Geometry File Format\n";
  mesh_stream << "Generated by \n";
  mesh_stream << "node id off\n";
  mesh_stream << "element id given\n";
  mesh_stream << "part\n";
  mesh_stream << std::setw(10) << 1 << "\n";
  mesh_stream << "uns-elements\n";
  mesh_stream << "coordinates\n";

  // mapping between nodal index and your coordinates
  typedef std::map<int, Point> mesh_nodes_map_t;
  typedef mesh_nodes_map_t::iterator mesh_nodes_iterator;
  mesh_nodes_map_t mesh_nodes_map;

  // Map for grouping elements of the same type
  typedef std::map<ElemType, std::vector<const Elem *> > ensight_parts_map_t;
  typedef ensight_parts_map_t::iterator ensight_parts_iterator;
  ensight_parts_map_t ensight_parts_map;

  const MeshBase & the_mesh = MeshOutput<MeshBase>::mesh();

  // Construct the various required maps
  {
    MeshBase::const_element_iterator       el     = the_mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = the_mesh.active_local_elements_end();

    for ( ; el != end_el ; ++el)
      {
        const Elem * elem = *el;
        ensight_parts_map[elem->type()].push_back(elem);

        for (unsigned int i = 0; i < elem->n_nodes(); i++)
          mesh_nodes_map[elem->node_id(i)] = elem->point(i);
      }
  }

  // Write number of local points
  mesh_stream << std::setw(10) << mesh_nodes_map.size() << "\n";

  // write x, y, and z node positions, build mapping between
  // ensight and libmesh node numbers.
  std::map <int, int> ensight_node_index;
  {
    mesh_nodes_iterator no_it = mesh_nodes_map.begin();
    const mesh_nodes_iterator no_end_it = mesh_nodes_map.end();

    for (unsigned direction=0; direction<3; ++direction)
      {
        for (int i = 1; no_it != no_end_it; ++no_it, i++)
          {
            mesh_stream << std::setw(12)
                        << std::setprecision(5)
                        << std::scientific
                        << no_it->second(direction)
                        << "\n";
            ensight_node_index[no_it->first] = i;
          }

        // Reset iterator to the beginning of the map
        no_it = mesh_nodes_map.begin();
      }
  }

  // Write parts
  {
    ensight_parts_iterator parts_it  =  ensight_parts_map.begin();
    const ensight_parts_iterator end_parts_it  =  ensight_parts_map.end();

    for (; parts_it != end_parts_it; ++parts_it)
      {
        // Look up this ElemType in the map, error if not present.
        std::map<ElemType, std::string>::iterator name_it = _element_map.find(parts_it->first);
        if (name_it == _element_map.end())
          libmesh_error_msg("Error: Unsupported ElemType " << parts_it->first << " for EnsightIO.");

        // Write element type
        mesh_stream << "\n" << name_it->second << "\n";

        std::vector<const Elem *> elem_ref  = parts_it->second;

        // Write number of element
        mesh_stream << std::setw(10) << elem_ref.size() << "\n";

        // Write element id
        for (std::size_t i = 0; i < elem_ref.size(); i++)
          mesh_stream << std::setw(10) << elem_ref[i]->id() << "\n";

        // Write connectivity
        for (std::size_t i = 0; i < elem_ref.size(); i++)
          {
            for (unsigned int j = 0; j < elem_ref[i]->n_nodes(); j++)
              {
                // tests!
                if (parts_it->first == QUAD9 && i==4)
                  continue;

                // tests!
                if (parts_it->first == HEX27 &&
                    (i==4    || i ==10 || i == 12 ||
                     i == 13 || i ==14 || i == 16 || i == 22))
                  continue;

                mesh_stream << std::setw(10) << ensight_node_index[elem_ref[i]->node_id(j)];
              }
            mesh_stream << "\n";
          }
      }
  }
}





void EnsightIO::write_case()
{
  std::ostringstream case_file;
  case_file << _ensight_file_name << ".case";

  // Open a stream for writing the case file.
  std::ofstream case_stream(case_file.str().c_str());

  case_stream << "FORMAT\n";
  case_stream << "type:  ensight gold\n\n";
  case_stream << "GEOMETRY\n";
  case_stream << "model:            1     " << _ensight_file_name << ".geo" << "***\n";

  system_vars_map_t::iterator sys_it = _system_vars_map.begin();
  const system_vars_map_t::iterator sys_end  = _system_vars_map.end();

  // Write Variable per node section
  if (sys_it != sys_end)
    case_stream << "\n\nVARIABLE\n";

  for (; sys_it != sys_end; ++sys_it)
    {
      for (std::size_t i=0; i < sys_it->second.EnsightScalars.size(); i++)
        {
          Scalars scalar = sys_it->second.EnsightScalars[i];
          case_stream << "scalar per node:   1  "
                      << scalar.description << " "
                      << _ensight_file_name << "_" << scalar.scalar_name << ".scl***\n";
        }

      for (std::size_t i=0; i < sys_it->second.EnsightVectors.size(); i++)
        {
          Vectors vec = sys_it->second.EnsightVectors[i];
          case_stream << "vector per node:      1    "
                      << vec.description << " "
                      << _ensight_file_name << "_" << vec.description << ".vec***\n";
        }

      // Write time step section
      if (_time_steps.size() != 0)
        {
          case_stream << "\n\nTIME\n";
          case_stream << "time set:             1\n";
          case_stream << "number of steps:   " << std::setw(10) << _time_steps.size() << "\n";
          case_stream << "filename start number:   " << std::setw(10) << 0 << "\n";
          case_stream << "filename increment:  " << std::setw(10) << 1 << "\n";
          case_stream << "time values:\n";
          for (std::size_t i = 0; i < _time_steps.size(); i++)
            case_stream << std::setw(12) << std::setprecision(5) << std::scientific << _time_steps[i] << "\n";
        }
    }
}


// Write scalar and vector solution
void EnsightIO::write_solution_ascii()
{
  system_vars_map_t::iterator sys_it = _system_vars_map.begin();
  const system_vars_map_t::iterator sys_end = _system_vars_map.end();

  for (; sys_it != sys_end; ++sys_it)
    {
      for (std::size_t i = 0; i < sys_it->second.EnsightScalars.size(); i++)
        this->write_scalar_ascii(sys_it->first,
                                 sys_it->second.EnsightScalars[i].scalar_name);

      for (std::size_t i = 0; i < sys_it->second.EnsightVectors.size(); i++)
        this->write_vector_ascii(sys_it->first,
                                 sys_it->second.EnsightVectors[i].components,
                                 sys_it->second.EnsightVectors[i].description);
    }
}


void EnsightIO::write_scalar_ascii(const std::string & sys,
                                   const std::string & var_name)
{
  // Construct scalar variable filename
  std::ostringstream scl_file;
  scl_file << _ensight_file_name
           << "_"
           << var_name
           << ".scl"
           << std::setw(3)
           << std::setprecision(0)
           << std::setfill('0')
           << std::right
           << _time_steps.size()-1;

  // Open a stream and start writing scalar variable info.
  std::ofstream scl_stream(scl_file.str().c_str());
  scl_stream << "Per node scalar value\n";
  scl_stream << "part\n";
  scl_stream << std::setw(10) << 1 << "\n";
  scl_stream << "coordinates\n";

  const MeshBase & the_mesh = MeshOutput<MeshBase>::mesh();
  const unsigned int dim = the_mesh.mesh_dimension();
  const System & system = _equation_systems.get_system(sys);
  const DofMap & dof_map = system.get_dof_map();
  int var = system.variable_number(var_name);

  std::vector<dof_id_type> dof_indices_scl;

  // Loop over active local elements, construct the nodal solution, and write it to file.
  MeshBase::const_element_iterator       el     = the_mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = the_mesh.active_local_elements_end();

  // Map from node id -> solution value.  We end up just writing this
  // map out in order, not sure what would happen if there were holes
  // in the numbering...
  typedef std::map<int, Real> map_local_soln;
  typedef map_local_soln::iterator local_soln_iterator;
  map_local_soln local_soln;

  std::vector<Number> elem_soln;
  std::vector<Number> nodal_soln;

  for ( ; el != end_el ; ++el)
    {
      const Elem * elem = *el;

      const FEType & fe_type = system.variable_type(var);

      dof_map.dof_indices (elem, dof_indices_scl, var);

      elem_soln.resize(dof_indices_scl.size());

      for (std::size_t i = 0; i < dof_indices_scl.size(); i++)
        elem_soln[i] = system.current_solution(dof_indices_scl[i]);

      FEInterface::nodal_soln (dim, fe_type, elem, elem_soln, nodal_soln);

      libmesh_assert_equal_to (nodal_soln.size(), elem->n_nodes());

#ifdef LIBMESH_USE_COMPLEX_NUMBERS
      libmesh_error_msg("Complex-valued Ensight output not yet supported");
#endif

      for (unsigned int n=0; n<elem->n_nodes(); n++)
        local_soln[elem->node_id(n)] = libmesh_real(nodal_soln[n]);
    }

  {
    local_soln_iterator it = local_soln.begin();
    const local_soln_iterator it_end = local_soln.end();
    for ( ; it != it_end; ++it)
      scl_stream << std::setw(12)
                 << std::setprecision(5)
                 << std::scientific
                 << it->second
                 << "\n";
  }
}


void EnsightIO::write_vector_ascii(const std::string & sys,
                                   const std::vector<std::string> & vec,
                                   const std::string & var_name)
{
  // Construct vector variable filename
  std::ostringstream vec_file;
  vec_file << _ensight_file_name
           << "_"
           << var_name
           << ".vec"
           << std::setw(3)
           << std::setprecision(0)
           << std::setfill('0')
           << std::right
           << _time_steps.size()-1;

  // Open a stream and start writing vector variable info.
  std::ofstream vec_stream(vec_file.str().c_str());
  vec_stream << "Per vector per value\n";
  vec_stream << "part\n";
  vec_stream << std::setw(10) << 1 << "\n";
  vec_stream << "coordinates\n";

  // Get a constant reference to the mesh object.
  const MeshBase & the_mesh = MeshOutput<MeshBase>::mesh();

  // The dimension that we are running
  const unsigned int dim = the_mesh.mesh_dimension();

  const System & system = _equation_systems.get_system(sys);

  const DofMap & dof_map = system.get_dof_map();

  const unsigned int u_var = system.variable_number(vec[0]);
  const unsigned int v_var = system.variable_number(vec[1]);
  const unsigned int w_var = (dim==3) ? system.variable_number(vec[2]) : 0;

  std::vector<dof_id_type> dof_indices_u;
  std::vector<dof_id_type> dof_indices_v;
  std::vector<dof_id_type> dof_indices_w;

  // Now we will loop over all the elements in the mesh.
  MeshBase::const_element_iterator       el     = the_mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = the_mesh.active_local_elements_end();

  // Map from node id -> solution value.  We end up just writing this
  // map out in order, not sure what would happen if there were holes
  // in the numbering...
  typedef std::map<int,std::vector<Real> > map_local_soln;
  typedef map_local_soln::iterator  local_soln_iterator;
  map_local_soln local_soln;

  for ( ; el != end_el ; ++el)
    {
      const Elem * elem = *el;

      const FEType & fe_type = system.variable_type(u_var);

      dof_map.dof_indices (elem, dof_indices_u, u_var);
      dof_map.dof_indices (elem, dof_indices_v, v_var);
      if (dim==3)
        dof_map.dof_indices (elem, dof_indices_w, w_var);

      std::vector<Number> elem_soln_u;
      std::vector<Number> elem_soln_v;
      std::vector<Number> elem_soln_w;

      std::vector<Number> nodal_soln_u;
      std::vector<Number> nodal_soln_v;
      std::vector<Number> nodal_soln_w;

      elem_soln_u.resize(dof_indices_u.size());
      elem_soln_v.resize(dof_indices_v.size());
      if (dim == 3)
        elem_soln_w.resize(dof_indices_w.size());

      for (std::size_t i = 0; i < dof_indices_u.size(); i++)
        {
          elem_soln_u[i] = system.current_solution(dof_indices_u[i]);
          elem_soln_v[i] = system.current_solution(dof_indices_v[i]);
          if (dim==3)
            elem_soln_w[i] = system.current_solution(dof_indices_w[i]);
        }

      FEInterface::nodal_soln (dim, fe_type, elem, elem_soln_u, nodal_soln_u);
      FEInterface::nodal_soln (dim, fe_type, elem, elem_soln_v, nodal_soln_v);
      if (dim == 3)
        FEInterface::nodal_soln (dim, fe_type, elem, elem_soln_w, nodal_soln_w);

      libmesh_assert_equal_to (nodal_soln_u.size(), elem->n_nodes());
      libmesh_assert_equal_to (nodal_soln_v.size(), elem->n_nodes());

#ifdef LIBMESH_ENABLE_COMPLEX
      libmesh_error_msg("Complex-valued Ensight output not yet supported");
#endif

      for (unsigned int n=0; n<elem->n_nodes(); n++)
        {
          std::vector<Real> node_vec(3);
          node_vec[0] = libmesh_real(nodal_soln_u[n]);
          node_vec[1] = libmesh_real(nodal_soln_v[n]);
          node_vec[2] = 0.0;
          if (dim==3)
            node_vec[2] = libmesh_real(nodal_soln_w[n]);
          local_soln[elem->node_id(n)] = node_vec;
        }
    }

  {
    local_soln_iterator it = local_soln.begin();
    const local_soln_iterator it_end = local_soln.end();

    for (unsigned dir=0; dir<3; ++dir)
      {
        for (; it != it_end; ++it)
          vec_stream << std::setw(12)
                     << std::scientific
                     << std::setprecision(5)
                     << it->second[dir]
                     << "\n";

        // Reset the iterator to the beginning of the map
        it = local_soln.begin();
      }
  }
}

} // namespace libMesh
