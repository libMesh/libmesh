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

// C++ includes
#include <iomanip>
#include <fstream>
#include <vector>
#include <ctype.h> // isspace

// Local includes
#include "libmesh/libmesh_config.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/gmv_io.h"
#include "libmesh/mesh_base.h"
#include "libmesh/elem.h"
#include "libmesh/equation_systems.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/string_to_enum.h"

// Wrap everything in a GMVLib namespace and
// use extern "C" to avoid name mangling.
#ifdef LIBMESH_HAVE_GMV
namespace GMVLib
{
extern "C"
{
#include "gmvread.h"
}
}
#endif

// anonymous namespace to hold local data
namespace
{
using namespace libMesh;

/**
 * Defines mapping from libMesh element types to GMV element types.
 * Note: Not all of the GMV element types have an identity mapping
 * to libmesh node numbering, but the node mappings do all happen to
 * be their own inverse, that is, pairs of nodes are simply swapped
 * between the two definitions.  Therefore we need only one node map
 * for both reading and writing.
 */
struct ElementDefinition {
  // GMV element name
  std::string label;

  // Used to map libmesh nodes to GMV for writing
  std::vector<unsigned> node_map;
};


// maps from a libMesh element type to the proper GMV
// ElementDefinition.  Placing the data structure here in this
// anonymous namespace gives us the benefits of a global variable
// without the nasty side-effects.
std::map<ElemType, ElementDefinition> eletypes;

// Helper function to fill up eletypes map
void add_eletype_entry(ElemType libmesh_elem_type,
                       const unsigned * node_map,
                       const std::string & gmv_label,
                       unsigned nodes_size )
{
  // If map entry does not exist, this will create it
  ElementDefinition & map_entry = eletypes[libmesh_elem_type];

  // Set the label
  map_entry.label = gmv_label;

  // Use the "swap trick" from Scott Meyer's "Effective STL" to swap
  // an unnamed temporary vector into the map_entry's vector.  Note:
  // the vector(iter, iter) constructor is used.
  std::vector<unsigned int>(node_map,
                            node_map+nodes_size).swap(map_entry.node_map);
}


// ------------------------------------------------------------
// helper function to initialize the eletypes map
void init_eletypes ()
{
  if (eletypes.empty())
    {
      // This should happen only once.  The first time this method
      // is called the eletypes data struture will be empty, and
      // we will fill it.  Any subsequent calls will find an initialized
      // eletypes map and will do nothing.

      // EDGE2
      {
        const unsigned int node_map[] = {0,1};
        add_eletype_entry(EDGE2, node_map, "line 2", 2);
      }

      // LINE3
      {
        const unsigned int node_map[] = {0,1,2};
        add_eletype_entry(EDGE3, node_map, "3line 3", 3);
      }

      // TRI3
      {
        const unsigned int node_map[] = {0,1,2};
        add_eletype_entry(TRI3, node_map, "tri3 3", 3);
      }

      // TRI6
      {
        const unsigned int node_map[] = {0,1,2,3,4,5};
        add_eletype_entry(TRI6, node_map, "6tri 6", 6);
      }

      // QUAD4
      {
        const unsigned int node_map[] = {0,1,2,3};
        add_eletype_entry(QUAD4, node_map, "quad 4", 4);
      }

      // QUAD8, QUAD9
      {
        const unsigned int node_map[] = {0,1,2,3,4,5,6,7};
        add_eletype_entry(QUAD8, node_map, "8quad 8", 8);

        // QUAD9 was not supported by GMV but it gets the same entry, even the label (is that correct?)
        eletypes[QUAD9] = eletypes[QUAD8];
      }

      // HEX8
      {
        const unsigned int node_map[] = {0,1,2,3,4,5,6,7};
        add_eletype_entry(HEX8, node_map, "phex8 8", 8);
      }

      // HEX20, HEX27
      {
        // Note: This map is its own inverse
        const unsigned int node_map[] = {0,1,2,3,4,5,6,7,8,9,10,11,16,17,18,19,12,13,14,15};
        add_eletype_entry(HEX20, node_map, "phex20 20", 20);

        // HEX27 was not supported by GMV but it gets the same entry, even the label (is that correct?)
        eletypes[HEX27] = eletypes[HEX20];
      }

      // TET4
      {
        // This is correct, see write_ascii_old_impl() to confirm.
        // This map is also its own inverse.
        const unsigned node_map[] = {0,2,1,3};
        add_eletype_entry(TET4, node_map, "tet 4", 4);
      }

      // TET10
      {
        const unsigned int node_map[] = {0,1,2,3,4,5,6,7,8,9};
        add_eletype_entry(TET10, node_map, "ptet10 10", 10);
      }

      // PRISM6
      {
        const unsigned int node_map[] = {0,1,2,3,4,5};
        add_eletype_entry(PRISM6, node_map, "pprism6 6", 6);
      }

      // PRISM15, PRISM18
      {
        // Note: This map is its own inverse
        const unsigned int node_map[] = {0,1,2,3,4,5,6,7,8,12,13,14, 9,10,11};
        add_eletype_entry(PRISM15, node_map, "pprism15 15", 15);

        // PRISM18 was not supported by GMV but it gets the same entry, even the label (is that correct?)
        eletypes[PRISM18] = eletypes[PRISM15];
      }
      //==============================
    }
}

} // end anonymous namespace


namespace libMesh
{

// Initialize the static data members by calling the static build functions.
std::map<std::string, ElemType> GMVIO::_reading_element_map = GMVIO::build_reading_element_map();



// Static function used to build the _reading_element_map.
std::map<std::string, ElemType> GMVIO::build_reading_element_map()
{
  std::map<std::string, ElemType> ret;

  ret["line"]     = EDGE2;
  ret["tri"]      = TRI3;
  ret["tri3"]     = TRI3;
  ret["quad"]     = QUAD4;
  ret["tet"]      = TET4;
  ret["ptet4"]    = TET4;
  ret["hex"]      = HEX8;
  ret["phex8"]    = HEX8;
  ret["prism"]    = PRISM6;
  ret["pprism6"]  = PRISM6;
  ret["phex20"]   = HEX20;
  ret["phex27"]   = HEX27;
  ret["pprism15"] = PRISM15;
  ret["ptet10"]   = TET10;
  ret["6tri"]     = TRI6;
  ret["8quad"]    = QUAD8;
  ret["3line"]    = EDGE3;

  // Unsupported/Unused types
  // ret["vface2d"]
  // ret["vface3d"]
  // ret["pyramid"]
  // ret["ppyrmd5"]
  // ret["ppyrmd13"]

  return ret;
}


GMVIO::GMVIO (const MeshBase & mesh) :
  MeshOutput<MeshBase>    (mesh),
  _binary                 (false),
  _discontinuous          (false),
  _partitioning           (true),
  _write_subdomain_id_as_material (false),
  _subdivide_second_order (true),
  _p_levels               (true),
  _next_elem_id           (0)
{
}



GMVIO::GMVIO (MeshBase & mesh) :
  MeshInput<MeshBase> (mesh),
  MeshOutput<MeshBase>(mesh),
  _binary (false),
  _discontinuous          (false),
  _partitioning           (true),
  _write_subdomain_id_as_material (false),
  _subdivide_second_order (true),
  _p_levels               (true),
  _next_elem_id           (0)
{
}



void GMVIO::write (const std::string & fname)
{
  if (this->binary())
    this->write_binary (fname);
  else
    this->write_ascii_old_impl  (fname);
}



void GMVIO::write_nodal_data (const std::string & fname,
                              const std::vector<Number> & soln,
                              const std::vector<std::string> & names)
{
  LOG_SCOPE("write_nodal_data()", "GMVIO");

  if (this->binary())
    this->write_binary (fname, &soln, &names);
  else
    this->write_ascii_old_impl  (fname, &soln, &names);
}



void GMVIO::write_ascii_new_impl (const std::string & fname,
                                  const std::vector<Number> * v,
                                  const std::vector<std::string> * solution_names)
{
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

  libMesh::err << "WARNING:  GMVIO::write_ascii_new_impl() not infinite-element aware!"
               << std::endl;
  libmesh_here();

  // Set it to our current precision
  this->write_ascii_old_impl (fname, v, solution_names);

#else

  // Get a reference to the mesh
  const MeshBase & mesh = MeshOutput<MeshBase>::mesh();

  // This is a parallel_only function
  const unsigned int n_active_elem = mesh.n_active_elem();

  if (MeshOutput<MeshBase>::mesh().processor_id() != 0)
    return;

  // Open the output file stream
  std::ofstream out_stream (fname.c_str());

  out_stream << std::setprecision(this->ascii_precision());

  // Make sure it opened correctly
  if (!out_stream.good())
    libmesh_file_error(fname.c_str());

  unsigned int mesh_max_p_level = 0;

  // Begin interfacing with the GMV data file
  {
    out_stream << "gmvinput ascii\n\n";

    // write the nodes
    out_stream << "nodes " << mesh.n_nodes() << "\n";
    for (unsigned int n=0; n<mesh.n_nodes(); n++)
      out_stream << mesh.point(n)(0) << " ";
    out_stream << "\n";

    for (unsigned int n=0; n<mesh.n_nodes(); n++)
#if LIBMESH_DIM > 1
      out_stream << mesh.point(n)(1) << " ";
#else
    out_stream << 0. << " ";
#endif
    out_stream << "\n";

    for (unsigned int n=0; n<mesh.n_nodes(); n++)
#if LIBMESH_DIM > 2
      out_stream << mesh.point(n)(2) << " ";
#else
    out_stream << 0. << " ";
#endif
    out_stream << "\n\n";
  }

  {
    // write the connectivity
    out_stream << "cells " << n_active_elem << "\n";

    // initialize the eletypes map (eletypes is a file-global variable)
    init_eletypes();

    MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
    const MeshBase::const_element_iterator end = mesh.active_elements_end();

    for ( ; it != end; ++it)
      {
        const Elem * elem = *it;

        mesh_max_p_level = std::max(mesh_max_p_level,
                                    elem->p_level());

        // Make sure we have a valid entry for
        // the current element type.
        libmesh_assert (eletypes.count(elem->type()));

        const ElementDefinition & ele = eletypes[elem->type()];

        // The element mapper better not require any more nodes
        // than are present in the current element!
        libmesh_assert_less_equal (ele.node_map.size(), elem->n_nodes());

        out_stream << ele.label << "\n";
        for (std::size_t i=0; i < ele.node_map.size(); i++)
          out_stream << elem->node_id(ele.node_map[i])+1 << " ";
        out_stream << "\n";
      }
    out_stream << "\n";
  }

  // optionally write the partition information
  if (this->partitioning())
    {
      if (this->write_subdomain_id_as_material())
        libmesh_error_msg("Not yet supported in GMVIO::write_ascii_new_impl");

      else // write processor IDs as materials.  This is the default
        {
          out_stream << "material "
                     << mesh.n_partitions()
            // Note: GMV may give you errors like
            // Error, material for cell 1 is greater than 1
            // Error, material for cell 2 is greater than 1
            // Error, material for cell 3 is greater than 1
            // ... because you put the wrong number of partitions here.
            // To ensure you write the correct number of materials, call
            // mesh.recalculate_n_partitions() before writing out the
            // mesh.
            // Note: we can't call it now because the Mesh is const here and
            // it is a non-const function.
                     << " 0\n";

          for (unsigned int proc=0; proc<mesh.n_partitions(); proc++)
            out_stream << "proc_" << proc << "\n";

          MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
          const MeshBase::const_element_iterator end = mesh.active_elements_end();

          // FIXME - don't we need to use an ElementDefinition here? - RHS
          for ( ; it != end; ++it)
            out_stream << (*it)->processor_id()+1 << "\n";
          out_stream << "\n";
        }
    }

  // If there are *any* variables at all in the system (including
  // p level, or arbitrary cell-based data)
  // to write, the gmv file needs to contain the word "variable"
  // on a line by itself.
  bool write_variable = false;

  // 1.) p-levels
  if (this->p_levels() && mesh_max_p_level)
    write_variable = true;

  // 2.) solution data
  if ((solution_names != libmesh_nullptr) && (v != libmesh_nullptr))
    write_variable = true;

  // 3.) cell-centered data
  if ( !(this->_cell_centered_data.empty()) )
    write_variable = true;

  if (write_variable)
    out_stream << "variable\n";

  //   if ((this->p_levels() && mesh_max_p_level) ||
  //     ((solution_names != libmesh_nullptr) && (v != libmesh_nullptr)))
  //     out_stream << "variable\n";

  // optionally write the polynomial degree information
  if (this->p_levels() && mesh_max_p_level)
    {
      out_stream << "p_level 0\n";

      MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
      const MeshBase::const_element_iterator end = mesh.active_elements_end();

      for ( ; it != end; ++it)
        {
          const Elem * elem = *it;

          const ElementDefinition & ele = eletypes[elem->type()];

          // The element mapper better not require any more nodes
          // than are present in the current element!
          libmesh_assert_less_equal (ele.node_map.size(), elem->n_nodes());

          for (std::size_t i=0; i < ele.node_map.size(); i++)
            out_stream << elem->p_level() << " ";
        }
      out_stream << "\n\n";
    }


  // optionally write cell-centered data
  if ( !(this->_cell_centered_data.empty()) )
    {
      std::map<std::string, const std::vector<Real> *>::iterator       it  = this->_cell_centered_data.begin();
      const std::map<std::string, const std::vector<Real> *>::iterator end = this->_cell_centered_data.end();

      for (; it != end; ++it)
        {
          // write out the variable name, followed by a zero.
          out_stream << (*it).first << " 0\n";

          const std::vector<Real> * the_array = (*it).second;

          // Loop over active elements, write out cell data.  If second-order cells
          // are split into sub-elements, the sub-elements inherit their parent's
          // cell-centered data.
          MeshBase::const_element_iterator       elem_it  = mesh.active_elements_begin();
          const MeshBase::const_element_iterator elem_end = mesh.active_elements_end();

          for (; elem_it != elem_end; ++elem_it)
            {
              const Elem * e = *elem_it;

              // Use the element's ID to find the value.
              libmesh_assert_less (e->id(), the_array->size());
              const Real the_value = the_array->operator[](e->id());

              if (this->subdivide_second_order())
                for (unsigned int se=0; se < e->n_sub_elem(); se++)
                  out_stream << the_value << " ";
              else
                out_stream << the_value << " ";
            }

          out_stream << "\n\n";
        }
    }


  // optionally write the data
  if ((solution_names != libmesh_nullptr) && (v != libmesh_nullptr))
    {
      const unsigned int n_vars = solution_names->size();

      if (!(v->size() == mesh.n_nodes()*n_vars))
        libMesh::err << "ERROR: v->size()=" << v->size()
                     << ", mesh.n_nodes()=" << mesh.n_nodes()
                     << ", n_vars=" << n_vars
                     << ", mesh.n_nodes()*n_vars=" << mesh.n_nodes()*n_vars
                     << "\n";

      libmesh_assert_equal_to (v->size(), mesh.n_nodes()*n_vars);

      for (unsigned int c=0; c<n_vars; c++)
        {

#ifdef LIBMESH_USE_COMPLEX_NUMBERS

          // in case of complex data, write _three_ data sets
          // for each component

          // this is the real part
          out_stream << "r_" << (*solution_names)[c] << " 1\n";

          for (unsigned int n=0; n<mesh.n_nodes(); n++)
            out_stream << (*v)[n*n_vars + c].real() << " ";

          out_stream << "\n\n";

          // this is the imaginary part
          out_stream << "i_" << (*solution_names)[c] << " 1\n";

          for (unsigned int n=0; n<mesh.n_nodes(); n++)
            out_stream << (*v)[n*n_vars + c].imag() << " ";

          out_stream << "\n\n";

          // this is the magnitude
          out_stream << "a_" << (*solution_names)[c] << " 1\n";
          for (unsigned int n=0; n<mesh.n_nodes(); n++)
            out_stream << std::abs((*v)[n*n_vars + c]) << " ";

          out_stream << "\n\n";

#else

          out_stream << (*solution_names)[c] << " 1\n";

          for (unsigned int n=0; n<mesh.n_nodes(); n++)
            out_stream << (*v)[n*n_vars + c] << " ";

          out_stream << "\n\n";

#endif
        }

    }

  // If we wrote any variables, we have to close the variable section now
  if (write_variable)
    out_stream << "endvars\n";


  // end of the file
  out_stream << "\nendgmv\n";

#endif
}






void GMVIO::write_ascii_old_impl (const std::string & fname,
                                  const std::vector<Number> * v,
                                  const std::vector<std::string> * solution_names)
{
  // Get a reference to the mesh
  const MeshBase & mesh = MeshOutput<MeshBase>::mesh();

  // Use a MeshSerializer object to gather a parallel mesh before outputting it.
  // Note that we cast away constness here (which is bad), but the destructor of
  // the MeshSerializer object reparallelizes the Mesh, hopefully keeping it
  // "logically const" outside the context of this function...
  MeshSerializer serialize(const_cast<MeshBase &>(mesh),
                           !MeshOutput<MeshBase>::_is_parallel_format);

  // These are parallel_only functions
  const dof_id_type n_active_elem = mesh.n_active_elem(),
    n_active_sub_elem = mesh.n_active_sub_elem();

  if (MeshOutput<MeshBase>::mesh().processor_id() != 0)
    return;

  // Open the output file stream
  std::ofstream out_stream (fname.c_str());

  // Set it to our current precision
  out_stream << std::setprecision(this->ascii_precision());

  // Make sure it opened correctly
  if (!out_stream.good())
    libmesh_file_error(fname.c_str());

  // Make sure our nodes are contiguous and serialized
  libmesh_assert_equal_to (mesh.n_nodes(), mesh.max_node_id());

  // libmesh_assert (mesh.is_serial());
  // if (!mesh.is_serial())
  //   {
  //     if (MeshOutput<MeshBase>::mesh().processor_id() == 0)
  //       libMesh::err << "Error: GMVIO cannot yet write a DistributedMesh solution"
  //                     << std::endl;
  //     return;
  //   }

  unsigned int mesh_max_p_level = 0;

  // Begin interfacing with the GMV data file

  // FIXME - if subdivide_second_order() is off,
  // we probably should only be writing the
  // vertex nodes - RHS
  {
    // write the nodes

    out_stream << "gmvinput ascii\n\n";
    out_stream << "nodes " << mesh.n_nodes() << '\n';
    for (unsigned int n=0; n<mesh.n_nodes(); n++)
      out_stream << mesh.point(n)(0) << " ";

    out_stream << '\n';

    for (unsigned int n=0; n<mesh.n_nodes(); n++)
#if LIBMESH_DIM > 1
      out_stream << mesh.point(n)(1) << " ";
#else
    out_stream << 0. << " ";
#endif

    out_stream << '\n';

    for (unsigned int n=0; n<mesh.n_nodes(); n++)
#if LIBMESH_DIM > 2
      out_stream << mesh.point(n)(2) << " ";
#else
    out_stream << 0. << " ";
#endif

    out_stream << '\n' << '\n';
  }



  {
    // write the connectivity

    out_stream << "cells ";
    if (this->subdivide_second_order())
      out_stream << n_active_sub_elem;
    else
      out_stream << n_active_elem;
    out_stream << '\n';

    // The same temporary storage will be used for each element
    std::vector<dof_id_type> conn;

    MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
    const MeshBase::const_element_iterator end = mesh.active_elements_end();

    for ( ; it != end; ++it)
      {
        const Elem * elem = *it;

        mesh_max_p_level = std::max(mesh_max_p_level,
                                    elem->p_level());

        switch (elem->dim())
          {
          case 1:
            {
              if (this->subdivide_second_order())
                for (unsigned int se=0; se<(*it)->n_sub_elem(); se++)
                  {
                    out_stream << "line 2\n";
                    (*it)->connectivity(se, TECPLOT, conn);
                    for (std::size_t i=0; i<conn.size(); i++)
                      out_stream << conn[i] << " ";

                    out_stream << '\n';
                  }
              else
                {
                  out_stream << "line 2\n";
                  if ((*it)->default_order() == FIRST)
                    (*it)->connectivity(0, TECPLOT, conn);
                  else
                    {
                      UniquePtr<Elem> lo_elem = Elem::build(Elem::first_order_equivalent_type((*it)->type()));
                      for (unsigned int i = 0; i != lo_elem->n_nodes(); ++i)
                        lo_elem->set_node(i) = (*it)->node_ptr(i);
                      lo_elem->connectivity(0, TECPLOT, conn);
                    }
                  for (std::size_t i=0; i<conn.size(); i++)
                    out_stream << conn[i] << " ";

                  out_stream << '\n';
                }

              break;
            }

          case 2:
            {
              if (this->subdivide_second_order())
                for (unsigned int se=0; se<(*it)->n_sub_elem(); se++)
                  {
                    // Quad elements
                    if (((*it)->type() == QUAD4) ||
                        ((*it)->type() == QUAD8) || // Note: QUAD8 will be output as one central quad and
                        // four surrounding triangles (though they will be written
                        // to GMV as QUAD4s).
                        ((*it)->type() == QUAD9)
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
                        || ((*it)->type() == INFQUAD4)
                        || ((*it)->type() == INFQUAD6)
#endif
                        )
                      {
                        out_stream << "quad 4\n";
                        (*it)->connectivity(se, TECPLOT, conn);
                        for (std::size_t i=0; i<conn.size(); i++)
                          out_stream << conn[i] << " ";
                      }

                    // Triangle elements
                    else if (((*it)->type() == TRI3) ||
                             ((*it)->type() == TRI6))
                      {
                        out_stream << "tri 3\n";
                        (*it)->connectivity(se, TECPLOT, conn);
                        for (unsigned int i=0; i<3; i++)
                          out_stream << conn[i] << " ";
                      }
                    else
                      libmesh_error_msg("Unsupported element type: " << Utility::enum_to_string((*it)->type()));
                  }
              else // !this->subdivide_second_order()
                {
                  // Quad elements
                  if (((*it)->type() == QUAD4)
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
                      || ((*it)->type() == INFQUAD4)
#endif
                      )
                    {
                      (*it)->connectivity(0, TECPLOT, conn);
                      out_stream << "quad 4\n";
                      for (std::size_t i=0; i<conn.size(); i++)
                        out_stream << conn[i] << " ";
                    }
                  else if (((*it)->type() == QUAD8) ||
                           ((*it)->type() == QUAD9)
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
                           || ((*it)->type() == INFQUAD6)
#endif
                           )
                    {
                      UniquePtr<Elem> lo_elem = Elem::build(Elem::first_order_equivalent_type((*it)->type()));
                      for (unsigned int i = 0; i != lo_elem->n_nodes(); ++i)
                        lo_elem->set_node(i) = (*it)->node_ptr(i);
                      lo_elem->connectivity(0, TECPLOT, conn);
                      out_stream << "quad 4\n";
                      for (std::size_t i=0; i<conn.size(); i++)
                        out_stream << conn[i] << " ";
                    }
                  else if ((*it)->type() == TRI3)
                    {
                      (*it)->connectivity(0, TECPLOT, conn);
                      out_stream << "tri 3\n";
                      for (unsigned int i=0; i<3; i++)
                        out_stream << conn[i] << " ";
                    }
                  else if ((*it)->type() == TRI6)
                    {
                      UniquePtr<Elem> lo_elem = Elem::build(Elem::first_order_equivalent_type((*it)->type()));
                      for (unsigned int i = 0; i != lo_elem->n_nodes(); ++i)
                        lo_elem->set_node(i) = (*it)->node_ptr(i);
                      lo_elem->connectivity(0, TECPLOT, conn);
                      out_stream << "tri 3\n";
                      for (unsigned int i=0; i<3; i++)
                        out_stream << conn[i] << " ";
                    }

                  out_stream << '\n';
                }

              break;
            }

          case 3:
            {
              if (this->subdivide_second_order())
                for (unsigned int se=0; se<(*it)->n_sub_elem(); se++)
                  {

#ifndef  LIBMESH_ENABLE_INFINITE_ELEMENTS
                    if (((*it)->type() == HEX8)   ||
                        ((*it)->type() == HEX27))
                      {
                        out_stream << "phex8 8\n";
                        (*it)->connectivity(se, TECPLOT, conn);
                        for (std::size_t i=0; i<conn.size(); i++)
                          out_stream << conn[i] << " ";
                      }

                    else if ((*it)->type() == HEX20)
                      {
                        out_stream << "phex20 20\n";
                        out_stream << (*it)->node_id(0)+1  << " "
                                   << (*it)->node_id(1)+1  << " "
                                   << (*it)->node_id(2)+1  << " "
                                   << (*it)->node_id(3)+1  << " "
                                   << (*it)->node_id(4)+1  << " "
                                   << (*it)->node_id(5)+1  << " "
                                   << (*it)->node_id(6)+1  << " "
                                   << (*it)->node_id(7)+1  << " "
                                   << (*it)->node_id(8)+1  << " "
                                   << (*it)->node_id(9)+1  << " "
                                   << (*it)->node_id(10)+1 << " "
                                   << (*it)->node_id(11)+1 << " "
                                   << (*it)->node_id(16)+1 << " "
                                   << (*it)->node_id(17)+1 << " "
                                   << (*it)->node_id(18)+1 << " "
                                   << (*it)->node_id(19)+1 << " "
                                   << (*it)->node_id(12)+1 << " "
                                   << (*it)->node_id(13)+1 << " "
                                   << (*it)->node_id(14)+1 << " "
                                   << (*it)->node_id(15)+1 << " ";
                      }

                    // According to our copy of gmvread.c, this is the
                    // mapping for the Hex27 element.  Unfortunately,
                    // I tried it and Paraview does not seem to be
                    // able to read Hex27 elements.  Since this is
                    // unlikely to change any time soon, we'll
                    // continue to write out Hex27 elements as 8 Hex8
                    // sub-elements.
                    //
                    // TODO:
                    // 1.) If we really wanted to use this code for
                    // something, we'd want to avoid repeating the
                    // hard-coded node ordering from the HEX20 case.
                    // These should both be able to use
                    // ElementDefinitions.
                    // 2.) You would also need to change
                    // Hex27::n_sub_elem() to return 1, not sure how
                    // much other code that would affect...

                    // else if ((*it)->type() == HEX27)
                    //   {
                    //     out_stream << "phex27 27\n";
                    //     out_stream << (*it)->node_id(0)+1  << " "
                    //                << (*it)->node_id(1)+1  << " "
                    //                << (*it)->node_id(2)+1  << " "
                    //                << (*it)->node_id(3)+1  << " "
                    //                << (*it)->node_id(4)+1  << " "
                    //                << (*it)->node_id(5)+1  << " "
                    //                << (*it)->node_id(6)+1  << " "
                    //                << (*it)->node_id(7)+1  << " "
                    //                << (*it)->node_id(8)+1  << " "
                    //                << (*it)->node_id(9)+1  << " "
                    //                << (*it)->node_id(10)+1 << " "
                    //                << (*it)->node_id(11)+1 << " "
                    //                << (*it)->node_id(16)+1 << " "
                    //                << (*it)->node_id(17)+1 << " "
                    //                << (*it)->node_id(18)+1 << " "
                    //                << (*it)->node_id(19)+1 << " "
                    //                << (*it)->node_id(12)+1 << " "
                    //                << (*it)->node_id(13)+1 << " "
                    //                << (*it)->node_id(14)+1 << " "
                    //                << (*it)->node_id(15)+1 << " "
                    //       // mid-face nodes
                    //                << (*it)->node_id(21)+1 << " " // GMV front
                    //                << (*it)->node_id(22)+1 << " " // GMV right
                    //                << (*it)->node_id(23)+1 << " " // GMV back
                    //                << (*it)->node_id(24)+1 << " " // GMV left
                    //                << (*it)->node_id(20)+1 << " " // GMV bottom
                    //                << (*it)->node_id(25)+1 << " " // GMV top
                    //       // center node
                    //                << (*it)->node_id(26)+1 << " ";
                    //   }

#else // LIBMESH_ENABLE_INFINITE_ELEMENTS

                    // In case of infinite elements, HEX20
                    // should be handled just like the
                    // INFHEX16, since these connect to each other
                    if (((*it)->type() == HEX8)     ||
                        ((*it)->type() == HEX27)    ||
                        ((*it)->type() == INFHEX8)  ||
                        ((*it)->type() == INFHEX16) ||
                        ((*it)->type() == INFHEX18) ||
                        ((*it)->type() == HEX20))
                      {
                        out_stream << "phex8 8\n";
                        (*it)->connectivity(se, TECPLOT, conn);
                        for (std::size_t i=0; i<conn.size(); i++)
                          out_stream << conn[i] << " ";
                      }
#endif

                    else if (((*it)->type() == TET4)  ||
                             ((*it)->type() == TET10))
                      {
                        out_stream << "tet 4\n";
                        // Tecplot connectivity returns 8 entries for
                        // the Tet, enough to store it as a degenerate Hex.
                        // For GMV we only pick out the four relevant node
                        // indices.
                        (*it)->connectivity(se, TECPLOT, conn);
                        out_stream << conn[0] << " "  // libmesh tet node 0
                                   << conn[2] << " "  // libmesh tet node 2
                                   << conn[1] << " "  // libmesh tet node 1
                                   << conn[4] << " "; // libmesh tet node 3
                      }
#ifndef  LIBMESH_ENABLE_INFINITE_ELEMENTS
                    else if (((*it)->type() == PRISM6)  ||
                             ((*it)->type() == PRISM15) ||
                             ((*it)->type() == PRISM18) ||
                             ((*it)->type() == PYRAMID5))
#else
                    else if (((*it)->type() == PRISM6)     ||
                             ((*it)->type() == PRISM15)    ||
                             ((*it)->type() == PRISM18)    ||
                             ((*it)->type() == PYRAMID5)   ||
                             ((*it)->type() == INFPRISM6)  ||
                             ((*it)->type() == INFPRISM12))
#endif
                      {
                        // Note that the prisms are treated as
                        // degenerated phex8's.
                        out_stream << "phex8 8\n";
                        (*it)->connectivity(se, TECPLOT, conn);
                        for (std::size_t i=0; i<conn.size(); i++)
                          out_stream << conn[i] << " ";
                      }

                    else
                      libmesh_error_msg("Encountered an unrecognized element " \
                                        << "type: " << (*it)->type()  \
                                        << "\nPossibly a dim-1 dimensional " \
                                        << "element?  Aborting...");

                    out_stream << '\n';
                  }
              else // !this->subdivide_second_order()
                {
                  UniquePtr<Elem> lo_elem = Elem::build(Elem::first_order_equivalent_type((*it)->type()));
                  for (unsigned int i = 0; i != lo_elem->n_nodes(); ++i)
                    lo_elem->set_node(i) = (*it)->node_ptr(i);
                  if ((lo_elem->type() == HEX8)
#ifdef  LIBMESH_ENABLE_INFINITE_ELEMENTS
                      || (lo_elem->type() == HEX27)
#endif
                      )
                    {
                      out_stream << "phex8 8\n";
                      lo_elem->connectivity(0, TECPLOT, conn);
                      for (std::size_t i=0; i<conn.size(); i++)
                        out_stream << conn[i] << " ";
                    }

                  else if (lo_elem->type() == TET4)
                    {
                      out_stream << "tet 4\n";
                      lo_elem->connectivity(0, TECPLOT, conn);
                      out_stream << conn[0] << " "
                                 << conn[2] << " "
                                 << conn[1] << " "
                                 << conn[4] << " ";
                    }
                  else if ((lo_elem->type() == PRISM6)
#ifdef  LIBMESH_ENABLE_INFINITE_ELEMENTS
                           || (lo_elem->type() == INFPRISM6)
#endif
                           )
                    {
                      // Note that the prisms are treated as
                      // degenerated phex8's.
                      out_stream << "phex8 8\n";
                      lo_elem->connectivity(0, TECPLOT, conn);
                      for (std::size_t i=0; i<conn.size(); i++)
                        out_stream << conn[i] << " ";
                    }

                  else
                    libmesh_error_msg("Encountered an unrecognized element " \
                                      << "type.  Possibly a dim-1 dimensional " \
                                      << "element?  Aborting...");

                  out_stream << '\n';
                }

              break;
            }

          default:
            libmesh_error_msg("Unsupported element dimension: " <<
                              elem->dim());
          }
      }

    out_stream << '\n';
  }



  // optionally write the partition information
  if (this->partitioning())
    {
      if (this->write_subdomain_id_as_material())
        {
          // Subdomain IDs can be non-contiguous and need not
          // necessarily start at 0.  Furthermore, since the user is
          // free to define subdomain_id_type to be a signed type, we
          // can't even assume max(subdomain_id) >= # unique subdomain ids.

          // We build a map<subdomain_id, unsigned> to associate to each
          // user-selected subdomain ID a unique, contiguous unsigned value
          // which we can write to file.
          std::map<subdomain_id_type, unsigned> sbdid_map;
          typedef std::map<subdomain_id_type, unsigned>::iterator sbdid_map_iter;
          {
            MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
            const MeshBase::const_element_iterator end = mesh.active_elements_end();

            for ( ; it != end; ++it)
              {
                // Try to insert with dummy value
                sbdid_map.insert( std::make_pair((*it)->subdomain_id(), 0) );
              }
          }

          // Map is created, iterate through it to set indices.  They will be
          // used repeatedly below.
          {
            unsigned ctr=0;
            for (sbdid_map_iter it=sbdid_map.begin(); it != sbdid_map.end(); ++it)
              (*it).second = ctr++;
          }

          out_stream << "material "
                     << sbdid_map.size()
                     << " 0\n";

          for (std::size_t sbdid=0; sbdid<sbdid_map.size(); sbdid++)
            out_stream << "proc_" << sbdid << "\n";

          MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
          const MeshBase::const_element_iterator end = mesh.active_elements_end();

          for ( ; it != end; ++it)
            {
              // Find the unique index for (*it)->subdomain_id(), print that to file
              sbdid_map_iter map_iter = sbdid_map.find( (*it)->subdomain_id() );
              unsigned gmv_mat_number = (*map_iter).second;

              if (this->subdivide_second_order())
                for (unsigned int se=0; se<(*it)->n_sub_elem(); se++)
                  out_stream << gmv_mat_number+1 << '\n';
              else
                out_stream << gmv_mat_number+1 << "\n";
            }
          out_stream << '\n';

        }
      else // write processor IDs as materials.  This is the default
        {
          out_stream << "material "
                     << mesh.n_partitions()
                     << " 0"<< '\n';

          for (unsigned int proc=0; proc<mesh.n_partitions(); proc++)
            out_stream << "proc_" << proc << '\n';

          MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
          const MeshBase::const_element_iterator end = mesh.active_elements_end();

          for ( ; it != end; ++it)
            if (this->subdivide_second_order())
              for (unsigned int se=0; se<(*it)->n_sub_elem(); se++)
                out_stream << (*it)->processor_id()+1 << '\n';
            else
              out_stream << (*it)->processor_id()+1 << '\n';

          out_stream << '\n';
        }
    }


  // If there are *any* variables at all in the system (including
  // p level, or arbitrary cell-based data)
  // to write, the gmv file needs to contain the word "variable"
  // on a line by itself.
  bool write_variable = false;

  // 1.) p-levels
  if (this->p_levels() && mesh_max_p_level)
    write_variable = true;

  // 2.) solution data
  if ((solution_names != libmesh_nullptr) && (v != libmesh_nullptr))
    write_variable = true;

  // 3.) cell-centered data
  if ( !(this->_cell_centered_data.empty()) )
    write_variable = true;

  if (write_variable)
    out_stream << "variable\n";


  // optionally write the p-level information
  if (this->p_levels() && mesh_max_p_level)
    {
      out_stream << "p_level 0\n";

      MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
      const MeshBase::const_element_iterator end = mesh.active_elements_end();

      for ( ; it != end; ++it)
        if (this->subdivide_second_order())
          for (unsigned int se=0; se<(*it)->n_sub_elem(); se++)
            out_stream << (*it)->p_level() << " ";
        else
          out_stream << (*it)->p_level() << " ";
      out_stream << "\n\n";
    }




  // optionally write cell-centered data
  if ( !(this->_cell_centered_data.empty()) )
    {
      std::map<std::string, const std::vector<Real> *>::iterator       it  = this->_cell_centered_data.begin();
      const std::map<std::string, const std::vector<Real> *>::iterator end = this->_cell_centered_data.end();

      for (; it != end; ++it)
        {
          // write out the variable name, followed by a zero.
          out_stream << (*it).first << " 0\n";

          const std::vector<Real> * the_array = (*it).second;

          // Loop over active elements, write out cell data.  If second-order cells
          // are split into sub-elements, the sub-elements inherit their parent's
          // cell-centered data.
          MeshBase::const_element_iterator       elem_it  = mesh.active_elements_begin();
          const MeshBase::const_element_iterator elem_end = mesh.active_elements_end();

          for (; elem_it != elem_end; ++elem_it)
            {
              const Elem * e = *elem_it;

              // Use the element's ID to find the value...
              libmesh_assert_less (e->id(), the_array->size());
              const Real the_value = the_array->operator[](e->id());

              if (this->subdivide_second_order())
                for (unsigned int se=0; se < e->n_sub_elem(); se++)
                  out_stream << the_value << " ";
              else
                out_stream << the_value << " ";
            }

          out_stream << "\n\n";
        }
    }




  // optionally write the data
  if ((solution_names != libmesh_nullptr) &&
      (v != libmesh_nullptr))
    {
      const unsigned int n_vars =
        cast_int<unsigned int>(solution_names->size());

      if (!(v->size() == mesh.n_nodes()*n_vars))
        libMesh::err << "ERROR: v->size()=" << v->size()
                     << ", mesh.n_nodes()=" << mesh.n_nodes()
                     << ", n_vars=" << n_vars
                     << ", mesh.n_nodes()*n_vars=" << mesh.n_nodes()*n_vars
                     << std::endl;

      libmesh_assert_equal_to (v->size(), mesh.n_nodes()*n_vars);

      for (unsigned int c=0; c<n_vars; c++)
        {

#ifdef LIBMESH_USE_COMPLEX_NUMBERS

          // in case of complex data, write _tree_ data sets
          // for each component

          // this is the real part
          out_stream << "r_" << (*solution_names)[c] << " 1\n";

          for (unsigned int n=0; n<mesh.n_nodes(); n++)
            out_stream << (*v)[n*n_vars + c].real() << " ";

          out_stream << '\n' << '\n';


          // this is the imaginary part
          out_stream << "i_" << (*solution_names)[c] << " 1\n";

          for (unsigned int n=0; n<mesh.n_nodes(); n++)
            out_stream << (*v)[n*n_vars + c].imag() << " ";

          out_stream << '\n' << '\n';

          // this is the magnitude
          out_stream << "a_" << (*solution_names)[c] << " 1\n";
          for (unsigned int n=0; n<mesh.n_nodes(); n++)
            out_stream << std::abs((*v)[n*n_vars + c]) << " ";

          out_stream << '\n' << '\n';

#else

          out_stream << (*solution_names)[c] << " 1\n";

          for (unsigned int n=0; n<mesh.n_nodes(); n++)
            out_stream << (*v)[n*n_vars + c] << " ";

          out_stream << '\n' << '\n';

#endif
        }

    }

  // If we wrote any variables, we have to close the variable section now
  if (write_variable)
    out_stream << "endvars\n";


  // end of the file
  out_stream << "\nendgmv\n";
}







void GMVIO::write_binary (const std::string & fname,
                          const std::vector<Number> * vec,
                          const std::vector<std::string> * solution_names)
{
  // We use the file-scope global variable eletypes for mapping nodes
  // from GMV to libmesh indices, so initialize that data now.
  init_eletypes();

  // Get a reference to the mesh
  const MeshBase & mesh = MeshOutput<MeshBase>::mesh();

  // This is a parallel_only function
  const dof_id_type n_active_elem = mesh.n_active_elem();

  if (MeshOutput<MeshBase>::mesh().processor_id() != 0)
    return;

  std::ofstream out_stream (fname.c_str());

  libmesh_assert (out_stream.good());

  unsigned int mesh_max_p_level = 0;

  std::string buffer;

  // Begin interfacing with the GMV data file
  {
    // write the nodes
    buffer = "gmvinput";
    out_stream.write(buffer.c_str(), buffer.size());

    buffer = "ieeei4r4";
    out_stream.write(buffer.c_str(), buffer.size());
  }



  // write the nodes
  {
    buffer = "nodes   ";
    out_stream.write(buffer.c_str(), buffer.size());

    unsigned int tempint = mesh.n_nodes();
    out_stream.write(reinterpret_cast<char *>(&tempint), sizeof(unsigned int));

    // write the x coordinate
    std::vector<float> temp(mesh.n_nodes());
    for (unsigned int v=0; v<mesh.n_nodes(); v++)
      temp[v] = static_cast<float>(mesh.point(v)(0));
    out_stream.write(reinterpret_cast<char *>(&temp[0]), sizeof(float)*mesh.n_nodes());

    // write the y coordinate
    for (unsigned int v=0; v<mesh.n_nodes(); v++)
      {
#if LIBMESH_DIM > 1
      temp[v] = static_cast<float>(mesh.point(v)(1));
#else
      temp[v] = 0.;
#endif
      }
    out_stream.write(reinterpret_cast<char *>(&temp[0]), sizeof(float)*mesh.n_nodes());

    // write the z coordinate
    for (unsigned int v=0; v<mesh.n_nodes(); v++)
      {
#if LIBMESH_DIM > 2
      temp[v] = static_cast<float>(mesh.point(v)(2));
#else
      temp[v] = 0.;
#endif
      }
    out_stream.write(reinterpret_cast<char *>(&temp[0]), sizeof(float)*mesh.n_nodes());
  }


  // write the connectivity
  {
    buffer = "cells   ";
    out_stream.write(buffer.c_str(), buffer.size());

    unsigned int tempint = n_active_elem;
    out_stream.write(reinterpret_cast<char *>(&tempint), sizeof(unsigned int));

    MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
    const MeshBase::const_element_iterator end = mesh.active_elements_end();

    for ( ; it != end; ++it)
      {
        const Elem * elem = *it;

        mesh_max_p_level = std::max(mesh_max_p_level,
                                    elem->p_level());

        // The ElementDefinition label contains the GMV name
        // and the number of nodes for the respective
        // element.
        const ElementDefinition & ed = eletypes[elem->type()];

        // The ElementDefinition labels all have a name followed by a
        // whitespace and then the number of nodes.  To write the
        // string in binary, we want to drop everything after that
        // space, and then pad the string out to 8 characters.
        buffer = ed.label;

        // Erase everything from the first whitespace character to the end of the string.
        buffer.erase(buffer.find_first_of(' '), std::string::npos);

        // Append whitespaces until the string is exactly 8 characters long.
        while (buffer.size() < 8)
          buffer.insert(buffer.end(), ' ');

        // Debugging:
        // libMesh::out << "Writing element with name = '" << buffer << "'" << std::endl;

        // Finally, write the 8 character stream to file.
        out_stream.write(buffer.c_str(), buffer.size());

        // Write the number of nodes that this type of element has.
        // Note: don't want to use the number of nodes of the element,
        // because certain elements, like QUAD9 and HEX27 only support
        // being written out as lower-order elements (QUAD8 and HEX20,
        // respectively).
        tempint = ed.node_map.size();
        out_stream.write(reinterpret_cast<char *>(&tempint), sizeof(unsigned int));

        // Write the element connectivity
        for (std::size_t i=0; i<ed.node_map.size(); i++)
          {
            dof_id_type id = elem->node_id(ed.node_map[i]) + 1;
            out_stream.write(reinterpret_cast<char *>(&id), sizeof(dof_id_type));
          }
      }
  }



  // optionally write the partition information
  if (this->partitioning())
    {
      if (this->write_subdomain_id_as_material())
        libmesh_error_msg("Not yet supported in GMVIO::write_binary");

      else
        {
          buffer = "material";
          out_stream.write(buffer.c_str(), buffer.size());

          unsigned int tmpint = mesh.n_processors();
          out_stream.write(reinterpret_cast<char *>(&tmpint), sizeof(unsigned int));

          tmpint = 0; // IDs are cell based
          out_stream.write(reinterpret_cast<char *>(&tmpint), sizeof(unsigned int));

          for (unsigned int proc=0; proc<mesh.n_processors(); proc++)
            {
              // We have to write exactly 8 bytes.  This means that
              // the processor id can be up to 3 characters long, and
              // we pad it with blank characters on the end.
              std::ostringstream oss;
              oss << "proc_" << std::setw(3) << std::left << proc;
              out_stream.write(oss.str().c_str(), oss.str().size());
            }

          std::vector<unsigned int> proc_id (n_active_elem);

          unsigned int n=0;

          MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
          const MeshBase::const_element_iterator end = mesh.active_elements_end();

          for ( ; it != end; ++it)
            {
              const Elem * elem = *it;

              // We no longer write sub-elems in GMV, so just assign a
              // processor id value to each element.
              proc_id[n++] = elem->processor_id() + 1;
            }


          out_stream.write(reinterpret_cast<char *>(&proc_id[0]),
                           sizeof(unsigned int)*proc_id.size());
        }
    }

  // If there are *any* variables at all in the system (including
  // p level, or arbitrary cell-based data)
  // to write, the gmv file needs to contain the word "variable"
  // on a line by itself.
  bool write_variable = false;

  // 1.) p-levels
  if (this->p_levels() && mesh_max_p_level)
    write_variable = true;

  // 2.) solution data
  if ((solution_names != libmesh_nullptr) && (vec != libmesh_nullptr))
    write_variable = true;

  //   // 3.) cell-centered data - unsupported
  //   if ( !(this->_cell_centered_data.empty()) )
  //     write_variable = true;

  if (write_variable)
    {
      buffer = "variable";
      out_stream.write(buffer.c_str(), buffer.size());
    }

  // optionally write the partition information
  if (this->p_levels() && mesh_max_p_level)
    {
      unsigned int n_floats = n_active_elem;
      for (unsigned int i=0; i != mesh.mesh_dimension(); ++i)
        n_floats *= 2;

      std::vector<float> temp(n_floats);

      buffer = "p_level";
      out_stream.write(buffer.c_str(), buffer.size());

      unsigned int tempint = 0; // p levels are cell data
      out_stream.write(reinterpret_cast<char *>(&tempint), sizeof(unsigned int));

      MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
      const MeshBase::const_element_iterator end = mesh.active_elements_end();
      unsigned int n=0;

      for (; it != end; ++it)
        for (unsigned int se=0; se<(*it)->n_sub_elem(); se++)
          temp[n++] = static_cast<float>( (*it)->p_level() );

      out_stream.write(reinterpret_cast<char *>(&temp[0]),
                       sizeof(float)*n_floats);
    }


  // optionally write cell-centered data
  if ( !(this->_cell_centered_data.empty()) )
    {
      libMesh::err << "Cell-centered data not (yet) supported in binary I/O mode!" << std::endl;

      //        std::map<std::string, const std::vector<Real> * >::iterator       it  = this->_cell_centered_data.begin();
      //        const std::map<std::string, const std::vector<Real> * >::iterator end = this->_cell_centered_data.end();

      //        for (; it != end; ++it)
      //  {
      //    // Write out the variable name ...
      //    out_stream.write(it->first.c_str(), it->first.size());

      //    // ... followed by a zero.
      //    unsigned int tempint = 0; // 0 signifies cell data
      //    out_stream.write(reinterpret_cast<char *>(&tempint), sizeof(unsigned int));

      //    // Get a pointer to the array of cell-centered data values
      //    const std::vector<Real> * the_array = (*it).second;

      //   // Since the_array might contain zeros (for inactive elements) we need to
      //   // make a copy of it containing just values for active elements.
      //   const unsigned int n_floats = n_active_elem * (1<<mesh.mesh_dimension());
      //   std::vector<float> temp(n_floats);

      //   MeshBase::const_element_iterator       elem_it  = mesh.active_elements_begin();
      //   const MeshBase::const_element_iterator elem_end = mesh.active_elements_end();
      //   unsigned int n=0;

      //   for (; elem_it != elem_end; ++elem_it)
      //     {
      //       // If there's a seg-fault, it will probably be here!
      //       const float the_value = static_cast<float>(the_array->operator[]((*elem_it)->id()));

      //       for (unsigned int se=0; se<(*elem_it)->n_sub_elem(); se++)
      // temp[n++] = the_value;
      //     }


      //    // Write "the_array" directly to the file
      //    out_stream.write(reinterpret_cast<char *>(&temp[0]),
      //      sizeof(float)*n_floats);
      //  }
    }




  // optionally write the data
  if ((solution_names != libmesh_nullptr) &&
      (vec != libmesh_nullptr))
    {
      std::vector<float> temp(mesh.n_nodes());

      const unsigned int n_vars =
        cast_int<unsigned int>(solution_names->size());

      for (unsigned int c=0; c<n_vars; c++)
        {

#ifdef LIBMESH_USE_COMPLEX_NUMBERS
          // for complex data, write three datasets


          // Real part
          buffer = "r_";
          out_stream.write(buffer.c_str(), buffer.size());

          buffer = (*solution_names)[c];
          out_stream.write(buffer.c_str(), buffer.size());

          unsigned int tempint = 1; // always do nodal data
          out_stream.write(reinterpret_cast<char *>(&tempint), sizeof(unsigned int));

          for (unsigned int n=0; n<mesh.n_nodes(); n++)
            temp[n] = static_cast<float>( (*vec)[n*n_vars + c].real() );

          out_stream.write(reinterpret_cast<char *>(&temp[0]), sizeof(float)*mesh.n_nodes());


          // imaginary part
          buffer = "i_";
          out_stream.write(buffer.c_str(), buffer.size());

          buffer = (*solution_names)[c];
          out_stream.write(buffer.c_str(), buffer.size());

          out_stream.write(reinterpret_cast<char *>(&tempint), sizeof(unsigned int));

          for (unsigned int n=0; n<mesh.n_nodes(); n++)
            temp[n] = static_cast<float>( (*vec)[n*n_vars + c].imag() );

          out_stream.write(reinterpret_cast<char *>(&temp[0]), sizeof(float)*mesh.n_nodes());

          // magnitude
          buffer = "a_";
          out_stream.write(buffer.c_str(), buffer.size());
          buffer = (*solution_names)[c];
          out_stream.write(buffer.c_str(), buffer.size());

          out_stream.write(reinterpret_cast<char *>(&tempint), sizeof(unsigned int));

          for (unsigned int n=0; n<mesh.n_nodes(); n++)
            temp[n] = static_cast<float>(std::abs((*vec)[n*n_vars + c]));

          out_stream.write(reinterpret_cast<char *>(&temp[0]), sizeof(float)*mesh.n_nodes());

#else

          buffer = (*solution_names)[c];
          out_stream.write(buffer.c_str(), buffer.size());

          unsigned int tempint = 1; // always do nodal data
          out_stream.write(reinterpret_cast<char *>(&tempint), sizeof(unsigned int));

          for (unsigned int n=0; n<mesh.n_nodes(); n++)
            temp[n] = static_cast<float>((*vec)[n*n_vars + c]);

          out_stream.write(reinterpret_cast<char *>(&temp[0]), sizeof(float)*mesh.n_nodes());

#endif
        }
    }

  // If we wrote any variables, we have to close the variable section now
  if (write_variable)
    {
      buffer = "endvars ";
      out_stream.write(buffer.c_str(), buffer.size());
    }

  // end the file
  buffer = "endgmv  ";
  out_stream.write(buffer.c_str(), buffer.size());
}









void GMVIO::write_discontinuous_gmv (const std::string & name,
                                     const EquationSystems & es,
                                     const bool write_partitioning,
                                     const std::set<std::string> * system_names) const
{
  std::vector<std::string> solution_names;
  std::vector<Number>      v;

  // Get a reference to the mesh
  const MeshBase & mesh = MeshOutput<MeshBase>::mesh();

  es.build_variable_names  (solution_names, libmesh_nullptr, system_names);
  es.build_discontinuous_solution_vector (v, system_names);

  // These are parallel_only functions
  const unsigned int n_active_elem = mesh.n_active_elem();

  if (mesh.processor_id() != 0)
    return;

  std::ofstream out_stream(name.c_str());

  libmesh_assert (out_stream.good());

  // Begin interfacing with the GMV data file
  {

    // write the nodes
    out_stream << "gmvinput ascii" << std::endl << std::endl;

    // Compute the total weight
    {
      MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
      const MeshBase::const_element_iterator end = mesh.active_elements_end();

      unsigned int tw=0;

      for ( ; it != end; ++it)
        tw += (*it)->n_nodes();

      out_stream << "nodes " << tw << std::endl;
    }



    // Write all the x values
    {
      MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
      const MeshBase::const_element_iterator end = mesh.active_elements_end();

      for ( ; it != end; ++it)
        for (unsigned int n=0; n<(*it)->n_nodes(); n++)
          out_stream << (*it)->point(n)(0) << " ";

      out_stream << std::endl;
    }


    // Write all the y values
    {
      MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
      const MeshBase::const_element_iterator end = mesh.active_elements_end();

      for ( ; it != end; ++it)
        for (unsigned int n=0; n<(*it)->n_nodes(); n++)
#if LIBMESH_DIM > 1
          out_stream << (*it)->point(n)(1) << " ";
#else
      out_stream << 0. << " ";
#endif

      out_stream << std::endl;
    }


    // Write all the z values
    {
      MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
      const MeshBase::const_element_iterator end = mesh.active_elements_end();

      for ( ; it != end; ++it)
        for (unsigned int n=0; n<(*it)->n_nodes(); n++)
#if LIBMESH_DIM > 2
          out_stream << (*it)->point(n)(2) << " ";
#else
      out_stream << 0. << " ";
#endif

      out_stream << std::endl << std::endl;
    }
  }



  {
    // write the connectivity

    out_stream << "cells " << n_active_elem << std::endl;

    MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
    const MeshBase::const_element_iterator end = mesh.active_elements_end();

    unsigned int nn=1;

    switch (mesh.mesh_dimension())
      {
      case 1:
        {
          for ( ; it != end; ++it)
            for (unsigned int se=0; se<(*it)->n_sub_elem(); se++)
              {
                if (((*it)->type() == EDGE2) ||
                    ((*it)->type() == EDGE3) ||
                    ((*it)->type() == EDGE4)
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
                    || ((*it)->type() == INFEDGE2)
#endif
                    )
                  {
                    out_stream << "line 2" << std::endl;
                    for (unsigned int i=0; i<(*it)->n_nodes(); i++)
                      out_stream << nn++ << " ";

                  }
                else
                  libmesh_error_msg("Unsupported 1D element type: " << Utility::enum_to_string((*it)->type()));

                out_stream << std::endl;
              }

          break;
        }

      case 2:
        {
          for ( ; it != end; ++it)
            for (unsigned int se=0; se<(*it)->n_sub_elem(); se++)
              {
                if (((*it)->type() == QUAD4) ||
                    ((*it)->type() == QUAD8) || // Note: QUAD8 will be output as one central quad and
                    // four surrounding triangles (though they will be written
                    // to GMV as QUAD4s).
                    ((*it)->type() == QUAD9)
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
                    || ((*it)->type() == INFQUAD4)
                    || ((*it)->type() == INFQUAD6)
#endif
                    )
                  {
                    out_stream << "quad 4" << std::endl;
                    for (unsigned int i=0; i<(*it)->n_nodes(); i++)
                      out_stream << nn++ << " ";

                  }
                else if (((*it)->type() == TRI3) ||
                         ((*it)->type() == TRI6))
                  {
                    out_stream << "tri 3" << std::endl;
                    for (unsigned int i=0; i<(*it)->n_nodes(); i++)
                      out_stream << nn++ << " ";

                  }
                else
                  libmesh_error_msg("Unsupported 2D element type: " << Utility::enum_to_string((*it)->type()));

                out_stream << std::endl;
              }

          break;
        }


      case 3:
        {
          for ( ; it != end; ++it)
            for (unsigned int se=0; se<(*it)->n_sub_elem(); se++)
              {
                if (((*it)->type() == HEX8) ||
                    ((*it)->type() == HEX20) ||
                    ((*it)->type() == HEX27)
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
                    || ((*it)->type() == INFHEX8)
                    || ((*it)->type() == INFHEX16)
                    || ((*it)->type() == INFHEX18)
#endif
                    )
                  {
                    out_stream << "phex8 8" << std::endl;
                    for (unsigned int i=0; i<(*it)->n_nodes(); i++)
                      out_stream << nn++ << " ";
                  }
                else if (((*it)->type() == PRISM6) ||
                         ((*it)->type() == PRISM15) ||
                         ((*it)->type() == PRISM18)
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
                         || ((*it)->type() == INFPRISM6)
                         || ((*it)->type() == INFPRISM12)
#endif
                         )
                  {
                    out_stream << "pprism6 6" << std::endl;
                    for (unsigned int i=0; i<(*it)->n_nodes(); i++)
                      out_stream << nn++ << " ";
                  }
                else if (((*it)->type() == TET4) ||
                         ((*it)->type() == TET10))
                  {
                    out_stream << "tet 4" << std::endl;
                    for (unsigned int i=0; i<(*it)->n_nodes(); i++)
                      out_stream << nn++ << " ";
                  }
                else
                  libmesh_error_msg("Unsupported 3D element type: " << Utility::enum_to_string((*it)->type()));

                out_stream << std::endl;
              }

          break;
        }

      default:
        libmesh_error_msg("Unsupported mesh dimension: " << mesh.mesh_dimension());
      }

    out_stream << std::endl;
  }



  // optionally write the partition information
  if (write_partitioning)
    {
      if (_write_subdomain_id_as_material)
        libmesh_error_msg("Not yet supported in GMVIO::write_discontinuous_gmv");

      else
        {
          out_stream << "material "
                     << mesh.n_processors()
                     << " 0"<< std::endl;

          for (unsigned int proc=0; proc<mesh.n_processors(); proc++)
            out_stream << "proc_" << proc << std::endl;

          MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
          const MeshBase::const_element_iterator end = mesh.active_elements_end();

          for ( ; it != end; ++it)
            out_stream << (*it)->processor_id()+1 << std::endl;

          out_stream << std::endl;
        }
    }


  // Writing cell-centered data is not yet supported in discontinuous GMV files.
  if ( !(this->_cell_centered_data.empty()) )
    {
      libMesh::err << "Cell-centered data not (yet) supported for discontinuous GMV files!" << std::endl;
    }



  // write the data
  {
    const unsigned int n_vars =
      cast_int<unsigned int>(solution_names.size());

    //    libmesh_assert_equal_to (v.size(), tw*n_vars);

    out_stream << "variable" << std::endl;


    for (unsigned int c=0; c<n_vars; c++)
      {

#ifdef LIBMESH_USE_COMPLEX_NUMBERS

        // in case of complex data, write _tree_ data sets
        // for each component

        // this is the real part
        out_stream << "r_" << solution_names[c] << " 1" << std::endl;
        {
          MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
          const MeshBase::const_element_iterator end = mesh.active_elements_end();

          for ( ; it != end; ++it)
            for (unsigned int n=0; n<(*it)->n_nodes(); n++)
              out_stream << v[(n++)*n_vars + c].real() << " ";
        }
        out_stream << std::endl << std::endl;


        // this is the imaginary part
        out_stream << "i_" << solution_names[c] << " 1" << std::endl;
        {
          MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
          const MeshBase::const_element_iterator end = mesh.active_elements_end();

          for ( ; it != end; ++it)
            for (unsigned int n=0; n<(*it)->n_nodes(); n++)
              out_stream << v[(n++)*n_vars + c].imag() << " ";
        }
        out_stream << std::endl << std::endl;

        // this is the magnitude
        out_stream << "a_" << solution_names[c] << " 1" << std::endl;
        {
          MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
          const MeshBase::const_element_iterator end = mesh.active_elements_end();

          for ( ; it != end; ++it)
            for (unsigned int n=0; n<(*it)->n_nodes(); n++)
              out_stream << std::abs(v[(n++)*n_vars + c]) << " ";
        }
        out_stream << std::endl << std::endl;

#else

        out_stream << solution_names[c] << " 1" << std::endl;
        {
          MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
          const MeshBase::const_element_iterator end = mesh.active_elements_end();

          unsigned int nn=0;

          for ( ; it != end; ++it)
            for (unsigned int n=0; n<(*it)->n_nodes(); n++)
              out_stream << v[(nn++)*n_vars + c] << " ";
        }
        out_stream << std::endl << std::endl;

#endif

      }

    out_stream << "endvars" << std::endl;
  }


  // end of the file
  out_stream << std::endl << "endgmv" << std::endl;
}





void GMVIO::add_cell_centered_data (const std::string &       cell_centered_data_name,
                                    const std::vector<Real> * cell_centered_data_vals)
{
  libmesh_assert(cell_centered_data_vals);

  // Make sure there are *at least* enough entries for all the active elements.
  // There can also be entries for inactive elements, they will be ignored.
  // libmesh_assert_greater_equal (cell_centered_data_vals->size(),
  //                           MeshOutput<MeshBase>::mesh().n_active_elem());
  this->_cell_centered_data[cell_centered_data_name] = cell_centered_data_vals;
}






void GMVIO::read (const std::string & name)
{
  // This is a serial-only process for now;
  // the Mesh should be read on processor 0 and
  // broadcast later
  libmesh_assert_equal_to (MeshOutput<MeshBase>::mesh().processor_id(), 0);

  _next_elem_id = 0;

  libmesh_experimental();

#ifndef LIBMESH_HAVE_GMV

  libmesh_error_msg("Cannot read GMV file " << name << " without the GMV API.");

#else
  // We use the file-scope global variable eletypes for mapping nodes
  // from GMV to libmesh indices, so initialize that data now.
  init_eletypes();

  // Clear the mesh so we are sure to start from a pristeen state.
  MeshBase & mesh = MeshInput<MeshBase>::mesh();
  mesh.clear();

  // Keep track of what kinds of elements this file contains
  elems_of_dimension.clear();
  elems_of_dimension.resize(4, false);

  // It is apparently possible for gmv files to contain
  // a "fromfile" directive (?) But we currently don't make
  // any use of this feature in LibMesh.  Nonzero return val
  // from any function usually means an error has occurred.
  int ierr = GMVLib::gmvread_open_fromfileskip(const_cast<char *>(name.c_str()));
  if (ierr != 0)
    libmesh_error_msg("GMVLib::gmvread_open_fromfileskip failed!");


  // Loop through file until GMVEND.
  int iend = 0;
  while (iend == 0)
    {
      GMVLib::gmvread_data();

      // Check for GMVEND.
      if (GMVLib::gmv_data.keyword == GMVEND)
        {
          iend = 1;
          GMVLib::gmvread_close();
          break;
        }

      // Check for GMVERROR.
      if (GMVLib::gmv_data.keyword == GMVERROR)
        libmesh_error_msg("Encountered GMVERROR while reading!");

      // Process the data.
      switch (GMVLib::gmv_data.keyword)
        {
        case NODES:
          {
            //libMesh::out << "Reading nodes." << std::endl;

            if (GMVLib::gmv_data.num2 == NODES)
              this->_read_nodes();

            else if (GMVLib::gmv_data.num2 == NODE_V)
              libmesh_error_msg("Unsupported GMV data type NODE_V!");

            break;
          }

        case CELLS:
          {
            // Read 1 cell at a time
            // libMesh::out << "\nReading one cell." << std::endl;
            this->_read_one_cell();
            break;
          }

        case MATERIAL:
          {
            // keyword == 6
            // These are the materials, which we use to specify the mesh
            // partitioning.
            this->_read_materials();
            break;
          }

        case VARIABLE:
          {
            // keyword == 8
            // This is a field variable.

            // Check to see if we're done reading variables and break out.
            if (GMVLib::gmv_data.datatype == ENDKEYWORD)
              {
                // libMesh::out << "Done reading GMV variables." << std::endl;
                break;
              }

            if (GMVLib::gmv_data.datatype == NODE)
              {
                // libMesh::out << "Reading node field data for variable "
                //   << GMVLib::gmv_data.name1 << std::endl;
                this->_read_var();
                break;
              }

            else
              {
                libMesh::err << "Warning: Skipping variable: "
                             << GMVLib::gmv_data.name1
                             << " which is of unsupported GMV datatype "
                             << GMVLib::gmv_data.datatype
                             << ".  Nodal field data is currently the only type currently supported."
                             << std::endl;
                break;
              }

          }

        default:
          libmesh_error_msg("Encountered unknown GMV keyword " << GMVLib::gmv_data.keyword);

        } // end switch
    } // end while

  // Set the mesh dimension to the largest encountered for an element
  for (unsigned char i=0; i!=4; ++i)
    if (elems_of_dimension[i])
      mesh.set_mesh_dimension(i);

#if LIBMESH_DIM < 3
  if (mesh.mesh_dimension() > LIBMESH_DIM)
    libmesh_error_msg("Cannot open dimension " \
                      << mesh.mesh_dimension()            \
                      << " mesh file when configured without "        \
                      << mesh.mesh_dimension()                        \
                      << "D support.");
#endif

  // Done reading in the mesh, now call find_neighbors, etc.
  // mesh.find_neighbors();

  // Don't allow renumbering of nodes and elements when calling
  // prepare_for_use().  It is usually confusing for users when the
  // Mesh gets renumbered automatically, since sometimes there are
  // other parts of the simulation which rely on the original
  // ordering.
  mesh.allow_renumbering(false);
  mesh.prepare_for_use();
#endif
}




void GMVIO::_read_var()
{
#ifdef LIBMESH_HAVE_GMV

  // Copy all the variable's values into a local storage vector.
  _nodal_data.insert ( std::make_pair(std::string(GMVLib::gmv_data.name1),
                                      std::vector<Number>(GMVLib::gmv_data.doubledata1, GMVLib::gmv_data.doubledata1+GMVLib::gmv_data.num) ) );
#endif
}



void GMVIO::_read_materials()
{
#ifdef LIBMESH_HAVE_GMV

  // LibMesh assigns materials on a per-cell basis
  libmesh_assert_equal_to (GMVLib::gmv_data.datatype, CELL);

  //   // Material names: LibMesh has no use for these currently...
  //   libMesh::out << "Number of material names="
  //     << GMVLib::gmv_data.num
  //     << std::endl;

  //   for (int i = 0; i < GMVLib::gmv_data.num; i++)
  //     {
  //       // Build a 32-char string from the appropriate entries
  //       std::string mat_string(&GMVLib::gmv_data.chardata1[i*33], 32);

  //       libMesh::out << "Material name " << i << ": " << mat_string << std::endl;
  //     }

  //   // Material labels: These correspond to (1-based) CPU IDs, and
  //   // there should be 1 of these for each element.
  //   libMesh::out << "Number of material labels = "
  //     << GMVLib::gmv_data.nlongdata1
  //     << std::endl;

  for (int i = 0; i < GMVLib::gmv_data.nlongdata1; i++)
    {
      // Debugging Info
      // libMesh::out << "Material ID " << i << ": "
      // << GMVLib::gmv_data.longdata1[i]
      // << std::endl;

      MeshInput<MeshBase>::mesh().elem_ref(i).processor_id() =
        cast_int<processor_id_type>(GMVLib::gmv_data.longdata1[i]-1);
    }

#endif
}




void GMVIO::_read_nodes()
{
#ifdef LIBMESH_HAVE_GMV
  // Debugging
  // libMesh::out << "gmv_data.datatype = " << GMVLib::gmv_data.datatype << std::endl;

  // LibMesh writes UNSTRUCT=100 node data
  libmesh_assert_equal_to (GMVLib::gmv_data.datatype, UNSTRUCT);

  // The nodal data is stored in gmv_data.doubledata{1,2,3}
  // and is nnodes long
  for (int i = 0; i < GMVLib::gmv_data.num; i++)
    {
      // Debugging
      // libMesh::out << "(x,y,z)="
      //              << "("
      //              << GMVLib::gmv_data.doubledata1[i] << ","
      //              << GMVLib::gmv_data.doubledata2[i] << ","
      //              << GMVLib::gmv_data.doubledata3[i]
      //              << ")"
      //              << std::endl;

      // Add the point to the Mesh
      MeshInput<MeshBase>::mesh().add_point(Point(GMVLib::gmv_data.doubledata1[i],
                                                  GMVLib::gmv_data.doubledata2[i],
                                                  GMVLib::gmv_data.doubledata3[i]), i);
    }
#endif
}


void GMVIO::_read_one_cell()
{
#ifdef LIBMESH_HAVE_GMV
  // Debugging
  // libMesh::out << "gmv_data.datatype=" << GMVLib::gmv_data.datatype << std::endl;

  // This is either a REGULAR=111 cell or
  // the ENDKEYWORD=207 of the cells
#ifndef NDEBUG
  bool recognized =
    (GMVLib::gmv_data.datatype==REGULAR) ||
    (GMVLib::gmv_data.datatype==ENDKEYWORD);
#endif
  libmesh_assert (recognized);

  MeshBase & mesh = MeshInput<MeshBase>::mesh();

  if (GMVLib::gmv_data.datatype == REGULAR)
    {
      // Debugging
      // libMesh::out << "Name of the cell is: " << GMVLib::gmv_data.name1 << std::endl;
      // libMesh::out << "Cell has " << GMVLib::gmv_data.num2 << " vertices." << std::endl;

      // We need a mapping from GMV element types to LibMesh
      // ElemTypes.  Basically the reverse of the eletypes
      // std::map above.
      //
      // FIXME: Since Quad9's apparently don't exist for GMV, and since
      // In general we write linear sub-elements to GMV files, we need
      // to be careful to read back in exactly what we wrote out...
      //
      // The gmv_data.name1 field is padded with blank characters when
      // it is read in by GMV, so the function that converts it to the
      // libmesh element type needs to take that into account.
      ElemType type = this->gmv_elem_to_libmesh_elem(GMVLib::gmv_data.name1);

      Elem * elem = Elem::build(type).release();
      elem->set_id(_next_elem_id++);

      // Get the ElementDefinition object for this element type
      const ElementDefinition & eledef = eletypes[type];

      // Print out the connectivity information for
      // this cell.
      for (int i=0; i<GMVLib::gmv_data.num2; i++)
        {
          // Debugging
          // libMesh::out << "Vertex " << i << " is node " << GMVLib::gmv_data.longdata1[i] << std::endl;

          // Map index i to GMV's numbering scheme
          unsigned mapped_i = eledef.node_map[i];

          // Note: Node numbers (as stored in libmesh) are 1-based
          elem->set_node(i) = mesh.node_ptr
            (cast_int<dof_id_type>(GMVLib::gmv_data.longdata1[mapped_i]-1));
        }

      elems_of_dimension[elem->dim()] = true;

      // Add the newly-created element to the mesh
      mesh.add_elem(elem);
    }


  if (GMVLib::gmv_data.datatype == ENDKEYWORD)
    {
      // There isn't a cell to read, so we just return
      return;
    }

#endif
}


ElemType GMVIO::gmv_elem_to_libmesh_elem(std::string elemname)
{
  // Erase any whitespace padding in name coming from gmv before performing comparison.
  elemname.erase(std::remove_if(elemname.begin(), elemname.end(), isspace), elemname.end());

  // Look up the string in our string->ElemType name.
  std::map<std::string, ElemType>::iterator it = _reading_element_map.find(elemname);

  if (it == _reading_element_map.end())
    libmesh_error_msg("Uknown/unsupported element: " << elemname << " was read.");

  return it->second;
}




void GMVIO::copy_nodal_solution(EquationSystems & es)
{
  // Check for easy return if there isn't any nodal data
  if (_nodal_data.empty())
    {
      libMesh::err << "Unable to copy nodal solution: No nodal "
                   << "solution has been read in from file." << std::endl;
      return;
    }

  // Be sure there is at least one system
  libmesh_assert (es.n_systems());

  // Keep track of variable names which have been found and
  // copied already.  This could be used to prevent us from
  // e.g. copying the same var into 2 different systems ...
  // but this seems unlikely.  Also, it is used to tell if
  // any variables which were read in were not successfully
  // copied to the EquationSystems.
  std::set<std::string> vars_copied;

  // For each entry in the nodal data map, try to find a system
  // that has the same variable key name.
  for (unsigned int sys=0; sys<es.n_systems(); ++sys)
    {
      // Get a generic refernence to the current System
      System & system = es.get_system(sys);

      // And a reference to that system's dof_map
      // const DofMap & dof_map = system.get_dof_map();

      // For each var entry in the _nodal_data map, try to find
      // that var in the system
      std::map<std::string, std::vector<Number> >::iterator it = _nodal_data.begin();
      const std::map<std::string, std::vector<Number> >::iterator end = _nodal_data.end();
      for (; it != end; ++it)
        {
          std::string var_name = it->first;
          // libMesh::out << "Searching for var " << var_name << " in system " << sys << std::endl;

          if (system.has_variable(var_name))
            {
              // Check if there are as many nodes in the mesh as there are entries
              // in the stored nodal data vector
              libmesh_assert_equal_to ( it->second.size(), MeshInput<MeshBase>::mesh().n_nodes() );

              const unsigned int var_num = system.variable_number(var_name);

              // libMesh::out << "Variable "
              // << var_name
              // << " is variable "
              // << var_num
              // << " in system " << sys << std::endl;

              // The only type of nodal data we can read in from GMV is for
              // linear LAGRANGE type elements.
              const FEType & fe_type = system.variable_type(var_num);
              if ((fe_type.order.get_order() != FIRST) || (fe_type.family != LAGRANGE))
                {
                  libMesh::err << "Only FIRST-order LAGRANGE variables can be read from GMV files. "
                               << "Skipping variable " << var_name << std::endl;
                  break;
                }


              // Loop over the stored vector's entries, inserting them into
              // the System's solution if appropriate.
              for (std::size_t i=0; i<it->second.size(); ++i)
                {
                  // Since this var came from a GMV file, the index i corresponds to
                  // the (single) DOF value of the current variable for node i.
                  const unsigned int dof_index =
                    MeshInput<MeshBase>::mesh().node_ptr(i)->dof_number(sys,      // system #
                                                                        var_num,  // var #
                                                                        0);       // component #, always zero for LAGRANGE

                  // libMesh::out << "Value " << i << ": "
                  //     << it->second [i]
                  //     << ", dof index="
                  //     << dof_index << std::endl;

                  // If the dof_index is local to this processor, set the value
                  if ((dof_index >= system.solution->first_local_index()) &&
                      (dof_index <  system.solution->last_local_index()))
                    system.solution->set (dof_index, it->second [i]);
                } // end loop over my GMVIO's copy of the solution

              // Add the most recently copied var to the set of copied vars
              vars_copied.insert (var_name);
            } // end if (system.has_variable)
        } // end for loop over _nodal_data

      // Communicate parallel values before going to the next system.
      system.solution->close();
      system.update();

    } // end loop over all systems



  // Warn the user if any GMV variables were not successfully copied over to the EquationSystems object
  {
    std::map<std::string, std::vector<Number> >::iterator it = _nodal_data.begin();
    const std::map<std::string, std::vector<Number> >::iterator end = _nodal_data.end();

    for (; it != end; ++it)
      {
        if (vars_copied.find( it->first ) == vars_copied.end())
          {
            libMesh::err << "Warning: Variable "
                         << it->first
                         << " was not copied to the EquationSystems object."
                         << std::endl;
          }
      }
  }

}

} // namespace libMesh
