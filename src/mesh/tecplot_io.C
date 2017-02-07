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
#include <fstream>
#include <iomanip>
#include <sstream>

// Local includes
#include "libmesh/libmesh_config.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/tecplot_io.h"
#include "libmesh/mesh_base.h"
#include "libmesh/elem.h"
#include "libmesh/parallel.h"

#ifdef LIBMESH_HAVE_TECPLOT_API
extern "C" {
# include <TECIO.h>
}
#endif


namespace libMesh
{


//--------------------------------------------------------
// Macros for handling Tecplot API data

#ifdef LIBMESH_HAVE_TECPLOT_API

namespace
{
class TecplotMacros
{
public:
  TecplotMacros(const dof_id_type n_nodes,
                const unsigned int n_vars,
                const dof_id_type n_cells,
                const unsigned int n_vert);
  float & nd(const std::size_t i, const std::size_t j);
  int   & cd(const std::size_t i, const std::size_t j);
  std::vector<float> nodalData;
  std::vector<int>   connData;
  //float* nodalData;
  //int *   connData;

  void set_n_cells (const dof_id_type nc);

  const dof_id_type n_nodes;
  const unsigned int n_vars;
  dof_id_type n_cells;
  const unsigned int n_vert;
};
}



inline
TecplotMacros::TecplotMacros(const dof_id_type nn,
                             const unsigned int nvar,
                             const dof_id_type nc,
                             const unsigned int nvrt) :
  n_nodes(nn),
  n_vars(nvar),
  n_cells(nc),
  n_vert(nvrt)
{
  nodalData.resize(n_nodes*n_vars);
  connData.resize(n_cells*n_vert);
}



inline
float & TecplotMacros::nd(const std::size_t i, const std::size_t j)
{
  return nodalData[(i)*(n_nodes) + (j)];
}



inline
int & TecplotMacros::cd(const std::size_t i, const std::size_t j)
{
  return connData[(i) + (j)*(n_vert)];
}


inline
void TecplotMacros::set_n_cells (const dof_id_type nc)
{
  n_cells = nc;
  connData.resize(n_cells*n_vert);
}
#endif
//--------------------------------------------------------



// ------------------------------------------------------------
// TecplotIO  members
TecplotIO::TecplotIO (const MeshBase & mesh_in,
                      const bool binary_in,
                      const double time_in,
                      const int strand_offset_in) :
  MeshOutput<MeshBase> (mesh_in),
  _binary (binary_in),
  _time (time_in),
  _strand_offset (strand_offset_in),
  _zone_title ("zone"),
  _ascii_append(false)
{
  // Gather a list of subdomain ids in the mesh.
  // We must do this now, while we have every
  // processor's attention
  // (some of the write methods only execute on processor 0).
  mesh_in.subdomain_ids (_subdomain_ids);
}



bool & TecplotIO::binary ()
{
  return _binary;
}



double & TecplotIO::time ()
{
  return _time;
}



int & TecplotIO::strand_offset ()
{
  return _strand_offset;
}



std::string & TecplotIO::zone_title ()
{
  return _zone_title;
}


bool & TecplotIO::ascii_append ()
{
  return _ascii_append;
}


void TecplotIO::write (const std::string & fname)
{
  if (this->mesh().processor_id() == 0)
    {
      if (this->binary())
        this->write_binary (fname);
      else
        this->write_ascii  (fname);
    }
}



void TecplotIO::write_nodal_data (const std::string & fname,
                                  const std::vector<Number> & soln,
                                  const std::vector<std::string> & names)
{
  LOG_SCOPE("write_nodal_data()", "TecplotIO");

  if (this->mesh().processor_id() == 0)
    {
      if (this->binary())
        this->write_binary (fname, &soln, &names);
      else
        this->write_ascii  (fname, &soln, &names);
    }
}



unsigned TecplotIO::elem_dimension()
{
  // Get a constant reference to the mesh.
  const MeshBase & the_mesh = MeshOutput<MeshBase>::mesh();

  std::vector<unsigned> elem_dims(3);

  // Loop over all the elements and mark the proper dimension entry in
  // the elem_dims vector.
  MeshBase::const_element_iterator       it  = the_mesh.active_elements_begin();
  const MeshBase::const_element_iterator end = the_mesh.active_elements_end();
  for ( ; it != end; ++it)
    elem_dims[(*it)->dim() - 1] = 1;

  // Detect and disallow (for now) the writing of mixed dimension meshes.
  if (std::count(elem_dims.begin(), elem_dims.end(), 1) > 1)
    libmesh_error_msg("Error, cannot write Mesh with mixed element dimensions to Tecplot file!");

  if (elem_dims[0])
    return 1;
  else if (elem_dims[1])
    return 2;
  else if (elem_dims[2])
    return 3;
  else
    libmesh_error_msg("No 1, 2, or 3D elements detected!");
}



void TecplotIO::write_ascii (const std::string & fname,
                             const std::vector<Number> * v,
                             const std::vector<std::string> * solution_names)
{
  // Should only do this on processor 0!
  libmesh_assert_equal_to (this->mesh().processor_id(), 0);

  // Create an output stream, possibly in append mode.
  std::ofstream out_stream(fname.c_str(), _ascii_append ? std::ofstream::app : std::ofstream::out);

  // Make sure it opened correctly
  if (!out_stream.good())
    libmesh_file_error(fname.c_str());

  // Get a constant reference to the mesh.
  const MeshBase & the_mesh = MeshOutput<MeshBase>::mesh();

  // Write header to stream
  {
    {
      // TODO: We used to print out the SVN revision here when we did keyword expansions...
      out_stream << "# For a description of the Tecplot format see the Tecplot User's guide.\n"
                 << "#\n";
    }

    out_stream << "Variables=x,y,z";

    if (solution_names != libmesh_nullptr)
      for (std::size_t n=0; n<solution_names->size(); n++)
        {
#ifdef LIBMESH_USE_REAL_NUMBERS

          // Write variable names for real variables
          out_stream << "," << (*solution_names)[n];

#else

          // Write variable names for complex variables
          out_stream << "," << "r_"   << (*solution_names)[n]
                     << "," << "i_"   << (*solution_names)[n]
                     << "," << "a_"   << (*solution_names)[n];

#endif
        }

    out_stream << '\n';

    out_stream << "Zone f=fepoint, n=" << the_mesh.n_nodes() << ", e=" << the_mesh.n_active_sub_elem();

    // We cannot choose the element type simply based on the mesh
    // dimension... there might be 1D elements living in a 3D mesh.
    // So look at the elements which are actually in the Mesh, and
    // choose either "lineseg", "quadrilateral", or "brick" depending
    // on if the elements are 1, 2, or 3D.

    // Write the element type we've determined to the header.
    out_stream << ", et=";

    switch (this->elem_dimension())
      {
      case 1:
        out_stream << "lineseg";
        break;
      case 2:
        out_stream << "quadrilateral";
        break;
      case 3:
        out_stream << "brick";
        break;
      default:
        libmesh_error_msg("Unsupported element dimension: " << this->elem_dimension());
      }

    // Output the time in the header
    out_stream << ", t=\"T " << _time << "\"";

    // Use default mesh color = black
    out_stream << ", c=black\n";

  } // finished writing header

  for (unsigned int i=0; i<the_mesh.n_nodes(); i++)
    {
      // Print the point without a newline
      the_mesh.point(i).write_unformatted(out_stream, false);

      if ((v != libmesh_nullptr) && (solution_names != libmesh_nullptr))
        {
          const std::size_t n_vars = solution_names->size();


          for (std::size_t c=0; c<n_vars; c++)
            {
#ifdef LIBMESH_USE_REAL_NUMBERS
              // Write real data
              out_stream << std::setprecision(this->ascii_precision())
                         << (*v)[i*n_vars + c] << " ";

#else
              // Write complex data
              out_stream << std::setprecision(this->ascii_precision())
                         << (*v)[i*n_vars + c].real() << " "
                         << (*v)[i*n_vars + c].imag() << " "
                         << std::abs((*v)[i*n_vars + c]) << " ";

#endif
            }
        }

      // Write a new line after the data for this node
      out_stream << '\n';
    }

  MeshBase::const_element_iterator       it  = the_mesh.active_elements_begin();
  const MeshBase::const_element_iterator end = the_mesh.active_elements_end();

  for ( ; it != end; ++it)
    (*it)->write_connectivity(out_stream, TECPLOT);
}



void TecplotIO::write_binary (const std::string & fname,
                              const std::vector<Number> * vec,
                              const std::vector<std::string> * solution_names)
{
  //-----------------------------------------------------------
  // Call the ASCII output function if configure did not detect
  // the Tecplot binary API
#ifndef LIBMESH_HAVE_TECPLOT_API

  libMesh::err << "WARNING: Tecplot Binary files require the Tecplot API." << std::endl
               << "Continuing with ASCII output."
               << std::endl;

  if (this->mesh().processor_id() == 0)
    this->write_ascii (fname, vec, solution_names);
  return;



  //------------------------------------------------------------
  // New binary formats, time aware and whatnot
#elif defined(LIBMESH_HAVE_TECPLOT_API_112)

  // Get a constant reference to the mesh.
  const MeshBase & the_mesh = MeshOutput<MeshBase>::mesh();

  // Required variables
  std::string tecplot_variable_names;
  int
    ierr      =  0,
    file_type =  0, // full
    is_double =  0,
#ifdef DEBUG
    tec_debug =  1,
#else
    tec_debug =  0,
#endif
    cell_type   = -1,
    nn_per_elem = -1;

  switch (this->elem_dimension())
    {
    case 1:
      cell_type   = 1;  // FELINESEG
      nn_per_elem = 2;
      break;

    case 2:
      cell_type   = 3; // FEQUADRILATERAL
      nn_per_elem = 4;
      break;

    case 3:
      cell_type   = 5; // FEBRICK
      nn_per_elem = 8;
      break;

    default:
      libmesh_error_msg("Unsupported element dimension: " << this->elem_dimension());
    }

  // Build a string containing all the variable names to pass to Tecplot
  {
    tecplot_variable_names += "x, y, z";

    if (solution_names != libmesh_nullptr)
      {
        for (std::size_t name=0; name<solution_names->size(); name++)
          {
#ifdef LIBMESH_USE_REAL_NUMBERS

            tecplot_variable_names += ", ";
            tecplot_variable_names += (*solution_names)[name];

#else

            tecplot_variable_names += ", ";
            tecplot_variable_names += "r_";
            tecplot_variable_names += (*solution_names)[name];
            tecplot_variable_names += ", ";
            tecplot_variable_names += "i_";
            tecplot_variable_names += (*solution_names)[name];
            tecplot_variable_names += ", ";
            tecplot_variable_names += "a_";
            tecplot_variable_names += (*solution_names)[name];

#endif
          }
      }
  }

  // Instantiate a TecplotMacros interface.  In 2D the most nodes per
  // face should be 4, in 3D it's 8.


  TecplotMacros tm(the_mesh.n_nodes(),
#ifdef LIBMESH_USE_REAL_NUMBERS
                   (3 + ((solution_names == libmesh_nullptr) ? 0 :
                         cast_int<unsigned int>(solution_names->size()))),
#else
                   (3 + 3*((solution_names == libmesh_nullptr) ? 0 :
                           cast_int<unsigned int>(solution_names->size()))),
#endif
                   the_mesh.n_active_sub_elem(),
                   nn_per_elem
                   );


  // Copy the nodes and data to the TecplotMacros class. Note that we store
  // everything as a float here since the eye doesn't require a double to
  // understand what is going on
  for (unsigned int v=0; v<the_mesh.n_nodes(); v++)
    {
      tm.nd(0,v) = static_cast<float>(the_mesh.point(v)(0));
      tm.nd(1,v) = static_cast<float>(the_mesh.point(v)(1));
      tm.nd(2,v) = static_cast<float>(the_mesh.point(v)(2));

      if ((vec != libmesh_nullptr) &&
          (solution_names != libmesh_nullptr))
        {
          const std::size_t n_vars = solution_names->size();

          for (std::size_t c=0; c<n_vars; c++)
            {
#ifdef LIBMESH_USE_REAL_NUMBERS

              tm.nd((3+c),v)     = static_cast<float>((*vec)[v*n_vars + c]);
#else
              tm.nd((3+3*c),v)   = static_cast<float>((*vec)[v*n_vars + c].real());
              tm.nd((3+3*c+1),v) = static_cast<float>((*vec)[v*n_vars + c].imag());
              tm.nd((3+3*c+2),v) = static_cast<float>(std::abs((*vec)[v*n_vars + c]));
#endif
            }
        }
    }


  // Initialize the file
  ierr = TECINI112 (libmesh_nullptr,
                    const_cast<char *>(tecplot_variable_names.c_str()),
                    const_cast<char *>(fname.c_str()),
                    const_cast<char *>("."),
                    &file_type,
                    &tec_debug,
                    &is_double);

  if (ierr)
    libmesh_file_error(fname);

  // A zone for each subdomain
  bool firstzone=true;
  for (std::set<subdomain_id_type>::const_iterator sbd_it=_subdomain_ids.begin();
       sbd_it!=_subdomain_ids.end(); ++sbd_it)
    {
      // Copy the connectivity for this subdomain
      {
        MeshBase::const_element_iterator       it  = the_mesh.active_subdomain_elements_begin (*sbd_it);
        const MeshBase::const_element_iterator end = the_mesh.active_subdomain_elements_end   (*sbd_it);

        unsigned int n_subcells_in_subdomain=0;

        for (; it != end; ++it)
          n_subcells_in_subdomain += (*it)->n_sub_elem();

        // update the connectivty array to include only the elements in this subdomain
        tm.set_n_cells (n_subcells_in_subdomain);

        unsigned int te = 0;

        for (it  = the_mesh.active_subdomain_elements_begin (*sbd_it);
             it != end; ++it)
          {
            std::vector<dof_id_type> conn;
            for (unsigned int se=0; se<(*it)->n_sub_elem(); se++)
              {
                (*it)->connectivity(se, TECPLOT, conn);

                for (std::size_t node=0; node<conn.size(); node++)
                  tm.cd(node,te) = conn[node];

                te++;
              }
          }
      }


      // Ready to call the Tecplot API for this subdomain
      {
        int
          num_nodes   = static_cast<int>(the_mesh.n_nodes()),
          num_cells   = static_cast<int>(tm.n_cells),
          num_faces   = 0,
          i_cell_max  = 0,
          j_cell_max  = 0,
          k_cell_max  = 0,
          strand_id   = std::max(*sbd_it,static_cast<subdomain_id_type>(1)) + this->strand_offset(),
          parent_zone = 0,
          is_block    = 1,
          num_face_connect   = 0,
          face_neighbor_mode = 0,
          tot_num_face_nodes = 0,
          num_connect_boundary_faces = 0,
          tot_num_boundary_connect   = 0,
          share_connect_from_zone=0;

        std::vector<int>
          passive_var_list    (tm.n_vars, 0),
          share_var_from_zone (tm.n_vars, 1); // We only write data for the first zone, all other
        // zones will share from this one.

        // get the subdomain name from libMesh, if there is one.
        std::string subdomain_name = the_mesh.subdomain_name(*sbd_it);
        std::ostringstream zone_name;
        zone_name << this->zone_title();

        // We will title this
        // "{zone_title()}_{subdomain_name}", or
        // "{zone_title()}_{subdomain_id}", or
        // "{zone_title()}"
        if (subdomain_name.size())
          {
            zone_name << "_";
            zone_name << subdomain_name;
          }
        else if (_subdomain_ids.size() > 1)
          {
            zone_name << "_";
            zone_name << *sbd_it;
          }

        ierr = TECZNE112 (const_cast<char *>(zone_name.str().c_str()),
                          &cell_type,
                          &num_nodes,
                          &num_cells,
                          &num_faces,
                          &i_cell_max,
                          &j_cell_max,
                          &k_cell_max,
                          &_time,
                          &strand_id,
                          &parent_zone,
                          &is_block,
                          &num_face_connect,
                          &face_neighbor_mode,
                          &tot_num_face_nodes,
                          &num_connect_boundary_faces,
                          &tot_num_boundary_connect,
                          &passive_var_list[0],
                          libmesh_nullptr, // = all are node centered
                          (firstzone) ? libmesh_nullptr : &share_var_from_zone[0],
                          &share_connect_from_zone);

        if (ierr)
          libmesh_file_error(fname);

        // Write *all* the data for the first zone, then share it with the others
        if (firstzone)
          {
            int total = cast_int<int>
#ifdef LIBMESH_USE_REAL_NUMBERS
              ((3 + ((solution_names == libmesh_nullptr) ? 0 : solution_names->size()))*num_nodes);
#else
            ((3 + 3*((solution_names == libmesh_nullptr) ? 0 : solution_names->size()))*num_nodes);
#endif


            ierr = TECDAT112 (&total,
                              &tm.nodalData[0],
                              &is_double);

            if (ierr)
              libmesh_file_error(fname);
          }

        // Write the connectivity
        ierr = TECNOD112 (&tm.connData[0]);

        if (ierr)
          libmesh_file_error(fname);
      }

      firstzone = false;
    }

  // Done, close the file.
  ierr = TECEND112 ();

  if (ierr)
    libmesh_file_error(fname);




  //------------------------------------------------------------
  // Legacy binary format
#else

  // Get a constant reference to the mesh.
  const MeshBase & the_mesh = MeshOutput<MeshBase>::mesh();

  // Tecplot binary output only good for dim=2,3
  if (the_mesh.mesh_dimension() == 1)
    {
      this->write_ascii (fname, vec, solution_names);

      return;
    }

  // Required variables
  std::string tecplot_variable_names;
  int is_double =  0,
    tec_debug =  0,
    cell_type = ((the_mesh.mesh_dimension()==2) ? (1) : (3));

  // Build a string containing all the variable names to pass to Tecplot
  {
    tecplot_variable_names += "x, y, z";

    if (solution_names != libmesh_nullptr)
      {
        for (std::size_t name=0; name<solution_names->size(); name++)
          {
#ifdef LIBMESH_USE_REAL_NUMBERS

            tecplot_variable_names += ", ";
            tecplot_variable_names += (*solution_names)[name];

#else

            tecplot_variable_names += ", ";
            tecplot_variable_names += "r_";
            tecplot_variable_names += (*solution_names)[name];
            tecplot_variable_names += ", ";
            tecplot_variable_names += "i_";
            tecplot_variable_names += (*solution_names)[name];
            tecplot_variable_names += ", ";
            tecplot_variable_names += "a_";
            tecplot_variable_names += (*solution_names)[name];

#endif
          }
      }
  }

  // Instantiate a TecplotMacros interface.  In 2D the most nodes per
  // face should be 4, in 3D it's 8.


  TecplotMacros tm(cast_int<unsigned int>(the_mesh.n_nodes()),
                   cast_int<unsigned int>
#ifdef LIBMESH_USE_REAL_NUMBERS
                   (3 + ((solution_names == libmesh_nullptr) ? 0 : solution_names->size())),
#else
                   (3 + 3*((solution_names == libmesh_nullptr) ? 0 : solution_names->size())),
#endif
                   cast_int<unsigned int>
                   (the_mesh.n_active_sub_elem()),
                   ((the_mesh.mesh_dimension() == 2) ? 4 : 8)
                   );


  // Copy the nodes and data to the TecplotMacros class. Note that we store
  // everything as a float here since the eye doesn't require a double to
  // understand what is going on
  for (unsigned int v=0; v<the_mesh.n_nodes(); v++)
    {
      tm.nd(0,v) = static_cast<float>(the_mesh.point(v)(0));
      tm.nd(1,v) = static_cast<float>(the_mesh.point(v)(1));
      tm.nd(2,v) = static_cast<float>(the_mesh.point(v)(2));

      if ((vec != libmesh_nullptr) &&
          (solution_names != libmesh_nullptr))
        {
          const std::size_t n_vars = solution_names->size();

          for (std::size_t c=0; c<n_vars; c++)
            {
#ifdef LIBMESH_USE_REAL_NUMBERS

              tm.nd((3+c),v)     = static_cast<float>((*vec)[v*n_vars + c]);
#else
              tm.nd((3+3*c),v)   = static_cast<float>((*vec)[v*n_vars + c].real());
              tm.nd((3+3*c+1),v) = static_cast<float>((*vec)[v*n_vars + c].imag());
              tm.nd((3+3*c+2),v) = static_cast<float>(std::abs((*vec)[v*n_vars + c]));
#endif
            }
        }
    }


  // Copy the connectivity
  {
    unsigned int te = 0;

    MeshBase::const_element_iterator       it  = the_mesh.active_elements_begin();
    const MeshBase::const_element_iterator end = the_mesh.active_elements_end();

    for ( ; it != end; ++it)
      {
        std::vector<dof_id_type> conn;
        for (unsigned int se=0; se<(*it)->n_sub_elem(); se++)
          {
            (*it)->connectivity(se, TECPLOT, conn);

            for (std::size_t node=0; node<conn.size(); node++)
              tm.cd(node,te) = conn[node];

            te++;
          }
      }
  }


  // Ready to call the Tecplot API
  {
    int ierr = 0,
      num_nodes = static_cast<int>(the_mesh.n_nodes()),
      num_cells = static_cast<int>(the_mesh.n_active_sub_elem());


    ierr = TECINI (libmesh_nullptr,
                   (char *) tecplot_variable_names.c_str(),
                   (char *) fname.c_str(),
                   (char *) ".",
                   &tec_debug,
                   &is_double);

    if (ierr)
      libmesh_file_error(fname);


    ierr = TECZNE (libmesh_nullptr,
                   &num_nodes,
                   &num_cells,
                   &cell_type,
                   (char *) "FEBLOCK",
                   libmesh_nullptr);

    if (ierr)
      libmesh_file_error(fname);


    int total =
#ifdef LIBMESH_USE_REAL_NUMBERS
      ((3 + ((solution_names == libmesh_nullptr) ? 0 : solution_names->size()))*num_nodes);
#else
    ((3 + 3*((solution_names == libmesh_nullptr) ? 0 : solution_names->size()))*num_nodes);
#endif


    ierr = TECDAT (&total,
                   &tm.nodalData[0],
                   &is_double);

    if (ierr)
      libmesh_file_error(fname);

    ierr = TECNOD (&tm.connData[0]);

    if (ierr)
      libmesh_file_error(fname);

    ierr = TECEND ();

    if (ierr)
      libmesh_file_error(fname);
  }

#endif
}

} // namespace libMesh
