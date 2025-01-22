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

// Local includes
#include "libmesh/stl_io.h"

#include "libmesh/distributed_mesh.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/face_tri3.h"
#include "libmesh/libmesh_config.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/mesh_base.h"
#include "libmesh/point.h"
#include "libmesh/utility.h"

// gzstream for reading compressed files as a stream
#ifdef LIBMESH_HAVE_GZSTREAM
# include "libmesh/ignore_warnings.h" // shadowing in gzstream.h
# include "gzstream.h"
# include "libmesh/restore_warnings.h"
#endif

// C++ includes
#include <iomanip>
#include <fstream>
#include <vector>
#include <cctype> // isspace

#ifdef LIBMESH_HAVE_CXX11_REGEX
#include <regex>
#endif

namespace {

void test_and_add(std::unique_ptr<libMesh::Elem> e,
                  libMesh::MeshBase & mesh)
{
  const libMesh::Real volume = e->volume();
  if (volume < libMesh::TOLERANCE * libMesh::TOLERANCE)
    libmesh_warning
      ("Warning: STL file contained sliver element with volume " <<
       volume << "!\n");
  mesh.add_elem(std::move(e));
}

}


namespace libMesh
{

STLIO::STLIO (const MeshBase & mesh) :
  MeshOutput<MeshBase>    (mesh),
  _subdivide_second_order (true)
{
}



STLIO::STLIO (MeshBase & mesh) :
  MeshInput<MeshBase> (mesh),
  MeshOutput<MeshBase>(mesh),
  _subdivide_second_order (true)
{
}



void STLIO::write (const std::string & fname)
{
  LOG_SCOPE("write()", "STLIO");

  // Get a reference to the mesh
  const MeshBase * mesh_ptr = &MeshOutput<MeshBase>::mesh();

  libmesh_parallel_only(mesh_ptr->comm());

  // If necessary, serialize to proc 0, which will do all the writing
  std::unique_ptr<DistributedMesh> mesh_copy;
  std::unique_ptr<MeshSerializer> serializer;
  if (!mesh_ptr->is_serial())
    {
      mesh_copy = std::make_unique<DistributedMesh>(*mesh_ptr);
      mesh_ptr = mesh_copy.get();
      serializer = std::make_unique<MeshSerializer>(*mesh_copy, true, true);
    }

  const MeshBase & mesh = *mesh_ptr;

  if (mesh.processor_id() != 0)
    return;

  // Open the output file stream
  std::ofstream out_stream (fname.c_str());

  out_stream << std::setprecision(this->ascii_precision());

  // Make sure it opened correctly
  if (!out_stream.good())
    libmesh_file_error(fname.c_str());

  // Write the header
  out_stream << "solid " << this->_name << '\n';

  // Write all our triangles
  // scream if we see a non-triangle
  for (auto elem : mesh.element_ptr_range())
    {
      // Get to this later
      if (elem->type() == TRI6 || elem->type() == TRI7)
        {
          if (this->_subdivide_second_order)
            libmesh_not_implemented();

          libmesh_error_msg("Tried to write a non-linear triangle to an STL file");
        }

      libmesh_error_msg_if(elem->type() != TRI3,
                           "Tried to write a non-triangle to an STL file");

      auto n = (elem->point(1)-elem->point(0)).cross(elem->point(2)-elem->point(0));

      // Other STL files have slivers, I guess ours can too
      if (auto length = n.norm())
        n /= length;

      out_stream << "facet normal " <<
        n(0) << ' ' <<
        n(1) << ' ' <<
        n(2) << "\n"
        "    outer loop\n"
        "        vertex " <<
        elem->point(0)(0) << ' ' <<
        elem->point(0)(1) << ' ' <<
        elem->point(0)(2) << "\n"
        "        vertex " <<
        elem->point(1)(0) << ' ' <<
        elem->point(1)(1) << ' ' <<
        elem->point(1)(2) << "\n"
        "        vertex " <<
        elem->point(2)(0) << ' ' <<
        elem->point(2)(1) << ' ' <<
        elem->point(2)(2) << "\n"
        "    endloop\n"
        "endfacet\n";
    }

  out_stream << "endsolid" << std::endl;
}



void STLIO::read (const std::string & filename)
{
  LOG_SCOPE("read()", "STLIO");

  // This is a serial-only process for now;
  // the Mesh should be read on processor 0 and
  // broadcast later
  libmesh_assert_equal_to (MeshOutput<MeshBase>::mesh().processor_id(), 0);

  // Clear the mesh so we are sure to start from a pristine state.
  MeshBase & mesh = MeshInput<MeshBase>::mesh();
  mesh.clear();
  mesh.set_mesh_dimension(2);

  std::unique_ptr<std::istream> fstream =
    this->open_file(filename);

  char c;
  const char * expected_header = "solid!";
  bool is_ascii_stl = false,
       is_binary_stl = false;
  while (fstream->get(c))
    {
      if (c == ' ')
        continue;
      if (c == *expected_header)
        {
          ++expected_header;
          if (*expected_header == '!')
          {
            is_ascii_stl = true;
            break;
          }
        }
      else
        {
          is_binary_stl = true; // probably
          break;
        }
    }

  if (is_ascii_stl)
    {
      fstream->seekg(0);
      this->read_ascii(*fstream);
    }
  else if (is_binary_stl)
    {
      fstream->seekg(0, std::ios_base::end);
      std::size_t length = fstream->tellg();
      fstream->seekg(0);
      this->read_binary(*fstream, length);
    }
  else
    libmesh_error_msg("Failed to read an STL header in " << filename);
}


std::unique_ptr<std::istream> STLIO::open_file (const std::string & filename)
{
  std::string_view basename = Utility::basename_of(filename);
  const bool gzipped_file = (basename.rfind(".gz") == basename.size() - 3);

  std::unique_ptr<std::istream> file;

  if (gzipped_file)
    {
#ifdef LIBMESH_HAVE_GZSTREAM
      auto inf = std::make_unique<igzstream>();
      libmesh_assert(inf);
      inf->open(filename.c_str(), std::ios::in);
      file = std::move(inf);
#else
      libmesh_error_msg("ERROR: need gzstream to handle .gz files!!!");
#endif
    }
  else
    {
      auto inf = std::make_unique<std::ifstream>();
      libmesh_assert(inf);

      std::string new_name = Utility::unzip_file(filename);

      inf->open(new_name.c_str(), std::ios::in);
      file = std::move(inf);
    }

  return file;
}



void STLIO::read_ascii (std::istream & file)
{
  LOG_SCOPE("read_ascii()", "STLIO");

#ifndef LIBMESH_HAVE_CXX11_REGEX
  libmesh_not_implemented();  // What is your compiler?!?  Email us!
  libmesh_ignore(file);
#else

  MeshBase & mesh = MeshInput<MeshBase>::mesh();

  const std::regex all_expected_chars
    ("^[\\w\\d\\-\\.\\^\\$]*$");

  const std::regex start_regex
    ("^\\s*solid");
  const std::regex start_with_name_regex
    ("^\\s*solid\\s+(\\w+)");

  const std::regex start_facet_regex
    ("^\\s*facet");

  // We'll ignore facets' normals for now
  /*
  const std::regex facet_with_normal_regex
    ("^\\s*facet\\s+normal"
     "\\s+([+-]?(\\d+([.]\\d*)?([eE][+-]?\\d+)?|[.]\\d+([eE][+-]?\\d+)?)"
     "\\s+[+-]?(\\d+([.]\\d*)?([eE][+-]?\\d+)?|[.]\\d+([eE][+-]?\\d+)?)"
     "\\s+[+-]?(\\d+([.]\\d*)?([eE][+-]?\\d+)?|[.]\\d+([eE][+-]?\\d+)?))");
  */

  const std::regex start_loop_regex
    ("^\\s*outer\\s+loop");

  const std::regex vertex_regex
    ("^\\s*vertex"
     "\\s+([+-]?(\\d+([.]\\d*)?([eE][+-]?\\d+)?|[.]\\d+([eE][+-]?\\d+)?)"
     "\\s+[+-]?(\\d+([.]\\d*)?([eE][+-]?\\d+)?|[.]\\d+([eE][+-]?\\d+)?)"
     "\\s+[+-]?(\\d+([.]\\d*)?([eE][+-]?\\d+)?|[.]\\d+([eE][+-]?\\d+)?))");

  const std::regex end_facet_regex
    ("^\\s*end\\s*facet");

  const std::regex end_loop_regex
    ("^\\s*end\\s*loop");

  bool have_started = false;
  bool in_facet = false;
  bool in_vertex_loop = false;
  std::unique_ptr<Tri3> triangle;
  int next_vertex = 0;
  int line_num = 0;

  std::unordered_map<Point, Node *> mesh_points;

  for (std::string line; std::getline(file, line);)
    {
      ++line_num;
      std::smatch sm;

      for (char & c : line)
        {
          libmesh_error_msg_if
            (c != '\n' && c != '\r' &&
             (c < ' ' || c > '~'),
             "Found non-ASCII character " << int(c) << " on line " <<
             line_num << "\nSTLIO does not yet support binary files.");
        }

      if (std::regex_search(line, sm, start_regex))
        {
          if (std::regex_search(line, sm, start_with_name_regex))
            this->_name = sm[1];
          libmesh_error_msg_if
            (have_started,
             "Found two 'solid' lines starting the same STL file?");
          have_started = true;
        }

      else if (std::regex_search(line, sm, start_facet_regex))
        {
          libmesh_error_msg_if
            (in_facet,
             "Found a repeated 'facet' line with no 'endfacet' between them.");
          in_facet = true;

          /*
          if (std::regex_search(line, sm, facet_with_normal_regex))
            {
              const std::string normalvec = sm[1];

              // Try to be compatible with higher-than-double T; maybe
              // *someone* out there is doing STL in 128 bits? ...
              std::stringstream ss(normalvec);
              ss >> normal(0);
              ss >> normal(1);
              ss >> normal(2);
              }
          */
        }

      else if (std::regex_search(line, end_facet_regex))
        {
          libmesh_error_msg_if
            (!in_facet,
             "Found an 'endfacet' line with no matching 'facet'.");
          in_facet = false;
        }

      else if (std::regex_search(line, start_loop_regex))
        {
          libmesh_error_msg_if
            (!in_facet,
             "Found an 'outer loop' line with no matching 'facet'.");
          libmesh_error_msg_if
            (in_vertex_loop,
             "Found a repeated 'outer loop' line with no 'endloop' between them.");
          in_vertex_loop = true;
          triangle = std::make_unique<Tri3>();
          next_vertex = 0;
        }

      else if (std::regex_search(line, sm, vertex_regex))
        {
          const std::string normalvec = sm[1];

          // Try to be compatible with higher-than-double-precision T;
          // maybe *someone* out there is doing STL in 128 bits? ...
          std::stringstream ss(normalvec);
          Point p;
          ss >> p(0);
          ss >> p(1);
          ss >> p(2);

          Node * node;
          if (auto it = mesh_points.find(p); it != mesh_points.end())
            {
              node = it->second;
            }
          else
            {
              node = mesh.add_point(p);
              mesh_points[p] = node;
            }

          libmesh_error_msg_if
            (next_vertex > 2,
             "Found more than 3 vertices in a loop; STLIO only supports Tri3.");
          triangle->set_node(next_vertex++) = node;
        }

      else if (std::regex_search(line, end_loop_regex))
        {
          libmesh_error_msg_if
            (next_vertex != 3,
             "Found an 'endloop' line after only seeing " << next_vertex << " vertices in 'loop'.");
          libmesh_error_msg_if
            (!in_vertex_loop,
             "Found an 'endloop' line with no matching 'loop'.");
          in_vertex_loop = false;

          test_and_add(std::move(triangle), mesh);
        }
    }

  libmesh_error_msg_if
    (in_facet,
     "File ended without ending a facet first.");

  libmesh_error_msg_if
    (in_vertex_loop,
     "File ended without ending an outer loop first.");
#endif // LIBMESH_HAVE_CXX11_REGEX
}



void STLIO::read_binary (std::istream & file,
                         std::size_t input_size)
{
  LOG_SCOPE("read_binary()", "STLIO");

  MeshBase & mesh = MeshInput<MeshBase>::mesh();

  char header_buffer[80];

  // 80-character header which is generally ignored - Wikipedia
  file.read(header_buffer, 80);

  // What's our endianness here?  Input binary files are specified to
  // be little endian, and we might need to do conversions.

  uint32_t test_int = 0x87654321;
  const bool big_endian = ((*reinterpret_cast<char *>(&test_int)) != 0x21);
  const Utility::ReverseBytes endian_fix{big_endian};

  // C++ doesn't specify a size for float, but we really need 4-byte
  // floats here.  Fortunately basically every implementation ever
  // uses exactly 4 bytes for float.
  if constexpr (sizeof(float) != 4)
    libmesh_error_msg("Trying to read 4 byte floats without 4-byte float?");

  uint32_t n_elem = 0;

  file.read(reinterpret_cast<char*>(&n_elem), 4);
  endian_fix(n_elem);

  libmesh_error_msg_if
    (input_size &&
     (input_size < 84+n_elem*50),
     "Not enough data for " << n_elem << " STL triangles in " <<
     input_size << " uncompressed bytes.");

  std::unique_ptr<Tri3> triangle;
  std::unordered_map<Point, Node *> mesh_points;

  for (unsigned int e : make_range(n_elem))
  {
    libmesh_ignore(e);

    triangle = std::make_unique<Tri3>();

    // We'll ignore facets' normals for now
    char ignored_buffer[12];
    file.read(ignored_buffer, 12);

    // Read vertex locations
    for (int i : make_range(3))
      {
        float point_buffer[3];
        file.read(reinterpret_cast<char*>(point_buffer), 12);
        endian_fix(point_buffer[0]);
        endian_fix(point_buffer[1]);
        endian_fix(point_buffer[2]);
        const Point p {point_buffer[0],
                       point_buffer[1],
                       point_buffer[2]};

        Node * node;
        if (auto it = mesh_points.find(p); it != mesh_points.end())
          {
            node = it->second;
          }
        else
          {
            node = mesh.add_point(p);
            mesh_points[p] = node;
          }

        triangle->set_node(i) = node;
      }

    // The 2-byte "attribute byte count" is unstandardized; typically
    // 0, or sometimes triangle color.  Ignore it.
    file.read(ignored_buffer, 2);

    test_and_add(std::move(triangle), mesh);
  }
}


} // namespace libMesh
