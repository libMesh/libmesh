// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_GMSH_IO_H
#define LIBMESH_GMSH_IO_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/mesh_input.h"
#include "libmesh/mesh_output.h"

// C++ includes
#include <cstddef>

namespace libMesh
{

// Forward declarations
template <typename> class MeshBaseTempl;
typedef MeshBaseTempl<Real> MeshBase;



/**
 * \brief Reading and writing meshes in the Gmsh format.
 *
 * For a full description of the Gmsh format and to obtain the
 * GMSH software see
 * <a href="http://http://www.geuz.org/gmsh/">the Gmsh home page</a>
 *
 * \author John W. Peterson
 * \date 2004, 2014
 * \author Martin Luthi
 * \date 2005
 */
class GmshIO : public MeshInput<MeshBase>,
               public MeshOutput<MeshBase>
{
public:

  /**
   * Constructor.  Takes a non-const Mesh reference which it
   * will fill up with elements via the read() command.
   */
  explicit
  GmshIO (MeshBase & mesh);

  /**
   * Constructor.  Takes a reference to a constant mesh object.
   * This constructor will only allow us to write the mesh.
   */
  explicit
  GmshIO (const MeshBase & mesh);

  /**
   * Reads in a mesh in the Gmsh *.msh format from the ASCII file
   * given by name.
   *
   * \note The user is responsible for calling Mesh::prepare_for_use()
   * after reading the mesh and before using it.
   *
   * The physical group names defined in the Gmsh-file are stored, depending
   * on the element dimension, as either subdomain name or as side name. The
   * IDs of the former can be retrieved by using the MeshBase::get_id_by_name
   * method; the IDs of the latter can be retrieved by using the
   * MeshBase::get_boundary_info and the BoundaryInfo::get_id_by_name methods.
   */
  virtual void read (const std::string & name) override;

  /**
   * This method implements writing a mesh to a specified file
   * in the Gmsh *.msh format.
   */
  virtual void write (const std::string & name) override;

  /**
   * Bring in base class functionality for name resolution and to
   * avoid warnings about hidden overloaded virtual functions.
   */
  using MeshOutput<MeshBase>::write_nodal_data;

  /**
   * This method implements writing a mesh with nodal data to a
   * specified file where the nodal data and variable names are provided.
   */
  virtual void write_nodal_data (const std::string &,
                                 const std::vector<Number> &,
                                 const std::vector<std::string> &) override;

  /**
   * Flag indicating whether or not to write a binary file.  While binary
   * files may end up being smaller than equivalent ASCII files, they will
   * almost certainly take longer to write.  The reason for this is that
   * the ostream::write() function which is used to write "binary" data to
   * streams, only takes a pointer to char as its first argument.  This means
   * if you want to write anything other than a buffer of chars, you first
   * have to use a strange memcpy hack to get the data into the desired format.
   * See the templated to_binary_stream() function below.
   */
  bool & binary ();

  /**
   * Access to the flag which controls whether boundary elements are
   * written to the Mesh file.
   */
  bool & write_lower_dimensional_elements ();

private:
  /**
   * Implementation of the read() function.  This function
   * is called by the public interface function and implements
   * reading the file.
   */
  void read_mesh (std::istream & in);

  /**
   * This method implements writing a mesh to a
   * specified file.  This will write an ASCII *.msh file.
   */
  void write_mesh (std::ostream & out);

  /**
   * This method implements writing a mesh with nodal data to a specified file
   * where the nodal data and variable names are optionally provided.  This
   * will write an ASCII or binary *.pos file, depending on the binary flag.
   */
  void write_post (const std::string &,
                   const std::vector<Number> * = nullptr,
                   const std::vector<std::string> * = nullptr);

  /**
   * Flag to write binary data.
   */
  bool _binary;

  /**
   * If true, lower-dimensional elements based on the boundary
   * conditions get written to the output file.
   */
  bool _write_lower_dimensional_elements;

  /**
   * Defines mapping from libMesh element types to Gmsh element types or vice-versa.
   */
  struct ElementDefinition
  {
    ElementDefinition(ElemType type_in,
                      unsigned gmsh_type_in,
                      unsigned dim_in,
                      unsigned nnodes_in) :
      type(type_in),
      gmsh_type(gmsh_type_in),
      dim(dim_in),
      nnodes(nnodes_in)
    {}

    ElemType type;
    unsigned int gmsh_type;
    unsigned int dim;
    unsigned int nnodes;
    std::vector<unsigned int> nodes;
  };

  /**
   * struct which holds a map from Gmsh to libMesh element numberings
   * and vice-versa.
   */
  struct ElementMaps
  {
    // Helper function to add a (key, value) pair to both maps
    void add_def(const ElementDefinition & eledef)
    {
      out.insert(std::make_pair(eledef.type, eledef));
      in.insert(std::make_pair(eledef.gmsh_type, eledef));
    }

    std::map<ElemType, ElementDefinition> out;
    std::map<unsigned int, ElementDefinition> in;
  };

  /**
   * A static ElementMaps object that is built statically and used by
   * all instances of this class.
   */
  static ElementMaps _element_maps;

  /**
   * A static function used to construct the _element_maps struct,
   * statically.
   */
  static ElementMaps build_element_maps();
};


} // namespace libMesh

#endif // LIBMESH_GMSH_IO_H
