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



#ifndef LIBMESH_DYNA_IO_H
#define LIBMESH_DYNA_IO_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/mesh_input.h"
#include "libmesh/mesh_output.h"

// C++ includes
#include <algorithm>
#include <cstddef>
#include <map>
#include <vector>

namespace libMesh
{

// Forward declarations
class MeshBase;



/**
 * \brief Reading and writing meshes in (a subset of) LS-DYNA format.
 *
 * The initial implementation only handles cards in the format
 * described in "Geometry import to LS-DYNA", for isogeometric
 * analysis.
 *
 * \author Roy H. Stogner
 * \date 2019
 */
class DynaIO : public MeshInput<MeshBase>
{
public:

  /**
   * Constructor.  Takes a non-const Mesh reference which it
   * will fill up with elements via the read() command.
   */
  explicit
  DynaIO (MeshBase & mesh);

  /**
   * Reads in a mesh in the Dyna format from the ASCII file given by
   * name.
   *
   * \note The user is responsible for calling Mesh::prepare_for_use()
   * after reading the mesh and before using it.
   *
   * The patch ids defined in the Dyna file are stored as subdomain
   * ids.
   */
  virtual void read (const std::string & name) override;

private:
  /**
   * Implementation of the read() function.  This function
   * is called by the public interface function and implements
   * reading the file.
   */
  void read_mesh (std::istream & in);

  /**
   * The integer type DYNA uses
   */
  typedef int32_t dyna_int_type;

  /**
   * How many can we find on a line?
   */
  static const int max_ints_per_line = 10;

  /**
   * The floating-point type DYNA uses
   */
  typedef double dyna_fp_type;

  /**
   * How many can we find on a line?
   */
  static const int max_fps_per_line = 5;

  /**
   * Defines mapping from libMesh element types to LS-DYNA element
   * types or vice-versa.
   *
   * For the foreseeable future only isotropic p elements, with the
   * same polynomial degree in every direction, are supported.
   */
  struct ElementDefinition
  {
    ElementDefinition(ElemType type_in,
                      unsigned dyna_type_in,
                      unsigned dim_in,
                      unsigned p_in) :
      type(type_in),
      dyna_type(dyna_type_in),
      dim(dim_in),
      p(p_in)
    {
      const unsigned int n_nodes = Elem::type_to_n_nodes_map[type_in];
      nodes.resize(n_nodes);
      std::iota(nodes.begin(), nodes.end(), 0);
    }

    ElementDefinition(ElemType type_in,
                      unsigned dyna_type_in,
                      unsigned dim_in,
                      unsigned p_in,
                      std::vector<unsigned int> && nodes_in) :
      type(type_in),
      dyna_type(dyna_type_in),
      dim(dim_in),
      p(p_in),
      nodes(nodes_in)
    {}

    ElemType type;
    unsigned int dyna_type;
    unsigned int dim;
    unsigned int p;
    std::vector<unsigned int> nodes;
  };

  /**
   * struct which holds a map from LS-DYNA to libMesh element numberings
   * and vice-versa.
   */
  struct ElementMaps
  {
    // Helper function to add a (key, value) pair to both maps
    void add_def(const ElementDefinition & eledef)
    {
      out.insert(std::make_pair(eledef.type, eledef));
      in.insert(std::make_pair(std::make_tuple(eledef.dyna_type, eledef.dim, eledef.p), eledef));
    }

    std::map<ElemType, ElementDefinition> out;
    std::map<std::tuple<unsigned int, unsigned int, unsigned int>, ElementDefinition> in;
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

#endif // LIBMESH_DYNA_IO_H
