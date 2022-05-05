// The libMesh Finite Element Library.
// Copyright (C) 2002-2022 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
   *
   * By default \p keep_spline_nodes is true, and when reading a BEXT
   * file we retain the spline nodes and their constraints on FE
   * nodes.  If \p keep_spline_nodes is false, we only keep the FE
   * elements, which are no longer constrained to have higher
   * continuity.
   */
  explicit
  DynaIO (MeshBase & mesh,
          bool keep_spline_nodes = true);

  /**
   * Reads in a mesh in the Dyna format from the ASCII file given by
   * name.
   *
   * \note The user is responsible for calling Mesh::prepare_for_use()
   * after reading the mesh and before using it.
   *
   * \note To safely use DynaIO::add_spline_constraints with a
   * DistributedMesh, currently the user must
   * allow_remote_element_removal(false) and allow_renumbering(false)
   * before the mesh is read.
   *
   * The patch ids defined in the Dyna file are stored as subdomain
   * ids.
   *
   * The spline nodes defined in the Dyna file are added to the mesh
   * with type NodeElem.  The only connection between spline nodes and
   * finite element nodes will be user constraint equations, so using
   * a space-filling-curve partitioner for these meshes might be a
   * good idea.
   */
  virtual void read (const std::string & name) override;

  /**
   * Constrains finite element degrees of freedom in terms of spline
   * degrees of freedom by adding user-defined constraint rows to \p
   * sys
   */
  void add_spline_constraints(DofMap & dof_map,
                              unsigned int sys_num,
                              unsigned int var_num);

  /**
   * Removes any spline nodes (both NodeElem and Node), leaving only
   * the FE mesh generated from those splines.  Also removes node
   * constraints to the now-missing nodes.
   */
  void clear_spline_nodes();

  /**
   * The integer type DYNA uses
   */
  typedef int32_t dyna_int_type;

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
                      dyna_int_type dyna_type_in,
                      dyna_int_type dim_in,
                      dyna_int_type p_in);

    ElementDefinition(ElemType type_in,
                      dyna_int_type dyna_type_in,
                      dyna_int_type dim_in,
                      dyna_int_type p_in,
                      std::vector<unsigned int> && nodes_in);

    ElemType type;
    dyna_int_type dyna_type;
    dyna_int_type dim;
    dyna_int_type p;
    std::vector<unsigned int> nodes;
  };


  /**
   * Finds the ElementDefinition corresponding to a particular element
   * type.
   */
  static const ElementDefinition &
  find_elem_definition(dyna_int_type dyna_elem, int dim, int p);

  static const ElementDefinition &
  find_elem_definition(ElemType libmesh_elem, int dim, int p);

private:
  // Keep track of spline node indexing, so as to enable adding
  // constraint rows easily later.
  std::vector<Node *> spline_node_ptrs;
  std::unordered_map<Node *, Elem *> spline_nodeelem_ptrs;

  /**
   * Whether to keep or eventually discard spline nodes
   */
  bool _keep_spline_nodes;

  /**
   * Implementation of the read() function.  This function
   * is called by the public interface function and implements
   * reading the file.
   */
  void read_mesh (std::istream & in);

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
   * struct which holds a map from LS-DYNA to libMesh element numberings
   * and vice-versa.
   */
  struct ElementMaps
  {
    // Helper function to add a (key, value) pair to both maps
    void add_def(const ElementDefinition & eledef)
    {
      out.emplace(eledef.type, eledef);
      in.emplace(std::make_tuple(eledef.dyna_type, eledef.dim, eledef.p), eledef);
    }

    std::map<ElemType, ElementDefinition> out;

    typedef std::tuple<dyna_int_type, dyna_int_type, dyna_int_type> Key;

    std::map<Key, ElementDefinition> in;
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
