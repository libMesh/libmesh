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



#ifndef LIBMESH_TREE_BASE_H
#define LIBMESH_TREE_BASE_H

// Local includes
#include "libmesh/reference_counted_object.h"
#include "libmesh/id_types.h" // subdomain_id_type
#include "libmesh/libmesh_common.h" // TOLERANCE

// C++ includes
#include <set>
#include <ostream>

namespace libMesh
{


// Forward Declarations
class TreeBase;
class MeshBase;
class Point;
class Elem;


namespace Trees
{
/**
 * \p enum defining how to build the tree.  \p NODES will populate the
 * tree with nodes and then replace the nodes with element
 * connectivity, \p ELEMENTS will populate the tree with the elements
 * directly.  LOCAL_ELEMENTS will populate the tree only with elements
 * from the current processor.  This experimental capability may be
 * useful if you do not wish to include off-processor elements in the
 * search for a Point.
 *
 * \author Daniel Dreyer
 * \date 2003
 * \brief Base class for different Tree types.
 */
enum BuildType : int {NODES=0,
                ELEMENTS,
                LOCAL_ELEMENTS,
                INVALID_BUILD_TYPE };
}

/**
 * This is the base class for trees, it allows pointer
 * usage of trees.
 */
class TreeBase : public ReferenceCountedObject<TreeBase>
{
protected:
  /**
   * Constructor.  Protected.
   */
  explicit
  TreeBase (const MeshBase & m) : mesh(m) {}

public:
  /**
   * Destructor.
   */
  virtual ~TreeBase() = default;

  /**
   * Prints the nodes.
   */
  virtual void print_nodes(std::ostream & out_stream=libMesh::out) const = 0;

  /**
   * Prints the nodes.
   */
  virtual void print_elements(std::ostream & out_stream=libMesh::out) const = 0;

  /**
   * \returns The number of active bins.
   */
  virtual unsigned int n_active_bins() const = 0;

  /**
   * \returns A pointer to the element containing point p,
   * optionally restricted to a set of allowed subdomains,
   * optionally using a non-zero relative tolerance for searches.
   */
  virtual const Elem * find_element(const Point & p,
                                    const std::set<subdomain_id_type> * allowed_subdomains = nullptr,
                                    Real relative_tol = TOLERANCE) const = 0;

  /**
   * Fills \p candidate_elements with any elements containing the
   * specified point \p p,
   * optionally restricted to a set of allowed subdomains,
   * optionally using a non-default relative tolerance for searches.
   */
  virtual void find_elements(const Point & p,
                             std::set<const Elem *> & candidate_elements,
                             const std::set<subdomain_id_type> * allowed_subdomains = nullptr,
                             Real relative_tol = TOLERANCE) const = 0;

protected:

  /**
   * Constant reference to a mesh.  Declared
   * at construction.
   */
  const MeshBase & mesh;
};

} // namespace libMesh


#endif // LIBMESH_TREE_BASE_H
