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



#ifndef LIBMESH_FACE_TRI_H
#define LIBMESH_FACE_TRI_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/face.h"

namespace libMesh
{

/**
 * The \p Tri is an element in 2D composed of 3 sides.
 * It looks like this:
 * \verbatim
 *          ^
 *         / \
 *        /   \
 *       /     \
 *      /       \
 *     /         \
 *     -----------
 * \endverbatim
 *
 * \author Benjamin S. Kirk
 * \date 2002
 * \brief The base class for all triangular element types.
 */
class Tri : public Face
{
public:

  /**
   * Default triangular element, takes number of nodes and
   * parent. Derived classes implement 'true' elements.
   */
  Tri (const unsigned int nn,
       Elem * p,
       Node ** nodelinkdata) :
    Face(nn, Tri::n_sides(), p, _elemlinks_data, nodelinkdata)
  {
    // Make sure the interior parent isn't undefined
    if (LIBMESH_DIM > 2)
      this->set_interior_parent(nullptr);
  }

  Tri (Tri &&) = delete;
  Tri (const Tri &) = delete;
  Tri & operator= (const Tri &) = delete;
  Tri & operator= (Tri &&) = delete;
  virtual ~Tri() = default;

  /**
   * \returns The \p Point associated with local \p Node \p i,
   * in master element rather than physical coordinates.
   */
  virtual Point master_point (const unsigned int i) const override final
  {
    libmesh_assert_less(i, this->n_nodes());
    return Point(_master_points[i][0],
                 _master_points[i][1],
                 _master_points[i][2]);
  }

  /**
   * \returns 3.  All tri-derivatives are guaranteed to have at
   * least 3 nodes.
   */
  virtual unsigned int n_nodes() const override { return 3; }

  /**
   * \returns 3.
   */
  virtual unsigned int n_sides() const override final { return 3; }

  /**
   * \returns 3.  All triangles have 3 vertices.
   */
  virtual unsigned int n_vertices() const override final { return 3; }

  /**
   * \returns 3.  All triangles have 3 edges.
   */
  virtual unsigned int n_edges() const override final { return 3; }

  /**
   * \returns 4.
   */
  virtual unsigned int n_children() const override final { return 4; }

  /**
   * \returns \p true if the specified child is on the
   * specified side.
   */
  virtual bool is_child_on_side(const unsigned int c,
                                const unsigned int s) const override final;

  /**
   * Don't hide Elem::key() defined in the base class.
   */
  using Elem::key;

  /**
   * \returns An id associated with the \p s side of this element.
   * The id is not necessarily unique, but should be close.  This is
   * particularly useful in the \p MeshBase::find_neighbors() routine.
   */
  virtual dof_id_type key (const unsigned int s) const override;

  /**
   * \returns An id associated with the global node ids of this
   * element.  The id is not necessarily unique, but should be
   * close.
   */
  virtual dof_id_type key () const override final;

  /**
   * \returns \p Tri3::side_nodes_map[side][side_node] after doing some range checking.
   */
  virtual unsigned int local_side_node(unsigned int side,
                                       unsigned int side_node) const override;

  /**
   * Calls local_side_node(edge, edge_node). For 2D elements, there is an implied
   * equivalence between edges and sides, e.g. n_edges() == n_sides(), so we treat
   * these two functions the same.
   */
  virtual unsigned int local_edge_node(unsigned int edge,
                                       unsigned int edge_node) const override;

  /**
   * \returns A primitive (2-noded) edge for edge i.
   */
  virtual std::unique_ptr<Elem> side_ptr (const unsigned int i) override final;

  /**
   * Rebuilds an EDGE2 coincident with face i.
   */
  virtual void side_ptr (std::unique_ptr<Elem> & elem,
                         const unsigned int i) override final;

  /**
   * \returns A quantitative assessment of element quality based on
   * the quality metric \p q specified by the user.
   */
  virtual Real quality (const ElemQuality q) const override;

  /**
   * \returns The suggested quality bounds for the hex based on quality
   * measure \p q.  These are the values suggested by the CUBIT User's
   * Manual.
   */
  virtual std::pair<Real, Real> qual_bounds (const ElemQuality q) const override;

  /**
   * Three sides, one orientation each.
   */
  virtual unsigned int n_permutations() const override final { return 3; }

protected:

  /**
   * Data for links to parent/neighbor/interior_parent elements.
   */
  Elem * _elemlinks_data[4+(LIBMESH_DIM>2)];

  /**
   * Master element node locations
   */
  static const Real _master_points[6][3];
};

} // namespace libMesh

#endif // LIBMESH_FACE_TRI_H
