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



#ifndef LIBMESH_CELL_TET_H
#define LIBMESH_CELL_TET_H

// Local includes
#include "libmesh/cell.h"

namespace libMesh
{

/**
 * The \p Tet is an element in 3D composed of 4 sides.
 *
 * \author Benjamin S. Kirk
 * \date 2002
 * \brief The base class for all tetrahedral element types.
 */
template <typename RealType = Real>
class TetTempl : public CellTempl<RealType>
{
public:
  typedef TetTempl<RealType> Tet;
  typedef CellTempl<RealType> Cell;
  typedef ElemTempl<RealType> Elem;
  typedef PointTempl<RealType> Point;
  typedef NodeTempl<RealType> Node;

  /**
   * Default tetrahedral element, takes number of nodes and
   * parent. Derived classes implement 'true' elements.
   */
  TetTempl (const unsigned int nn, Elem * p, Node ** nodelinkdata) :
    Cell(nn, Tet::n_sides(), p, _elemlinks_data, nodelinkdata),
    _diagonal_selection(INVALID_DIAG)
  {
    // Make sure the interior parent isn't undefined
    if (LIBMESH_DIM > 3)
      this->set_interior_parent(nullptr);
  }

  TetTempl (Tet &&) = delete;
  TetTempl (const Tet &) = delete;
  Tet & operator= (const Tet &) = delete;
  Tet & operator= (Tet &&) = delete;
  virtual ~TetTempl() = default;

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
   * \returns 4.
   */
  virtual unsigned int n_sides() const override final { return 4; }

  /**
   * \returns 4.  All tetrahedra have 4 vertices.
   */
  virtual unsigned int n_vertices() const override final { return 4; }

  /**
   * \returns 6.  All tetrahedra have 6 edges.
   */
  virtual unsigned int n_edges() const override final { return 6; }

  /**
   * \returns 4.  All tetrahedra have 4 faces.
   */
  virtual unsigned int n_faces() const override final { return 4; }

  /**
   * \returns 8.
   */
  virtual unsigned int n_children() const override final { return 8; }

  /**
   * \returns \p true if the specified edge is on the specified side.
   */
  virtual bool is_edge_on_side(const unsigned int e,
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
   * \returns \p Tet4::side_nodes_map[side][side_node] after doing some range checking.
   */
  virtual unsigned int which_node_am_i(unsigned int side,
                                       unsigned int side_node) const override;

  /**
   * \returns A primitive (3-noded) triangle for face i.
   */
  virtual std::unique_ptr<Elem> side_ptr (const unsigned int i) override final;

  /**
   * Rebuilds a primitive (3-noded) triangle for face i.
   */
  virtual void side_ptr (std::unique_ptr<Elem> & side, const unsigned int i) override final;

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
   * This enumeration keeps track of which diagonal is selected during
   * refinement.  In general there are three possible diagonals to
   * choose when splitting the octahedron, and by choosing the shortest
   * one we obtain the best element shape.
   */
  enum Diagonal
    {
      DIAG_02_13=0,    // diagonal between edges (0,2) and (1,3)
      DIAG_03_12=1,    // diagonal between edges (0,3) and (1,2)
      DIAG_01_23=2,    // diagonal between edges (0,1) and (2,3)
      INVALID_DIAG=99  // diagonal not yet selected
    };

  /**
   * \returns The diagonal that has been selected during refinement.
   */
  Diagonal diagonal_selection () const { return _diagonal_selection; }

  /**
   * Allows the user to select the diagonal for the refinement.  This
   * function may only be called before the element is ever refined.
   */
  void select_diagonal (const Diagonal diag) const;



#ifdef LIBMESH_ENABLE_AMR


  /**
   * Tetrahedral elements permute the embedding matrix depending on which
   * interior diagonal is used to subdivide into child elements.
   * But we want to cache topology data based on that matrix.  So we return a
   * "version number" based on the diagonal selection.
   */
  virtual unsigned int embedding_matrix_version () const override final
  {
    this->choose_diagonal();
    return this->diagonal_selection();
  }

#endif // LIBMESH_ENABLE_AMR



protected:

  /**
   * Data for links to parent/neighbor/interior_parent elements.
   */
  Elem * _elemlinks_data[5+(LIBMESH_DIM>3)];

  /**
   * Master element node locations
   */
  static const Real _master_points[10][3];

  /**
   * Called by descendant classes with appropriate data to determine
   * if child c is on side s.  Only works if LIBMESH_ENABLE_AMR.
   */
  bool is_child_on_side_helper(const unsigned int c,
                               const unsigned int s,
                               const unsigned int checked_nodes[][3] ) const;

  /**
   * The currently-selected diagonal used during refinement.
   * Initialized to INVALID_DIAG.
   */
  mutable Diagonal _diagonal_selection;

  /**
   * Derived classes use this function to select an initial
   * diagonal during refinement. The optimal choice is the shortest
   * of the three.
   */
  void choose_diagonal() const;
};

// ------------------------------------------------------------
// Tet class static member initializations


// We need to require C++11...
template <typename RealType>
const Real TetTempl<RealType>::_master_points[10][3] =
  {
    {0, 0, 0},
    {1, 0, 0},
    {0, 1, 0},
    {0, 0, 1},
    {0.5, 0, 0},
    {0.5, 0.5, 0},
    {0, 0.5, 0},
    {0, 0, 0.5},
    {0.5, 0, 0.5},
    {0, 0.5, 0.5}
  };

typedef TetTempl<Real> Tet;

} // namespace libMesh

#endif // LIBMESH_CELL_TET_H
