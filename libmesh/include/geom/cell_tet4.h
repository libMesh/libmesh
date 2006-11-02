// $Id: cell_tet4.h,v 1.14 2006-11-02 13:50:36 jwpeterson Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __cell_tet4_h__
#define __cell_tet4_h__

// C++ includes

// Local includes
#include "cell_tet.h"




/**
 * The \p Tet4 is an element in 3D composed of 4 nodes.
 * It is numbered like this:
   \verbatim
  TET4:
        3
        o
       /|\
      / | \
     /  |  \
  0 o...|...o 2
     \  |  /
      \ | /
       \|/
        o
        1
  \endverbatim
 */

// ------------------------------------------------------------
// Tet4 class definition
class Tet4 : public Tet
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  Tet4  (Elem* p=NULL);
  
  /**
   * @returns \p TET4
   */
  ElemType type () const { return TET4; }
  
  /**
   * @returns 4
   */
  unsigned int n_nodes() const { return 4; }

  /**
   * @returns 1
   */
  unsigned int n_sub_elem() const { return 1; }
    
  /**
   * @returns true iff the specified (local) node number is a vertex.
   */
  virtual bool is_vertex(const unsigned int i) const;

  /**
   * @returns true iff the specified (local) node number is an edge.
   */
  virtual bool is_edge(const unsigned int i) const;

  /**
   * @returns true iff the specified (local) node number is a face.
   */
  virtual bool is_face(const unsigned int i) const;
  
  /*
   * @returns true iff the specified (local) node number is on the
   * specified side
   */
  virtual bool is_node_on_side(const unsigned int n,
			       const unsigned int s) const;
  
  /*
   * @returns true iff the specified (local) node number is on the
   * specified edge
   */
  virtual bool is_node_on_edge(const unsigned int n,
			       const unsigned int e) const;
  
  /*
   * @returns true iff the element map is definitely affine within
   * numerical tolerances
   */
  virtual bool has_affine_map () const { return true; }

  /**
   * @returns FIRST
   */
  Order default_order() const { return FIRST; }
  
  /**
   * Builds a \p TRI3 built coincident with face i.  
   * The \p AutoPtr<Elem> handles the memory aspect.
   */
  AutoPtr<Elem> build_side (const unsigned int i) const;

  /**
   * Builds a \p EDGE2 built coincident with face i.  
   * The \p AutoPtr<Elem> handles the memory aspect.
   */
  AutoPtr<Elem> build_edge (const unsigned int i) const;

  virtual void connectivity(const unsigned int sc,
			    const IOPackage iop,
			    std::vector<unsigned int>& conn) const;

  /**
   * This maps the \f$ j^{th} \f$ node of the \f$ i^{th} \f$ side to
   * element node numbers.
   */
  static const unsigned int side_nodes_map[4][3];

  /**
   * This maps the \f$ j^{th} \f$ node of the \f$ i^{th} \f$ edge to
   * element node numbers.
   */
  static const unsigned int edge_nodes_map[6][2];
  
  /**
   * An optimized method for computing the area of a
   * 4-node tetrahedron.
   */
  virtual Real volume () const;

  /**
   * Returns the min and max *dihedral* angles for the tetrahedron.
   * Note there are 6 dihedral angles (angles between the planar
   * faces) for the Tet4.  Dihedral angles near 180 deg. are generally
   * bad for interplation.  Small dihedral angles are not necessarily
   * bad for interplation, but they can effect the stiffness matrix
   * condition number.
   */
   std::pair<Real, Real> min_and_max_angle() const;

protected:

  
#ifdef ENABLE_AMR
  
  /**
   * Matrix used to create the elements children.
   */
  float embedding_matrix (const unsigned int i,
			  const unsigned int j,
			  const unsigned int k) const;
  //  { return _embedding_matrix[i][j][k]; }

  /**
   * Matrix that computes new nodal locations/solution values
   * from current nodes/solution.  
   */
  static const float _embedding_matrix[8][4][4];

  /**
   * This enumeration keeps track of which diagonal is selected during
   * refinement.  In general there are three possible diagonals to
   * choose when splitting the octahedron, and by choosing the shortest
   * one we obtain the best element shape.
   */
  enum Diagonal
    {DIAG_02_13=0,    // diagonal between edges (0,2) and (1,3)
     DIAG_03_12=1,    // diagonal between edges (0,3) and (1,2)
     DIAG_01_23=2,    // diagonal between edges (0,1) and (2,3)
     INVALID_DIAG=99  // diagonal not yet selected
    };

  mutable Diagonal _diagonal_selection;
  
#endif
  
};



// ------------------------------------------------------------
// Tet4 class member functions
inline
Tet4::Tet4(Elem* p) :
  Tet(Tet4::n_nodes(), p),
  _diagonal_selection(INVALID_DIAG)
{
}


#endif
