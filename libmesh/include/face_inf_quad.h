// $Id: face_inf_quad.h,v 1.2 2003-03-26 01:08:15 ddreyer Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __face_inf_quad_h__
#define __face_inf_quad_h__

// C++ includes

// Local includes
#include "mesh_config.h"
#ifdef ENABLE_INFINITE_ELEMENTS

#include "elem.h"


// Forward declarations



/**
 * The \p InfQuad is an abstract element type that lives in
 * two dimensions.  Here, an infinite face is @e always a quadrilateral,
 * so this class is directly derived from \p Elem, without an intermediate
 * \p InfFace class or so.
 * It looks like this:
   \verbatim
                    
                                   closer to infinity
          |           |
          |           |
   side 2 |           | side 1
          |           |
          |           |
           -----------             base side

             side 0
                    
  \endverbatim
 */
// ------------------------------------------------------------
// InfQuad class definition
class InfQuad : public Elem
{
public:

  /**
   * Constructor.  Derived classes implement 'true' elements.
   */
  InfQuad (const unsigned int nn, const Elem* p);

/*   /\** */
/*    * @returns 4.  All quad-derivatives are guaranteed to have at */
/*    * least 4 nodes. */
/*    *\/ */
/*   unsigned int n_nodes() const { return 4; } */

  /**
   * @returns 2, the dimensionality of the object.
   */
  unsigned int dim() const { return 2; }

  /**
   * @returns 2 for the base, 1 otherwise
   */
  unsigned int n_children_per_side(const unsigned int s) const;

  /**
   * @returns 3.  Infinite faces have one side less
   * than their conventional counterparts, since one
   * side is supposed to be located at infinity.
   */
  unsigned int n_sides() const { return 3; }

  /**
   * @returns 4.  All infinite quads (in our setting) have 4 vertices.
   */
  unsigned int n_vertices() const { return 4; }

  /**
   * @returns 3.  All infinite quads have 1 edge in the
   * base, and 2 perpendicular to the base.
   */
  unsigned int n_edges() const { return 6; }

  /**
   * @returns 0.  All 2D elements have no faces, just
   * edges.
   */
  unsigned int n_faces() const { return 0; }
  
  /**
   * @returns 2
   */
  unsigned int n_children() const { return 2; }

  /**
   * @returns a primitive (2-noded) edge or infedge for 
   * edge \p i.
   */
  AutoPtr<Elem> side (const unsigned int i) const;
  
  /**
   * Based on the quality metric \p q specified by the user,
   * returns a quantitative assessment of element quality.
   */
  Real quality (const ElemQuality q) const;

  /**
   * Returns the suggested quality bounds for
   * the hex based on quality measure q.  These are
   * the values suggested by the CUBIT User's Manual.
   */
  std::pair<Real, Real> qual_bounds (const ElemQuality q) const;


#ifdef ENABLE_INFINITE_ELEMENTS

  /**
   * @returns \p true.  All classes derived from \p InfQuad
   * are infinite elements. 
   */
  bool infinite () const { return true; }

#endif


protected:
  
#ifdef ENABLE_AMR
  
  /**
   * Matrix that allows children to inherit boundary conditions.
   */
  unsigned int side_children_matrix (const unsigned int i,
				     const unsigned int j) const
  { return _side_children_matrix[i][j]; }

#endif  



  
private:
  
#ifdef ENABLE_AMR
  
  /**
   * Matrix that tells which children share which of
   * my sides.
   */
  static const unsigned int _side_children_matrix[3][2];
  
#endif    

};



// ------------------------------------------------------------
// InfQuad class member functions
inline
InfQuad::InfQuad(const unsigned int nn, const Elem* p) :
  Elem(nn, InfQuad::n_sides(), p) 
{
}



inline
unsigned int InfQuad::n_children_per_side(const unsigned int s) const
{
  assert (s < this->n_sides());

  switch (s)
  {
    case 0:
      // every infinite face has 2 children in the base edge
      return 2;

    default:
      // on infinite edges only 1 child
      return 1;
  }
}


#endif // ifdef ENABLE_INFINITE_ELEMENTS

#endif




