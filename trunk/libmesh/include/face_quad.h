// $Id: face_quad.h,v 1.4 2003-01-21 19:24:34 benkirk Exp $

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



#ifndef __quad_h__
#define __quad_h__


// C++ includes


// Local includes
#include "mesh_common.h"
#include "face.h"


// Forward declarations


/**
 * The \p QUAD is an element in 2D composed of 4 sides.
 * It looks like this:
 * \verbatim
 *                   
 *         ----------- 
 *        |           |
 *        |           |
 *        |           | 
 *        |           |
 *        |           |
 *         ----------- 
 *                   
 * \endverbatim
 */

// ------------------------------------------------------------
// Quad class definition
class Quad : public Face
{
public:

  /**
   * Default quadrilateral element, takes number of nodes and 
   * parent. Derived classes implement 'true' elements.
   */
  Quad  (const unsigned int nn, Face* p);
  
  /**
   * @returns 4.  All quad-derivatives are guaranteed to have at
   * least 4 nodes.
   */
  unsigned int n_nodes() const { return 4; };

  /**
   * @returns 4
   */
  unsigned int n_sides() const { return 4; };  

  /**
   * @returns 4.  All quadrilaterals have 4 vertices.
   */
  unsigned int n_vertices() const { return 4; };

  /**
   * @returns 4.  All quadrilaterals have 4 edges.
   */
  unsigned int n_edges() const { return 4; };

  /**
   * @returns 4
   */
  unsigned int n_children() const { return 4; };

  /**
   * @returns a primitive (2-noded) edge for 
   * edge i.
   */
  AutoPtr<Elem> side (const unsigned int i) const;
  
  /**
   * Based on the quality metric q specified by the user,
   * returns a quantitative assessment of element quality.
   */
  real quality (const ElemQuality q) const;

  /**
   * Returns the suggested quality bounds for
   * the hex based on quality measure q.  These are
   * the values suggested by the CUBIT User's Manual.
   */
  std::pair<real, real> qual_bounds (const ElemQuality q) const;

private:
    
};


// ------------------------------------------------------------
// Quad class member functions
inline
Quad::Quad(const unsigned int nn, Face* p) :
  Face(nn, Quad::n_sides(), p)
{
};

#endif
