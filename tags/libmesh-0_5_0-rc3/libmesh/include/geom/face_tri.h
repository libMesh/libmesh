// $Id: face_tri.h,v 1.6 2005-06-06 16:23:56 knezed01 Exp $

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



#ifndef __tri_h__
#define __tri_h__


// C++ includes


// Local includes
#include "libmesh_common.h"
#include "face.h"


// Forward declarations





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
 *               
 * \endverbatim
 */

// ------------------------------------------------------------
// Tri class definition
class Tri : public Face
{
public:

  /**
   * Default triangular element, takes number of nodes and 
   * parent. Derived classes implement 'true' elements.
   */
  Tri (const unsigned int nn,
       const Elem* p) :
    Face(nn, Tri::n_sides(), p) {}
 
  /**
   * Constructor.  Explicitly specifies the number of
   * nodes and neighbors for which storage will be allocated.
   */
  Tri (const unsigned int nn,
       const unsigned int ns,
       const Elem* p) :
    Face(nn, ns, p) {}
    
  /**
   * @returns 3.  All tri-derivatives are guaranteed to have at
   * least 3 nodes.
   */
  unsigned int n_nodes() const { return 3; }

  /**
   * @returns 3
   */
  unsigned int n_sides() const { return 3; }

  /**
   * @returns 3.  All triangles have 3 vertices.
   */
  unsigned int n_vertices() const { return 3; }

  /**
   * @returns 3.  All triangles have 3 edges.
   */
  unsigned int n_edges() const { return 3; }

  /**
   * @returns 4
   */
  unsigned int n_children() const { return 4; }
  
  /**
   * @returns an id associated with the \p s side of this element.
   * The id is not necessariy unique, but should be close.  This is
   * particularly useful in the \p MeshBase::find_neighbors() routine.
   */
  unsigned int key (const unsigned int s) const;

  /**
   * @returns a primitive (2-noded) edge for 
   * edge i.
   */
  AutoPtr<DofObject> side (const unsigned int i) const;

  /**
   * Based on the quality metric q specified by the user,
   * returns a quantitative assessment of element quality.
   */
  Real quality (const ElemQuality q) const;

  /**
   * Returns the suggested quality bounds for
   * the hex based on quality measure q.  These are
   * the values suggested by the CUBIT User's Manual.
   */
  std::pair<Real, Real> qual_bounds (const ElemQuality q) const;

};

// ------------------------------------------------------------
// Tri class member functions

#endif


