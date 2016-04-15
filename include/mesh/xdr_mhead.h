// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_XDR_MHEAD_H
#define LIBMESH_XDR_MHEAD_H

// Local includes
#include "libmesh/xdr_head.h" // for base class
#include "libmesh/enum_elem_type.h" // for ElemType

// C++ includes
#include <vector>

namespace libMesh
{

// Forward declarations
class XdrMESH;

/**
 * The \p XdrMHEAD class.
 * This class is responsible
 * for reading/writing \p xdr mesh file headers.
 *
 * \author Bill Barth
 * \author Robert McLay
 * \author John W. Peterson
 * \date 2000
 */
class XdrMHEAD : public XdrHEAD
{
  friend class XdrMESH;
public:
  /**
   * Constructor.  Initializes the number of blocks in the mesh to 1
   * and the number of levels to zero.
   */
  XdrMHEAD() : _n_blocks(1) {}

  /**
   * Destructor.
   */
  ~XdrMHEAD()                          {}

  /**
   * Set the number of
   * elements in the mesh.
   */
  void setNumEl(int numel)             { m_numel = numel; }

  /**
   * Get the number of
   * elements in the mesh.
   */
  int  getNumEl() const                { return m_numel; }

  /**
   * Set the mesh weighting.
   * You probably shouldn't
   * set this yourself ...
   */
  void setSumWghts(int sumWghts)       { m_sumWghts = sumWghts; }

  /**
   * Get the mesh weighting.
   */
  int  getSumWghts() const             { return m_sumWghts; }

  /**
   * A mesh block by definition contains
   * only a single type of element.
   *
   * @return The number of mesh blocks.
   */
  unsigned int get_n_blocks() const { return _n_blocks; }

  /**
   * Sets the number of mesh blocks.
   */
  void set_n_blocks(const unsigned int nb) { this->_n_blocks = nb; }

  /**
   * Element block types are defined in enum_elem_type.h.
   * They may be for example TRI3, TRI6, QUAD4, etc.
   *
   * @return A writeable reference to the vector of element block types.
   */
  void get_block_elt_types(std::vector<ElemType> & bet) const { bet = block_elt_types; }

  /**
   * Set the vector of element block types
   */
  void set_block_elt_types(const std::vector<ElemType> & bet) { block_elt_types = bet; }

  /**
   * The size of each element block is
   * the total number of a given type of
   * element in the mesh.
   *
   * @return The vector of block sizes
   */
  void get_num_elem_each_block(std::vector<unsigned int> & neeb) const { neeb = num_elem_each_block; }

  /**
   * Set the vector of block sizes
   */
  void set_num_elem_each_block(const std::vector<unsigned int> & neeb) { num_elem_each_block = neeb; }


private:

  /**
   * DEAL mesh specific variables:
   *
   *
   * Tells the total number of element
   * blocks.  An element block is
   * contains only a single type of
   * element.
   */
  unsigned int _n_blocks;


  /**
   * A vector of length n_blocks
   * which describes the elemnt type
   * in each block e.g. TRI, QUAD, etc.
   * Note: The element type uniquely
   * defines the number of nodes for
   * that element.
   * @see enum_elem_type.h for more
   */
  std::vector<ElemType> block_elt_types;

  /**
   * A vector of length n_blocks
   * containing the number of elements
   * in each block.
   */
  std::vector<unsigned int> num_elem_each_block;

};

} // namespace libMesh

#endif // LIBMESH_XDR_MHEAD_H
