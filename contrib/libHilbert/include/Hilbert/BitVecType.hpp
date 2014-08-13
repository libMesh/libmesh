/*
 * Copyright (C) 2007-2014 Benjamin S. Kirk
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
 */

#ifndef _BITVECTYPE_HPP_
#define _BITVECTYPE_HPP_

#include "BigBitVec.hpp"
#include <ostream>

// Specific extensions to libHilbert-0.2.
namespace Hilbert
{
  typedef FBV_UINT /**/ inttype;

  // Forward declarations
  class BitVecType;
  struct HilbertIndices;

  /**
   * Define a simple struct to facilitate communicating
   * Hilbert indices using MPI-derived datatypes.
   */
  struct HilbertIndices
  {
    inttype rack0;
    inttype rack1;
    inttype rack2;

    // Default constructor
    HilbertIndices () :
      rack0(0),
      rack1(0),
      rack2(0)
    {}

    // Constructor from a BitVecType
    HilbertIndices (const CBigBitVec &bv)
    {
      assert (bv.rackCount() == 3);

      this->rack0 = bv.racks()[0];
      this->rack1 = bv.racks()[1];
      this->rack2 = bv.racks()[2];
    }

    // Assignment operator to a HilbertIndices from a BitVecType
    HilbertIndices&
    operator=(const CBigBitVec &bv)
    {
      assert (bv.rackCount() == 3);

      this->rack0 = bv.racks()[0];
      this->rack1 = bv.racks()[1];
      this->rack2 = bv.racks()[2];

      return *this;
    }

    // Comparison operator
    bool
    operator<(const HilbertIndices& other) const
    {
      if (other.rack2 != this->rack2) return this->rack2 < other.rack2;
      if (other.rack1 != this->rack1) return this->rack1 < other.rack1;

      return this->rack0 < other.rack0;
    }

    // Equality comparison
    bool
    operator==(const HilbertIndices& other) const
    {
      if (other.rack2 != this->rack2) return false;
      if (other.rack1 != this->rack1) return false;

      return other.rack0 == this->rack0;
    }

    // <= operator
    bool
    operator<=(const HilbertIndices& other) const
    {
      if (*this == other) return true;
      return *this < other;
    }

    friend std::ostream&
    operator << (std::ostream& os, const Hilbert::HilbertIndices& t)
    {
      os << t.rack2 << "_"
	 << t.rack1 << "_"
	 << t.rack0;

      return os;
    }
  };


  /**
   * Augment the libHilbert CBigBitVec class
   * with a convenient assignment operator from
   * a HilbertIndices object.
   */
  class BitVecType : public CBigBitVec
  {
  public:

    BitVecType(unsigned int size = 3*sizeof(double)*sizeof(inttype)) :
      CBigBitVec(size) {};

    BitVecType&
      operator=(const HilbertIndices &hi)
      {
	assert (this->rackCount() == 3);

	this->racks()[0] = hi.rack0;
	this->racks()[1] = hi.rack1;
	this->racks()[2] = hi.rack2;

	return *this;
      }

    friend std::ostream&
      operator << (std::ostream& os, const Hilbert::BitVecType& t)
      {
	assert (t.rackCount() == 3);

	os << t.racks()[2] << "_"
	   << t.racks()[1] << "_"
	   << t.racks()[0];

	return os;
      }
  };
}

#endif
