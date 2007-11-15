/*
 * Copyright (C) 2007 Benjamin S. Kirk
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef __hilbert_h__
#define __hilbert_h__

#include <Hilbert.hpp>

// Specific extensions to libHilbert-0.2.
namespace Hilbert 
{
  typedef FBV_UINT   /**/ inttype;

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
  };
}


void __hilbert_max_op (Hilbert::HilbertIndices *in, Hilbert::HilbertIndices *inout, int *len, void *);
void __hilbert_min_op (Hilbert::HilbertIndices *in, Hilbert::HilbertIndices *inout, int *len, void *);

#endif
