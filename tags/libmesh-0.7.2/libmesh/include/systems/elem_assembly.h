// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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

#ifndef __elem_assembly_h__
#define __elem_assembly_h__

// Local includes
#include "reference_counted_object.h"

namespace libMesh
{

class FEMContext;

/**
 * ElemAssembly provides a per-element (interior and boundary) assembly
 * functionality.
 *
 * @author David J. Knezevic, 2011
 */
class ElemAssembly : public ReferenceCountedObject<ElemAssembly>
{
public:

  /**
   * Constructor.  Initializes required
   * data structures.
   */
  ElemAssembly () {};
  
  /**
   * Destructor.
   */
  virtual ~ElemAssembly () {};
  
  /**
   * Perform the element interior assembly.
   */
  virtual void interior_assembly(FEMContext& ) { }

  /**
   * Perform the element boundary assembly.
   */  
  virtual void boundary_assembly(FEMContext& ) { }

};
 
}

#endif
