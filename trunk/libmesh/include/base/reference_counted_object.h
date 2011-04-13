 // $Id$

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



#ifndef __reference_counted_object_h__
#define __reference_counted_object_h__

// C++ includes
#include <typeinfo>

// Local includes
#include "reference_counter.h"

namespace libMesh
{



/**
 * This class implements reference counting. Any class that
 * is properly derived from this class will get reference counted, provided
 * that the library is configured with \p --enable-reference-counting
 * and you are compiling with \p DEBUG defined.
 * For example, the following is sufficient to define the class \p Foo
 * as a reference counted class:
 *
 * \code
 * class Foo : public ReferenceCountedObject<Foo>
 * {
 *  public:
 *
 *    Foo  () {}
 *
 *    ~Foo () {}
 *
 *    void bar ();
 * 
 *  private:
 * };
 *
 * \endcode
 *
 * \par
 * If the library is configured with \p --disable-reference-counting
 * or \p DEBUG is not defined then this class does nothing.
 * All members are inlined and empty, so they should effectively disappear.
 *
 * \author Benjamin S. Kirk, 2002-2007
 */

// ------------------------------------------------------------
// ReferenceCounter class definition
template <typename T>
class ReferenceCountedObject : public ReferenceCounter
{
protected:

  /**
   * Constructor. Protected so that you cannont
   * instantiate a \p ReferenceCountedObject, only derive
   * from it.
   */
  ReferenceCountedObject ()
  {
#if defined(LIBMESH_ENABLE_REFERENCE_COUNTING) && defined(DEBUG)
    
    increment_constructor_count(typeid(T).name());

#endif
  }

  /**
   * Also, increment the counter if the copy-constructor is called.
   */
  ReferenceCountedObject (const ReferenceCountedObject&)
  {
#if defined(LIBMESH_ENABLE_REFERENCE_COUNTING) && defined(DEBUG)
    
    increment_constructor_count(typeid(T).name());

#endif
  }

public:
  
  /**
   * Destructor.
   */
  virtual ~ReferenceCountedObject ()
  {
#if defined(LIBMESH_ENABLE_REFERENCE_COUNTING) && defined(DEBUG)
    
    increment_destructor_count(typeid(T).name());

#endif
  }
  
private:

};


} // namespace libMesh


#endif // end #ifndef __reference_counted_object_h__




