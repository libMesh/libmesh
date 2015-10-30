// The libMesh Finite Element Library.
// Copyright (C) 2002-2013 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_LIBMESH_SINGLETON_H
#define LIBMESH_LIBMESH_SINGLETON_H

#include "libmesh/libmesh_common.h"

namespace libMesh {

/**
 * Base class for all library singleton objects.
 */
class Singleton
{
protected:

  /**
   * Constructor.  Adds the derived object to the singleton cache list.
   */
  Singleton();

  /**
   * Destructor.
   */
  virtual ~Singleton() {}

public:

  /**
   * Abstract base class for runtime singleton setup.
   * This will be called from the \p LibMeshInit constructor.
   */
  class Setup
  {
  protected:
    /**
     * Constructor.  Adds the derived object to the setup cache list.
     */
    Setup ();

  public:
    /**
     * Destructor.
     */
    virtual ~Setup() {}

    /**
     * Setup method.  Importantly, this is called *after main()* from the
     * \p LibMeshInit constructor.
     */
    virtual void setup () = 0;
  };

  /**
   * Setup function.  Initializes any derived \p Singleton::Setup objects.
   * objects.
   */
  static void setup();

  /**
   * Cleanup function.  Removes all dynamically created \p Singleton
   * objects.
   */
  static void cleanup();
};

} // namespace libMesh

#endif // LIBMESH_LIBMESH_SINGLETON_H
