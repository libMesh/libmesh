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



#ifndef LIBMESH_FACTORY_H
#define LIBMESH_FACTORY_H


// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/auto_ptr.h"

// C++ includes
#include <cstddef>
#include <map>
#include <string>

namespace libMesh
{



/**
 * Factory class defintion.
 */
template <class Base>
class Factory
{
protected:

  /**
   * Constructor. Takes the name to be mapped.
   */
  Factory (const std::string & name);

public:

  /**
   * Destructor. (Empty.)
   */
  virtual ~Factory () {}

  /**
   * Builds an object of type Base identified by name.
   */
  static UniquePtr<Base> build (const std::string & name);

  /**
   * Create a Base class.  Force this to be implemented
   * later.
   */
  virtual UniquePtr<Base> create () = 0;


protected:

  /**
   * Map from a name to a Factory<Base> * pointer.
   */
  static std::map<std::string, Factory<Base> *> & factory_map();
};



/**
 * Factory implementation class.
 */
template <class Derived, class Base>
class FactoryImp: public Factory<Base>
{
public:

  /**
   * Constructor.  Takes a name as input.
   */
  FactoryImp (const std::string & name) : Factory<Base>(name) { }

  /**
   * Destructor.  Empty.
   */
  ~FactoryImp () {}

private:

  /**
   * @returns a new object of type Derived.
   */
  virtual UniquePtr<Base> create () libmesh_override;
};



// -----------------------------------------------------
// Factory members
template <class Base>
inline
Factory<Base>::Factory (const std::string & name)
{
  // Make sure we haven't already added this name
  // to the map
  libmesh_assert (!factory_map().count(name));

  factory_map()[name] = this;
}



template <class Base>
inline
UniquePtr<Base> Factory<Base>::build (const std::string & name)
{
  // name not found in the map
  if (!factory_map().count(name))
    {
      libMesh::err << "Tried to build an unknown type: " << name << std::endl;

      libMesh::err << "valid options are:" << std::endl;

      for (typename std::map<std::string,Factory<Base> *>::const_iterator
             it = factory_map().begin(); it != factory_map().end(); ++it)
        libMesh::err << "  " << it->first << std::endl;

      libmesh_error_msg("Exiting...");

      // We'll never get here
      return UniquePtr<Base>();
    }

  Factory<Base> * f = factory_map()[name];
  return UniquePtr<Base>(f->create());
}



// Note - this cannot be inlined!
// template <class Base>
// std::map<std::string, Factory<Base> *> & Factory<Base>::factory_map()
// {
//   static std::map<std::string, Factory<Base> *> _factory_map;

//   return _factory_map;
// }



template <class Derived, class Base>
inline
UniquePtr<Base> FactoryImp<Derived,Base>::create ()
{
  return UniquePtr<Base>(new Derived);
}


} // namespace libMesh

#endif // LIBMESH_FACTORY_H
