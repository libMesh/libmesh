// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef __factory_h__
#define __factory_h__


// System & C++ includes
#include <string>
#include <map>

// Local includes
#include "libmesh_common.h"
#include "auto_ptr.h"

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
  Factory (const std::string& name);

public:

  /**
   * Destructor. (Empty.)
   */
  virtual ~Factory () {}

  /**
   * Builds an object of type Base identified by name.
   */
  static AutoPtr<Base> build (const std::string& name);

  /**
   * Create a Base class.  Force this to be implemented
   * later.
   */
  virtual AutoPtr<Base> create () = 0;


protected:

  /**
   * Map from a name to a Factory<Base>* pointer.
   */
  static std::map<std::string, Factory<Base>*>& factory_map();
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
  FactoryImp (const std::string& name) : Factory<Base>(name) { }

  /**
   * Destructor.  Empty.
   */
  ~FactoryImp () {}

private:

  /**
   * @returns a new object of type Derived.
   */
  AutoPtr<Base> create ();

};



// -----------------------------------------------------
// Factory members
template <class Base>
inline
Factory<Base>::Factory (const std::string& name)
{
  // Make sure we haven't already added this name
  // to the map
  libmesh_assert (!factory_map().count(name));

  factory_map()[name] = this;
}



template <class Base>
inline
AutoPtr<Base> Factory<Base>::build (const std::string& name)
{
  // name not found in the map
  if (!factory_map().count(name))
    {
      libMesh::err << "Tried to build an unknown type: " << name << std::endl;

      libMesh::err << "valid options are:" << std::endl;

      for (typename std::map<std::string,Factory<Base>*>::const_iterator
	     it = factory_map().begin(); it != factory_map().end(); ++it)
        libMesh::err << "  " << it->first << std::endl;

      libmesh_error();

      // Do this the stoopid way for IBM xlC
      AutoPtr<Base> ret_val (NULL);

      return ret_val;
    }

  // Do this the stoopid way for IBM xlC
  Factory<Base> *f = factory_map()[name];

  AutoPtr<Base> ret_val (f->create());

  return ret_val;
}



// Note - this cannot be inlined!
// template <class Base>
// std::map<std::string, Factory<Base>*>& Factory<Base>::factory_map()
// {
//   static std::map<std::string, Factory<Base>*> _factory_map;

//   return _factory_map;
// }



template <class Derived, class Base>
inline
AutoPtr<Base> FactoryImp<Derived,Base>::create ()
{
  // Do this the stoopid way for IBM xlC
  AutoPtr<Base> ret_val (new Derived);

  return ret_val;
}


} // namespace libMesh

#endif
