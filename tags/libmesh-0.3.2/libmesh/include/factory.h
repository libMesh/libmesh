// $Id: factory.h,v 1.6 2003-02-13 22:56:07 benkirk Exp $

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



#ifndef __factory_h__
#define __factory_h__


// System & C++ includes
#include <string>
#include <map>

// Local includes
#include "mesh_common.h"
#include "auto_ptr.h"





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
  Factory(const std::string& name);

public:

  /**
   * Destructor. (Empty.)
   */
  virtual ~Factory() {}

  /**
   * Builds an object of type Base identified by name.
   */
  static AutoPtr<Base> build(const std::string& name);

  /**
   * Create a Base class.  Force this to be implemented
   * later.
   */
  virtual AutoPtr<Base> create() = 0;

private:

  /**
   * Map from a name to a Factory<Base>* pointer.
   */
  static std::map<std::string, Factory<Base>*> factory_map;
};



/**
 * Factory class implementation.
 */
template<class Derived, class Base>
class FactoryImp: public Factory<Base>
{
public:

  /**
   * Constructor.  Takes a name as input.
   */
  FactoryImp(const std::string& name) : Factory<Base>(name) { }

  /**
   * Destructor.  Empty.
   */
  ~FactoryImp() {}

private:

  /**
   * @returns a new object of type Derived. 
   */
  AutoPtr<Base> create() { return AutoPtr<Base>(new Derived); }
};




// -----------------------------------------------------
// Factory members
template<class Base>
inline
Factory<Base>::Factory(const std::string& name)
{
  if(!factory_map.empty())
    assert(!factory_map.count(name));
  
  factory_map[name] = this;
}



template<class Base>
inline
AutoPtr<Base> Factory<Base>::build(const std::string& name)
{  
  if(!factory_map.count(name))
    {
      std::cerr << "Tried to build an unknown type: " << name << std::endl;

      std::cerr << "valid options are:" << std::endl;
      
      for (typename std::map<std::string,Factory<Base>*>::const_iterator
	     it = factory_map.begin(); it != factory_map.end(); ++it)
        std::cerr << "  " << it->first << std::endl;
      
      return AutoPtr<Base>(NULL);
    }
  return AutoPtr<Base>(factory_map[name]->create());  
}



#endif
