// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson
  
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

#ifndef __single_predicates_h__
#define __single_predicates_h__

#include <vector>
#include "enum_elem_type.h"

/**
 * This file declares several predicates in the Predicates namespace.  They
 * are called "single predicates" since the purpose of each one is to act
 * as a single functor which returns true or false depending on the result
 * of the operator() function.  The single predicates are used together
 * as building blocks to create the "multi predicates" which can be found
 * in the multi_predicates.h header file.
 *
 * @author John W. Peterson, 2004
 */
namespace Predicates
{
  // Forward declaration
  template <typename T> class abstract_multi_predicate;
  
  // abstract single predicate.  Derived classes must implement the clone()
  // function.  Be careful using it since it allocates memory!  The clone()
  // function is necessary since the predicate class has pure virtual
  // functions.
  template <typename T>
  struct predicate
  {
    virtual ~predicate() {}
    virtual bool operator()(const T& it) const = 0;


  protected:
    friend class abstract_multi_predicate<T>;
    virtual predicate* clone() const = 0;

  };


  // The is_null predicate returns true if the underlying pointer is NULL.
  template <typename T>
  struct is_null : predicate<T>
  {
    virtual ~is_null() {}
    virtual bool operator()(const T& it) const { return *it == NULL; }

  protected:
    virtual predicate<T>* clone() const { return new is_null<T>(*this); }
  };
  
  // The not_null predicate simply returns true if the pointer is not NULL.
  template <typename T>
  struct not_null : is_null<T>
  {
    virtual bool operator()(const T& it) const { return !is_null<T>::operator()(it); }

  protected:
    virtual predicate<T>* clone() const { return new not_null<T>(*this); }
  };


  // The active predicate returns true if the pointer is active.
  template <typename T>
  struct active : predicate<T>
  {
    virtual ~active() {}
    virtual bool operator()(const T& it) const { return (*it)->active(); }
    
  protected:
    virtual predicate<T>* clone() const { return new active<T>(*this); }
  };

  // The not active predicate returns true when the pointer is inactive
  template <typename T>
  struct not_active : active<T>
  {
    virtual bool operator()(const T& it) const { return !active<T>::operator()(it); }

  protected:
    virtual predicate<T>* clone() const { return new not_active<T>(*this); }
  };


  // The subactive predicate returns true if the pointer is subactive.
  template <typename T>
  struct subactive : predicate<T>
  {
    virtual ~subactive() {}
    virtual bool operator()(const T& it) const { return (*it)->subactive(); }
    
  protected:
    virtual predicate<T>* clone() const { return new subactive<T>(*this); }
  };

  // The not_subactive predicate returns true when the pointer is not subactive
  template <typename T>
  struct not_subactive : subactive<T>
  {
    virtual bool operator()(const T& it) const { return !subactive<T>::operator()(it); }

  protected:
    virtual predicate<T>* clone() const { return new not_subactive<T>(*this); }
  };



  // The pid predicate returns true if the pointers
  // processor id matches a given processor id.
  template <typename T>
  struct pid : predicate<T>
  {
    // Constructor
    pid(const unsigned int p) : _pid(p) {}
    virtual ~pid() {}

    // op()
    virtual bool operator()(const T& it) const { return (*it)->processor_id() == _pid; }
    
  protected:
    virtual predicate<T>* clone() const { return new pid<T>(*this); }
    const unsigned int _pid;
  };


  
  // The not_pid predicate returns ture if the pointers
  // processor id does _not_ match p.
  template <typename T>
  struct not_pid : pid<T>
  {
    not_pid(const unsigned int p) : pid<T>(p) {}

    virtual bool operator()(const T& it) const { return !pid<T>::operator()(it); }

  protected:
    virtual predicate<T>* clone() const { return new not_pid<T>(*this); }
  };


  // The elem_type predicate returns true if the pointers
  // type matches the given type.  Of course, this one can only
  // be instantiated for objects which return Elem*s when dereferened.
  template <typename T>
  struct elem_type : predicate<T>
  {
    // Constructor
    elem_type (const ElemType t) : _elem_type(t) {}
    virtual ~elem_type() {}
    
    virtual bool operator()(const T& it) const { return (*it)->type() == _elem_type; }

  protected:
    virtual predicate<T>* clone() const { return new elem_type<T>(*this); }
    const ElemType _elem_type;
  };




  
  // The level predicate returns true if the pointers level
  // matches the given level.
  template <typename T>
  struct level : predicate<T>
  {
    // Constructor
    level (const unsigned int l) : _level(l) {}
    virtual ~level() {}

    virtual bool operator()(const T& it) const { return (*it)->level() == _level; }
    
  protected:
    virtual predicate<T>* clone() const { return new level<T>(*this); }
    const unsigned int _level;
  };


  
  // The not_level predicate returns true if the pointers level
  // _does not_ match the given level.
  template <typename T>
  struct not_level : level<T>
  {
    // Constructor
    not_level(const unsigned int l) : level<T>(l) {}
  
    virtual bool operator()(const T& it) const { return !level<T>::operator()(it); }

  protected:
    virtual predicate<T>* clone() const { return new not_level<T>(*this); }
  };


  

  // The null_neighbor predicate returns true if the pointer has any
  // NULL neigbors.
  template <typename T>
  struct null_neighbor : predicate<T>
  {
    virtual ~null_neighbor() {}
    virtual bool operator()(const T& it) const
    {
      return (*it)->on_boundary();
    }
    
  protected:
    virtual predicate<T>* clone() const { return new null_neighbor<T>(*this); }
  };



  // This predicate simply forwards the work of determining whether
  // a particular side is on the boundary to the iterator itself, which
  // has more information.
  template <typename T>
  struct boundary_side : predicate<T>
  {
    virtual ~boundary_side() {}
    virtual bool operator()(const T& it) const
    {
      return it.side_on_boundary();
    }

  protected:
    virtual predicate<T>* clone() const { return new boundary_side<T>(*this); }
  };
}


#endif
