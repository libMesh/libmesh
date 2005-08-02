// $Id: multi_predicates.h,v 1.6 2005-08-02 21:44:02 jwpeterson Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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

#ifndef __multi_predicates_h__
#define __multi_predicates_h__

#include <vector>
#include "single_predicates.h"

/**
 * This namespace defines several multi_predicates which are used by
 * the element and node iterators.  These classes are not in general
 * used by the user, although they could be.
 *
 * @author John W. Peterson, 2004
 */
namespace Predicates
{
  
  // Empty place-holder base class for multi_predicates
  struct multi_predicate {};
  

  // This class represents a generic combination of more than one predicate.
  // It is meant to be derived from to actually be used.
  template <typename T>
  struct abstract_multi_predicate : multi_predicate
  {
    // virtual destructor.
    virtual ~abstract_multi_predicate()
    {
      // Clean-up vector
      for (unsigned int i=0; i<_predicates.size(); ++i)
	delete _predicates[i];
    }

    // operator= (perform deep copy of entries in _predicates vector
    abstract_multi_predicate& operator=(const abstract_multi_predicate& rhs)
    {
      // First clear out the predicates vector
      for (unsigned int i=0; i<_predicates.size(); ++i)
	delete _predicates[i];

      // Now copy over the information from the rhs.
      this->deep_copy(rhs);
      
      return *this;
    }
    
    // operator() checks all the predicates in the vector.
    virtual bool operator()(const T& it) const
    {
      for (unsigned int i=0; i<_predicates.size(); ++i)
	{
	  const predicate<T>* pred = _predicates[i];

	  assert (pred != NULL);
	  
	  if ( ! (*pred)(it) )
	    return false;
	}
      
      return true;
    }
    
  protected:
    // Do not instantiate the base class.
    abstract_multi_predicate() {}

    // Copy constructor.
    abstract_multi_predicate(const abstract_multi_predicate& rhs)
    {
      this->deep_copy(rhs);
    }

    // The deep_copy function is used by both the op= and
    // copy constructors.  This function uses the default (empty)
    // copy constructor for the predicate class.
    void deep_copy(const abstract_multi_predicate& rhs)
    {
      for (unsigned int i=0; i<rhs._predicates.size(); ++i)
	_predicates.push_back(rhs._predicates[i]->clone());
    }
    
    // Predicates to be evaluated.
    std::vector<predicate<T>*> _predicates;
  };
  
  

  // Instantiation of the IsNull abstract_multi_predicate.
  // This would be used to iterate over NULL entries in a container.
  template <typename T>
  struct IsNull : abstract_multi_predicate<T>
  {
    // Constructor, pushes back a single predicate
    IsNull()
    {
      this->_predicates.push_back(new is_null<T>);
    }
  };
  




  
  // Instantiation for the NotNull abstract_multi_predicate
  template <typename T>
  struct NotNull : abstract_multi_predicate<T>
  {
    // Constructor, pushes back a single predicate
    NotNull()
    {
      this->_predicates.push_back(new not_null<T>);
    }
  };





  // Instantiation for the Active abstract_multi_predicate
  template <typename T>
  struct Active : abstract_multi_predicate<T>
  {
    // Constructor, pushes back two single predicates
    Active()
    {
      this->_predicates.push_back(new not_null<T>);
      this->_predicates.push_back(new active<T>);
    }
  };



  // Instantiation for the NotActive abstract_multi_predicate
  template <typename T>
  struct NotActive : abstract_multi_predicate<T>
  {
    // Constructor, pushes back two single predicates
    NotActive()
    {
      this->_predicates.push_back(new not_null<T>);
      this->_predicates.push_back(new not_active<T>);
    }
  };




  
  // Instantiation for the Local abstract_multi_predicate
  template <typename T>
  struct Local : abstract_multi_predicate<T>
  {
    // Constructor, pushes back two single predicates
    Local()
    {
      this->_predicates.push_back(new not_null<T>);
      this->_predicates.push_back(new pid<T>(libMesh::processor_id()));
    }

  };


  // Instantiation for the Type abstract_multi_predicate
  template <typename T>
  struct Type : abstract_multi_predicate<T>
  {
    Type(const ElemType type)
    {
      this->_predicates.push_back(new not_null<T>);
      this->_predicates.push_back(new elem_type<T>(type));
    }
  };

  
  
  // Instantiation for the ActiveType abstract_multi_predicate
  template <typename T>
  struct ActiveType : abstract_multi_predicate<T>
  {
    ActiveType(const ElemType type)
    {
      this->_predicates.push_back(new not_null<T>);
      this->_predicates.push_back(new active<T>);
      this->_predicates.push_back(new elem_type<T>(type));
    }
  };



  // Instantiation for the ActivePID abstract_multi_predicate
  template <typename T>
  struct ActivePID : abstract_multi_predicate<T>
  {
    ActivePID(const unsigned int proc_id)
    {
      this->_predicates.push_back(new not_null<T>);
      this->_predicates.push_back(new active<T>);
      this->_predicates.push_back(new pid<T>(proc_id));
    }
  };




  
  // Instantiation for the ActiveLocal abstract_multi_predicate
  template <typename T>
  struct ActiveLocal : abstract_multi_predicate<T>
  {
    ActiveLocal()
    {
      this->_predicates.push_back(new not_null<T>);
      this->_predicates.push_back(new active<T>);
      this->_predicates.push_back(new pid<T>(libMesh::processor_id()));
    }
  };





  // Instantiation for the PID abstract_multi_predicate
  template <typename T>
  struct PID : abstract_multi_predicate<T>
  {
    PID(const unsigned int proc_id)
    {
      this->_predicates.push_back(new not_null<T>);
      this->_predicates.push_back(new pid<T>(proc_id));
    }
  };




  // Instantiation for the Level abstract_multi_predicate
  template <typename T>
  struct Level : abstract_multi_predicate<T>
  {
    Level(const unsigned int l)
    {
      this->_predicates.push_back(new not_null<T>);
      this->_predicates.push_back(new level<T>(l));
    }
  };




  // Instantiation for the NotLevel abstract_multi_predicate
  template <typename T>
  struct NotLevel : abstract_multi_predicate<T>
  {
    NotLevel(const unsigned int l)
    {
      this->_predicates.push_back(new not_null<T>);
      this->_predicates.push_back(new not_level<T>(l));
    }
  };



  // Instantiation for the ActiveOnBoundary abstract_multi_predicate
  template <typename T>
  struct ActiveOnBoundary : abstract_multi_predicate<T>
  {
    ActiveOnBoundary()
    {
      this->_predicates.push_back(new not_null<T>);
      this->_predicates.push_back(new active<T>);
      this->_predicates.push_back(new null_neighbor<T>);
    }
  };


  
  // Instantiation for the BoundarySide abstract_multi_predicate
  template <typename T>
  struct BoundarySide : abstract_multi_predicate<T>
  {
    BoundarySide()
    {
      this->_predicates.push_back(new boundary_side<T>);
    }
  };
  
}


#endif
