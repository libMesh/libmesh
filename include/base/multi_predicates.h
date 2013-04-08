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

#ifndef LIBMESH_MULTI_PREDICATES_H
#define LIBMESH_MULTI_PREDICATES_H

// Local includes
#include "libmesh/single_predicates.h"

// C++ includes
#include <vector>

namespace libMesh
{

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

	  libmesh_assert (pred);

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




  // Instantiation for the Ancestor abstract_multi_predicate
  template <typename T>
  struct Ancestor : abstract_multi_predicate<T>
  {
    // Constructor, pushes back two single predicates
    Ancestor()
    {
      this->_predicates.push_back(new not_null<T>);
      this->_predicates.push_back(new ancestor<T>);
    }
  };




  // Instantiation for the NotAncestor abstract_multi_predicate
  template <typename T>
  struct NotAncestor : abstract_multi_predicate<T>
  {
    // Constructor, pushes back two single predicates
    NotAncestor()
    {
      this->_predicates.push_back(new not_null<T>);
      this->_predicates.push_back(new not_ancestor<T>);
    }
  };




  // Instantiation for the SubActive abstract_multi_predicate
  template <typename T>
  struct SubActive : abstract_multi_predicate<T>
  {
    // Constructor, pushes back two single predicates
    SubActive()
    {
      this->_predicates.push_back(new not_null<T>);
      this->_predicates.push_back(new subactive<T>);
    }
  };




  // Instantiation for the NotSubActive abstract_multi_predicate
  template <typename T>
  struct NotSubActive : abstract_multi_predicate<T>
  {
    // Constructor, pushes back two single predicates
    NotSubActive()
    {
      this->_predicates.push_back(new not_null<T>);
      this->_predicates.push_back(new not_subactive<T>);
    }
  };



  // Instantiation for the Local abstract_multi_predicate
  template <typename T>
  struct Local : abstract_multi_predicate<T>
  {
    // Constructor, pushes back two single predicates
    Local(const processor_id_type my_pid)
    {
      this->_predicates.push_back(new not_null<T>);
      this->_predicates.push_back(new pid<T>(my_pid));
    }

  };


  // Instantiation for the NotLocal abstract_multi_predicate
  template <typename T>
  struct NotLocal : abstract_multi_predicate<T>
  {
    // Constructor, pushes back two single predicates
    NotLocal(const processor_id_type my_pid)
    {
      this->_predicates.push_back(new not_null<T>);
      this->_predicates.push_back(new not_pid<T>(my_pid));
    }

  };


  // Instantiation for the ActiveNotLocal abstract_multi_predicate
  template <typename T>
  struct ActiveNotLocal : abstract_multi_predicate<T>
  {
    // Constructor, pushes back two single predicates
    ActiveNotLocal(const processor_id_type my_pid)
    {
      this->_predicates.push_back(new not_null<T>);
      this->_predicates.push_back(new active<T>);
      this->_predicates.push_back(new not_pid<T>(my_pid));
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
    ActivePID(const processor_id_type proc_id)
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
    ActiveLocal(const processor_id_type my_pid)
    {
      this->_predicates.push_back(new not_null<T>);
      this->_predicates.push_back(new active<T>);
      this->_predicates.push_back(new pid<T>(my_pid));
    }
  };





  // Instantiation for the PID abstract_multi_predicate
  template <typename T>
  struct PID : abstract_multi_predicate<T>
  {
    PID(const processor_id_type proc_id)
    {
      this->_predicates.push_back(new not_null<T>);
      this->_predicates.push_back(new pid<T>(proc_id));
    }
  };



  // Instantiation for the NotPID abstract_multi_predicate
  template <typename T>
  struct NotPID : abstract_multi_predicate<T>
  {
    NotPID(const processor_id_type proc_id)
    {
      this->_predicates.push_back(new not_null<T>);
      this->_predicates.push_back(new not_pid<T>(proc_id));
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




  // Instantiation for the LocalLevel abstract_multi_predicate
  template <typename T>
  struct LocalLevel : abstract_multi_predicate<T>
  {
    LocalLevel(const processor_id_type my_pid,
	       const unsigned int l)
    {
      this->_predicates.push_back(new not_null<T>);
      this->_predicates.push_back(new pid<T>(my_pid));
      this->_predicates.push_back(new level<T>(l));
    }
  };




  // Instantiation for the LocalNotLevel abstract_multi_predicate
  template <typename T>
  struct LocalNotLevel : abstract_multi_predicate<T>
  {
    LocalNotLevel(const processor_id_type my_pid,
		  const unsigned int l)
    {
      this->_predicates.push_back(new not_null<T>);
      this->_predicates.push_back(new pid<T>(my_pid));
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



  // Instantiation for the ActiveLocalSubdomain abstract_multi_predicate
  template <typename T>
  struct ActiveLocalSubdomain : abstract_multi_predicate<T>
  {
    ActiveLocalSubdomain(const processor_id_type my_pid,
			 const subdomain_id_type subdomain_id)
    {
      this->_predicates.push_back(new not_null<T>);
      this->_predicates.push_back(new active<T>);
      this->_predicates.push_back(new pid<T>(my_pid));
      this->_predicates.push_back(new subdomain<T>(subdomain_id));
    }
  };



  // Instantiation for the ActiveSubdomain abstract_multi_predicate
  template <typename T>
  struct ActiveSubdomain : abstract_multi_predicate<T>
  {
    ActiveSubdomain(const subdomain_id_type subdomain_id)
    {
      this->_predicates.push_back(new not_null<T>);
      this->_predicates.push_back(new active<T>);
      this->_predicates.push_back(new subdomain<T>(subdomain_id));
    }
  };

}


} // namespace libMesh

#endif // LIBMESH_MULTI_PREDICATES_H
