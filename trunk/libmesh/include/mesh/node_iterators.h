// $Id: node_iterators.h,v 1.4 2005-02-22 22:17:33 jwpeterson Exp $

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

#ifndef __node_iterators_h__
#define __node_iterators_h__

// C++ includes
#include <vector>

// Local includes
#include "predicated_iterator.h"
#include "node.h"




/**
 * The basic_node_iterator class simply redefines
 * the predicate to return whether or not what is currently
 * iterated to is NULL. The typedefs node_iterator and const_nodeiterator
 * should be used to instantiate actual iterators.  Since other
 * classes will be derived from this one, we must be able to
 * specify whether or not to advance() in the constructor.
 * This is the third (default=true) parameter in the constructor.
 */
template <class T>
class basic_node_iterator : public PredicatedIterator<T>
{
public:
  /**
   * Constructor.
   */
  basic_node_iterator(const std::pair<T,T>& p,
		      const bool b=true)
    : PredicatedIterator<T>(p.first, p.second)
  {
    if (b) this->advance();
  }

protected:
  /**
   * Redefinition of the predicate
   */
  virtual bool predicate() const { return (*this->_current != NULL); }
};


/**
 * Specialization of the basic_node_iterator class
 * for \p std::vector<Node*>::iterator.  This is what
 * users will create if they want to iterate over
 * all the nodes in the mesh.
 */
typedef basic_node_iterator<std::vector<Node*>::iterator>
node_iterator;


/**
 * Specialization of the basic_node_iterator class
 * for \p std::vector<const Node*>::const_iterator.  This is what
 * users will create if they want to iterate over
 * all the elements in the mesh in const functions.
 */
typedef basic_node_iterator<std::vector<const Node*>::const_iterator>
const_node_iterator; 






/**
 * The basic_active_node_iterator class is an un-specialized templated
 * class which iterates over active nodes.  Active nodes are
 * those for which the member function active() returns true.  Use
 * the typedefs active_node_iterator and const_active_node_iterator
 * to instantiate non-const and const versions of the iterator,
 * respectively.  Note the "false" argument to the basic_node_iterator
 * constructor.  This basically tells that constructor to NOT advance();
 */
template <class T>
class basic_active_node_iterator : public basic_node_iterator<T>
{
public:

  /**
   * Constructor
   */
  basic_active_node_iterator(const std::pair<T,T>& p,
			     const bool b=true)
    : basic_node_iterator<T>(p, false)
  {
    if (b) this->advance();
  }

  
protected:
 
  /**
   * Re-Definition of the predicate.  Test the base class predicate
   * first, this will avoid calling the \p active member on \p NULL
   * nodes.
   */
  virtual bool predicate() const
  {
    return (basic_node_iterator<T>::predicate() &&
					   (*this->_current)->active());
  }
};


/**
 * Specialization of the basic_active_node_iterator class for
 * for \p std::vector<Node*>::iterator.  This is what users
 * will create if they want to iterate over the all the active
 * nodes in the mesh.
 */
typedef basic_active_node_iterator<std::vector<Node*>::iterator>
active_node_iterator;


/**
 * Specialization of the basic_active_node_iterator class for
 * for \p std::vector<const Node*>::const_iterator.  This is what users
 * will create if they want to iterate over the all the active
 * nodes in the mesh in const functions.
 */
typedef basic_active_node_iterator<std::vector<const Node*>::const_iterator>
const_active_node_iterator; 







/**
 * The basic_active_node_iterator class is an un-specialized templated
 * class which iterates over inactive nodes.  Inactive nodes are
 * those for which the member function active() returns false.  Use
 * the typedefs not_active_node_iterator and const_not_active_node_iterator
 * to instantiate non-const and const versions of the iterator,
 * respectively.  Note the "false" argument to the basic_node_iterator
 * constructor.  This basically tells that constructor to NOT advance();
 */
template <class T>
class basic_not_active_node_iterator : public basic_active_node_iterator<T>
{
public:

  /**
   * Constructor
   */
  basic_not_active_node_iterator(const std::pair<T,T>& p,
				 const bool b=true)
    : basic_active_node_iterator<T>(p, false)
  {
    if (b) this->advance();
  }
  
  
protected:
 
  /**
   * Re-Definition of the predicate. 
   */
  virtual bool predicate() const { return !basic_active_node_iterator<T>::predicate(); }
};


/**
 * Specialization of the basic_not_active_node_iterator class for
 * for \p std::vector<Node*>::iterator.  This is what users
 * will create if they want to iterate over the all the inactive
 * nodes in the mesh.
 */
typedef basic_not_active_node_iterator<std::vector<Node*>::iterator>
not_active_node_iterator;


/**
 * Specialization of the basic_not_active_node_iterator class for
 * for \p std::vector<const Node*>::const_iterator.  This is what users
 * will create if they want to iterate over the all the inactive
 * nodes in the mesh in const functions.
 */
typedef basic_not_active_node_iterator<std::vector<const Node*>::const_iterator>
const_not_active_node_iterator; 



#endif // #ifndef __node_iterators_h__


