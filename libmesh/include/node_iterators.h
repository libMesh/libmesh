// $Id: node_iterators.h,v 1.4 2003-02-25 16:26:46 ddreyer Exp $

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

#ifndef __node_iterators_h__
#define __node_iterators_h__

// C++ includes

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
typedef basic_node_iterator<std::vector<Node*>::iterator> node_iterator;




/**
 * Specialization of the basic_node_iterator class
 * for \p std::vector<Node*>::const_iterator.  This is what
 * users will create if they want to iterate over
 * all the elements in the mesh in const functions.
 */
typedef basic_node_iterator<std::vector<Node*>::const_iterator> const_node_iterator; 



#endif


