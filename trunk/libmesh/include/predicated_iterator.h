// $Id: predicated_iterator.h,v 1.1 2003-02-12 05:41:29 jwpeterson Exp $

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

#ifndef __predicated_iterator_h__
#define __predicated_iterator_h__

// C++ includes
#include <iterator>

// Local includes


/**
 * This class defines an iterator whose op++
 * may be redefined based on any sort of predicate.
 * To make the syntax as simple as possible, 
 * iterators derived from this class need only
 * reimplement the pure virtual function called
 * predicate.  Thanks to Michael, Bill, and Ben
 * for fine tuning.
 *
 * \author John W. Peterson, 2002-2003
 */
template <class T>
class PredicatedIterator : public std::iterator<std::forward_iterator_tag, T>
{
  /**
   * 
   */
private:
  PredicatedIterator() {}
  
protected:
  /**
   * Constructor.  Initializes the current pointer
   * and the end pointer.  The end pointer is necessary
   * so that we do not iterate past-the-end while
   * evaluating the predicate.  Note that in the derived
   * classes, we must immediately advance to the first element in
   * the range which satisfies the predicate.  This constructor
   * is protected to prevent accidental attempts at instantiating
   * this abstract base class.
   */
  PredicatedIterator(T const &c, T const &e) : _current(c), _end(e) {}

public:

  /**
   * Copy constructor.
   */
  PredicatedIterator(const PredicatedIterator& p) :
    _current(p._current),
    _end(p._end)
  {
    std::cout << "In PredicatedIterator(const PredicatedIterator& p)"
	      << std::endl;
  };

  
  /**
   * We need a typedef which defines the type of values
   * that we're iterating over.
   */
  typedef typename std::iterator_traits<T>::value_type value_type; // =Elem* for example?

  /**
   * The predicate.  Must be redefined in classes
   * derived from this class.
   */
  virtual bool predicate() const = 0;
    
  /**
   * Prefix op++, i.e. ++i.
   */
  PredicatedIterator& operator++()
  {
    if (_current != _end)
      {
	++_current;
	advance();
      }
    
    return *this;
  }
  
  /**
   * Postfix op++ i.e. i++.
   * Less efficient in that
   * it requires the creation of a temporary.  Therefore
   * you should avoid using it.
   */
  PredicatedIterator& operator++(int)
  {
    std::cerr << "Don't use this! It creates an unnecessary temporary!" << std::endl;
    PredicatedIterator& tmp = *this;
    this->operator++();
    return tmp;
  }
  
  /**
   * A dereferencing operator.  Returns a reference
   * to what is currently being iterated to.
   */
  //value_type& operator* () { return *_current; }

  /**
   * A const dereferencing operator.  Returns a const
   * reference to what is currently being iterated to.
   */
  const value_type& operator* () const { return *_current; }
  
  /**
   * A dereferencing operator.  Returns a pointer to
   * what is currently being iterated to.
   */
  value_type* operator-> () { return _current.operator->(); }
  //value_type* operator-> () { return &(this->operator*()); } // Deal II
  
  /**
   * A const dereferencing operator.  Returns a const pointer
   * to what is currently being iterated to.
   */
  //const value_type* operator-> () const { return _current.operator->(); }
  //const value_type* operator-> () const { return &(this->operator*()); } // Deal II
  
  
  /**
   * Assignment operator.
   */
  PredicatedIterator& operator= (const PredicatedIterator& p)
  {
    _current = p._current;
    _end     = p._end;
    return *this;
  }
  /**
   * The equivalency operator.
   */
  bool operator== (PredicatedIterator const &rhs) const
  {
    return _current == rhs._current;
  }
  
  /**
   * The not-equal operator.
   */
  bool operator!= (PredicatedIterator const &rhs) const
  {
    return _current != rhs._current;
  }
  
protected:

  /**
   * This function automatically advances _current using
   * prefix op++ to the next element which satisfies the
   * predicate.  It's not callable by the user, but maybe
   * it could be?
   */
  bool advance()
  {
    while (_current != _end)
      {
	if (predicate())
	  return true;
	
	++_current;
      }
    
    return false;
  }

  /**
   * We need to know what we're currently pointing to. 
   */
  T _current;

  /**
   * We also need to know the end of the range so
   * we don't iterate past it while using op++.
   */
  T _end;
};



#endif




