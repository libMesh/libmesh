// $Id: elem_iterators.h,v 1.5 2003-02-21 22:40:59 jwpeterson Exp $

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

#ifndef __elem_iterators_h__
#define __elem_iterators_h__

// C++ includes

// Local includes
#include "predicated_iterator.h"
#include "elem.h"


/**
 * The basic_elem_iterator class simply redefines
 * the predicate to return whether or not what is currently
 * iterated to is NULL. The typedefs elem_iterator and const_elem_iterator
 * should be used to instantiate actual iterators.  Since other
 * classes will be derived from this one, we must be able to
 * specify whether or not to advance() in the constructor.
 * This is the third (default=true) parameter in the constructor.
 */
template <class T>
class basic_elem_iterator : public PredicatedIterator<T>
{
public:
  /**
   * Constructor.
   */
  basic_elem_iterator(const std::pair<T,T>& p,
		      const bool b=true)
    : PredicatedIterator<T>(p.first, p.second)
  {
    if (b) advance();
  }

protected:
  /**
   * Redefinition of the predicate
   */
  virtual bool predicate() const { return *_current != NULL; }
};




/**
 * Specialization of the basic_elem_iterator class
 * for \p std::vector<Elem*>::iterator.  This is what
 * users will create if they want to iterate over
 * all the elements in the mesh.  The neighbor iterator
 * may be used to access the _neighbors array in the
 * Elem class.
 */
typedef basic_elem_iterator<std::vector<Elem*>::iterator> elem_iterator;


/**
 * Specialization of the basic_elem_iterator class
 * for \p Elem**.  This is what
 * users will create if they want to iterate over
 * all the neighbors of an element in non-const functions.
 */
typedef basic_elem_iterator<Elem**> neighbor_iterator;




/**
 * Specialization of the basic_elem_iterator class
 * for \p std::vector<Elem*>::const_iterator.  This is what
 * users will create if they want to iterate over
 * all the elements in the mesh in const functions.
 * The const_neighbor_iterator may be used to access the
 * _neighbors array in the Elem class in const functions.
 */
typedef basic_elem_iterator<std::vector<Elem*>::const_iterator> const_elem_iterator; 


/**
 * Specialization of the basic_elem_iterator class
 * for \p const Elem**.  This is what
 * users will create if they want to iterate over
 * all the neighbors of an element in const functions.
 */
typedef basic_elem_iterator<const Elem**> const_neighbor_iterator;




/**
 * The basic_active_elem_iterator class is an un-specialized templated
 * class which iterates over active elements.  Active elements are
 * those for which the member function active() returns true.  Use
 * the typedefs active_elem_iterator and const_active_elem_iterator
 * to instantiate non-const and const versions of the iterator,
 * respectively.  Note the "false" argument to the basic_elem_iterator
 * constructor.  This basically tells that constructor to NOT advance();
 * This class uses virtual inheritance because it will eventually
 * be used in a "diamond" inheritance hierarchy.
 */
template <class T>
class basic_active_elem_iterator : virtual public basic_elem_iterator<T>
{
public:

  /**
   * Constructor
   */
  basic_active_elem_iterator(const std::pair<T,T>& p,
			     const bool b=true)
    : basic_elem_iterator<T>(p, false)
  {
    if (b) advance();
  }

  
protected:
  /**
   * Re-Definition of the predicate.
   */
  virtual bool predicate() const { return (*_current)->active(); }
};







/**
 * Specialization of the basic_active_elem_iterator class for
 * for \p std::vector<Elem*>::iterator.  This is what users
 * will create if they want to iterate over the all the active
 * elements in the mesh.
 */
typedef basic_active_elem_iterator<std::vector<Elem*>::iterator> active_elem_iterator;




/**
 * Specialization of the basic_active_elem_iterator class for
 * for \p std::vector<Elem*>::const_iterator.  This is what users
 * will create if they want to iterate over the all the active
 * elements in the mesh in const functions.
 */
typedef basic_active_elem_iterator<std::vector<Elem*>::const_iterator> const_active_elem_iterator; 





/**
 * The basic_not_active_elem_iterator
 * iterates over inactive elements.  Inactive elements are
 * those for which the member function active() returns false.  Use
 * the typedefs not_active_elem_iterator and const_not_active_elem_iterator
 * to instantiate non-const and const versions of the iterator,
 * respectively.  Note the "false" argument to the basic_active_elem_iterator
 * constructor.  This basically tells that constructor to NOT advance();
 * Note that this class needs to explicilty instantiate the basic_elem_iterator
 * in its initialization list since its parent is a virtual public basic_elem_iterator.
 */
template <class T>
class basic_not_active_elem_iterator : public basic_active_elem_iterator<T>
{
public:
  /**
   * Constructor.  Explicilty initializes the base class constructor.
   * Note that it passes the false parameter so that the base class
   * does not call the advance() function.
   */
  basic_not_active_elem_iterator(const std::pair<T,T>& p)
    : basic_elem_iterator<T>(p, false),
      basic_active_elem_iterator<T>(p, false)
  {
    advance();
  }
  
protected:
  /**
   * Re-definition of the predicate.  Returns the negation of the
   * basic_active_elem_iterator predicate.
   */
  virtual bool predicate() const { return !(basic_active_elem_iterator<T>::predicate()); }
};




/**
 * Specialization of the basic_not_active_elem_iterator class for
 * for \p std::vector<Elem*>::iterator.  This is what users
 * will create if they want to iterate over the all the inactive
 * elements in the mesh.
 */
typedef basic_not_active_elem_iterator<std::vector<Elem*>::iterator> not_active_elem_iterator;




/**
 * Specialization of the basic_not_active_elem_iterator class for
 * for \p std::vector<Elem*>::const_iterator.  This is what users
 * will create if they want to iterate over the all the inactive
 * elements in the mesh in const functions.
 */
typedef basic_not_active_elem_iterator<std::vector<Elem*>::const_iterator> const_not_active_elem_iterator; 




/**
 * The const_type_elem_iterator class is an unspecialized template
 * which only iterates over a certain type of element.
 * The typedefs type_elem_iterator and const_type_elem_iterator
 * should be used to instantiate actual iterators.  This class
 * adds the additional internal state which describes the type
 * of element.  Note the "false" argument to the basic_elem_iterator
 * constructor.  This basically tells that constructor to NOT advance();
 * This class uses virtual inheritance 
 */
template <class T>
class basic_type_elem_iterator : virtual public basic_elem_iterator<T>
{
public:
  /**
   * Constructor.
   */
  basic_type_elem_iterator(const std::pair<T,T>& p,
			   const ElemType t,
			   const bool b=true)
    : basic_elem_iterator<T>(p, false),
      _type(t)
  {
    if (b) advance();
  }
  
protected:
  /**
   * Definition of the predicate.
   */
  virtual bool predicate() const { return (*_current)->type() == _type; }

  const ElemType _type;
};



/**
 * Specialization of the basic_type_elem_iterator for
 * \p std::vector<Elem*>::iterator.  This is what users
 * create when they want to iterate over all elements
 * of a specific type.
 */
typedef basic_type_elem_iterator<std::vector<Elem*>::iterator> type_elem_iterator;


/**
 * Specialization of the basic_type_elem_iterator for
 * \p std::vector<Elem*>::const_iterator.  This is what users
 * create when they want to iterate over all elements
 * of a specific type in const functions.
 */
typedef basic_type_elem_iterator<std::vector<Elem*>::const_iterator> const_type_elem_iterator;







/**
 * The const_type_elem_iterator class is an unspecialized template
 * which only iterates over a certain type of element.
 * The typedefs type_elem_iterator and const_type_elem_iterator
 * should be used to instantiate actual iterators.  This class
 * adds the additional internal state which describes the type
 * of element.  Note the "false" argument to the basic_elem_iterator
 * constructor.  This basically tells that constructor to NOT advance();
 * Note that this class needs to explicitly instantiate the basic_elem_iterator
 * in its initialization list since its parent is a virtual public
 * basic_elem_iterator.
 */
template <class T>
class basic_not_type_elem_iterator : public basic_type_elem_iterator<T>
{
public:
  /**
   * Constructor.
   */
  basic_not_type_elem_iterator(const std::pair<T,T>& p,
			       const ElemType t,
			       const bool b=true)
    : basic_elem_iterator<T>(p, false),
      basic_type_elem_iterator<T>(p, t, false)
  {
    if (b) advance();
  }
  
protected:
  /**
   * Definition of the predicate.
   */
  virtual bool predicate() const { return !basic_type_elem_iterator<T>::predicate(); }
};


/**
 * Specialization of the basic_not_type_elem_iterator for
 * \p std::vector<Elem*>::iterator.  This is what users
 * create when they want to iterate over all elements
 * which are NOT of a specific type in non-const functions.
 */
typedef basic_not_type_elem_iterator<std::vector<Elem*>::iterator> not_type_elem_iterator;


/**
 * Specialization of the basic_not_type_elem_iterator for
 * \p std::vector<Elem*>::const_iterator.  This is what users
 * create when they want to iterate over all elements
 * which are NOT of a specific type in const functions.
 */
typedef basic_not_type_elem_iterator<std::vector<Elem*>::const_iterator> const_not_type_elem_iterator;




/**
 * The basic_active_type_elem_iterator is a templated unspecialized
 * iterator which can be used to iterate only over active elements
 * of a specific type.  The typedefs active_type_elem_iterator
 * and const_active_type_elem_iterator should be used to instantiate
 * actual iterators.  This class uses multiple inheritance (in a diamond
 * pattern no less) so that it can mimic the behavior of both the
 * active_elem_iterator and type_elem_iterator classes without
 * reimplementing those functions.  Because its base classes use virtual
 * inheritance, we must explicitly call the constructor for their
 * common parent, the basic_elem_iterator, here.
 */
template <class T>
class basic_active_type_elem_iterator : public basic_type_elem_iterator<T>, public basic_active_elem_iterator<T>
{
public:
  /**
   * Constructor.
   */
  basic_active_type_elem_iterator(const std::pair<T,T>& p,
				  const ElemType t)
    : basic_elem_iterator<T>(p,false),
      basic_type_elem_iterator<T>(p, t, false),
      basic_active_elem_iterator<T>(p, false) 
  {
    advance();
  }
  
protected:
  /**
   * Definition of the predicate
   */
  virtual bool predicate() const
  {
    return (basic_active_elem_iterator<T>::predicate() && basic_type_elem_iterator<T>::predicate());
  }
};







/**
 * Specialization of the basic_active_type_elem_iterator for
 * \p std::vector<Elem*>::iterator.  This is what users
 * create when they want to iterate over all active elements
 * of a specific type.
 */
typedef basic_active_type_elem_iterator<std::vector<Elem*>::iterator> active_type_elem_iterator;



/**
 * Specialization of the basic_active_type_elem_iterator for
 * \p std::vector<Elem*>::const_iterator.  This is what users
 * create when they want to iterate over all active elements
 * of a specific type.
 */
typedef basic_active_type_elem_iterator<std::vector<Elem*>::const_iterator> const_active_type_elem_iterator; 





/**
 * The basic_level_elem_iterator is a templated unspecialized
 * iterator which can be used to iterate only over elements of a certain level.
 * The typedefs level_elem_iterator
 * and const_level_elem_iterator should be used to instantiate
 * actual iterators.  The third and final constructor argument,
 * which has a default true value, tells whether or not to
 * actually call advance() in the constructor.  Derived classes
 * will want to set this to false so that advance() is not called
 * under the wrong predicate!
 */
template <class T>
class basic_level_elem_iterator : public basic_elem_iterator<T>
{
public:
  /**
   * Constructor.  Explicitly calls the base class constructor.
   */
  basic_level_elem_iterator(const std::pair<T,T>& p,
			    const unsigned int l,
			    const bool b=true)
    : basic_elem_iterator<T>(p, false),
      _level(l)
  {
    if (b) advance();
  }

protected:
  /**
   * Definition of the predicate
   */
  virtual bool predicate() const { return (*_current)->level() == _level; }
  
  const unsigned int _level;
};





/**
 * Specialization of the basic_level_elem_iterator for
 * \p std::vector<Elem*>::iterator.  This is what users
 * create when they want to iterate over all active elements
 * of a specific level.
 */
typedef basic_level_elem_iterator<std::vector<Elem*>::iterator> level_elem_iterator;




/**
 * Specialization of the basic_level_elem_iterator for
 * \p std::vector<Elem*>::const_iterator.  This is what users
 * create when they want to iterate over all active elements
 * of a specific type.
 */
typedef basic_level_elem_iterator<std::vector<Elem*>::const_iterator> const_level_elem_iterator;





/**
 * This class is similar to the basic_level_elem_iterator, but
 * it only iterates over elements whose level is not equal to
 * _level.  Note that it calls its parent's constructor, with
 * a false for the third parameter
 */
template <class T>
class basic_not_level_elem_iterator : public basic_level_elem_iterator<T>
{
public:
  /**
   * Constructor.
   */
  basic_not_level_elem_iterator(const std::pair<T,T>& p,
				const unsigned int l)
    : basic_level_elem_iterator<T>(p, l, false)
  {
    advance();
  }
  
protected:
  /**
   * Redefinition of the predicate.  Simply negate the
   * predicate in the parent class.
   */
  virtual bool predicate() const { return !(basic_level_elem_iterator<T>::predicate()); }
};



/**
 * Specialization of the basic_not_level_elem_iterator for
 * \p std::vector<Elem*>::iterator.  This is what users
 * create when they want to iterate over all active elements
 * NOT of a specific level.
 */
typedef basic_not_level_elem_iterator<std::vector<Elem*>::iterator> not_level_elem_iterator;


/**
 * Specialization of the basic_not_level_elem_iterator for
 * \p std::vector<Elem*>::const_iterator.  This is what users
 * create when they want to iterate over all active elements
 * NOT of a specific type.
 */
typedef basic_not_level_elem_iterator<std::vector<Elem*>::const_iterator> const_not_level_elem_iterator;









/**
 * The basic_pid_elem_iterator is a templated unspecialized
 * iterator which can be used to iterate only over elements 
 * with a certain processor id.
 * The typedefs pid_elem_iterator
 * and const_pid_elem_iterator should be used to instantiate
 * actual iterators.  The third and final constructor argument,
 * which has a default true value, tells whether or not to
 * actually call advance() in the constructor.  Derived classes
 * will want to set this to false so that advance() is not called
 * under the wrong predicate!
 */
template <class T>
class basic_pid_elem_iterator : virtual public basic_elem_iterator<T>
{
public:
  basic_pid_elem_iterator(const std::pair<T,T>& p,
			  const unsigned int pid,
			  const bool b=true)
    : basic_elem_iterator<T>(p, false),
      _pid(pid)
  {
    if (b) advance();
  }

protected:
  /**
   * Definition of the predicate
   */
  virtual bool predicate() const { return (*_current)->processor_id() == _pid; }

  const unsigned int _pid;
};

/**
 * Specialization of the basic_pid_elem_iterator for
 * \p std::vector<Elem*>::iterator.  This is what users
 * create when they want to iterate over all elements
 * with a specific processor id.
 */
typedef basic_pid_elem_iterator<std::vector<Elem*>::iterator> pid_elem_iterator;


/**
 * Specialization of the basic_pid_elem_iterator for
 * \p std::vector<Elem*>::const_iterator.  This is what users
 * create when they want to iterate over all elements
 * with a specific processor id.
 */
typedef basic_pid_elem_iterator<std::vector<Elem*>::const_iterator> const_pid_elem_iterator;





/**
 * The basic_not_pid_elem_iterator is a templated unspecialized
 * iterator which can be used to iterate only over elements 
 * which do NOT have a certain processor id.
 * The typedefs pid_not_elem_iterator
 * and const_not_pid_elem_iterator should be used to instantiate
 * actual iterators.  The third and final constructor argument,
 * which has a default true value, tells whether or not to
 * actually call advance() in the constructor.  Derived classes
 * will want to set this to false so that advance() is not called
 * under the wrong predicate!  Note that we need to explicitly
 * instantiate the base class since the basic_pid_iterator is
 * a virtual public basic_elem_iterator.
 */
template <class T>
class basic_not_pid_elem_iterator : public basic_pid_elem_iterator<T>
{
public:
  basic_not_pid_elem_iterator(const std::pair<T,T>& p,
			      const unsigned int pid)
    : basic_elem_iterator<T>(p, false),
      basic_pid_elem_iterator<T>(p, pid, false)
  {
    advance();
  }

protected:
  /**
   * Definition of the predicate - negates the predicate of the parent class.
   */
  virtual bool predicate() const { return !(basic_pid_elem_iterator<T>::predicate()); }
};


/**
 * Specialization of the basic_not_pid_elem_iterator for
 * \p std::vector<Elem*>::iterator.  This is what users
 * create when they want to iterate over all elements
 * that DO NOT have a specific processor id.
 */
typedef basic_not_pid_elem_iterator<std::vector<Elem*>::iterator> not_pid_elem_iterator;


/**
 * Specialization of the basic_not_pid_elem_iterator for
 * \p std::vector<Elem*>::const_iterator.  This is what users
 * create when they want to iterate over all active elements
 * with a specific processor id.
 */
typedef basic_not_pid_elem_iterator<std::vector<Elem*>::const_iterator> const_not_pid_elem_iterator;



/**
 * The basic_active_pid_elem_iterator is a templated unspecialized
 * iterator which can be used to iterate only over active elements
 * on a specific processor.  The typedefs active_pid_elem_iterator
 * and const_active_pid_elem_iterator should be used to instantiate
 * actual iterators.  This class uses multiple inheritance (in a diamond
 * pattern no less) so that it can mimic the behavior of both the
 * active_elem_iterator and pid_elem_iterator classes without
 * reimplementing those functions.  Because its base classes use virtual
 * inheritance, we must explicitly call the constructor for their
 * common parent, the basic_elem_iterator, here.
 */
template <class T>
class basic_active_pid_elem_iterator : public basic_pid_elem_iterator<T>, public basic_active_elem_iterator<T>
{
public:
  /**
   * Constructor.
   */
  basic_active_pid_elem_iterator(const std::pair<T,T>& p,
				  const unsigned int pid)
    : basic_elem_iterator<T>(p,false),
      basic_pid_elem_iterator<T>(p, pid, false),
      basic_active_elem_iterator<T>(p, false)
  {
    advance();
  }

protected:
  /**
   * Definition of the predicate
   */
  virtual bool predicate() const
  {
    return (basic_active_elem_iterator<T>::predicate() && basic_pid_elem_iterator<T>::predicate());
  }
};


/**
 * Specialization of the basic_active_pid_elem_iterator for
 * \p std::vector<Elem*>::iterator.  This is what users
 * create when they want to iterate over all active elements
 * that have a specific processor id in non-const functions.
 */
typedef basic_active_pid_elem_iterator<std::vector<Elem*>::iterator> active_pid_elem_iterator;


/**
 * Specialization of the basic_active_pid_elem_iterator for
 * \p std::vector<Elem*>::const_iterator.  This is what users
 * create when they want to iterate over all active elements
 * that have a specific processor id in const functions.
 */
typedef basic_active_pid_elem_iterator<std::vector<Elem*>::const_iterator> const_active_pid_elem_iterator;

#endif


