// $Id: variant_filter_iterator.h,v 1.6 2004-11-09 22:51:35 benkirk Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2004  Benjamin S. Kirk, John W. Peterson
  
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

#ifndef __variant_filter_iterator_h__
#define __variant_filter_iterator_h__



#include <iterator>
#include <iostream>
#include <algorithm> // for std::swap
#include <stdlib.h>  // for abort()

/**
 * Original Authors: Corwin Joy          * Michael Gradman
 *                   cjoy@houston.rr.com * Michael.Gradman@caminus.com
 * Caminus, Suite 1150, Two Allen Center, 1200 Smith Street, Houston, TX 77002
 * This class is an extension of variant_bidirectional_iterator to a
 * filter_iterator similar to boost's.  The filter iterator is modeled
 * after a forward_iterator since to go backward and forward requires
 * the storage of both a "begin" and "end" iterator to avoid stepping
 * off the end or the beginning.  To reduce complexity, we only allow
 * traversal in one direction.
 *
 * @author John W. Peterson, 2004.
 */
template<class Type, class Predicate>
#if defined(__GNUC__) && (__GNUC__ < 3)  && !defined(__INTEL_COMPILER)
class variant_filter_iterator : public std::forward_iterator<std::forward_iterator_tag, Type>
#else
class variant_filter_iterator : public std::iterator<std::forward_iterator_tag,  Type>
#endif
{
#ifdef __SUNPRO_CC // make these public in this case
 public:           // otherwise we'll need *lots* of friends
#else
 private:
#endif

  /**
   * Abstract base class for the iterator type.
   */
  struct IterBase
  {
    virtual ~IterBase() {}
    virtual  IterBase* clone() const = 0 ;
    virtual Type& operator*() const = 0;    // <-- CUSTOM INTERFACE METHOD
    virtual void operator++() = 0;          // <-- CUSTOM INTERFACE METHOD
    virtual bool equal(const IterBase *other) const = 0;
  };


  
  /**
   * Abstract base class for the predicate.
   */
  struct PredBase
  {
    virtual ~PredBase() {}
    virtual PredBase* clone() const = 0 ;
    virtual bool operator()(const IterBase* in) const = 0;
  };


  
  /**
   * The actual iterator object is held as a template parameter here.
   */
  template<typename IterType>
  struct Iter : IterBase
  {
    /**
     * Constructor
     */
    Iter (const IterType& v) :
      iter_data (v) {}

    /**
     * Destructor
     */
    virtual ~Iter () {}
    
    /**
     * @returns a copy of this object as a pointer to
     * the base (non-templated) class.
     */
    virtual IterBase* clone() const
    {
#ifdef __SUNPRO_CC      
      variant_filter_iterator::Iter<IterType> *copy = 
	new variant_filter_iterator::Iter<IterType>(iter_data);
#else
      Iter<IterType> *copy = 
	new Iter<IterType>(iter_data);      
#endif
      
      return copy;
    }

    /**
     * Custom interface method.
     */
    virtual Type& operator*() const   // <-- CUSTOM INTERFACE METHOD
    {
      return *iter_data;
    }

    /**
     * Custom interface method.
     */
    virtual void operator++()         // <-- CUSTOM INTERFACE METHOD
    {
      ++iter_data;
    }

    /**
     * Use dynamic_cast to convert the base pointer
     * passed in to the derived type.  If the cast
     * fails it means you compared two different derived
     * classes.
     */
    virtual bool equal(const IterBase *other) const
    {
#if defined(__SUNPRO_CC) || (defined(__GNUC__) && (__GNUC__ < 3)  && !defined(__INTEL_COMPILER))
      const variant_filter_iterator::Iter<IterType>* p = 
	dynamic_cast<const variant_filter_iterator::Iter<IterType>*>(other);
#else      
      const Iter<IterType>* p = 
	dynamic_cast<const Iter<IterType>*>(other);      
#endif
      // Check for failed cast
      if (p == NULL)
	{
	  std::cerr << "Dynamic cast failed in Iter::equal(...)" << std::endl;
	  abort();
	}
      
      return (iter_data == p->iter_data);
    }

    /**
     * This is the iterator passed by the user.
     */
    IterType iter_data;
  };


  
  /**
   * The actual predicate is held as a template parameter here.
   * There are two template arguments here, one for the actual type
   * of the predicate and one for the iterator type.
   */
  template <typename IterType, typename PredType>
  struct Pred : PredBase
  {
    /**
     * Constructor
     */
    Pred (const PredType& v) :
      pred_data (v) {}

    /**
     * Destructor
     */
    virtual ~Pred () {}

    /**
     * Returns a copy of this object as a pointer to the base class.
     */
    virtual PredBase* clone() const
    {
#ifdef __SUNPRO_CC
      variant_filter_iterator::Pred<IterType,PredType> *copy = 
	new variant_filter_iterator::Pred<IterType,PredType>(pred_data);
#else
      Pred<IterType,PredType> *copy = 
	new Pred<IterType,PredType>(pred_data);
#endif
      
      return copy;
    }
    
    /**
     * Re-implementation of op()
     */
    virtual bool operator() (const IterBase* in) const
    {
      assert (in != NULL);
      
      // Attempt downcast
#if defined(__SUNPRO_CC) || (defined(__GNUC__) && (__GNUC__ < 3)  && !defined(__INTEL_COMPILER))
      const variant_filter_iterator::Iter<IterType>* p =
	dynamic_cast<const variant_filter_iterator::Iter<IterType>* >(in);
#else
      const Iter<IterType>* p =
	dynamic_cast<const Iter<IterType>* >(in);
#endif
      
      // Check for failure
      if ( p == NULL )
	{
	  std::cerr << "Dynamic cast failed in Pred::op()" << std::endl;
	  std::cerr << "typeid(IterType).name()=" << typeid(IterType).name() << std::endl;
	  abort();
	}
      
      // Return result of op() for the user's predicate.
      return pred_data(p->iter_data);
    }
    
    /**
     * This is the predicate passed in by the user.
     */
    PredType pred_data;
  };


  
  /**
   * Polymorphic pointer to the object.  Don't confuse
   * with the data pointer located in the \p Iter!
   */
  IterBase* data;

  /**
   * Also have a polymorphic pointer to the end object,
   * this prevents iterating past the end.
   */
  IterBase* end;

  /**
   * The predicate object.  Must have op() capable of
   * operating on IterBase* pointers.  Therefore it has
   * to follow the same paradigm as \p IterBase.
   */
  PredBase* pred;


  
public:

  /**
   * The iterator this class provides.
   */
  typedef variant_filter_iterator<Type, Predicate> Iterator;

  /**
   * Templated Constructor.  Allows you to construct the iterator
   * and predicate from any types.  Also advances the data pointer
   * to the first entry which satisfies the predicate.
   */
  template<typename IterType, class PredType>
  variant_filter_iterator (const IterType& d,
			   const IterType& e,
			   const PredType& p ):
    data ( new Iter<IterType>(d) ),
    end  ( new Iter<IterType>(e) ),
    pred ( new Pred<IterType,PredType>(p) )
  {
    this->satisfy_predicate();
  }
  
  /**
   * Default Constructor.
   */
  variant_filter_iterator () :
    data(NULL),
    end(NULL),
    pred(NULL) {}

  /**
   * Copy Constructor.
   * Copy the internal data instead of sharing it.
   */
  variant_filter_iterator (const Iterator& rhs) :
    data (rhs.data != NULL ? rhs.data->clone() : NULL),
    end  (rhs.end  != NULL ? rhs.end->clone()  : NULL),
    pred (rhs.pred != NULL ? rhs.pred->clone() : NULL) {}
  
  /**
   * Destructor
   */
  ~variant_filter_iterator()
  {
    delete data; data = NULL;
    delete end;  end  = NULL;
    delete pred; pred = NULL;
  }
  
  /**
   * unary op*() forwards on to \p Iter::op*()
   */
  Type& operator*() const
  {
    return **data;
  }


  /**
   * op->()
   */
  Type* operator->() const
  {
    return (&**this);
  }
  
  /**
   * op++() forwards on to \p Iter::op++()
   */
  Iterator& operator++()
  {
    ++*data;
    this->satisfy_predicate();
    return (*this);
  }

  /**
   * postfix op++(), creates a temporary!
   */
  const Iterator operator++(int) // const here to prevent iterator++++ type operations
  {
    Iterator oldValue(*this); // standard is to return old value
    ++*data;
    this->satisfy_predicate();
    return oldValue;
  }

  /**
   * forwards on the the equal function defined for the
   * IterBase pointer.  Possibly also compare the end pointers,
   * but this is usually not important and would require an
   * additional dynamic_cast.
   */
  bool equal(const variant_filter_iterator& other) const
  {
    return data->equal(other.data);
  }

  /**
   * swap, used to implement op=
   */
  void swap(Iterator& lhs, Iterator& rhs)
  {
    // Swap the data pointers
    std::swap (lhs.data, rhs.data);

    // Swap the end pointers
    std::swap (lhs.end, rhs.end);

    // Also swap the predicate objects.
    std::swap (lhs.pred, rhs.pred);
  }

  /**
   * Assignment operator.
   */
  Iterator& operator=(const Iterator& rhs)
  {
    Iterator temp(rhs);
    swap(temp, *this);
    return *this;
  }


  
private:
  
  /**
   * Advances the data pointer until it reaches
   * the end or the predicate is satisfied.
   */
  void satisfy_predicate()
  {
    while ( !data->equal(end) && !(*pred)(data) )
      ++(*data);
  }  
};






//---------------------------------------------------------------------------
// op==
template<class Type, class Predicate>
inline
bool operator==(const variant_filter_iterator<Type, Predicate>& x,
		const variant_filter_iterator<Type, Predicate>& y)
{
  return x.equal(y);
}



// op!=
template<class Type, class Predicate>
inline
bool operator!=(const variant_filter_iterator<Type, Predicate>& x,
		const variant_filter_iterator<Type, Predicate>& y)
{
  return !(x == y);
}



#endif
