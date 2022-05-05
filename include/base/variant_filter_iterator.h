// The libMesh Finite Element Library.
// Copyright (C) 2002-2022 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_VARIANT_FILTER_ITERATOR_H
#define LIBMESH_VARIANT_FILTER_ITERATOR_H


// C++ includes
#include <algorithm> // for std::swap
#include <cstddef>
#include <cstdlib>   // for std::abort()
#include <iterator>

#if defined(__GNUC__) && (__GNUC__ < 3)  && !defined(__INTEL_COMPILER)
#include <typeinfo>
#endif

// Local includes
#include "libmesh/libmesh_common.h" // for cast_ptr()

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
 * \author John W. Peterson
 * \date 2004
 */
template<class Predicate, class Type, class ReferenceType = Type &, class PointerType = Type *,
                          class ConstType = const Type, class ConstReferenceType = const Type &,
                          class ConstPointerType = const Type *>
class variant_filter_iterator
{
public:
  using iterator_category = std::forward_iterator_tag;
  using value_type = Type;
  using difference_type = std::ptrdiff_t;
  using pointer = PointerType;
  using reference = ReferenceType;

  /**
   * Shortcut name for the fully-qualified typename.
   */
  typedef variant_filter_iterator<Predicate, Type, ReferenceType, PointerType,
                                  ConstType, ConstReferenceType, ConstPointerType> Iterator;



public:
  /**
   * Abstract base class for the iterator type.  Ideally these mixin classes would be protected,
   * but due to the fact that different templated versions of the same class (which are not related
   * by inheritance) need to be able to see each other's IterBase and PredBase members.  Thus, the
   * mixin classes are in the public interface.
   */
  struct IterBase
  {
    virtual ~IterBase() = default;
    virtual  IterBase * clone() const = 0;

    /**
     * Dereferences the iterator.
     */
    virtual ReferenceType operator*() const = 0;

    /**
     * Pre-increments the iterator.
     */
    virtual IterBase & operator++() = 0;

    virtual bool equal(const IterBase * other) const = 0;

    // typedef typename variant_filter_iterator<Predicate, Type, const Type &, const Type *>::IterBase const_IterBase;
    typedef typename variant_filter_iterator<Predicate, ConstType, ConstReferenceType, ConstPointerType>::IterBase const_IterBase;

    /**
     * Similar to the \p clone() function.
     *
     * \returns A pointer to a copy of a different type.
     */
    virtual const_IterBase * const_clone() const = 0;
  };







  /**
   * Abstract base class for the predicate.
   */
  struct PredBase
  {
    virtual ~PredBase() = default;
    virtual PredBase * clone() const = 0;
    virtual bool operator()(const IterBase * in) const = 0;

    // typedef typename variant_filter_iterator<Predicate, Type, const Type &, const Type *>::PredBase const_PredBase;
    typedef typename variant_filter_iterator<Predicate, ConstType, ConstReferenceType, ConstPointerType>::PredBase const_PredBase;

    /**
     * Similar to the \p clone() function.
     *
     * \returns A pointer to a copy of a different type.
     */
    virtual const_PredBase * const_clone() const = 0;
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
    Iter (const IterType & v) :
      iter_data (v)
    {
      // libMesh::out << "In Iter<IterType>::Iter(const IterType & v)" << std::endl;
    }


    /**
     * Copy Constructor.
     */
    Iter (const Iter & other) :
      iter_data(other.iter_data)
    {}


    /**
     * Destructor
     */
    virtual ~Iter () = default;

    /**
     * \returns A copy of this object as a pointer to
     * the base (non-templated) class.
     */
    virtual IterBase * clone() const override
    {
#ifdef __SUNPRO_CC
      variant_filter_iterator::Iter<IterType> * copy =
        new variant_filter_iterator::Iter<IterType>(iter_data);
#else
      Iter<IterType> * copy =
        new Iter<IterType>(iter_data);
#endif

      return copy;
    }

    /**
     * \returns A copy of this object as a pointer to a
     * different type of object.
     */
    virtual typename IterBase::const_IterBase * const_clone() const override
    {
      /**
       * Important typedef for const_iterators.  Notice the weird syntax!  Does it compile everywhere?
       */
      // typedef typename variant_filter_iterator<Predicate, Type, const Type &, const Type *>::template Iter<IterType> const_Iter;
      typedef typename variant_filter_iterator<Predicate, ConstType, ConstReferenceType, ConstPointerType>::template Iter<IterType> const_Iter;

      typename IterBase::const_IterBase * copy =
        new const_Iter(iter_data);

      return copy;
    }

    /**
     * Dereferences the iterator.
     */
    virtual ReferenceType operator*() const override
    {
      return * this->iter_ptr();
    }

    /**
     * Pre-increments the iterator.
     */
    virtual Iter & operator++() override
    {
      ++iter_data;
      return *this;
    }

    /**
     * Use a dynamic cast to convert the base pointer
     * passed in to the derived type.  If the cast
     * fails it means you compared two different derived
     * classes.
     */
    virtual bool equal(const IterBase * other) const override
    {
#if defined(__SUNPRO_CC) || (defined(__GNUC__) && (__GNUC__ < 3)  && !defined(__INTEL_COMPILER))
      const variant_filter_iterator::Iter<IterType> * p =
        libMesh::cast_ptr<const variant_filter_iterator::Iter<IterType> *>(other);
#else
      const Iter<IterType> * p =
        libMesh::cast_ptr<const Iter<IterType> *>(other);
#endif

      return (iter_data == p->iter_data);
    }

    /**
     * This is the iterator passed by the user.
     */
    IterType iter_data;
  private:
    /**
     * This seems to work around a bug in g++ prior to version 10 and
     * clang++ prior to version 10
     */
    PointerType iter_ptr () const
    {
      return &*iter_data;
    }
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
    Pred (const PredType & v) :
      pred_data (v) {}

    /**
     * Destructor
     */
    virtual ~Pred () = default;

    /**
     * \returns A copy of this object as a pointer to the base class.
     */
    virtual PredBase * clone() const override
    {
#ifdef __SUNPRO_CC
      variant_filter_iterator::Pred<IterType,PredType> * copy =
        new variant_filter_iterator::Pred<IterType,PredType>(pred_data);
#else
      Pred<IterType,PredType> * copy =
        new Pred<IterType,PredType>(pred_data);
#endif

      return copy;
    }


    /**
     * The redefinition of the const_clone function for the Pred class.
     */
    virtual typename PredBase::const_PredBase * const_clone() const override
    {
      /**
       * Important typedef for const_iterators.
       */
      //      typedef typename variant_filter_iterator<Predicate, Type, const Type &, const Type *>::template Pred<IterType, PredType> const_Pred;
      typedef typename variant_filter_iterator<Predicate, ConstType, ConstReferenceType, ConstPointerType>::template Pred<IterType, PredType> const_Pred;


      typename PredBase::const_PredBase * copy =
        new const_Pred(pred_data);

      return copy;
    }




    /**
     * Re-implementation of op()
     */
    virtual bool operator() (const IterBase * in) const override
    {
      libmesh_assert(in);

      // Attempt downcast
#if defined(__SUNPRO_CC) || (defined(__GNUC__) && (__GNUC__ < 3)  && !defined(__INTEL_COMPILER))
      const variant_filter_iterator::Iter<IterType> * p =
        libMesh::cast_ptr<const variant_filter_iterator::Iter<IterType> * >(in);
#else
      const Iter<IterType> * p =
        libMesh::cast_ptr<const Iter<IterType> *>(in);
#endif

      // Return result of op() for the user's predicate.
      return pred_data(p->iter_data);
    }

    /**
     * This is the predicate passed in by the user.
     */
    PredType pred_data;
  };



public:
  /**
   * Ideally this private member data should have protected access.  However, if we want
   * a const_iterator to be constructable from an non-const one, templated versions of the
   * same class (not related by inheritance) will need to know about these private members.
   * Thus, they have public access.
   *
   * Polymorphic pointer to the object.  Don't confuse
   * with the data pointer located in the \p Iter!
   */
  IterBase * data;

  /**
   * Also have a polymorphic pointer to the end object,
   * this prevents iterating past the end.
   */
  IterBase * end;

  /**
   * The predicate object.  Must have op() capable of
   * operating on IterBase * pointers.  Therefore it has
   * to follow the same paradigm as \p IterBase.
   */
  PredBase * pred;



public:
  /**
   * Templated Constructor.  Allows you to construct the iterator
   * and predicate from any types.  Also advances the data pointer
   * to the first entry which satisfies the predicate.
   *
   * \note The initialization list uses the default IterBase copy
   * constructor.
   */
  template<typename PredType, typename IterType>
  variant_filter_iterator (const IterType & d,
                           const IterType & e,
                           const PredType & p ):
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
    data(nullptr),
    end(nullptr),
    pred(nullptr) {}

  /**
   * Copy Constructor.
   * Copy the internal data instead of sharing it.
   */
  variant_filter_iterator (const Iterator & rhs) :
    data (rhs.data != nullptr ? rhs.data->clone() : nullptr),
    end  (rhs.end  != nullptr ? rhs.end->clone()  : nullptr),
    pred (rhs.pred != nullptr ? rhs.pred->clone() : nullptr) {}



  /**
   * Copy construct from another (similar) variant_filter_iterator.
   * The Predicate is the same, but the Type, ReferenceType and
   * PointerType are different.  Example:
   * You are iterating over a std::vector<int *> with std::vector<int *>::iterator
   * Then, you have:
   * Type=int * ,  ReferenceType=int *& , PointerType=int **
   * On the other hand, when you iterate using std::vector<int *>::const_iterator
   * you have:
   * Type=int * const, ReferenceType=int * const & , PointerType=int * const *
   */
  template <class OtherType, class OtherReferenceType, class OtherPointerType,
            class OtherConstType, class OtherConstReferenceType, class OtherConstPointerType>
  variant_filter_iterator (const variant_filter_iterator<Predicate, OtherType, OtherReferenceType, OtherPointerType,
                                                         OtherConstType, OtherConstReferenceType, OtherConstPointerType> & rhs)
    : data (rhs.data != nullptr ? rhs.data->const_clone() : nullptr),
      end  (rhs.end  != nullptr ? rhs.end->const_clone()  : nullptr),
      pred (rhs.pred != nullptr ? rhs.pred->const_clone() : nullptr)
  {
    // libMesh::out << "Called templated copy constructor for variant_filter_iterator" << std::endl;
  }






  /**
   * Destructor
   */
  virtual ~variant_filter_iterator()
  {
    delete data; data = nullptr;
    delete end;  end  = nullptr;
    delete pred; pred = nullptr;
  }

  /**
   * unary op*() forwards on to \p Iter::op*()
   */
  ReferenceType operator*() const
  {
    return **data;
  }


  /**
   * op->()
   */
  PointerType operator->() const
  {
    return (&**this);
  }

  /**
   * op++() forwards on to \p Iter::op++()
   */
  Iterator & operator++()
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
   * Forwards to the \p equal() function defined for the
   * IterBase pointer.
   */
  bool equal(const variant_filter_iterator & other) const
  {
    return data->equal(other.data);
  }

  /**
   * swap, used to implement op=
   */
  void swap(Iterator & lhs, Iterator & rhs)
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
  Iterator & operator=(const Iterator & rhs)
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






// op==
template<class Predicate, class Type, class ReferenceType, class PointerType,
         class OtherConstType, class OtherConstReferenceType, class OtherConstPointerType>
inline
bool operator==(const variant_filter_iterator<Predicate, Type, ReferenceType, PointerType,
                                              OtherConstType, OtherConstReferenceType, OtherConstPointerType> & x,
                const variant_filter_iterator<Predicate, Type, ReferenceType, PointerType,
                                              OtherConstType, OtherConstReferenceType, OtherConstPointerType> & y)
{
  return x.equal(y);
}



// op!=
template<class Predicate, class Type, class ReferenceType, class PointerType,
         class OtherConstType, class OtherConstReferenceType, class OtherConstPointerType>
inline
bool operator!=(const variant_filter_iterator<Predicate, Type, ReferenceType, PointerType,
                                              OtherConstType, OtherConstReferenceType, OtherConstPointerType> & x,
                const variant_filter_iterator<Predicate, Type, ReferenceType, PointerType,
                                              OtherConstType, OtherConstReferenceType, OtherConstPointerType> & y)
{
  return !(x == y);
}



#endif // LIBMESH_VARIANT_FILTER_ITERATOR_H
