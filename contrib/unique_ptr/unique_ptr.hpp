///////////////////////////////////////////////////////////////////////////////
// UniquePtr.hpp header file
//
// Copyright 2009 Howard Hinnant, Ion Gazta&ntilde;aga.
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
// See http://www.boost.org/libs/foreach for documentation

// This is a C++03 emulation of std::unique_ptr placed in namespace boost.
// Reference http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2008/n2800.pdf
//   for the latest unique_ptr specification, and
//   reference http://www.open-std.org/jtc1/sc22/wg21/docs/lwg-active.html
//   for any pending issues against this specification.

#ifndef HINNANT_UNIQUE_PTR_HPP
#define HINNANT_UNIQUE_PTR_HPP

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits.hpp>
#include <boost/static_assert.hpp>
#include <boost/mpl/if.hpp>

namespace libMesh
{

namespace detail_UniquePtr
{

typedef char one;
struct two {one _[2];};

// An is_convertible<From, To> that considers From an rvalue (consistent with C++0X).
//   This is a simplified version neglecting the types function, array, void and abstract types
//   I had to make a special case out of is_convertible<T,T> to make move-only
//   types happy.

namespace is_conv_imp
{
template <class T> one test1(const T&);
template <class T> two test1(...);
template <class T> one test2(T);
template <class T> two test2(...);
template <class T> T source();
}

template <class T1, class T2>
struct is_convertible
{
    static const bool value = sizeof(is_conv_imp::test1<T2>(is_conv_imp::source<T1>())) == 1;
};

template <class T>
struct is_convertible<T, T>
{
    static const bool value = sizeof(is_conv_imp::test2<T>(is_conv_imp::source<T>())) == 1;
};

template <class T>
class rv
{
    T& r_;

public:
    explicit rv(T& r) : r_(r) {}
    T* operator->() {return &r_;}
    T& operator*() {return r_;}
};

template <class T>
struct identity
{
    typedef T type;
};

}  // detail_UniquePtr

template <class T>
inline
typename enable_if_c
<
    !detail_UniquePtr::is_convertible<T, detail_UniquePtr::rv<T> >::value,
    T&
>::type
move(T& t)
{
    return t;
}

template <class T>
inline
typename enable_if_c
<
    !detail_UniquePtr::is_convertible<T, detail_UniquePtr::rv<T> >::value,
    const T&
>::type
move(const T& t)
{
    return t;
}

template <class T>
inline
typename enable_if_c
<
    detail_UniquePtr::is_convertible<T, detail_UniquePtr::rv<T> >::value,
    T
>::type
move(T& t)
{
    return T(detail_UniquePtr::rv<T>(t));
}

template <class T>
inline
typename enable_if_c
<
    is_reference<T>::value,
    T
>::type
forward(typename detail_UniquePtr::identity<T>::type t)
{
    return t;
}

template <class T>
inline
typename enable_if_c
<
    !is_reference<T>::value,
    T
>::type
forward(typename detail_UniquePtr::identity<T>::type& t)
{
    return move(t);
}

template <class T>
inline
typename enable_if_c
<
    !is_reference<T>::value,
    T
>::type
forward(const typename detail_UniquePtr::identity<T>::type& t)
{
    return move(const_cast<T&>(t));
}

namespace detail_UniquePtr {

// A move-aware but stripped-down compressed_pair which only optimizes storage for T2
template <class T1, class T2, bool = is_empty<T2>::value>
class UniquePtr_storage
{
    T1 t1_;
    T2 t2_;

    typedef typename add_reference<T2>::type T2_reference;
    typedef typename add_reference<const T2>::type T2_const_reference;

    UniquePtr_storage(const UniquePtr_storage&);
    UniquePtr_storage& operator=(const UniquePtr_storage&);
public:
    operator rv<UniquePtr_storage>() {return rv<UniquePtr_storage>(*this);}

    UniquePtr_storage() : t1_(), t2_() {}

    explicit UniquePtr_storage(T1 t1)
        : t1_(move(t1)), t2_() {}

    UniquePtr_storage(T1 t1, T2 t2)
        : t1_(move(t1)), t2_(forward<T2>(t2)) {}

          T1& first()       {return t1_;}
    const T1& first() const {return t1_;}

          T2_reference second()       {return t2_;}
    T2_const_reference second() const {return t2_;}
};

template <class T1, class T2>
class UniquePtr_storage<T1, T2, true>
    : private T2
{
    T1 t1_;
    typedef T2 t2_;

    UniquePtr_storage(const UniquePtr_storage&);
    UniquePtr_storage& operator=(const UniquePtr_storage&);
public:
    operator rv<UniquePtr_storage>() {return rv<UniquePtr_storage>(*this);}

    UniquePtr_storage() : t1_() {}

    explicit UniquePtr_storage(T1 t1)
        : t1_(move(t1)) {}

    UniquePtr_storage(T1 t1, T2 t2)
        : t2_(move(t2)), t1_(move(t1)) {}

          T1& first()       {return t1_;}
    const T1& first() const {return t1_;}

          T2& second()       {return *this;}
    const T2& second() const {return *this;}
};

template <class T1, class T2, bool b>
inline
void
swap(UniquePtr_storage<T1, T2, b>& x, UniquePtr_storage<T1, T2, b>& y)
{
    using std::swap;
    swap(x.first(), y.first());
    swap(x.second(), y.second());
}

}  // detail_UniquePtr

template <class T>
struct default_delete
{
    default_delete() {}
    template <class U>
        default_delete(const default_delete<U>&,
            typename enable_if_c<detail_UniquePtr::is_convertible<U*, T*>::value>::type* = 0)
        {}

    void operator()(T* ptr) const
    {
        BOOST_STATIC_ASSERT(sizeof(T) > 0);
        delete ptr;
    }
};

template <class T>
struct default_delete<T[]>
{
    void operator()(T* ptr) const
    {
        BOOST_STATIC_ASSERT(sizeof(T) > 0);
        delete [] ptr;
    }

private:

    template <class U> void operator()(U*) const;
};

namespace detail_UniquePtr
{

namespace pointer_type_imp
{

template <class U> static two test(...);
template <class U> static one test(typename U::pointer* = 0);

}  // pointer_type_imp

template <class T>
struct has_pointer_type
{
    static const bool value = sizeof(pointer_type_imp::test<T>(0)) == 1;
};

namespace pointer_type_imp
{

template <class T, class D, bool = has_pointer_type<D>::value>
struct pointer_type
{
    typedef typename D::pointer type;
};

template <class T, class D>
struct pointer_type<T, D, false>
{
    typedef T* type;
};

}  // pointer_type_imp

template <class T, class D>
struct pointer_type
{
    typedef typename pointer_type_imp::pointer_type<T,
        typename boost::remove_reference<D>::type>::type type;
};

}  // detail_UniquePtr

template <class T, class D = default_delete<T> >
class UniquePtr
{
public:
    typedef T element_type;
    typedef D deleter_type;
    typedef typename detail_UniquePtr::pointer_type<element_type, deleter_type>::type pointer;

private:
    detail_UniquePtr::UniquePtr_storage<pointer, deleter_type> ptr_;

    typedef typename add_reference<deleter_type>::type deleter_reference;
    typedef typename add_reference<const deleter_type>::type deleter_const_reference;

    struct nat {int for_bool_;};

    UniquePtr(UniquePtr&);
    UniquePtr& operator=(UniquePtr&);

public:
    operator detail_UniquePtr::rv<UniquePtr>() {return detail_UniquePtr::rv<UniquePtr>(*this);}
    UniquePtr(detail_UniquePtr::rv<UniquePtr> r) : ptr_(r->release(), forward<deleter_type>(r->get_deleter())) {}
    UniquePtr& operator=(detail_UniquePtr::rv<UniquePtr> r)
    {
        reset(r->release());
        ptr_.second() = move(r->get_deleter());
        return *this;
    }

    UniquePtr()
        {
            BOOST_STATIC_ASSERT(!is_reference<deleter_type>::value);
            BOOST_STATIC_ASSERT(!is_pointer<deleter_type>::value);
        }

    explicit UniquePtr(pointer p)
        : ptr_(p)
        {
            BOOST_STATIC_ASSERT(!is_reference<deleter_type>::value);
            BOOST_STATIC_ASSERT(!is_pointer<deleter_type>::value);
        }

    UniquePtr(pointer p, typename mpl::if_<is_reference<D>,
                          volatile typename remove_reference<D>::type&, D>::type d)
        : ptr_(move(p), forward<D>(const_cast<typename add_reference<D>::type>(d))) {}

    template <class U, class E>
        UniquePtr(UniquePtr<U, E> u,
            typename enable_if_c
                <
                !boost::is_array<U>::value &&
                detail_UniquePtr::is_convertible<typename UniquePtr<U>::pointer, pointer>::value &&
                detail_UniquePtr::is_convertible<E, deleter_type>::value &&
                (
                    !is_reference<deleter_type>::value ||
                     is_same<deleter_type, E>::value
                )
                >::type* = 0)
            : ptr_(u.release(), forward<D>(forward<E>(u.get_deleter()))) {}

    ~UniquePtr() {reset();}

    UniquePtr& operator=(int nat::*)
    {
        reset();
        return *this;
    }

    template <class U, class E>
        UniquePtr&
        operator=(UniquePtr<U, E> u)
        {
            reset(u.release());
            ptr_.second() = move(u.get_deleter());
            return *this;
        }

    typename add_reference<T>::type operator*() const {return *get();}
    pointer operator->() const {return get();}
    pointer get() const {return ptr_.first();}
    deleter_reference       get_deleter()       {return ptr_.second();}
    deleter_const_reference get_deleter() const {return ptr_.second();}
    operator int nat::*() const {return get() ? &nat::for_bool_ : 0;}

    void reset(pointer p = pointer())
    {
        pointer t = get();
        if (t != pointer())
            get_deleter()(t);
        ptr_.first() = p;
    }

    pointer release()
    {
        pointer tmp = get();
        ptr_.first() = pointer();
        return tmp;
    }

    void swap(UniquePtr& u) {detail_UniquePtr::swap(ptr_, u.ptr_);}
};

template <class T, class D>
class UniquePtr<T[], D>
{
public:
    typedef T element_type;
    typedef D deleter_type;
    typedef typename detail_UniquePtr::pointer_type<element_type, deleter_type>::type pointer;

private:
    detail_UniquePtr::UniquePtr_storage<pointer, deleter_type> ptr_;

    typedef typename add_reference<deleter_type>::type deleter_reference;
    typedef typename add_reference<const deleter_type>::type deleter_const_reference;

    struct nat {int for_bool_;};

    UniquePtr(UniquePtr&);
    UniquePtr& operator=(UniquePtr&);

public:
    operator detail_UniquePtr::rv<UniquePtr>() {return detail_UniquePtr::rv<UniquePtr>(*this);}
    UniquePtr(detail_UniquePtr::rv<UniquePtr> r) : ptr_(r->release(), forward<deleter_type>(r->get_deleter())) {}
    UniquePtr& operator=(detail_UniquePtr::rv<UniquePtr> r)
    {
        reset(r->release());
        ptr_.second() = move(r->get_deleter());
        return *this;
    }

    UniquePtr()
        {
            BOOST_STATIC_ASSERT(!is_reference<deleter_type>::value);
            BOOST_STATIC_ASSERT(!is_pointer<deleter_type>::value);
        }

    explicit UniquePtr(pointer p)
        : ptr_(p)
        {
            BOOST_STATIC_ASSERT(!is_reference<deleter_type>::value);
            BOOST_STATIC_ASSERT(!is_pointer<deleter_type>::value);
        }

    UniquePtr(pointer p, typename mpl::if_<is_reference<D>,
                          volatile typename remove_reference<D>::type&, D>::type d)
        : ptr_(move(p), forward<D>(const_cast<typename add_reference<D>::type>(d))) {}

    ~UniquePtr() {reset();}

    T& operator[](size_t i) const {return get()[i];}
    pointer get() const {return ptr_.first();}
    deleter_reference       get_deleter()       {return ptr_.second();}
    deleter_const_reference get_deleter() const {return ptr_.second();}
    operator int nat::*() const {return get() ? &nat::for_bool_ : 0;}

    void reset(pointer p = pointer())
    {
        pointer t = get();
        if (t != pointer())
            get_deleter()(t);
        ptr_.first() = p;
    }

    pointer release()
    {
        pointer tmp = get();
        ptr_.first() = pointer();
        return tmp;
    }

    void swap(UniquePtr& u) {detail_UniquePtr::swap(ptr_, u.ptr_);}
private:
    template <class U>
        explicit UniquePtr(U,
            typename enable_if_c<detail_UniquePtr::is_convertible<U, pointer>::value>::type* = 0);

    template <class U>
        UniquePtr(U, typename mpl::if_<is_reference<D>,
                          volatile typename remove_reference<D>::type&, D>::type,
                          typename enable_if_c<detail_UniquePtr::is_convertible<U, pointer>::value>::type* = 0);
};

template<class T, class D>
inline
void
swap(UniquePtr<T, D>& x, UniquePtr<T, D>& y)
{
    x.swap(y);
}

template<class T1, class D1, class T2, class D2>
inline
bool
operator==(const UniquePtr<T1, D1>& x, const UniquePtr<T2, D2>& y)
{
    return x.get() == y.get();
}

template<class T1, class D1, class T2, class D2>
inline
bool
operator!=(const UniquePtr<T1, D1>& x, const UniquePtr<T2, D2>& y)
{
    return !(x == y);
}

template<class T1, class D1, class T2, class D2>
inline
bool
operator<(const UniquePtr<T1, D1>& x, const UniquePtr<T2, D2>& y)
{
    return x.get() < y.get();
}

template<class T1, class D1, class T2, class D2>
inline
bool
operator<=(const UniquePtr<T1, D1>& x, const UniquePtr<T2, D2>& y)
{
    return !(y < x);
}

template<class T1, class D1, class T2, class D2>
inline
bool
operator>(const UniquePtr<T1, D1>& x, const UniquePtr<T2, D2>& y)
{
    return y < x;
}

template<class T1, class D1, class T2, class D2>
inline
bool
operator>=(const UniquePtr<T1, D1>& x, const UniquePtr<T2, D2>& y)
{
    return !(x < y);
}

}  // libMesh

#endif  // HINNANT_UNIQUE_PTR_HPP
