// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA

#ifndef LIBMESH_SAFE_BOOL_H
#define LIBMESH_SAFE_BOOL_H

namespace libMesh
{

/**
 * This is a helper class which can be used to make pre-C++11 operator
 * bool() comparisons safer by making them behave a bit more like they
 * have the "explicit" keyword attached.  The idea is to define your
 * class using the CRTP idiom and implement the non-virtual boolean_test()
 * function, e.g.:
 *
 * class Foo : public safe_bool<Foo>
 * {
 * public:
 *   bool boolean_test() const
 *   {
 *     // Perform actual logic here to determine true/false return value.
 *     return false;
 *   }
 * };
 *
 * The idea is from:
 * https://en.wikibooks.org/wiki/More_C%2B%2B_Idioms/Safe_bool
 */
class safe_bool_base
{
public:
  typedef void (safe_bool_base::*bool_type)() const;
  void this_type_does_not_support_comparisons() const {}
protected:

  safe_bool_base() {}
  safe_bool_base(const safe_bool_base &) {}
  safe_bool_base & operator=(const safe_bool_base &) {return *this;}
  ~safe_bool_base() {}
};



template <typename T>
class safe_bool : private safe_bool_base
{
  // private or protected inheritance is very important here as it triggers the
  // access control violation in main.
public:
  operator bool_type() const
  {
    return (static_cast<const T *>(this))->boolean_test()
      ? &safe_bool_base::this_type_does_not_support_comparisons : 0;
  }
protected:
  ~safe_bool() {}
};



// Equality comparison operators between safe_bool<T> and regular bool
template <typename T>
bool operator==(const safe_bool<T> & lhs, bool b)
{
  return b == static_cast<bool>(lhs);
}

template <typename T>
bool operator==(bool b, const safe_bool<T> & rhs)
{
  return b == static_cast<bool>(rhs);
}



// Disallow equality comparison operators between safe_bool<T> and safe_bool<U>
template <typename T, typename U>
bool operator==(const safe_bool<T> & lhs,
                const safe_bool<U> & /*rhs*/)
{
  lhs.this_type_does_not_support_comparisons();
  return false;
}

template <typename T,typename U>
bool operator!=(const safe_bool<T> & lhs,
                const safe_bool<U> & /*rhs*/)
{
  lhs.this_type_does_not_support_comparisons();
  return false;
}

} // namespace libMesh

#endif // LIBMESH_SAFE_BOOL_H
