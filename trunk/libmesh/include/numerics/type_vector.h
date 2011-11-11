// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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



#ifndef __type_vector_h__
#define __type_vector_h__

// C++ includes
#include <cmath>

// Local includes
#include "libmesh_common.h"
#include "compare_types.h"

namespace libMesh
{

// Forward declaration for friend class
template <typename T>
class TypeTensor;


/**
 * This class defines a vector in \p LIBMESH_DIM dimensional space of type T.
 * T may either be Real or Complex.  The default constructor for
 * this class is protected, suggesting that you should not instantiate
 * one of these directly.  Instead use one of the derived types: \p Point
 * for a real-valued point in LIBMESH_DIM-space, or \p SpaceVector for a real
 * or complex-valued vector in LIBMESH_DIM-space.
 *
 * \author Benjamin S. Kirk, 2003. 
 */

template <typename T>
class TypeVector
{
template <typename T2>
friend class TypeVector;

friend class TypeTensor<T>;

protected:

  /**
   * Constructor.  By default sets all entries to 0.  Gives the vector 0 in
   * \p LIBMESH_DIM dimensions.
   */
  TypeVector  (const T x=0.,
	       const T y=0.,
	       const T z=0.);
  
public:

  /**
   * Copy-constructor.
   */
  template <typename T2>
  TypeVector (const TypeVector<T2>& p);
  
  /**
   * Destructor.
   */ 
  ~TypeVector ();

  /**
   * Assign to a vector without creating a temporary.
   */
  template <typename T2>
  void assign (const TypeVector<T2> &);

  /**
   * Return the \f$ i^{th} \f$ element of the vector.
   */
  const T & operator () (const unsigned int i) const;

  /**
   * Return a writeable reference to the \f$ i^{th} \f$ element of the vector.
   */
  T & operator () (const unsigned int i);
  
  /**
   * Add two vectors.
   */
  template <typename T2>
  TypeVector<typename CompareTypes<T, T2>::supertype> 
  operator + (const TypeVector<T2> &) const;

  /**
   * Add to this vector.
   */
  template <typename T2>
  const TypeVector<T> & operator += (const TypeVector<T2> &);
  
  /**
   * Add to this vector without creating a temporary.
   */
  template <typename T2>
  void add (const TypeVector<T2> &); 
  
  /**
   * Add a scaled value to this vector without
   * creating a temporary.
   */
  template <typename T2>
  void add_scaled (const TypeVector<T2> &, const T); 
  
  /**
   * Subtract two vectors.
   */
  template <typename T2>
  TypeVector<typename CompareTypes<T, T2>::supertype> 
  operator - (const TypeVector<T2> &) const;

  /**
   * Subtract from this vector.
   */
  template <typename T2>
  const TypeVector<T> & operator -= (const TypeVector<T2> &);

  /**
   * Subtract from this vector without creating a temporary.
   */
  template <typename T2>
  void subtract (const TypeVector<T2> &); 
  
  /**
   * Subtract a scaled value from this vector without
   * creating a temporary.
   */
  template <typename T2>
  void subtract_scaled (const TypeVector<T2> &, const T); 
  
  /**
   * Return the opposite of a vector
   */
  TypeVector<T> operator - () const;
  
  /**
   * Multiply a vector by a number, i.e. scale.
   */
  template <typename Scalar>
  typename boostcopy::enable_if_c<
    ScalarTraits<Scalar>::value,
    TypeVector<typename CompareTypes<T, Scalar>::supertype> >::type
  operator * (const Scalar) const;

  /**
   * Multiply this vector by a number, i.e. scale.
   */
  const TypeVector<T> & operator *= (const T);
  
  /**
   * Divide a vector by a number, i.e. scale.
   */
  template <typename Scalar>
  typename boostcopy::enable_if_c<
    ScalarTraits<Scalar>::value,
    TypeVector<typename CompareTypes<T, Scalar>::supertype> >::type
  operator / (const Scalar) const;

  /**
   * Divide this vector by a number, i.e. scale.
   */
  const TypeVector<T> & operator /= (const T);

  /**
   * Multiply 2 vectors together, i.e. dot-product.
   * The vectors may be of different types.
   */
  template <typename T2>
  typename CompareTypes<T, T2>::supertype
  operator * (const TypeVector<T2> &) const;

  /**
   * Cross 2 vectors together, i.e. cross-product.
   */
  template <typename T2>
  TypeVector<typename CompareTypes<T, T2>::supertype>
  cross(const TypeVector<T2> &) const;

  /**
   * Think of a vector as a \p dim dimensional vector.  This
   * will return a unit vector aligned in that direction.
   */
  TypeVector<T> unit() const;

  /**
   * Returns the magnitude of the vector, i.e. the square-root of the
   * sum of the elements squared.
   */
  Real size() const;

  /**
   * Returns the magnitude of the vector squared, i.e. the sum of the element
   * magnitudes squared.
   */
  Real size_sq() const;

  /**
   * Zero the vector in any dimension.
   */
  void zero();

  /**
   * @returns \p true iff two vectors occupy approximately the same
   * physical location in space, to within a relative tolerance of \p tol.
   */
  bool relative_fuzzy_equals(const TypeVector<T>& rhs, Real tol = TOLERANCE) const;
  
  /**
   * @returns \p true iff two vectors occupy approximately the same
   * physical location in space, to within an absolute tolerance of \p tol.
   */
  bool absolute_fuzzy_equals(const TypeVector<T>& rhs, Real tol = TOLERANCE) const;
  
  /**
   * @returns \p true iff two vectors occupy approximately the same
   * physical location in space, to within an absolute tolerance of \p TOLERANCE.
   */
  bool operator == (const TypeVector<T>& rhs) const;
  
  /**
   * @returns \p true iff two vectors do not occupy approximately the same
   * physical location in space.
   */
  bool operator != (const TypeVector<T>& rhs) const;
  
  /**
   * @returns \p true if this vector is "less"
   * than another.  Useful for sorting.
   * Also used for choosing some arbitrary basis function
   * orientations
   */
  bool operator < (const TypeVector<T>& rhs) const;
  
  /**
   * @returns \p true if this vector is "greater"
   * than another.  Useful for sorting.
   * Also used for choosing some arbitrary basis function
   * orientations
   */
  bool operator > (const TypeVector<T>& rhs) const;

  /**
   * Formatted print, by default to \p libMesh::out.
   */
  void print(std::ostream& os = libMesh::out) const;

  /**
   * Formatted print as above but allows you to do
   * Point p(1,2,3);
   * std::cout << p << std::endl;
   */
  friend std::ostream& operator << (std::ostream& os, const TypeVector<T>& t)
  {
    t.print(os);
    return os;
  }
  
  /**
   * Unformatted print to the stream \p out.  Simply prints the elements
   * of the vector separated by spaces.  Optionally prints a newline,
   * which it does by default.
   */ 
  void write_unformatted (std::ostream &out, const bool newline = true) const;
    
 protected:

  /**
   * The coordinates of the \p TypeVector
   */
  T _coords[LIBMESH_DIM];
};



//------------------------------------------------------
// Inline functions
template <typename T>
inline
TypeVector<T>::TypeVector (const T x,
			   const T y,
			   const T z)
{
  _coords[0] = x;

  if (LIBMESH_DIM > 1)
    {
      _coords[1] = y;

      if (LIBMESH_DIM == 3)
	_coords[2] = z;
      else
        libmesh_assert(z == 0.);
    }
  else
    libmesh_assert(y == 0.);
}



template <typename T>
template <typename T2>
inline
TypeVector<T>::TypeVector (const TypeVector<T2> &p)
{
  // copy the nodes from vector p to me
  for (unsigned int i=0; i<LIBMESH_DIM; i++)
    _coords[i] = p._coords[i];
}



template <typename T>
inline
TypeVector<T>::~TypeVector ()
{
}



template <typename T>
template <typename T2>
inline
void TypeVector<T>::assign (const TypeVector<T2> &p)
{
  for (unsigned int i=0; i<LIBMESH_DIM; i++)
    _coords[i] = p._coords[i];
}



template <typename T>
inline
const T & TypeVector<T>::operator () (const unsigned int i) const
{
  libmesh_assert (i<LIBMESH_DIM);

  return _coords[i];
}



template <typename T>
inline
T & TypeVector<T>::operator () (const unsigned int i)
{
  libmesh_assert (i<LIBMESH_DIM);
  
  return _coords[i];
}



template <typename T>
template <typename T2>
inline
TypeVector<typename CompareTypes<T, T2>::supertype> 
TypeVector<T>::operator + (const TypeVector<T2> &p) const
{
  typedef typename CompareTypes<T, T2>::supertype TS;
#if LIBMESH_DIM == 1
  return TypeVector<TS> (_coords[0] + p._coords[0]);
#endif

#if LIBMESH_DIM == 2 
  return TypeVector<TS> (_coords[0] + p._coords[0],
                         _coords[1] + p._coords[1]);
#endif

#if LIBMESH_DIM == 3
  return TypeVector<TS> (_coords[0] + p._coords[0],
                         _coords[1] + p._coords[1],
                         _coords[2] + p._coords[2]);
#endif
	       
}



template <typename T>
template <typename T2>
inline
const TypeVector<T> & TypeVector<T>::operator += (const TypeVector<T2> &p)
{
  this->add (p);

  return *this;
}



template <typename T>
template <typename T2>
inline
void TypeVector<T>::add (const TypeVector<T2> &p)
{
#if LIBMESH_DIM == 1
  _coords[0] += p._coords[0];
#endif
  
#if LIBMESH_DIM == 2
  _coords[0] += p._coords[0];
  _coords[1] += p._coords[1];
#endif
  
#if LIBMESH_DIM == 3
  _coords[0] += p._coords[0];
  _coords[1] += p._coords[1];
  _coords[2] += p._coords[2];
#endif

}



template <typename T>
template <typename T2>
inline
void TypeVector<T>::add_scaled (const TypeVector<T2> &p, const T factor)
{
#if LIBMESH_DIM == 1
  _coords[0] += factor*p(0);
#endif
  
#if LIBMESH_DIM == 2
  _coords[0] += factor*p(0);
  _coords[1] += factor*p(1);
#endif
  
#if LIBMESH_DIM == 3
  _coords[0] += factor*p(0);
  _coords[1] += factor*p(1);
  _coords[2] += factor*p(2);
#endif

}



template <typename T>
template <typename T2>
inline
TypeVector<typename CompareTypes<T, T2>::supertype> 
TypeVector<T>::operator - (const TypeVector<T2> &p) const
{
  typedef typename CompareTypes<T, T2>::supertype TS;

#if LIBMESH_DIM == 1
  return TypeVector<TS>(_coords[0] - p._coords[0]);
#endif

#if LIBMESH_DIM == 2 
  return TypeVector<TS>(_coords[0] - p._coords[0],
		        _coords[1] - p._coords[1]);
#endif

#if LIBMESH_DIM == 3
  return TypeVector<TS>(_coords[0] - p._coords[0],
		        _coords[1] - p._coords[1],
		        _coords[2] - p._coords[2]);
#endif

}



template <typename T>
template <typename T2>
inline
const TypeVector<T> & TypeVector<T>::operator -= (const TypeVector<T2> &p)
{
  this->subtract (p);

  return *this;
}



template <typename T>
template <typename T2>
inline
void TypeVector<T>::subtract (const TypeVector<T2>& p)
{
  for (unsigned int i=0; i<LIBMESH_DIM; i++)
    _coords[i] -= p._coords[i];
}



template <typename T>
template <typename T2>
inline
void TypeVector<T>::subtract_scaled (const TypeVector<T2> &p, const T factor)
{
  for (unsigned int i=0; i<LIBMESH_DIM; i++)
    _coords[i] -= factor*p(i);
}



template <typename T>
inline
TypeVector<T> TypeVector<T>::operator - () const
{
  
#if LIBMESH_DIM == 1
  return TypeVector(-_coords[0]);
#endif

#if LIBMESH_DIM == 2 
  return TypeVector(-_coords[0],
		    -_coords[1]);
#endif

#if LIBMESH_DIM == 3
  return TypeVector(-_coords[0],
		    -_coords[1], 
		    -_coords[2]);
#endif
  
}



template <typename T>
template <typename Scalar>
inline
typename boostcopy::enable_if_c<
  ScalarTraits<Scalar>::value,
  TypeVector<typename CompareTypes<T, Scalar>::supertype> >::type
TypeVector<T>::operator * (const Scalar factor) const
{
  typedef typename CompareTypes<T, Scalar>::supertype TS;

#if LIBMESH_DIM == 1
  return TypeVector<TS>(_coords[0]*factor);
#endif
  
#if LIBMESH_DIM == 2 
  return TypeVector<TS>(_coords[0]*factor,
		        _coords[1]*factor);
#endif
  
#if LIBMESH_DIM == 3
  return TypeVector<TS>(_coords[0]*factor,
		        _coords[1]*factor, 
		        _coords[2]*factor);
#endif  
}



template <typename T, typename Scalar>
inline
typename boostcopy::enable_if_c<
  ScalarTraits<Scalar>::value,
  TypeVector<typename CompareTypes<T, Scalar>::supertype> >::type
operator * (const Scalar factor,
            const TypeVector<T> &v)
{
  return v * factor;
}



template <typename T>
inline
const TypeVector<T> & TypeVector<T>::operator *= (const T factor)
{
#if LIBMESH_DIM == 1
  _coords[0] *= factor;
#endif
  
#if LIBMESH_DIM == 2
  _coords[0] *= factor;
  _coords[1] *= factor;
#endif
  
#if LIBMESH_DIM == 3
  _coords[0] *= factor;
  _coords[1] *= factor;
  _coords[2] *= factor;
#endif

  return *this;
}



template <typename T>
template <typename Scalar>
inline
typename boostcopy::enable_if_c<
  ScalarTraits<Scalar>::value,
  TypeVector<typename CompareTypes<T, Scalar>::supertype> >::type
TypeVector<T>::operator / (const Scalar factor) const
{
  libmesh_assert (factor != static_cast<T>(0.));

  typedef typename CompareTypes<T, Scalar>::supertype TS;
  
#if LIBMESH_DIM == 1
  return TypeVector<TS>(_coords[0]/factor);
#endif
  
#if LIBMESH_DIM == 2 
  return TypeVector<TS>(_coords[0]/factor,
		        _coords[1]/factor);
#endif
  
#if LIBMESH_DIM == 3
  return TypeVector<TS>(_coords[0]/factor,
		        _coords[1]/factor, 
		        _coords[2]/factor);
#endif
  
}




template <typename T>
inline
const TypeVector<T> &
TypeVector<T>::operator /= (const T factor)
{
  libmesh_assert (factor != static_cast<T>(0.));
  
  for (unsigned int i=0; i<LIBMESH_DIM; i++)
    _coords[i] /= factor;

  return *this;
}




template <typename T>
template <typename T2>
inline
typename CompareTypes<T, T2>::supertype
TypeVector<T>::operator * (const TypeVector<T2> &p) const
{
#if LIBMESH_DIM == 1
  return _coords[0]*p._coords[0];
#endif
  
#if LIBMESH_DIM == 2
  return (_coords[0]*p._coords[0] +
	  _coords[1]*p._coords[1]);
#endif
  
#if LIBMESH_DIM == 3
  return (_coords[0]*p(0) +
	  _coords[1]*p(1) +
	  _coords[2]*p(2));
#endif
}



template <typename T>
template <typename T2>
TypeVector<typename CompareTypes<T, T2>::supertype>
TypeVector<T>::cross(const TypeVector<T2>& p) const
{
  typedef typename CompareTypes<T, T2>::supertype TS;
  libmesh_assert (LIBMESH_DIM == 3);

  // |     i          j          k    |
  // |(*this)(0) (*this)(1) (*this)(2)|
  // |   p(0)       p(1)       p(2)   |

  return TypeVector<TS>( _coords[1]*p._coords[2] - _coords[2]*p._coords[1],
                        -_coords[0]*p._coords[2] + _coords[2]*p._coords[0],
                         _coords[0]*p._coords[1] - _coords[1]*p._coords[0]);
}



template <typename T>
inline
Real TypeVector<T>::size() const
{
  return std::sqrt(this->size_sq());  
}



template <typename T>
inline
void TypeVector<T>::zero()
{
  for (unsigned int i=0; i<LIBMESH_DIM; i++)
    _coords[i] = 0.;
}



template <typename T>
inline
Real TypeVector<T>::size_sq() const
{
#if LIBMESH_DIM == 1
  return (libmesh_norm(_coords[0]));
#endif
  
#if LIBMESH_DIM == 2
  return (libmesh_norm(_coords[0]) +
	  libmesh_norm(_coords[1]));
#endif
  
#if LIBMESH_DIM == 3
  return (libmesh_norm(_coords[0]) +
	  libmesh_norm(_coords[1]) + 
	  libmesh_norm(_coords[2]));
#endif
}



template <typename T>
inline
bool TypeVector<T>::absolute_fuzzy_equals(const TypeVector<T>& rhs, Real tol) const
{
#if LIBMESH_DIM == 1
  return (std::abs(_coords[0] - rhs._coords[0])
	  <= tol);
#endif

#if LIBMESH_DIM == 2
  return (std::abs(_coords[0] - rhs._coords[0]) +
	  std::abs(_coords[1] - rhs._coords[1])
	  <= tol);
#endif

#if LIBMESH_DIM == 3
  return (std::abs(_coords[0] - rhs._coords[0]) +
	  std::abs(_coords[1] - rhs._coords[1]) +
	  std::abs(_coords[2] - rhs._coords[2])
	  <= tol);
#endif
}



template <typename T>
inline
bool TypeVector<T>::relative_fuzzy_equals(const TypeVector<T>& rhs, Real tol) const
{
#if LIBMESH_DIM == 1
  return this->absolute_fuzzy_equals(rhs, tol *
				     (std::abs(_coords[0]) + std::abs(rhs._coords[0])));
#endif

#if LIBMESH_DIM == 2
  return this->absolute_fuzzy_equals(rhs, tol * 
				     (std::abs(_coords[0]) + std::abs(rhs._coords[0]) +
				      std::abs(_coords[1]) + std::abs(rhs._coords[1])));
#endif

#if LIBMESH_DIM == 3
  return this->absolute_fuzzy_equals(rhs, tol *
				     (std::abs(_coords[0]) + std::abs(rhs._coords[0]) +
				      std::abs(_coords[1]) + std::abs(rhs._coords[1]) +
				      std::abs(_coords[2]) + std::abs(rhs._coords[2])));
#endif
}



template <typename T>
inline
bool TypeVector<T>::operator == (const TypeVector<T>& rhs) const
{
#if LIBMESH_DIM == 1
  return (_coords[0] == rhs._coords[0]);
#endif

#if LIBMESH_DIM == 2
  return (_coords[0] == rhs._coords[0] &&
	  _coords[1] == rhs._coords[1]);
#endif

#if LIBMESH_DIM == 3
  return (_coords[0] == rhs._coords[0] &&
	  _coords[1] == rhs._coords[1] &&
	  _coords[2] == rhs._coords[2]);
#endif
}



template <>
inline
bool TypeVector<Real>::operator != (const TypeVector<Real>& rhs) const
{
  return (!(*this == rhs));
}

} // namespace libMesh

#endif // #define __type_vector_h__
