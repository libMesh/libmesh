// $Id: laspack_vector.h,v 1.10 2003-03-14 09:56:40 ddreyer Exp $

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




#ifndef __laspack_vector_h__
#define __laspack_vector_h__



#include "mesh_common.h"

#ifdef HAVE_LASPACK




// C++ includes
#include <stdio.h> // Needed by gcc 2.96

// Local includes
#include "numeric_vector.h"
#include "dense_vector.h"




#include <qvector.h>
#include <operats.h>


// Forward declarations
template <typename T> class LaspackInterface;



/**
 * Laspack vector.  Provides a nice interface to the
 * Laspack C-based data structures for serial vectors.
 *
 * @author Benjamin S. Kirk, 2002
 */

template <typename T> 
class LaspackVector : public NumericVector<T>
{
 public:

  /**
   *  Dummy-Constructor. Dimension=0
   */
  LaspackVector ();
  
  /**
   * Constructor. Set dimension to \p n and initialize all elements with zero.
   */
  LaspackVector (const unsigned int n);
    
  /**
   * Constructor. Set local dimension to \p n_local, the global dimension
   * to \p n, and initialize all elements with zero.
   */
  LaspackVector (const unsigned n,
		 const unsigned int n_local);
  
  /**
   * Destructor, deallocates memory. Made virtual to allow
   * for derived classes to behave properly.
   */
  ~LaspackVector ();

  /**
   * Call the assemble functions
   */
  void close (); 

  /**
   * @returns the \p LaspackVector to a pristine state.
   */
  void clear ();
  
  /**
   * Set all entries to zero. Equivalent to \p v = 0, but more obvious and
   * faster. 
   */
  void zero ();    

  /**
   * Change the dimension of the vector to \p N. The reserved memory for
   * this vector remains unchanged if possible, to make things faster, but
   * this may waste some memory, so take this in the back of your head.
   * However, if \p N==0 all memory is freed, i.e. if you want to resize
   * the vector and release the memory not needed, you have to first call
   * \p init(0) and then \p init(N). This cited behaviour is analogous
   * to that of the STL containers.
   *
   * On \p fast==false, the vector is filled by
   * zeros.
   */
    
  void init (const unsigned int N,
	     const unsigned int n_local,
	     const bool         fast=false);
    
  /**
   * call init with n_local = N,
   */
  void init (const unsigned int N,
	     const bool         fast=false);
    
//   /**
//    * Change the dimension to that of the
//    * vector \p V. The same applies as for
//    * the other \p init function.
//    *
//    * The elements of \p V are not copied, i.e.
//    * this function is the same as calling
//    * \p init(V.size(),fast).
//    */
//   void init (const NumericVector<T>& V,
// 	     const bool fast=false);

  /**
   * $U(0-N) = s$: fill all components.
   */
  NumericVector<T> & operator= (const T s);
    
  /**
   *  $U = V$: copy all components.
   */
  NumericVector<T> & operator= (const NumericVector<T> &V);
  
  /**
   *  $U = V$: copy all components.
   */
  LaspackVector<T> & operator= (const LaspackVector<T> &V);

  /**
   *  $U = V$: copy all components.
   */
  NumericVector<T> & operator= (const std::vector<T> &v);

  /**
   * @returns the minimum element in the vector.
   * In case of complex numbers, this returns the minimum
   * Real part.
   */
  Real min () const;
  
  /**
   * @returns the maximum element in the vector.
   * In case of complex numbers, this returns the maximum
   * Real part.
   */
  Real max () const;
  
  /**
   * @returns the $l_1$-norm of the vector, i.e.
   * the sum of the absolute values.
   */
  Real l1_norm () const;

  /**
   * @returns the $l_2$-norm of the vector, i.e.
   * the square root of the sum of the
   * squares of the elements.
   */
  Real l2_norm () const;

  /**
   * @returns the maximum absolute value of the
   * elements of this vector, which is the
   * $l_\infty$-norm of a vector.
   */
  Real linfty_norm () const;

  /**
   * @returns dimension of the vector. This
   * function was formerly called \p n(), but
   * was renamed to get the \p LaspackVector class
   * closer to the C++ standard library's
   * \p std::vector container.
   */
  unsigned int size () const;

  /**
   * @returns the local size of the vector
   * (index_stop-index_start)
   */
  unsigned int local_size() const;

  /**
   * @returns the index of the first vector element
   * actually stored on this processor
   */
  unsigned int first_local_index() const;

  /**
   * @returns the index of the last vector element
   * actually stored on this processor
   */
  unsigned int last_local_index() const;
    
  /**
   * Access components, returns \p U(i).
   */
  T operator() (const unsigned int i) const;
    
  /**
   * Addition operator.
   * Fast equivalent to \p U.add(1, V).
   */
  NumericVector<T> & operator += (const NumericVector<T> &V);

  /**
   * Subtraction operator.
   * Fast equivalent to \p U.add(-1, V).
   */
  NumericVector<T> & operator -= (const NumericVector<T> &V);
    
  /**
   * v(i) = value
   */
  void set (const unsigned int i, const T value);
    
  /**
   * v(i) += value
   */
  void add (const unsigned int i, const T value);
    
  /**
   * $U(0-DIM)+=s$.
   * Addition of \p s to all components. Note
   * that \p s is a scalar and not a vector.
   */
  void add (const T s);
    
  /**
   * U+=V.
   * Simple vector addition, equal to the
   * \p operator +=.
   */
  void add (const NumericVector<T>& V);

  /**
   * U+=a*V.
   * Simple vector addition, equal to the
   * \p operator +=.
   */
  void add (const T a, const NumericVector<T>& v);
  
  /**
   * U+=v where v is a std::vector<T> 
   * and you
   * want to specify WHERE to add it
   */
  void add_vector (const std::vector<T>& v,
		   const std::vector<unsigned int>& dof_indices);

  /**
   * U+=V where U and V are type 
   * NumericVector<T> and you
   * want to specify WHERE to add
   * the NumericVector<T> V 
   */
  void add_vector (const NumericVector<T>& V,
		   const std::vector<unsigned int>& dof_indices);
  
  /**
   * U+=V where U and V are type 
   * DenseVector<T> and you
   * want to specify WHERE to add
   * the DenseVector<T> V 
   */
  void add_vector (const DenseVector<T>& V,
		   const std::vector<unsigned int>& dof_indices);
  
  /**
   * Scale each element of the
   * vector by the given factor.
   */
  void scale (const T factor);
    
  /**
   * Creates a copy of the global vector in the
   * local vector \p v_local.
   */
  void localize (std::vector<T>& v_local) const;

  /**
   * Same, but fills a \p NumericVector<T> instead of
   * a \p std::vector.
   */
  void localize (NumericVector<T>& v_local) const;

  /**
   * Creates a local vector \p v_local containing
   * only information relevant to this processor, as
   * defined by the \p send_list.
   */
  void localize (NumericVector<T>& v_local,
		 const std::vector<unsigned int>& send_list) const;
  
  /**
   * Creates a local copy of the global vector in
   * \p v_local only on processor \p proc_id.  By
   * default the data is sent to processor 0.  This method
   * is useful for outputting data from one processor.
   */
  void localize_to_one (std::vector<T>& v_local,
			const unsigned int proc_id=0) const;
    
 private:

  /**
   * Actual Laspack vector datatype
   * to hold vector entries
   */
  QVector _vec;

  /**
   * Make other Laspack datatypes friends
   */
  friend class LaspackInterface<T>;
};



//----------------------- ----------------------------------
// LaspackVector inline methods
template <typename T>
inline
LaspackVector<T>::LaspackVector ()
{
}



template <typename T> 
inline
LaspackVector<T>::LaspackVector (const unsigned int n)
{
  init(n, n, false);
}



template <typename T> 
inline
LaspackVector<T>::LaspackVector (const unsigned int n,
				 const unsigned int n_local)
{
  init(n, n_local, false);
}



template <typename T> 
inline
LaspackVector<T>::~LaspackVector ()
{
  clear ();
}



template <typename T> 
inline
void LaspackVector<T>::init (const unsigned int n,
			     const unsigned int n_local,
			     const bool fast)
{
  // Laspack vectors only for serial cases.
  assert (n == n_local);

  // Only for uninitialized vectors
  if (initialized())
    {
      std::cerr << "ERROR: Vector already initialized!"
		<< std::endl;      
      error();
    }

  // create a sequential vector

  static int cnt = 0;
  char foo[80];
  sprintf(foo,  "Vec-%d", cnt++); 

  V_Constr(&_vec, const_cast<char*>(foo), n, Normal, _LPTrue);
    
  _is_initialized = true;
  
  // Optionally zero out all components
  if (fast == false)
    zero ();
  
  return;
}



template <typename T> 
inline
void LaspackVector<T>::init (const unsigned int n,
			  const bool fast)
{
  init(n,n,fast);
}



template <typename T> 
inline
void LaspackVector<T>::close ()
{
  assert (initialized());
  
  _is_closed = true;
  
  return;
}



template <typename T> 
inline
void LaspackVector<T>::clear ()
{
  if (initialized())
    {
      V_Destr (&_vec);
    }

  _is_closed = _is_initialized = false;
}



template <typename T> inline
void LaspackVector<T>::zero ()
{
  assert (initialized());

  V_SetAllCmp (&_vec, 0.);
}



template <typename T> 
inline
unsigned int LaspackVector<T>::size () const
{
  assert (initialized());

  return static_cast<unsigned int>(V_GetDim(const_cast<QVector*>(&_vec)));
}



template <typename T> 
inline
unsigned int LaspackVector<T>::local_size () const
{
  assert (initialized());
  
  return size();
}



template <typename T> 
inline
unsigned int LaspackVector<T>::first_local_index () const
{
  assert (initialized());
  
  return 0;
}



template <typename T> 
inline
unsigned int LaspackVector<T>::last_local_index () const
{
  assert (initialized());
  
  return size();
}



template <typename T> 
inline
void LaspackVector<T>::set (const unsigned int i, const T value)
{
  assert(initialized());
  assert(i<size());
  
  V_SetCmp (&_vec, i+1, value);
}



template <typename T> 
inline
void LaspackVector<T>::add (const unsigned int i, const T value)
{
  assert(initialized());
  assert(i<size());
  
  V_AddCmp (&_vec, i+1, value);
}



template <typename T> 
inline
T LaspackVector<T>::operator() (const unsigned int i) const
{
  assert (initialized());
  assert ( ((i >= first_local_index()) &&
	    (i <  last_local_index())) );

  
  return static_cast<T>(V_GetCmp(const_cast<QVector*>(&_vec), i+1));
}



template <typename T> 
inline
void LaspackVector<T>::add_vector (const std::vector<T>& v,
				   const std::vector<unsigned int>& dof_indices)
{
  assert (!v.empty());
  assert (v.size() == dof_indices.size());
  
  for (unsigned int i=0; i<v.size(); i++)
    add (dof_indices[i], v[i]);
}



template <typename T> 
inline
void LaspackVector<T>::add_vector (const NumericVector<T>& V,
				   const std::vector<unsigned int>& dof_indices)
{
  assert (V.size() == dof_indices.size());

  for (unsigned int i=0; i<V.size(); i++)
    add (dof_indices[i], V(i));
}



template <typename T> 
inline
void LaspackVector<T>::add_vector (const DenseVector<T>& V,
				   const std::vector<unsigned int>& dof_indices)
{
  assert (V.size() == dof_indices.size());

  for (unsigned int i=0; i<V.size(); i++)
    add (dof_indices[i], V(i));
}



template <typename T> 
inline
Real LaspackVector<T>::max() const
{
  assert (initialized());

  return static_cast<Real>(MaxNorm_V(const_cast<QVector*>(&_vec)));
}



#endif // #ifdef HAVE_LASPACK
#endif // #ifdef __laspack_vector_h__
