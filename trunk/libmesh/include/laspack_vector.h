// $Id: laspack_vector.h,v 1.3 2003-02-10 22:12:11 benkirk Exp $

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

#if  defined(HAVE_LASPACK) && !defined(USE_COMPLEX_NUMBERS)


// C++ includes

// Local includes
#include "numeric_vector.h"




namespace Laspack {
#include <qvector.h>
#include <operats.h>
}  


// Forward declarations
class LaspackInterface;



/**
 * Laspack vector. Provides a nice interface to the
 * Laspack C-based data structures for parallel vectors.
 *
 * @author Benjamin S. Kirk, 2002
 */

class LaspackVector : public NumericVector
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

//   /**
//    * @returns true if the vector has been initialized,
//    * false otherwise.
//    */
//   bool initialized() const { return (_vec != NULL); };
  
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
//   void init (const NumericVector& V,
// 	     const bool fast=false);

  /**
   * $U(0-N) = s$: fill all components.
   */
  NumericVector & operator= (const Complex s);
    
  /**
   *  $U = V$: copy all components.
   */
  NumericVector & operator= (const NumericVector &V);

  /**
   *  $U = V$: copy all components.
   */
  NumericVector & operator= (const std::vector<Complex> &v);

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
  Complex operator() (const unsigned int i) const;
    
  /**
   * Addition operator.
   * Fast equivalent to \p U.add(1, V).
   */
  NumericVector & operator += (const NumericVector &V);

  /**
   * Subtraction operator.
   * Fast equivalent to \p U.add(-1, V).
   */
  NumericVector & operator -= (const NumericVector &V);
    
  /**
   * v(i) = value
   */
  void set (const unsigned int i, const Complex value);
    
  /**
   * v(i) += value
   */
  void add (const unsigned int i, const Complex value);
    
  /**
   * $U(0-DIM)+=s$.
   * Addition of \p s to all components. Note
   * that \p s is a scalar and not a vector.
   */
  void add (const Complex s);
    
  /**
   * U+=V.
   * Simple vector addition, equal to the
   * \p operator +=.
   */
  void add (const NumericVector& V);

  /**
   * U+=a*V.
   * Simple vector addition, equal to the
   * \p operator +=.
   */
  void add (const Complex a, const NumericVector& v);
  
  /**
   * U+=v where v is a std::vector<Complex> 
   * and you
   * want to specify WHERE to add it
   */
  void add_vector (const std::vector<Complex>& v,
		   const std::vector<unsigned int>& dof_indices);

  /**
   * U+=V where U and V are type 
   * NumericVector and you
   * want to specify WHERE to add
   * the NumericVector V 
   */
  void add_vector (const NumericVector& V,
		   const std::vector<unsigned int>& dof_indices);
  
  /**
   * Scale each element of the
   * vector by the given factor.
   */
  void scale (const Complex factor);
    
  /**
   * Creates a copy of the global vector in the
   * local vector \p v_local.
   */
  void localize (std::vector<Complex>& v_local) const;

  /**
   * Same, but fills a \p NumericVector instead of
   * a \p std::vector.
   */
  void localize (NumericVector& v_local) const;

  /**
   * Creates a local vector \p v_local containing
   * only information relevant to this processor, as
   * defined by the \p send_list.
   */
  void localize (NumericVector& v_local,
		 const std::vector<unsigned int>& send_list) const;
  
  /**
   * Creates a local copy of the global vector in
   * \p v_local only on processor \p proc_id.  By
   * default the data is sent to processor 0.  This method
   * is useful for outputting data from one processor.
   */
  void localize_to_one (std::vector<Complex>& v_local,
			const unsigned int proc_id=0) const;
    
 private:

  /**
   * Actual Laspack vector datatype
   * to hold vector entries
   */
  Laspack::QVector _vec;

  /**
   * Make other Laspack datatypes friends
   */
  friend class LaspackInterface;
};



//----------------------- ----------------------------------
// LaspackVector inline methods
inline
LaspackVector::LaspackVector ()
{};



inline
LaspackVector::LaspackVector (const unsigned int n)
{
  init(n, n, false);
};



inline
LaspackVector::LaspackVector (const unsigned int n,
			      const unsigned int n_local)
{
  init(n, n_local, false);
};



inline
LaspackVector::~LaspackVector ()
{
  clear ();
};



inline
void LaspackVector::init (const unsigned int n,
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
    };

  // create a sequential vector
  using namespace Laspack;

  static int cnt = 0;
  char foo[80];
  sprintf(foo,  "Vec-%d", cnt++); 

  V_Constr(&_vec, const_cast<char*>(foo), n, Normal, True);
    
  _is_initialized = true;
  
  // Optionally zero out all components
  if (fast == false)
    zero ();
  
  return;
};



inline
void LaspackVector::init (const unsigned int n,
			  const bool fast)
{
  init(n,n,fast);
};



inline
void LaspackVector::close ()
{
  assert (initialized());
  
  _is_closed = true;
  
  return;
};



inline
void LaspackVector::clear ()
{
  if (initialized())
    {
      Laspack::V_Destr (&_vec);
    };

  _is_closed = _is_initialized = false;
};



inline
void LaspackVector::zero ()
{
  assert (initialized());

  Laspack::V_SetAllCmp (&_vec, 0.);
};



inline
unsigned int LaspackVector::size () const
{
  assert (initialized());

  return static_cast<unsigned int>(Laspack::V_GetDim(const_cast<Laspack::QVector*>(&_vec)));
};



inline
unsigned int LaspackVector::local_size () const
{
  assert (initialized());
  
  return size();
};



inline
unsigned int LaspackVector::first_local_index () const
{
  assert (initialized());
  
  return 0;
};



inline
unsigned int LaspackVector::last_local_index () const
{
  assert (initialized());
  
  return size();
};



inline
void LaspackVector::set (const unsigned int i, const Complex value)
{
  assert(initialized());
  assert(i<size());
  
  Laspack::V_SetCmp (&_vec, i+1, value);
};



inline
void LaspackVector::add (const unsigned int i, const Complex value)
{
  assert(initialized());
  assert(i<size());
  
  Laspack::V_AddCmp (&_vec, i+1, value);
};



inline
Complex LaspackVector::operator() (const unsigned int i) const
{
  assert (initialized());
  assert ( ((i >= first_local_index()) &&
	    (i <  last_local_index())) );

  
  return static_cast<Complex>(Laspack::V_GetCmp(const_cast<Laspack::QVector*>(&_vec), i+1));
};



inline
void LaspackVector::add_vector (const std::vector<Complex>& v,
				const std::vector<unsigned int>& dof_indices)
{
  assert (!v.empty());
  assert (v.size() == dof_indices.size());
  
  for (unsigned int i=0; i<v.size(); i++)
    add (dof_indices[i], v[i]);
};


inline
void LaspackVector::add_vector (const NumericVector& V,
				const std::vector<unsigned int>& dof_indices)
{
  assert (V.size() == dof_indices.size());

  for (unsigned int i=0; i<V.size(); i++)
    add (dof_indices[i], V(i));
};



inline
Real LaspackVector::min () const
{
  assert (initialized());

  Real min = 1.e30;

  for (unsigned int i=0; i<size(); i++)
    min = std::min (min, (*this)(i));

  return min;
};



inline
Real LaspackVector::max() const
{
  assert (initialized());

  return static_cast<Real>(Laspack::MaxNorm_V(const_cast<Laspack::QVector*>(&_vec)));
};



#endif // #ifdef HAVE_LASPACK
#endif // #ifdef __laspack_vector_h__
