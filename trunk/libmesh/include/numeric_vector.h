// $Id: numeric_vector.h,v 1.15 2003-05-15 23:34:34 benkirk Exp $

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



#ifndef __numeric_vector_h__
#define __numeric_vector_h__


// C++ includes
#include <vector>

// Local includes
#include "mesh_common.h"
#include "auto_ptr.h"
#include "enum_solver_package.h"
#include "reference_counted_object.h"
#include "libmesh.h"


// forward declarations
template <typename T> class NumericVector;
template <typename T> class DenseVector;
template <typename T> class SparseMatrix;


/**
 * Numeric vector. Provides a uniform interface
 * to vector storage schemes for different linear
 * algebra libraries.
 *
 * @author Benjamin S. Kirk, 2003
 */
template <typename T>
class NumericVector : public ReferenceCountedObject<NumericVector<T> >
{
public:

  /**
   *  Dummy-Constructor. Dimension=0
   */
  NumericVector ();
    
  /**
   * Constructor. Set dimension to \p n and initialize all elements with zero.
   */
  NumericVector (const unsigned int n);
    
  /**
   * Constructor. Set local dimension to \p n_local, the global dimension
   * to \p n, and initialize all elements with zero.
   */
  NumericVector (const unsigned n,
		 const unsigned int n_local);
    
public:

  /**
   * Destructor, deallocates memory. Made virtual to allow
   * for derived classes to behave properly.
   */
  virtual ~NumericVector ();

  /**
   * Builds a \p NumericVector using the linear solver package specified by
   * \p solver_package
   */
  static AutoPtr<NumericVector<T> > build(const SolverPackage solver_package = libMesh::default_solver_package());
  
  /**
   * @returns true if the vector has been initialized,
   * false otherwise.
   */
  virtual bool initialized() const { return _is_initialized; }

  /**
   * @returns true if the vector is closed and ready for
   * computation, false otherwise.
   */
  virtual bool closed() const { return _is_closed; }
  
  /**
   * Call the assemble functions
   */
  virtual void close () = 0; 

  /**
   * @returns the \p NumericVector<T> to a pristine state.
   */
  virtual void clear ();

  /**
   * Set all entries to zero. Equivalent to \p v = 0, but more obvious and
   * faster. 
   */
  virtual void zero () = 0;    

  /**
   * Creates a copy of this vector and returns it in an \p AutoPtr.
   * This must be overloaded in the derived classes.
   */
  virtual AutoPtr<NumericVector<T> > clone () const = 0;
  
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
    
  virtual void init (const unsigned int,
		     const unsigned int,
		     const bool = false) {}
  
  /**
   * call init with n_local = N,
   */
  virtual void init (const unsigned int,
		     const bool = false) {}
    
  //   /**
  //    * Change the dimension to that of the
  //    * vector \p V. The same applies as for
  //    * the other \p init function.
  //    *
  //    * The elements of \p V are not copied, i.e.
  //    * this function is the same as calling
  //    * \p init(V.size(),fast).
  //    */
  //   virtual void init (const NumericVector<T>&,
  // 		     const bool = false) {}

  /**
   * $U(0-N) = s$: fill all components.
   */
  virtual NumericVector<T> & operator= (const T s);
  
  /**
   *  $U = V$: copy all components.
   */
  virtual NumericVector<T> & operator= (const NumericVector<T> &V);

  /**
   *  $U = V$: copy all components.
   */
  virtual NumericVector<T> & operator= (const std::vector<T> &v);

  /**
   * @returns the minimum element in the vector.
   * In case of complex numbers, this returns the minimum
   * Real part.
   */
  virtual Real min () const = 0;
  
  /**
   * @returns the maximum element in the vector.
   * In case of complex numbers, this returns the maximum
   * Real part.
   */
  virtual Real max () const = 0;
  
  /**
   * @returns the $l_1$-norm of the vector, i.e.
   * the sum of the absolute values.
   */
  virtual Real l1_norm () const = 0;

  /**
   * @returns the $l_2$-norm of the vector, i.e.
   * the square root of the sum of the
   * squares of the elements.
   */
  virtual Real l2_norm () const = 0;

  /**
   * @returns the maximum absolute value of the
   * elements of this vector, which is the
   * $l_\infty$-norm of a vector.
   */
  virtual Real linfty_norm () const = 0;

  /**
   * @returns dimension of the vector. This
   * function was formerly called \p n(), but
   * was renamed to get the \p NumericVector<T> class
   * closer to the C++ standard library's
   * \p std::vector container.
   */
  virtual unsigned int size () const = 0;

  /**
   * @returns the local size of the vector
   * (index_stop-index_start)
   */
  virtual unsigned int local_size() const = 0;

  /**
   * @returns the index of the first vector element
   * actually stored on this processor.  Hint: the
   * minimum for this index is \p 0.
   */
  virtual unsigned int first_local_index() const = 0;

  /**
   * @returns the index+1 of the last vector element
   * actually stored on this processor.  Hint: the
   * maximum for this index is \p size().
   */
  virtual unsigned int last_local_index() const = 0;
    
  /**
   * Access components, returns \p U(i).
   */
  virtual T operator() (const unsigned int i) const = 0;
    
  /**
   * Addition operator.
   * Fast equivalent to \p U.add(1, V).
   */
  virtual NumericVector<T> & operator += (const NumericVector<T> &V) = 0;

  /**
   * Subtraction operator.
   * Fast equivalent to \p U.add(-1, V).
   */
  virtual NumericVector<T> & operator -= (const NumericVector<T> &V) = 0;
    
  /**
   * v(i) = value
   */
  virtual void set (const unsigned int i, const T value) = 0;
    
  /**
   * v(i) += value
   */
  virtual void add (const unsigned int i, const T value) = 0;
    
  /**
   * $U(0-DIM)+=s$.
   * Addition of \p s to all components. Note
   * that \p s is a scalar and not a vector.
   */
  virtual void add (const T s) = 0;
    
  /**
   * U+=V.
   * Simple vector addition, equal to the
   * \p operator +=.
   */
  virtual void add (const NumericVector<T>& V) = 0;

  /**
   * U+=a*V.
   * Simple vector addition, equal to the
   * \p operator +=.
   */
  virtual void add (const T a, const NumericVector<T>& v) = 0;

  /**
   * U+=v where v is a std::vector<T> 
   * and you
   * want to specify WHERE to add it
   */
  virtual void add_vector (const std::vector<T>& v,
			   const std::vector<unsigned int>& dof_indices) = 0;

  /**
   * U+=V where U and V are type 
   * NumericVector<T> and you
   * want to specify WHERE to add
   * the NumericVector<T> V 
   */
  virtual void add_vector (const NumericVector<T>& V,
			   const std::vector<unsigned int>& dof_indices) = 0;

  /**
   * U+=U+A*V.
   * Add the product of a Sparse matrix \p A
   * and a Numeric vector \p V to this Numeric vector.
   */
  virtual void add_vector (const NumericVector<T> &,
			   SparseMatrix<T> &) = 0;
      
  /**
   * U+=V where U and V are type 
   * DenseVector<T> and you
   * want to specify WHERE to add
   * the DenseVector<T> V 
   */
  virtual void add_vector (const DenseVector<T>& V,
			   const std::vector<unsigned int>& dof_indices) = 0;
    
  /**
   * Scale each element of the
   * vector by the given factor.
   */
  virtual void scale (const T factor) = 0;

  /**
   * Creates a copy of the global vector in the
   * local vector \p v_local.
   */
  virtual void localize (std::vector<T>& v_local) const = 0;

  /**
   * Same, but fills a \p NumericVector<T> instead of
   * a \p std::vector.
   */
  virtual void localize (NumericVector<T>& v_local) const = 0;

  /**
   * Creates a local vector \p v_local containing
   * only information relevant to this processor, as
   * defined by the \p send_list.
   */
  virtual void localize (NumericVector<T>& v_local,
			 const std::vector<unsigned int>& send_list) const = 0;

  /**
   * Creates a local copy of the global vector in
   * \p v_local only on processor \p proc_id.  By
   * default the data is sent to processor 0.  This method
   * is useful for outputting data from one processor.
   */
  virtual void localize_to_one (std::vector<T>& v_local,
				const unsigned int proc_id=0) const = 0;
    
  /**
   * @returns \p -1 when \p this is equivalent to \p other_vector,
   * up to the given \p threshold.  When differences occur,
   * the return value contains the first index where
   * the difference exceeded the threshold.  When
   * no threshold is given, the \p libMesh \p TOLERANCE
   * is used.
   */
  virtual int compare (const NumericVector<T> &other_vector,
		       const Real threshold = TOLERANCE) const;

  /**
   * Prints the contents of the vector to the screen.
   */
  virtual void print() const;
    
protected:
  
  /**
   * Flag to see if the Numeric
   * assemble routines have been called yet
   */
  bool _is_closed;
  
  /**
   * Flag to tell if init 
   * has been called yet
   */
  bool _is_initialized;
};


/*----------------------- Inline functions ----------------------------------*/



template <typename T>
inline
NumericVector<T>::NumericVector () :
  _is_closed(false),
  _is_initialized(false)
{}



template <typename T>
inline
NumericVector<T>::NumericVector (const unsigned int n) :
  _is_closed(false),
  _is_initialized(false)
{
  init(n, n, false);
}



template <typename T>
inline
NumericVector<T>::NumericVector (const unsigned int n,
				 const unsigned int n_local) :
  _is_closed(false),
  _is_initialized(false)
{
  init(n, n_local, false);
}



template <typename T>
inline
NumericVector<T>::~NumericVector ()
{
  clear ();
}



template <typename T>
inline
NumericVector<T> & NumericVector<T>::operator= (const T) 
{
  //  error();

  return *this;
}



template <typename T>
inline
NumericVector<T> & NumericVector<T>::operator= (const NumericVector<T>&) 
{
  //  error();

  return *this;
}



template <typename T>
inline
NumericVector<T> & NumericVector<T>::operator= (const std::vector<T>&) 
{
  //  error();

  return *this;
}



template <typename T>
inline
void NumericVector<T>::clear ()
{
  _is_closed      = false;
  _is_initialized = false;
}



// Full specialization of the print() member for complex
// variables.  This must precede the non-specialized
// version, at least according to icc v7.1
template <>
inline
void NumericVector<Complex>::print() const
{
  assert (initialized());
  std::cout << "Size\tglobal =  " << size()
	    << "\t\tlocal =  " << local_size() << std::endl;
  
  // std::complex<>::operator<<() is defined, but use this form
  std::cout << "#\tReal part\t\tImaginary part" << std::endl;
  for (unsigned int i=first_local_index(); i<last_local_index(); i++)
    std::cout << i << "\t" 
	      << (*this)(i).real() << "\t\t" 
	      << (*this)(i).imag() << std::endl;
}



template <typename T>
inline
void NumericVector<T>::print() const
{
  assert (initialized());
  std::cout << "Size\tglobal =  " << size()
	    << "\t\tlocal =  " << local_size() << std::endl;

  std::cout << "#\tValue" << std::endl;
  for (unsigned int i=first_local_index(); i<last_local_index(); i++)
    std::cout << i << "\t" << (*this)(i) << std::endl;
}


#endif  // #ifdef __numeric_vector_h__
