// $Id: petsc_vector.h,v 1.19 2003-03-26 13:55:24 benkirk Exp $

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




#ifndef __petsc_vector_h__
#define __petsc_vector_h__


#include "mesh_common.h"


#ifdef HAVE_PETSC



// TODO:[BSK} This seems necessary to use petsc on IBM Power3 at NERSC, but only there?  This will need to be wrapped in an ifdef with a variable set by configure
//#include <cmath>


// C++ includes
#include <vector>

// Local includes
#include "numeric_vector.h"

/*
 * Petsc include files.  PETSc with complex numbers 
 * is actually C++.
 */
# ifndef USE_COMPLEX_NUMBERS

extern "C" {
#include <petscvec.h>
}
// for easy switching between Petsc 2.1.0/2.1.1
// typedef Scalar PetscScalar;

#else

#include <petscvec.h>

#endif



// forward declarations
template <typename T> class PetscInterface;
template <typename T> class SparseMatrix;
template <typename T> class PetscMatrix;

/**
 * Petsc vector. Provides a nice interface to the
 * Petsc C-based data structures for parallel vectors.
 *
 * @author Benjamin S. Kirk, 2002
 */

template <typename T>
class PetscVector : public NumericVector<T>
{
public:

  /**
   *  Dummy-Constructor. Dimension=0
   */
  PetscVector ();
  
  /**
   * Constructor. Set dimension to \p n and initialize all elements with zero.
   */
  PetscVector (const unsigned int n);
    
  /**
   * Constructor. Set local dimension to \p n_local, the global dimension
   * to \p n, and initialize all elements with zero.
   */
  PetscVector (const unsigned n,
	       const unsigned int n_local);
    
  /**
   * Destructor, deallocates memory. Made virtual to allow
   * for derived classes to behave properly.
   */
  ~PetscVector ();

  /**
   * Call the assemble functions
   */
  void close (); 

  /**
   * @returns the \p PetscVector<T> to a pristine state.
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
  PetscVector<T> & operator= (const PetscVector<T> &V);

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
   * was renamed to get the \p PetscVector<T> class
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
   * U+=U+A*V.
   * Add the product of a Sparse matrix \p A
   * and a Numeric vector \p V to this Numeric vector.
   */
  void add_vector (const NumericVector<T> &V,
		   SparseMatrix<T> &A);
     
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
   * Actual Petsc vector datatype
   * to hold vector entries
   */
  Vec vec;
  
  /**
   * Make other Petsc datatypes friends
   */
  friend class PetscInterface<T>;
  friend class PetscMatrix<T>;
};


/*----------------------- Inline functions ----------------------------------*/



template <typename T>
inline
PetscVector<T>::PetscVector ()
{}



template <typename T>
inline
PetscVector<T>::PetscVector (const unsigned int n)
{
  init(n, n, false);
}



template <typename T>
inline
PetscVector<T>::PetscVector (const unsigned int n,
			     const unsigned int n_local)
{
  init(n, n_local, false);
}



template <typename T>
inline
PetscVector<T>::~PetscVector ()
{
  clear ();
}



template <typename T>
inline
void PetscVector<T>::init (const unsigned int n,
			   const unsigned int n_local,
			   const bool fast)
{
  int ierr=0;
  int petsc_n=static_cast<int>(n);
  int petsc_n_local=static_cast<int>(n_local);


  
  if (this->initialized())
    {
      std::cerr << "ERROR: Vector already initialized!"
		<< std::endl;
      
      error();
    }

  
  // create a sequential vector if on only 1 processor 
  if (n_local == n)
    {
      ierr = VecCreateSeq (PETSC_COMM_SELF, petsc_n, &vec);      CHKERRQ(ierr);
      ierr = VecSetFromOptions (vec);                            CHKERRQ(ierr);
    }
  // otherwise create an MPI-enabled vector
  else
    {
      assert (n_local < n);
      
      ierr = VecCreateMPI (PETSC_COMM_WORLD, petsc_n_local, petsc_n,
			   &vec);                       CHKERRQ(ierr);
      
      ierr = VecSetFromOptions (vec);                   CHKERRQ(ierr);
    }  
  
  this->_is_initialized = true;

  
  if (fast == false)
    zero ();

  
  
  return;
}



template <typename T>
inline
void PetscVector<T>::init (const unsigned int n,
			   const bool fast)
{
  init(n,n,fast);
}



template <typename T>
inline
void PetscVector<T>::close ()
{
  assert (this->initialized());
  
  int ierr=0;
  
  ierr = VecAssemblyBegin(vec); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(vec);   CHKERRQ(ierr);

  this->_is_closed = true;
  
  return;
}



template <typename T>
inline
void PetscVector<T>::clear ()
{
  if (this->initialized())
    {
      int ierr=0;

      ierr = VecDestroy(vec); CHKERRQ(ierr);
    }

  this->_is_closed = this->_is_initialized = false;
}



template <typename T>
inline
void PetscVector<T>::zero ()
{
  assert (this->initialized());
  
  int ierr=0;

  PetscScalar z=0.;
  
  ierr = VecSet (&z, vec); CHKERRQ(ierr);
}



template <typename T>
inline
unsigned int PetscVector<T>::size () const
{
  assert (this->initialized());
  
  int ierr=0, petsc_size=0;
  
  if (!this->initialized())
    return 0;
  
  ierr = VecGetSize(vec, &petsc_size); CHKERRQ(ierr);

  return static_cast<unsigned int>(petsc_size);
}



template <typename T>
inline
unsigned int PetscVector<T>::local_size () const
{
  assert (this->initialized());
  
  int ierr=0, petsc_size=0;
  
  ierr = VecGetLocalSize(vec, &petsc_size); CHKERRQ(ierr);
  
  return static_cast<unsigned int>(petsc_size);
}



template <typename T>
inline
unsigned int PetscVector<T>::first_local_index () const
{
  assert (this->initialized());
  
  int ierr=0, petsc_first=0, petsc_last=0;
  
  ierr = VecGetOwnershipRange (vec, &petsc_first, &petsc_last); CHKERRQ(ierr);
  
  return static_cast<unsigned int>(petsc_first);
}



template <typename T>
inline
unsigned int PetscVector<T>::last_local_index () const
{
  assert (this->initialized());
  
  int ierr=0, petsc_first=0, petsc_last=0;
  
  ierr = VecGetOwnershipRange (vec, &petsc_first, &petsc_last); CHKERRQ(ierr);
  
  return static_cast<unsigned int>(petsc_last);
}



template <typename T>
inline
T PetscVector<T>::operator() (const unsigned int i) const
{
  assert (this->initialized());
  assert ( ((i >= first_local_index()) &&
	    (i <  last_local_index())) );

  int ierr=0;
  PetscScalar *values, value=0.;
  

  ierr = VecGetArray(vec, &values); CHKERRQ(ierr);
  
  value = values[i];
  
  ierr = VecRestoreArray (vec, &values); CHKERRQ(ierr);
  
  return static_cast<T>(value);
}



template <typename T>
inline
Real PetscVector<T>::min () const
{
  assert (this->initialized());

  int index=0, ierr=0;
  PetscReal min=0.;

  ierr = VecMin (vec, &index, &min); CHKERRQ(ierr);

  // this return value is correct: VecMin returns a PetscReal
  return static_cast<Real>(min);
}



template <typename T>
inline
Real PetscVector<T>::max() const
{
  assert (this->initialized());

  int index=0, ierr=0;
  PetscReal max=0.;

  ierr = VecMax (vec, &index, &max); CHKERRQ(ierr);

  // this return value is correct: VecMax returns a PetscReal
  return static_cast<Real>(max);
}


#endif // #ifdef HAVE_PETSC
#endif // #ifdef __petsc_vector_h__








