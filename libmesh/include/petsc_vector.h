//    $Id: petsc_vector.h,v 1.8 2003-02-09 22:47:17 ddreyer Exp $

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



#include "mesh_common.h"


#ifdef HAVE_PETSC

#ifndef __petsc_vector_h__
#define __petsc_vector_h__


// TODO:[BSK} This seems necessary to use petsc on IBM Power3 at NERSC, but only there?  This will need to be wrapped in an ifdef with a variable set by configure
#include <cmath>


// C++ includes
#include <vector>



/*
 * Petsc include files.  PETSc with complex numbers 
 * is actually C++.
 */
# ifndef USE_COMPLEX_NUMBERS

namespace Petsc {
extern "C" {
#include "petscvec.h"
}
// for easy switching between Petsc 2.1.0/2.1.1
// typedef Scalar PetscScalar;
} 
using namespace Petsc;

#else

#include "petscvec.h"

#endif



// forward declarations
class PetscMatrix;
class PetscInterface;


/**
 * Petsc vector. Provides a nice interface to the
 * Petsc C-based data structures for parallel vectors.
 *
 * @author Benjamin S. Kirk, 2002
 */

class PetscVector
{
 public:

  /**
   *  Dummy-Constructor. Dimension=0
   */
  PetscVector ();
    
  /**
   * Copy-Constructor. Dimension set to that of V, all components are copied
   * from V
   */
  //PetscVector (const PetscVector& V);
    
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
   * @returns the \p PetscVector to a pristine state.
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
    
  /**
   * Change the dimension to that of the
   * vector \p V. The same applies as for
   * the other \p init function.
   *
   * The elements of \p V are not copied, i.e.
   * this function is the same as calling
   * \p init(V.size(),fast).
   */
  void init (const PetscVector& V,
	       const bool fast=false);

  /**
   * $U(0-N) = s$: fill all components.
   */
  PetscVector & operator= (const Complex s);
    
  /**
   *  $U = V$: copy all components.
   */
  PetscVector & operator= (const PetscVector &V);

  /**
   *  $U = V$: copy all components.
   */
  PetscVector & operator= (const std::vector<Complex> &v);

  /**
   * @returns the scalar product of
   * two vectors.  The return type
   * is the underlying type of
   * \p this vector, so the return
   * type and the accuracy with
   * which it the result is
   * computed depend on the order
   * of the arguments of this
   * vector.
   */
  Complex operator* (const PetscVector &V) const;

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
   * was renamed to get the \p PetscVector class
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
  PetscVector & operator += (const PetscVector &V);

  /**
   * Subtraction operator.
   * Fast equivalent to \p U.add(-1, V).
   */
  PetscVector & operator -= (const PetscVector &V);
    
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
  void add (const PetscVector& V);

  /**
   * U+=a*V.
   * Simple vector addition, equal to the
   * \p operator +=.
   */
  void add (const Complex a, const PetscVector& v);
  
  /**
   * U+=v where v is a std::vector<Complex> 
   * and you
   * want to specify WHERE to add it
   */
  void add_vector (const std::vector<Complex>& v,
		   const std::vector<unsigned int>& dof_indices);

  /**
   * U+=V where U and V are type 
   * PetscVector and you
   * want to specify WHERE to add
   * the PetscVector V 
   */
  void add_petsc_vector (const PetscVector& V,
			 const std::vector<unsigned int>& dof_indices);
    

  
  /**
   * Scale each element of the
   * vector by the given factor.
   */
  void scale (const Complex factor);

  /**
   * Scale each element of the
   * vector by a constant
   * value. This operator is an
   * alias to the @ref{scale}
   * function, except that it
   * returns a reference to itself.
   */
  PetscVector & operator *= (const Complex factor);
    
  /**
   * Creates a copy of the global vector in the
   * local vector \p v_local.
   */
  void localize (std::vector<Complex>& v_local) const;

  /**
   * Same, but fills a \p PetscVector instead of
   * a \p std::vector.
   */
  void localize (PetscVector& v_local) const;

  /**
   * Creates a local vector \p v_local containing
   * only information relevant to this processor, as
   * defined by the \p send_list.
   */
  void localize (PetscVector& v_local,
		 const std::vector<unsigned int>& send_list) const;

  /**
   * Creates a local copy of the global vector in
   * \p v_local only on processor \p proc_id.  By
   * default the data is sent to processor 0.  This method
   * is useful for outputting data from one processor.
   */
  void localize_to_one (std::vector<Complex>& v_local,
			const unsigned int proc_id=0) const;
    
  /**
   * Prints the contents of the vector to the screen.
   */
  void print() const;
    
 private:

  /**
   * Actual Petsc vector datatype
   * to hold vector entries
   */
  Vec vec;
  
  /**
   * Flag to see if the Petsc
   * assemble routines have been called yet
   */
  bool is_closed;
  
  /**
   * Flag to tell if init 
   * has been called yet
   */
  bool initialized;

  friend class PetscInterface;
};


/*----------------------- Inline functions ----------------------------------*/



inline
PetscVector::PetscVector () :
  is_closed(false),
  initialized(false)
{};



inline
PetscVector::PetscVector (const unsigned int n) :
  is_closed(false),
  initialized(false)
{
  init(n, n, false);
};



inline
PetscVector::PetscVector (const unsigned int n,
			  const unsigned int n_local) :
  is_closed(false),
  initialized(false)
{
  init(n, n_local, false);
};



inline
PetscVector::~PetscVector ()
{
  clear ();
};



inline
void PetscVector::init (const unsigned int n,
			const unsigned int n_local,
			const bool fast)
{
  int ierr=0;
  int petsc_n=static_cast<int>(n);
  int petsc_n_local=static_cast<int>(n_local);


  
  if (initialized)
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
  
  initialized = true;

  
  if (fast == false)
    zero ();

  
  
  return;
};



inline
void PetscVector::init (const unsigned int n,
			const bool fast)
{
  init(n,n,fast);
};



inline
void PetscVector::close ()
{
  assert (initialized);
  
  int ierr=0;
  
  ierr = VecAssemblyBegin(vec); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(vec);   CHKERRQ(ierr);

  is_closed = true;
  
  return;
};



inline
void PetscVector::clear ()
{
  if (initialized)
    {
      int ierr=0;

      ierr = VecDestroy(vec); CHKERRQ(ierr);
    };

  is_closed = initialized = false;
};



inline
void PetscVector::zero ()
{
  assert (initialized);
  
  int ierr=0;

  PetscScalar z=0.;
  
  ierr = VecSet (&z, vec); CHKERRQ(ierr);
};



inline
unsigned int PetscVector::size () const
{
  assert (initialized);
  
  int ierr=0, petsc_size=0;
  
  if (!initialized)
    return 0;
  
  ierr = VecGetSize(vec, &petsc_size); CHKERRQ(ierr);

  return static_cast<unsigned int>(petsc_size);
};



inline
unsigned int PetscVector::local_size () const
{
  assert (initialized);
  
  int ierr=0, petsc_size=0;
  
  ierr = VecGetLocalSize(vec, &petsc_size); CHKERRQ(ierr);
  
  return static_cast<unsigned int>(petsc_size);
};



inline
unsigned int PetscVector::first_local_index () const
{
  assert (initialized);
  
  int ierr=0, petsc_first=0, petsc_last=0;
  
  ierr = VecGetOwnershipRange (vec, &petsc_first, &petsc_last); CHKERRQ(ierr);
  
  return static_cast<unsigned int>(petsc_first);
};



inline
unsigned int PetscVector::last_local_index () const
{
  assert (initialized);
  
  int ierr=0, petsc_first=0, petsc_last=0;
  
  ierr = VecGetOwnershipRange (vec, &petsc_first, &petsc_last); CHKERRQ(ierr);
  
  return static_cast<unsigned int>(petsc_last);
};



inline
Complex PetscVector::operator() (const unsigned int i) const
{
  assert (initialized);
  assert ( ((i >= first_local_index()) &&
	    (i <  last_local_index())) );

  int ierr=0;
  PetscScalar *values, value=0.;
  

  ierr = VecGetArray(vec, &values); CHKERRQ(ierr);
  
  value = values[i];
  
  ierr = VecRestoreArray (vec, &values); CHKERRQ(ierr);
  
  return static_cast<Complex>(value);
};



inline
Real PetscVector::min () const
{
  assert (initialized);

  int index=0, ierr=0;
  PetscReal min=0.;

  ierr = VecMin (vec, &index, &min); CHKERRQ(ierr);

  // this return value is correct: VecMin returns a PetscReal
  return static_cast<Real>(min);
};



inline
Real PetscVector::max() const
{
  assert (initialized);

  int index=0, ierr=0;
  PetscReal max=0.;

  ierr = VecMax (vec, &index, &max); CHKERRQ(ierr);

  // this return value is correct: VecMax returns a PetscReal
  return static_cast<Real>(max);
};



inline
void PetscVector::print() const
{
  assert (initialized);
  
#ifndef USE_COMPLEX_NUMBERS

  for (unsigned int i=0; i<size(); i++)
    std::cout << (*this)(i) << std::endl;

#else
  // std::complex<>::operator<<() is defined, but use this form

  std::cout << "Real part:" << std::endl;
  for (unsigned int i=0; i<size(); i++)
    std::cout << (*this)(i).real() << std::endl;

  std::cout << std::endl << "Imaginary part:" << std::endl;
  for (unsigned int i=0; i<size(); i++)
    std::cout << (*this)(i).imag() << std::endl;


#endif
  
  return;
};


#endif
#endif
