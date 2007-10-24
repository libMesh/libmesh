// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson
  
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


#include "libmesh_common.h"


#ifdef HAVE_PETSC



// TODO:[BSK} This seems necessary to use petsc on IBM Power3 at NERSC, but only there?  This will need to be wrapped in an ifdef with a variable set by configure
//#include <cmath>


// C++ includes
#include <vector>

// Local includes
#include "numeric_vector.h"
#include "petsc_macro.h"

/**
 * Petsc include files.
 */
#ifndef USE_COMPLEX_NUMBERS
extern "C" {
# include <petscvec.h>
}
#else
# include <petscvec.h>
#endif



// forward declarations
template <typename T> class SparseMatrix;

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
  PetscVector (const unsigned int n,
	       const unsigned int n_local);

  /**
   * Constructor.  Creates a PetscVector assuming you already have a
   * valid PETSc Vec object.  In this case, v is NOT destroyed by the
   * PetscVector constructor when this object goes out of scope.
   * This allows ownership of v to remain with the original creator,
   * and to simply provide additional functionality with the PetscVector.
   */
  PetscVector(Vec v);
  
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
   * Creates a copy of this vector and returns it in an \p AutoPtr.
   */
  AutoPtr<NumericVector<T> > clone () const;
  
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
   * \f$U(0-N) = s\f$: fill all components.
   */
  NumericVector<T> & operator= (const T s);
    
  /**
   *  \f$U = V\f$: copy all components.
   */
  NumericVector<T> & operator= (const NumericVector<T> &V);

  /**
   *  \f$U = V\f$: copy all components.
   */
  PetscVector<T> & operator= (const PetscVector<T> &V);

  /**
   *  \f$U = V\f$: copy all components.
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
   * @returns the sum of values in a vector
   */
  T sum () const;

  /**
   * @returns the \f$l_1\f$-norm of the vector, i.e.
   * the sum of the absolute values.
   */
  Real l1_norm () const;

  /**
   * @returns the \f$l_2\f$-norm of the vector, i.e.
   * the square root of the sum of the
   * squares of the elements.
   */
  Real l2_norm () const;

  /**
   * @returns the maximum absolute value of the
   * elements of this vector, which is the
   * \f$l_\infty\f$-norm of a vector.
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
   * \f$U(0-DIM)+=s\f$.
   * Addition of \p s to all components. Note
   * that \p s is a scalar and not a vector.
   */
  void add (const T s);
    
  /**
   * \f$ U+=V \f$ .
   * Simple vector addition, equal to the
   * \p operator +=.
   */
  void add (const NumericVector<T>& V);

  /**
   * \f$ U+=a*V \f$ .
   * Simple vector addition, equal to the
   * \p operator +=.
   */
  void add (const T a, const NumericVector<T>& v);
  
  /**
   * \f$ U+=v \f$ where \p v is a std::vector<T> 
   * and you
   * want to specify WHERE to add it
   */
  void add_vector (const std::vector<T>& v,
		   const std::vector<unsigned int>& dof_indices);

  /**
   * \f$ U+=V \f$ where U and V are type 
   * \p NumericVector<T> and you
   * want to specify WHERE to add
   * the \p NumericVector<T> V 
   */
  void add_vector (const NumericVector<T>& V,
		   const std::vector<unsigned int>& dof_indices);


  /**
   * \f$U+=A*V\f$, add the product of a \p SparseMatrix \p A
   * and a \p NumericVector \p V to \p this, where \p this=U.
   */
  void add_vector (const NumericVector<T> &V,
		   const SparseMatrix<T> &A);
     
  /**
   * \f$U+=V \f$ where U and V are type 
   * DenseVector<T> and you
   * want to specify WHERE to add
   * the DenseVector<T> V 
   */
  void add_vector (const DenseVector<T>& V,
		   const std::vector<unsigned int>& dof_indices);
  
  /**
   * \f$ U=v \f$ where v is a DenseVector<T> 
   * and you want to specify WHERE to insert it
   */
  virtual void insert (const std::vector<T>& v,
		       const std::vector<unsigned int>& dof_indices);

  /**
   * \f$U=V\f$, where U and V are type 
   * NumericVector<T> and you
   * want to specify WHERE to insert
   * the NumericVector<T> V 
   */
  virtual void insert (const NumericVector<T>& V,
		       const std::vector<unsigned int>& dof_indices);
      
  /**
   * \f$ U+=V \f$ where U and V are type 
   * DenseVector<T> and you
   * want to specify WHERE to insert
   * the DenseVector<T> V 
   */
  virtual void insert (const DenseVector<T>& V,
		       const std::vector<unsigned int>& dof_indices);
    
  
  /**
   * Scale each element of the
   * vector by the given factor.
   */
  void scale (const T factor);

  /**
   * Computes the dot product, p = U.V
   */
  virtual Number dot(const NumericVector<T>& V) const;
  
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
   * Updates a local vector with selected values from neighboring
   * processors, as defined by \p send_list.
   */
  void localize (const unsigned int first_local_idx,
		 const unsigned int last_local_idx,
		 const std::vector<unsigned int>& send_list);
  
  /**
   * Creates a local copy of the global vector in
   * \p v_local only on processor \p proc_id.  By
   * default the data is sent to processor 0.  This method
   * is useful for outputting data from one processor.
   */
  void localize_to_one (std::vector<T>& v_local,
			const unsigned int proc_id=0) const;
  
  /**
   * Print the contents of the vector in Matlab
   * format. Optionally prints the
   * matrix to the file named \p name.  If \p name
   * is not specified it is dumped to the screen.
   */
  void print_matlab(const std::string name="NULL") const;

  /**
   * Creates a "subvector" from this vector using the rows indices
   * of the "rows" array.
   */
  virtual void create_subvector(NumericVector<T>& subvector,
				const std::vector<unsigned int>& rows) const;

  /**
   * Returns the raw PETSc vector context pointer.  Note this is generally
   * not required in user-level code. Just don't do anything crazy like
   * calling VecDestroy()!
   */
  Vec vec () { assert (_vec != NULL); return _vec; }


  
private:

  /**
   * Actual Petsc vector datatype
   * to hold vector entries
   */
  Vec _vec;

  /**
   * This boolean value should only be set to false
   * for the constructor which takes a PETSc Vec object. 
   */
  const bool _destroy_vec_on_exit;
};


/*----------------------- Inline functions ----------------------------------*/



template <typename T>
inline
PetscVector<T>::PetscVector ()
  : _destroy_vec_on_exit(true)
{}



template <typename T>
inline
PetscVector<T>::PetscVector (const unsigned int n)
  : _destroy_vec_on_exit(true)
{
  this->init(n, n, false);
}



template <typename T>
inline
PetscVector<T>::PetscVector (const unsigned int n,
			     const unsigned int n_local)
  : _destroy_vec_on_exit(true)
{
  this->init(n, n_local, false);
}





template <typename T>
inline
PetscVector<T>::PetscVector (Vec v)
  : _destroy_vec_on_exit(false)
{
  this->_vec = v;
  this->_is_initialized = true;
}




template <typename T>
inline
PetscVector<T>::~PetscVector ()
{
  this->clear ();
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


  // Clear initialized vectors 
  if (this->initialized())
    this->clear();

  
  // create a sequential vector if on only 1 processor 
  if (n_local == n)
    {
      ierr = VecCreateSeq (PETSC_COMM_SELF, petsc_n, &_vec);
             CHKERRABORT(PETSC_COMM_SELF,ierr);
      
      ierr = VecSetFromOptions (_vec);
             CHKERRABORT(PETSC_COMM_SELF,ierr);
    }
  // otherwise create an MPI-enabled vector
  else
    {
      assert (n_local < n);
      
      ierr = VecCreateMPI (libMesh::COMM_WORLD, petsc_n_local, petsc_n,
			   &_vec);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);
      
      ierr = VecSetFromOptions (_vec);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);
    }  
  
  this->_is_initialized = true;
  
  
  if (fast == false)
    this->zero ();
}



template <typename T>
inline
void PetscVector<T>::init (const unsigned int n,
			   const bool fast)
{
  this->init(n,n,fast);
}



template <typename T>
inline
void PetscVector<T>::close ()
{
  assert (this->initialized());
  
  int ierr=0;
  
  ierr = VecAssemblyBegin(_vec);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);
  ierr = VecAssemblyEnd(_vec);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  this->_is_closed = true;
}



template <typename T>
inline
void PetscVector<T>::clear ()
{
  if ((this->initialized()) && (this->_destroy_vec_on_exit))
    {
      int ierr=0;

      ierr = VecDestroy(_vec);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);
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

#if PETSC_VERSION_LESS_THAN(2,3,0)  
  
  // 2.2.x & earlier style
  ierr = VecSet (&z, _vec);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

#else
  
  // 2.3.x & newer
  ierr = VecSet (_vec, z);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

#endif
}



template <typename T>
inline
AutoPtr<NumericVector<T> > PetscVector<T>::clone () const
{
  AutoPtr<NumericVector<T> > cloned_vector (new PetscVector<T>);

  *cloned_vector = *this;

  return cloned_vector;
}



template <typename T>
inline
unsigned int PetscVector<T>::size () const
{
  assert (this->initialized());
  
  int ierr=0, petsc_size=0;
  
  if (!this->initialized())
    return 0;
  
  ierr = VecGetSize(_vec, &petsc_size);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  return static_cast<unsigned int>(petsc_size);
}



template <typename T>
inline
unsigned int PetscVector<T>::local_size () const
{
  assert (this->initialized());
  
  int ierr=0, petsc_size=0;
  
  ierr = VecGetLocalSize(_vec, &petsc_size);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);
  
  return static_cast<unsigned int>(petsc_size);
}



template <typename T>
inline
unsigned int PetscVector<T>::first_local_index () const
{
  assert (this->initialized());
  
  int ierr=0, petsc_first=0, petsc_last=0;
  
  ierr = VecGetOwnershipRange (_vec, &petsc_first, &petsc_last);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);
  
  return static_cast<unsigned int>(petsc_first);
}



template <typename T>
inline
unsigned int PetscVector<T>::last_local_index () const
{
  assert (this->initialized());
  
  int ierr=0, petsc_first=0, petsc_last=0;
  
  ierr = VecGetOwnershipRange (_vec, &petsc_first, &petsc_last);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);
  
  return static_cast<unsigned int>(petsc_last);
}



template <typename T>
inline
T PetscVector<T>::operator() (const unsigned int i) const
{
  assert (this->initialized());
  assert ( ((i >= this->first_local_index()) &&
	    (i <  this->last_local_index())) );

  int ierr=0;
  PetscScalar *values, value=0.;
  

  ierr = VecGetArray(_vec, &values);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);
  
  value = values[i - this->first_local_index()];
  
  ierr = VecRestoreArray (_vec, &values);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);
  
  return static_cast<T>(value);
}



template <typename T>
inline
Real PetscVector<T>::min () const
{
  assert (this->initialized());

  int index=0, ierr=0;
  PetscReal min=0.;

  ierr = VecMin (_vec, &index, &min);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

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

  ierr = VecMax (_vec, &index, &max);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  // this return value is correct: VecMax returns a PetscReal
  return static_cast<Real>(max);
}


#endif // #ifdef HAVE_PETSC
#endif // #ifdef __petsc_vector_h__








