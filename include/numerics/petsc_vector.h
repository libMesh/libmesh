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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA




#ifndef LIBMESH_PETSC_VECTOR_H
#define LIBMESH_PETSC_VECTOR_H


#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_HAVE_PETSC

// Local includes
#include "libmesh/numeric_vector.h"
#include "libmesh/petsc_macro.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/petsc_solver_exception.h"
#include LIBMESH_INCLUDE_UNORDERED_MAP

// Petsc include files.
#include <petscvec.h>

// C++ includes
#include <cstddef>
#include <cstring>
#include <vector>

#ifdef LIBMESH_HAVE_CXX11_THREAD
#include <atomic>
#include <mutex>
#endif

namespace libMesh
{

// forward declarations
template <typename T> class SparseMatrix;

/**
 * Petsc vector. Provides a nice interface to the
 * Petsc C-based data structures for parallel vectors.
 *
 * \author Benjamin S. Kirk
 * \date 2002
 */
template <typename T>
class PetscVector libmesh_final : public NumericVector<T>
{
public:

  /**
   *  Dummy-Constructor. Dimension=0
   */
  explicit
  PetscVector (const Parallel::Communicator & comm_in,
               const ParallelType type = AUTOMATIC);

  /**
   * Constructor. Set dimension to \p n and initialize all elements with zero.
   */
  explicit
  PetscVector (const Parallel::Communicator & comm_in,
               const numeric_index_type n,
               const ParallelType type = AUTOMATIC);

  /**
   * Constructor. Set local dimension to \p n_local, the global dimension
   * to \p n, and initialize all elements with zero.
   */
  PetscVector (const Parallel::Communicator & comm_in,
               const numeric_index_type n,
               const numeric_index_type n_local,
               const ParallelType type = AUTOMATIC);

  /**
   * Constructor. Set local dimension to \p n_local, the global
   * dimension to \p n, but additionally reserve memory for the
   * indices specified by the \p ghost argument.
   */
  PetscVector (const Parallel::Communicator & comm_in,
               const numeric_index_type N,
               const numeric_index_type n_local,
               const std::vector<numeric_index_type> & ghost,
               const ParallelType type = AUTOMATIC);

  /**
   * Constructor.  Creates a PetscVector assuming you already have a
   * valid PETSc Vec object.  In this case, \p v is NOT destroyed by the
   * PetscVector constructor when this object goes out of scope.
   * This allows ownership of \p v to remain with the original creator,
   * and to simply provide additional functionality with the PetscVector.
   */
  PetscVector(Vec v,
              const Parallel::Communicator & comm_in
              LIBMESH_CAN_DEFAULT_TO_COMMWORLD);

  /**
   * Destructor, deallocates memory. Made virtual to allow
   * for derived classes to behave properly.
   */
  ~PetscVector ();

  /**
   * Call the assemble functions
   */
  virtual void close () libmesh_override;

  /**
   * @returns the \p PetscVector<T> to a pristine state.
   */
  virtual void clear () libmesh_override;

  /**
   * Set all entries to zero. Equivalent to \p v = 0, but more obvious and
   * faster.
   */
  virtual void zero () libmesh_override;

  /**
   * Creates a vector which has the same type, size and partitioning
   * as this vector, but whose data is all zero.  Returns it in an \p
   * UniquePtr.
   */
  virtual UniquePtr<NumericVector<T> > zero_clone () const libmesh_override;

  /**
   * Creates a copy of this vector and returns it in an \p UniquePtr.
   */
  virtual UniquePtr<NumericVector<T> > clone () const libmesh_override;

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
  virtual void init (const numeric_index_type N,
                     const numeric_index_type n_local,
                     const bool         fast=false,
                     const ParallelType type=AUTOMATIC) libmesh_override;

  /**
   * call init with n_local = N,
   */
  virtual void init (const numeric_index_type N,
                     const bool         fast=false,
                     const ParallelType type=AUTOMATIC) libmesh_override;

  /**
   * Create a vector that holds tha local indices plus those specified
   * in the \p ghost argument.
   */
  virtual void init (const numeric_index_type /*N*/,
                     const numeric_index_type /*n_local*/,
                     const std::vector<numeric_index_type> & /*ghost*/,
                     const bool /*fast*/ = false,
                     const ParallelType = AUTOMATIC) libmesh_override;

  /**
   * Creates a vector that has the same dimension and storage type as
   * \p other, including ghost dofs.
   */
  virtual void init (const NumericVector<T> & other,
                     const bool fast = false) libmesh_override;

  /**
   * \f$U(0-N) = s\f$: fill all components.
   */
  virtual NumericVector<T> & operator= (const T s) libmesh_override;

  /**
   *  \f$U = V\f$: copy all components.
   */
  virtual NumericVector<T> & operator= (const NumericVector<T> & V) libmesh_override;

  /**
   *  \f$U = V\f$: copy all components.
   */
  PetscVector<T> & operator= (const PetscVector<T> & V);

  /**
   *  \f$U = V\f$: copy all components.
   */
  virtual NumericVector<T> & operator= (const std::vector<T> & v) libmesh_override;

  /**
   * @returns the minimum element in the vector.
   * In case of complex numbers, this returns the minimum
   * Real part.
   */
  virtual Real min () const libmesh_override;

  /**
   * @returns the maximum element in the vector.
   * In case of complex numbers, this returns the maximum
   * Real part.
   */
  virtual Real max () const libmesh_override;

  /**
   * @returns the sum of values in a vector
   */
  virtual T sum () const libmesh_override;

  /**
   * @returns the \f$l_1\f$-norm of the vector, i.e.
   * the sum of the absolute values.
   */
  virtual Real l1_norm () const libmesh_override;

  /**
   * @returns the \f$l_2\f$-norm of the vector, i.e.
   * the square root of the sum of the
   * squares of the elements.
   */
  virtual Real l2_norm () const libmesh_override;

  /**
   * @returns the maximum absolute value of the
   * elements of this vector, which is the
   * \f$l_\infty\f$-norm of a vector.
   */
  virtual Real linfty_norm () const libmesh_override;

  /**
   * @returns dimension of the vector. This
   * function was formerly called \p n(), but
   * was renamed to get the \p PetscVector<T> class
   * closer to the C++ standard library's
   * \p std::vector container.
   */
  virtual numeric_index_type size () const libmesh_override;

  /**
   * @returns the local size of the vector
   * (index_stop-index_start)
   */
  virtual numeric_index_type local_size() const libmesh_override;

  /**
   * @returns the index of the first vector element
   * actually stored on this processor
   */
  virtual numeric_index_type first_local_index() const libmesh_override;

  /**
   * @returns the index of the last vector element
   * actually stored on this processor
   */
  virtual numeric_index_type last_local_index() const libmesh_override;

  /**
   * Maps the global index \p i to the corresponding global index. If
   * the index is not a ghost cell, this is done by subtraction the
   * number of the first local index.  If it is a ghost cell, it has
   * to be looked up in the map.
   */
  numeric_index_type map_global_to_local_index(const numeric_index_type i) const;

  /**
   * Access components, returns \p U(i).
   */
  virtual T operator() (const numeric_index_type i) const libmesh_override;

  /**
   * Access multiple components at once.  \p values will *not* be
   * reallocated; it should already have enough space.  Overloaded
   * method that should be faster (probably much faster) than calling
   * \p operator() individually for each index.
   */
  virtual void get(const std::vector<numeric_index_type> & index,
                   T * values) const libmesh_override;

  /**
   * Addition operator.
   * Fast equivalent to \p U.add(1, V).
   */
  virtual NumericVector<T> & operator += (const NumericVector<T> & V) libmesh_override;

  /**
   * Subtraction operator.
   * Fast equivalent to \p U.add(-1, V).
   */
  virtual NumericVector<T> & operator -= (const NumericVector<T> & V) libmesh_override;

  /**
   * Replace each entry v_i of this vector by its reciprocal, 1/v_i.
   */
  virtual void reciprocal() libmesh_override;

  /**
   * Replace each entry v_i = real(v_i) + imag(v_i)
   * of this vector by its complex conjugate, real(v_i) - imag(v_i)
   */
  virtual void conjugate() libmesh_override;

  /**
   * v(i) = value
   */
  virtual void set (const numeric_index_type i,
                    const T value) libmesh_override;

  /**
   * v(i) += value
   */
  virtual void add (const numeric_index_type i,
                    const T value) libmesh_override;

  /**
   * \f$U(0-LIBMESH_DIM)+=s\f$.
   * Addition of \p s to all components. Note
   * that \p s is a scalar and not a vector.
   */
  virtual void add (const T s) libmesh_override;

  /**
   * \f$ U+=V \f$ .
   * Simple vector addition, equal to the
   * \p operator +=.
   */
  virtual void add (const NumericVector<T> & V) libmesh_override;

  /**
   * \f$ U+=a*V \f$ .
   * Simple vector addition, equal to the
   * \p operator +=.
   */
  virtual void add (const T a, const NumericVector<T> & v) libmesh_override;

  /**
   * We override two NumericVector<T>::add_vector() methods but don't
   * want to hide the other defaults.
   */
  using NumericVector<T>::add_vector;

  /**
   * \f$ U+=v \f$ where v is a pointer and each \p dof_indices[i]
   * specifies where to add value \p v[i]
   */
  virtual void add_vector (const T * v,
                           const std::vector<numeric_index_type> & dof_indices) libmesh_override;

  /**
   * \f$U+=A*V\f$, add the product of a \p SparseMatrix \p A
   * and a \p NumericVector \p V to \p this, where \p this=U.
   */
  virtual void add_vector (const NumericVector<T> & V,
                           const SparseMatrix<T> & A) libmesh_override;

  /**
   * \f$U+=A^T*V\f$, add the product of the transpose
   * of \p SparseMatrix \p A_trans and a \p NumericVector \p V to
   * \p this, where \p this=U.
   */
  virtual void add_vector_transpose (const NumericVector<T> & V,
                                     const SparseMatrix<T> & A_trans) libmesh_override;

  /**
   * \f$U+=A^H*V\f$, add the product of the conjugate-transpose
   * of \p SparseMatrix \p A_trans and a \p NumericVector \p V to
   * \p this, where \p this=U.
   */
  void add_vector_conjugate_transpose (const NumericVector<T> & V,
                                       const SparseMatrix<T> & A_trans);

  /**
   * We override one NumericVector<T>::insert() method but don't want
   * to hide the other defaults
   */
  using NumericVector<T>::insert;

  /**
   * \f$ U=v \f$ where v is a \p T[] or T *
   * and you want to specify WHERE to insert it
   */
  virtual void insert (const T * v,
                       const std::vector<numeric_index_type> & dof_indices) libmesh_override;

  /**
   * Scale each element of the
   * vector by the given factor.
   */
  virtual void scale (const T factor) libmesh_override;

  /**
   * Pointwise Division operator. ie divide every entry in this vector by the entry in v
   */
  virtual NumericVector<T> & operator /= (NumericVector<T> & v) libmesh_override;

  /**
   * v = abs(v)... that is, each entry in v is replaced
   * by its absolute value.
   */
  virtual void abs() libmesh_override;

  /**
   * Computes the dot product, p = U.V. Use complex-conjugate of V
   * in the complex-valued case.
   */
  virtual T dot(const NumericVector<T> &) const libmesh_override;

  /**
   * Computes the dot product, p = U.V. Do not use complex-conjugate of V
   * in the complex-valued case.
   */
  T indefinite_dot(const NumericVector<T> &) const;

  /**
   * Creates a copy of the global vector in the
   * local vector \p v_local.
   */
  virtual void localize (std::vector<T> & v_local) const libmesh_override;

  /**
   * Same, but fills a \p NumericVector<T> instead of
   * a \p std::vector.
   */
  virtual void localize (NumericVector<T> & v_local) const libmesh_override;

  /**
   * Creates a local vector \p v_local containing
   * only information relevant to this processor, as
   * defined by the \p send_list.
   */
  virtual void localize (NumericVector<T> & v_local,
                         const std::vector<numeric_index_type> & send_list) const libmesh_override;

  /**
   * Updates a local vector with selected values from neighboring
   * processors, as defined by \p send_list.
   */
  virtual void localize (const numeric_index_type first_local_idx,
                         const numeric_index_type last_local_idx,
                         const std::vector<numeric_index_type> & send_list) libmesh_override;

  /**
   * Creates a local copy of the global vector in
   * \p v_local only on processor \p proc_id.  By
   * default the data is sent to processor 0.  This method
   * is useful for outputting data from one processor.
   */
  virtual void localize_to_one (std::vector<T> & v_local,
                                const processor_id_type proc_id=0) const libmesh_override;

  /**
   * Computes the pointwise (i.e. component-wise) product of \p vec1
   * and \p vec2 and stores the result in \p *this.
   */
  virtual void pointwise_mult (const NumericVector<T> & vec1,
                               const NumericVector<T> & vec2) libmesh_override;

  /**
   * Print the contents of the vector in Matlab
   * format. Optionally prints the
   * matrix to the file named \p name.  If \p name
   * is not specified it is dumped to the screen.
   */
  virtual void print_matlab(const std::string & name = "") const libmesh_override;

  /**
   * Creates a "subvector" from this vector using the rows indices
   * of the "rows" array.
   */
  virtual void create_subvector(NumericVector<T> & subvector,
                                const std::vector<numeric_index_type> & rows) const libmesh_override;

  /**
   * Swaps the raw PETSc vector context pointers.
   */
  virtual void swap (NumericVector<T> & v) libmesh_override;

  /**
   * Returns the raw PETSc vector context pointer.  Note this is generally
   * not required in user-level code. Just don't do anything crazy like
   * calling LibMeshVecDestroy()!
   */
  Vec vec () { libmesh_assert (_vec); return _vec; }


private:

  /**
   * Actual Petsc vector datatype
   * to hold vector entries
   */
  Vec _vec;

  /**
   * If \p true, the actual Petsc array of the values of the vector is
   * currently accessible.  That means that the members \p _local_form
   * and \p _values are valid.
   */
#ifdef LIBMESH_HAVE_CXX11_THREAD
  // Note: we can't use std::atomic_flag here because we need load and store operations
  mutable std::atomic<bool> _array_is_present;
#else
  mutable bool _array_is_present;
#endif

  /**
   * First local index.
   *
   * Only valid when _array_is_present
   */
  mutable numeric_index_type _first;

  /**
   * Last local index.
   *
   * Only valid when _array_is_present
   */
  mutable numeric_index_type _last;

#ifndef NDEBUG
  /**
   * Size of the local form, for being used in assertations.  The
   * contents of this field are only valid if the vector is ghosted
   * and \p _array_is_present is \p true.
   */
  mutable numeric_index_type _local_size;
#endif

  /**
   * Petsc vector datatype to hold the local form of a ghosted vector.
   * The contents of this field are only valid if the vector is
   * ghosted and \p _array_is_present is \p true.
   */
  mutable Vec _local_form;

  /**
   * Pointer to the actual Petsc array of the values of the vector.
   * This pointer is only valid if \p _array_is_present is \p true.
   * We're using Petsc's VecGetArrayRead() function, which requires a
   * constant PetscScalar *, but _get_array and _restore_array are
   * const member functions, so _values also needs to be mutable
   * (otherwise it is a "const PetscScalar * const" in that context).
   */
  mutable const PetscScalar * _values;

  /**
   * Mutex for _get_array and _restore_array.  This is part of the
   * object to keep down thread contention when reading frmo multiple
   * PetscVectors simultaneously
   */
#ifdef LIBMESH_HAVE_CXX11_THREAD
  mutable std::mutex _petsc_vector_mutex;
#else
  mutable Threads::spin_mutex _petsc_vector_mutex;
#endif

  /**
   * Queries the array (and the local form if the vector is ghosted)
   * from Petsc.
   */
  void _get_array() const;

  /**
   * Restores the array (and the local form if the vector is ghosted)
   * to Petsc.
   */
  void _restore_array() const;

  /**
   * Type for map that maps global to local ghost cells.
   */
  typedef LIBMESH_BEST_UNORDERED_MAP<numeric_index_type,numeric_index_type> GlobalToLocalMap;

  /**
   * Map that maps global to local ghost cells (will be empty if not
   * in ghost cell mode).
   */
  GlobalToLocalMap _global_to_local_map;

  /**
   * This boolean value should only be set to false
   * for the constructor which takes a PETSc Vec object.
   */
  bool _destroy_vec_on_exit;
};


/*----------------------- Inline functions ----------------------------------*/



template <typename T>
inline
PetscVector<T>::PetscVector (const Parallel::Communicator & comm_in, const ParallelType ptype) :
  NumericVector<T>(comm_in, ptype),
  _array_is_present(false),
  _first(0),
  _last(0),
  _local_form(libmesh_nullptr),
  _values(libmesh_nullptr),
  _global_to_local_map(),
  _destroy_vec_on_exit(true)
{
  this->_type = ptype;
}



template <typename T>
inline
PetscVector<T>::PetscVector (const Parallel::Communicator & comm_in,
                             const numeric_index_type n,
                             const ParallelType ptype) :
  NumericVector<T>(comm_in, ptype),
  _array_is_present(false),
  _local_form(libmesh_nullptr),
  _values(libmesh_nullptr),
  _global_to_local_map(),
  _destroy_vec_on_exit(true)
{
  this->init(n, n, false, ptype);
}



template <typename T>
inline
PetscVector<T>::PetscVector (const Parallel::Communicator & comm_in,
                             const numeric_index_type n,
                             const numeric_index_type n_local,
                             const ParallelType ptype) :
  NumericVector<T>(comm_in, ptype),
  _array_is_present(false),
  _local_form(libmesh_nullptr),
  _values(libmesh_nullptr),
  _global_to_local_map(),
  _destroy_vec_on_exit(true)
{
  this->init(n, n_local, false, ptype);
}



template <typename T>
inline
PetscVector<T>::PetscVector (const Parallel::Communicator & comm_in,
                             const numeric_index_type n,
                             const numeric_index_type n_local,
                             const std::vector<numeric_index_type> & ghost,
                             const ParallelType ptype) :
  NumericVector<T>(comm_in, ptype),
  _array_is_present(false),
  _local_form(libmesh_nullptr),
  _values(libmesh_nullptr),
  _global_to_local_map(),
  _destroy_vec_on_exit(true)
{
  this->init(n, n_local, ghost, false, ptype);
}





template <typename T>
inline
PetscVector<T>::PetscVector (Vec v,
                             const Parallel::Communicator & comm_in) :
  NumericVector<T>(comm_in, AUTOMATIC),
  _array_is_present(false),
  _local_form(libmesh_nullptr),
  _values(libmesh_nullptr),
  _global_to_local_map(),
  _destroy_vec_on_exit(false)
{
  this->_vec = v;
  this->_is_closed = true;
  this->_is_initialized = true;

  /* We need to ask PETSc about the (local to global) ghost value
     mapping and create the inverse mapping out of it.  */
  PetscErrorCode ierr=0;
  PetscInt petsc_local_size=0;
  ierr = VecGetLocalSize(_vec, &petsc_local_size);
  LIBMESH_CHKERR(ierr);

  // Get the vector type from PETSc.
  // As of Petsc 3.0.0, the VecType #define lost its const-ness, so we
  // need to have it in the code
#if PETSC_VERSION_LESS_THAN(3,0,0) || !PETSC_VERSION_LESS_THAN(3,4,0)
  // Pre-3.0 and petsc-dev (as of October 2012) use non-const versions
  VecType ptype;
#else
  const VecType ptype;
#endif
  ierr = VecGetType(_vec, &ptype);
  LIBMESH_CHKERR(ierr);

  if((std::strcmp(ptype,VECSHARED) == 0) || (std::strcmp(ptype,VECMPI) == 0))
    {
#if PETSC_RELEASE_LESS_THAN(3,1,1)
      ISLocalToGlobalMapping mapping = _vec->mapping;
#else
      ISLocalToGlobalMapping mapping;
      ierr = VecGetLocalToGlobalMapping(_vec, &mapping);
      LIBMESH_CHKERR(ierr);
#endif

      // If is a sparsely stored vector, set up our new mapping
      if (mapping)
        {
          const numeric_index_type my_local_size = static_cast<numeric_index_type>(petsc_local_size);
          const numeric_index_type ghost_begin = static_cast<numeric_index_type>(petsc_local_size);
#if PETSC_RELEASE_LESS_THAN(3,4,0)
          const numeric_index_type ghost_end = static_cast<numeric_index_type>(mapping->n);
#else
          PetscInt n;
          ierr = ISLocalToGlobalMappingGetSize(mapping, &n);
          LIBMESH_CHKERR(ierr);
          const numeric_index_type ghost_end = static_cast<numeric_index_type>(n);
#endif
#if PETSC_RELEASE_LESS_THAN(3,1,1)
          const PetscInt * indices = mapping->indices;
#else
          const PetscInt * indices;
          ierr = ISLocalToGlobalMappingGetIndices(mapping,&indices);
          LIBMESH_CHKERR(ierr);
#endif
          for(numeric_index_type i=ghost_begin; i<ghost_end; i++)
            _global_to_local_map[indices[i]] = i-my_local_size;
          this->_type = GHOSTED;
#if !PETSC_RELEASE_LESS_THAN(3,1,1)
          ierr = ISLocalToGlobalMappingRestoreIndices(mapping, &indices);
          LIBMESH_CHKERR(ierr);
#endif
        }
      else
        this->_type = PARALLEL;
    }
  else
    this->_type = SERIAL;
}




template <typename T>
inline
PetscVector<T>::~PetscVector ()
{
  this->clear ();
}



template <typename T>
inline
void PetscVector<T>::init (const numeric_index_type n,
                           const numeric_index_type n_local,
                           const bool fast,
                           const ParallelType ptype)
{
  parallel_object_only();

  PetscErrorCode ierr=0;
  PetscInt petsc_n=static_cast<PetscInt>(n);
  PetscInt petsc_n_local=static_cast<PetscInt>(n_local);


  // Clear initialized vectors
  if (this->initialized())
    this->clear();

  if (ptype == AUTOMATIC)
    {
      if (n == n_local)
        this->_type = SERIAL;
      else
        this->_type = PARALLEL;
    }
  else
    this->_type = ptype;

  libmesh_assert ((this->_type==SERIAL && n==n_local) ||
                  this->_type==PARALLEL);

  // create a sequential vector if on only 1 processor
  if (this->_type == SERIAL)
    {
      ierr = VecCreateSeq (PETSC_COMM_SELF, petsc_n, &_vec);
      CHKERRABORT(PETSC_COMM_SELF,ierr);

      ierr = VecSetFromOptions (_vec);
      CHKERRABORT(PETSC_COMM_SELF,ierr);
    }
  // otherwise create an MPI-enabled vector
  else if (this->_type == PARALLEL)
    {
#ifdef LIBMESH_HAVE_MPI
      libmesh_assert_less_equal (n_local, n);
      ierr = VecCreateMPI (this->comm().get(), petsc_n_local, petsc_n,
                           &_vec);
      LIBMESH_CHKERR(ierr);
#else
      libmesh_assert_equal_to (n_local, n);
      ierr = VecCreateSeq (PETSC_COMM_SELF, petsc_n, &_vec);
      CHKERRABORT(PETSC_COMM_SELF,ierr);
#endif

      ierr = VecSetFromOptions (_vec);
      LIBMESH_CHKERR(ierr);
    }
  else
    libmesh_error_msg("Unsupported type " << this->_type);

  this->_is_initialized = true;
  this->_is_closed = true;


  if (fast == false)
    this->zero ();
}



template <typename T>
inline
void PetscVector<T>::init (const numeric_index_type n,
                           const bool fast,
                           const ParallelType ptype)
{
  this->init(n,n,fast,ptype);
}



template <typename T>
inline
void PetscVector<T>::init (const numeric_index_type n,
                           const numeric_index_type n_local,
                           const std::vector<numeric_index_type> & ghost,
                           const bool fast,
                           const ParallelType libmesh_dbg_var(ptype))
{
  parallel_object_only();

  PetscErrorCode ierr=0;
  PetscInt petsc_n=static_cast<PetscInt>(n);
  PetscInt petsc_n_local=static_cast<PetscInt>(n_local);
  PetscInt petsc_n_ghost=static_cast<PetscInt>(ghost.size());

  // If the mesh is not disjoint, every processor will either have
  // all the dofs, none of the dofs, or some non-zero dofs at the
  // boundary between processors.
  //
  // However we can't assert this, because someone might want to
  // construct a GHOSTED vector which doesn't include neighbor element
  // dofs.  Boyce tried to do so in user code, and we're going to want
  // to do so in System::project_vector().
  //
  // libmesh_assert(n_local == 0 || n_local == n || !ghost.empty());

  libmesh_assert(sizeof(PetscInt) == sizeof(numeric_index_type));

  PetscInt * petsc_ghost = ghost.empty() ? PETSC_NULL :
    const_cast<PetscInt *>(reinterpret_cast<const PetscInt *>(&ghost[0]));

  // Clear initialized vectors
  if (this->initialized())
    this->clear();

  libmesh_assert(ptype == AUTOMATIC || ptype == GHOSTED);
  this->_type = GHOSTED;

  /* Make the global-to-local ghost cell map.  */
  for(numeric_index_type i=0; i<ghost.size(); i++)
    {
      _global_to_local_map[ghost[i]] = i;
    }

  /* Create vector.  */
  ierr = VecCreateGhost (this->comm().get(), petsc_n_local, petsc_n,
                         petsc_n_ghost, petsc_ghost,
                         &_vec);
  LIBMESH_CHKERR(ierr);

  ierr = VecSetFromOptions (_vec);
  LIBMESH_CHKERR(ierr);

  this->_is_initialized = true;
  this->_is_closed = true;

  if (fast == false)
    this->zero ();
}



template <typename T>
inline
void PetscVector<T>::init (const NumericVector<T> & other,
                           const bool fast)
{
  parallel_object_only();

  // Clear initialized vectors
  if (this->initialized())
    this->clear();

  const PetscVector<T> & v = cast_ref<const PetscVector<T> &>(other);

  // Other vector should restore array.
  if(v.initialized())
    {
      v._restore_array();
    }

  this->_global_to_local_map = v._global_to_local_map;

  // Even if we're initializeing sizes based on an uninitialized or
  // unclosed vector, *this* vector is being initialized now and is
  // initially closed.
  this->_is_closed      = true; // v._is_closed;
  this->_is_initialized = true; // v._is_initialized;

  this->_type = v._type;

  // We want to have a valid Vec, even if it's initially of size zero
  PetscErrorCode ierr = VecDuplicate (v._vec, &this->_vec);
  LIBMESH_CHKERR(ierr);

  if (fast == false)
    this->zero ();
}



template <typename T>
inline
void PetscVector<T>::close ()
{
  parallel_object_only();

  this->_restore_array();

  PetscErrorCode ierr=0;

  ierr = VecAssemblyBegin(_vec);
  LIBMESH_CHKERR(ierr);
  ierr = VecAssemblyEnd(_vec);
  LIBMESH_CHKERR(ierr);

  if(this->type() == GHOSTED)
    {
      ierr = VecGhostUpdateBegin(_vec,INSERT_VALUES,SCATTER_FORWARD);
      LIBMESH_CHKERR(ierr);
      ierr = VecGhostUpdateEnd(_vec,INSERT_VALUES,SCATTER_FORWARD);
      LIBMESH_CHKERR(ierr);
    }

  this->_is_closed = true;
}



template <typename T>
inline
void PetscVector<T>::clear ()
{
  parallel_object_only();

  if (this->initialized())
    this->_restore_array();

  if ((this->initialized()) && (this->_destroy_vec_on_exit))
    {
      PetscErrorCode ierr=0;

      ierr = LibMeshVecDestroy(&_vec);
      LIBMESH_CHKERR(ierr);
    }

  this->_is_closed = this->_is_initialized = false;

  _global_to_local_map.clear();
}



template <typename T>
inline
void PetscVector<T>::zero ()
{
  parallel_object_only();

  libmesh_assert(this->closed());

  this->_restore_array();

  PetscErrorCode ierr=0;

  PetscScalar z=0.;

  if(this->type() != GHOSTED)
    {
      ierr = VecSet (_vec, z);
      LIBMESH_CHKERR(ierr);
    }
  else
    {
      /* Vectors that include ghost values require a special
         handling.  */
      Vec loc_vec;
      ierr = VecGhostGetLocalForm (_vec,&loc_vec);
      LIBMESH_CHKERR(ierr);

      ierr = VecSet (loc_vec, z);
      LIBMESH_CHKERR(ierr);

      ierr = VecGhostRestoreLocalForm (_vec,&loc_vec);
      LIBMESH_CHKERR(ierr);
    }
}



template <typename T>
inline
UniquePtr<NumericVector<T> > PetscVector<T>::zero_clone () const
{
  NumericVector<T> * cloned_vector = new PetscVector<T>(this->comm(), this->type());
  cloned_vector->init(*this);
  return UniquePtr<NumericVector<T> >(cloned_vector);
}



template <typename T>
inline
UniquePtr<NumericVector<T> > PetscVector<T>::clone () const
{
  NumericVector<T> * cloned_vector = new PetscVector<T>(this->comm(), this->type());
  cloned_vector->init(*this, true);
  *cloned_vector = *this;
  return UniquePtr<NumericVector<T> >(cloned_vector);
}



template <typename T>
inline
numeric_index_type PetscVector<T>::size () const
{
  libmesh_assert (this->initialized());

  PetscErrorCode ierr=0;
  PetscInt petsc_size=0;

  if (!this->initialized())
    return 0;

  ierr = VecGetSize(_vec, &petsc_size);
  LIBMESH_CHKERR(ierr);

  return static_cast<numeric_index_type>(petsc_size);
}



template <typename T>
inline
numeric_index_type PetscVector<T>::local_size () const
{
  libmesh_assert (this->initialized());

  PetscErrorCode ierr=0;
  PetscInt petsc_size=0;

  ierr = VecGetLocalSize(_vec, &petsc_size);
  LIBMESH_CHKERR(ierr);

  return static_cast<numeric_index_type>(petsc_size);
}



template <typename T>
inline
numeric_index_type PetscVector<T>::first_local_index () const
{
  libmesh_assert (this->initialized());

  numeric_index_type first = 0;

  if(_array_is_present) // Can we use cached values?
    first = _first;
  else
    {
      PetscErrorCode ierr=0;
      PetscInt petsc_first=0, petsc_last=0;
      ierr = VecGetOwnershipRange (_vec, &petsc_first, &petsc_last);
      LIBMESH_CHKERR(ierr);
      first = static_cast<numeric_index_type>(petsc_first);
    }

  return first;
}



template <typename T>
inline
numeric_index_type PetscVector<T>::last_local_index () const
{
  libmesh_assert (this->initialized());

  numeric_index_type last = 0;

  if(_array_is_present) // Can we use cached values?
    last = _last;
  else
    {
      PetscErrorCode ierr=0;
      PetscInt petsc_first=0, petsc_last=0;
      ierr = VecGetOwnershipRange (_vec, &petsc_first, &petsc_last);
      LIBMESH_CHKERR(ierr);
      last = static_cast<numeric_index_type>(petsc_last);
    }

  return last;
}



template <typename T>
inline
numeric_index_type PetscVector<T>::map_global_to_local_index (const numeric_index_type i) const
{
  libmesh_assert (this->initialized());

  numeric_index_type first=0;
  numeric_index_type last=0;

  if(_array_is_present) // Can we use cached values?
    {
      first = _first;
      last = _last;
    }
  else
    {
      PetscErrorCode ierr=0;
      PetscInt petsc_first=0, petsc_last=0;
      ierr = VecGetOwnershipRange (_vec, &petsc_first, &petsc_last);
      LIBMESH_CHKERR(ierr);
      first = static_cast<numeric_index_type>(petsc_first);
      last = static_cast<numeric_index_type>(petsc_last);
    }


  if((i>=first) && (i<last))
    {
      return i-first;
    }

  GlobalToLocalMap::const_iterator it = _global_to_local_map.find(i);
#ifndef NDEBUG
  const GlobalToLocalMap::const_iterator end = _global_to_local_map.end();
  if (it == end)
    {
      std::ostringstream error_message;
      error_message << "No index " << i << " in ghosted vector.\n"
                    << "Vector contains [" << first << ',' << last << ")\n";
      GlobalToLocalMap::const_iterator b = _global_to_local_map.begin();
      if (b == end)
        {
          error_message << "And empty ghost array.\n";
        }
      else
        {
          error_message << "And ghost array {" << b->first;
          for (++b; b != end; ++b)
            error_message << ',' << b->first;
          error_message << "}\n";
        }

      libmesh_error_msg(error_message.str());
    }
  libmesh_assert (it != _global_to_local_map.end());
#endif
  return it->second+last-first;
}



template <typename T>
inline
T PetscVector<T>::operator() (const numeric_index_type i) const
{
  this->_get_array();

  const numeric_index_type local_index = this->map_global_to_local_index(i);

#ifndef NDEBUG
  if(this->type() == GHOSTED)
    {
      libmesh_assert_less (local_index, _local_size);
    }
#endif

  return static_cast<T>(_values[local_index]);
}



template <typename T>
inline
void PetscVector<T>::get(const std::vector<numeric_index_type> & index,
                         T * values) const
{
  this->_get_array();

  const std::size_t num = index.size();

  for(std::size_t i=0; i<num; i++)
    {
      const numeric_index_type local_index = this->map_global_to_local_index(index[i]);
#ifndef NDEBUG
      if(this->type() == GHOSTED)
        {
          libmesh_assert_less (local_index, _local_size);
        }
#endif
      values[i] = static_cast<T>(_values[local_index]);
    }
}



template <typename T>
inline
Real PetscVector<T>::min () const
{
  parallel_object_only();

  this->_restore_array();

  PetscErrorCode ierr=0;
  PetscInt index=0;
  PetscReal returnval=0.;

  ierr = VecMin (_vec, &index, &returnval);
  LIBMESH_CHKERR(ierr);

  // this return value is correct: VecMin returns a PetscReal
  return static_cast<Real>(returnval);
}



template <typename T>
inline
Real PetscVector<T>::max() const
{
  parallel_object_only();

  this->_restore_array();

  PetscErrorCode ierr=0;
  PetscInt index=0;
  PetscReal returnval=0.;

  ierr = VecMax (_vec, &index, &returnval);
  LIBMESH_CHKERR(ierr);

  // this return value is correct: VecMax returns a PetscReal
  return static_cast<Real>(returnval);
}



template <typename T>
inline
void PetscVector<T>::swap (NumericVector<T> & other)
{
  parallel_object_only();

  NumericVector<T>::swap(other);

  PetscVector<T> & v = cast_ref<PetscVector<T> &>(other);

  std::swap(_vec, v._vec);
  std::swap(_destroy_vec_on_exit, v._destroy_vec_on_exit);
  std::swap(_global_to_local_map, v._global_to_local_map);

#ifdef LIBMESH_HAVE_CXX11_THREAD
  // Only truly atomic for v... but swap() doesn't really need to be thread safe!
  _array_is_present = v._array_is_present.exchange(_array_is_present);
#else
  std::swap(_array_is_present, v._array_is_present);
#endif

  std::swap(_local_form, v._local_form);
  std::swap(_values, v._values);
}





#ifdef LIBMESH_HAVE_CXX11
static_assert(sizeof(PetscInt) == sizeof(numeric_index_type),
              "PETSc and libMesh integer sizes must match!");
#endif


inline
PetscInt * numeric_petsc_cast(const numeric_index_type * p)
{
  return reinterpret_cast<PetscInt *>(const_cast<numeric_index_type *>(p));
}

} // namespace libMesh


#endif // #ifdef LIBMESH_HAVE_PETSC
#endif // LIBMESH_PETSC_VECTOR_H
