// The libMesh Finite Element Library.
// Copyright (C) 2002-2022 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/int_range.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/petsc_solver_exception.h"
#include "libmesh/parallel_only.h"
#include "libmesh/enum_to_string.h"

// PETSc include files.
#ifdef I
# define LIBMESH_SAW_I
#endif
#include <petscvec.h>
#ifndef LIBMESH_SAW_I
# undef I // Avoid complex.h contamination
#endif

// C++ includes
#include <cstddef>
#include <cstring>
#include <vector>
#include <unordered_map>
#include <limits>

#include <atomic>
#include <mutex>
#include <condition_variable>

namespace libMesh
{

// forward declarations
template <typename T> class SparseMatrix;

/**
 * This class provides a nice interface to PETSc's Vec object.  All
 * overridden virtual functions are documented in numeric_vector.h.
 *
 * \author Benjamin S. Kirk
 * \date 2002
 * \brief NumericVector interface to PETSc Vec.
 */
template <typename T>
class PetscVector final : public NumericVector<T>
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
              const Parallel::Communicator & comm_in);

  /**
   * Copy assignment operator.
   * Calls VecCopy after performing various checks.
   * \returns A reference to *this as the derived type.
   */
  PetscVector<T> & operator= (const PetscVector<T> & v);

  /**
   * This class manages a C-style struct (Vec) manually, so we
   * don't want to allow any automatic copy/move functions to be
   * generated, and we can't default the destructor.
   */
  PetscVector (PetscVector &&) = delete;
  PetscVector (const PetscVector &) = delete;
  PetscVector & operator= (PetscVector &&) = delete;
  virtual ~PetscVector ();

  virtual void close () override;

  /**
   * clear() is called from the destructor, so it should not throw.
   */
  virtual void clear () noexcept override;

  virtual void zero () override;

  virtual std::unique_ptr<NumericVector<T>> zero_clone () const override;

  virtual std::unique_ptr<NumericVector<T>> clone () const override;

  virtual void init (const numeric_index_type N,
                     const numeric_index_type n_local,
                     const bool fast=false,
                     const ParallelType type=AUTOMATIC) override;

  virtual void init (const numeric_index_type N,
                     const bool fast=false,
                     const ParallelType type=AUTOMATIC) override;

  virtual void init (const numeric_index_type N,
                     const numeric_index_type n_local,
                     const std::vector<numeric_index_type> & ghost,
                     const bool fast = false,
                     const ParallelType = AUTOMATIC) override;

  virtual void init (const NumericVector<T> & other,
                     const bool fast = false) override;

  virtual NumericVector<T> & operator= (const T s) override;

  virtual NumericVector<T> & operator= (const NumericVector<T> & v) override;

  virtual NumericVector<T> & operator= (const std::vector<T> & v) override;

  virtual Real min () const override;

  virtual Real max () const override;

  virtual T sum () const override;

  virtual Real l1_norm () const override;

  virtual Real l2_norm () const override;

  virtual Real linfty_norm () const override;

  virtual numeric_index_type size () const override;

  virtual numeric_index_type local_size() const override;

  virtual numeric_index_type first_local_index() const override;

  virtual numeric_index_type last_local_index() const override;

  /**
   * \returns The local index corresponding to global index \p i.
   *
   * If the index is not a ghost cell, this is done by subtraction the
   * number of the first local index.  If it is a ghost cell, it has
   * to be looked up in the map.
   */
  numeric_index_type map_global_to_local_index(const numeric_index_type i) const;

  virtual T operator() (const numeric_index_type i) const override;

  virtual void get(const std::vector<numeric_index_type> & index,
                   T * values) const override;

  /**
   * Get read/write access to the raw PETSc Vector data array.
   *
   * \note This is an advanced interface. In general you should avoid
   * using it unless you know what you are doing!
   */
  PetscScalar * get_array();

  /**
   * Get read only access to the raw PETSc Vector data array.
   *
   * \note This is an advanced interface. In general you should avoid
   * using it unless you know what you are doing!
   */
  const PetscScalar * get_array_read() const;

  /**
   * Restore the data array.
   *
   * \note This MUST be called after get_array() or get_array_read()
   * and before using any other interface functions on PetscVector.
   */
  void restore_array();

  virtual NumericVector<T> & operator += (const NumericVector<T> & v) override;

  virtual NumericVector<T> & operator -= (const NumericVector<T> & v) override;

  virtual void reciprocal() override;

  virtual void conjugate() override;

  virtual void set (const numeric_index_type i,
                    const T value) override;

  virtual void add (const numeric_index_type i,
                    const T value) override;

  virtual void add (const T s) override;

  virtual void add (const NumericVector<T> & v) override;

  virtual void add (const T a, const NumericVector<T> & v) override;

  /**
   * We override two NumericVector<T>::add_vector() methods but don't
   * want to hide the other defaults.
   */
  using NumericVector<T>::add_vector;

  virtual void add_vector (const T * v,
                           const std::vector<numeric_index_type> & dof_indices) override;

  virtual void add_vector (const NumericVector<T> & v,
                           const SparseMatrix<T> & A) override;

  virtual void add_vector_transpose (const NumericVector<T> & v,
                                     const SparseMatrix<T> & A) override;

  /**
   * \f$ U \leftarrow U + A^H v \f$.
   *
   * Adds the product of the conjugate-transpose of \p SparseMatrix \p
   * A and a \p NumericVector \p v to \p this.
   */
  void add_vector_conjugate_transpose (const NumericVector<T> & v,
                                       const SparseMatrix<T> & A);

  /**
   * We override one NumericVector<T>::insert() method but don't want
   * to hide the other defaults
   */
  using NumericVector<T>::insert;

  virtual void insert (const T * v,
                       const std::vector<numeric_index_type> & dof_indices) override;

  virtual void scale (const T factor) override;

  virtual NumericVector<T> & operator *= (const NumericVector<T> & v) override;

  virtual NumericVector<T> & operator /= (const NumericVector<T> & v) override;

  virtual void abs() override;

  virtual T dot(const NumericVector<T> & v) const override;

  /**
   * \returns The dot product of (*this) with the vector \p v.
   *
   * \note Does *not* use the complex-conjugate of v in the complex-valued case.
   */
  T indefinite_dot(const NumericVector<T> & v) const;

  virtual void localize (std::vector<T> & v_local) const override;

  virtual void localize (NumericVector<T> & v_local) const override;

  virtual void localize (NumericVector<T> & v_local,
                         const std::vector<numeric_index_type> & send_list) const override;

  virtual void localize (std::vector<T> & v_local,
                         const std::vector<numeric_index_type> & indices) const override;

  virtual void localize (const numeric_index_type first_local_idx,
                         const numeric_index_type last_local_idx,
                         const std::vector<numeric_index_type> & send_list) override;

  virtual void localize_to_one (std::vector<T> & v_local,
                                const processor_id_type proc_id=0) const override;

  virtual void pointwise_mult (const NumericVector<T> & vec1,
                               const NumericVector<T> & vec2) override;

  virtual void print_matlab(const std::string & name = "") const override;

  virtual void create_subvector(NumericVector<T> & subvector,
                                const std::vector<numeric_index_type> & rows) const override;

  virtual void swap (NumericVector<T> & v) override;

  virtual std::size_t max_allowed_id() const override;

  /**
   * \returns The raw PETSc Vec pointer.
   *
   * \note This is generally not required in user-level code.
   *
   * \note Don't do anything crazy like calling VecDestroy() on
   * it, or very bad things will likely happen!
   */
  Vec vec () { libmesh_assert (_vec); return _vec; }

  Vec vec () const { libmesh_assert (_vec); return _vec; }


private:

  /**
   * Actual PETSc vector datatype to hold vector entries.
   */
  Vec _vec;

  /**
   * If \p true, the actual PETSc array of the values of the vector is
   * currently accessible.  That means that the members \p _local_form
   * and \p _values are valid.
   */
#ifdef LIBMESH_HAVE_CXX11_THREAD
  // We can't use std::atomic_flag here because we need load and store operations.
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

  /**
   * Size of the local values from _get_array()
   */
  mutable numeric_index_type _local_size;

  /**
   * PETSc vector datatype to hold the local form of a ghosted vector.
   * The contents of this field are only valid if the vector is
   * ghosted and \p _array_is_present is \p true.
   */
  mutable Vec _local_form;

  /**
   * Pointer to the actual PETSc array of the values of the vector.
   * This pointer is only valid if \p _array_is_present is \p true.
   * We're using PETSc's VecGetArrayRead() function, which requires a
   * constant PetscScalar *, but _get_array and _restore_array are
   * const member functions, so _values also needs to be mutable
   * (otherwise it is a "const PetscScalar * const" in that context).
   */
  mutable const PetscScalar * _read_only_values;

  /**
   * Pointer to the actual PETSc array of the values of the vector.
   * This pointer is only valid if \p _array_is_present is \p true.
   * We're using PETSc's VecGetArrayRead() function, which requires a
   * constant PetscScalar *, but _get_array and _restore_array are
   * const member functions, so _values also needs to be mutable
   * (otherwise it is a "const PetscScalar * const" in that context).
   */
  mutable PetscScalar * _values;

  /**
   * Mutex for _get_array and _restore_array.  This is part of the
   * object to keep down thread contention when reading from multiple
   * PetscVectors simultaneously
   */
  mutable std::mutex _petsc_get_restore_array_mutex;

  /**
   * Queries the array (and the local form if the vector is ghosted)
   * from PETSc.
   *
   * \param read_only Whether or not a read only copy of the raw data
   */
  void _get_array(bool read_only) const;

  /**
   * Restores the array (and the local form if the vector is ghosted)
   * to PETSc.
   */
  void _restore_array() const;

  /**
   * Type for map that maps global to local ghost cells.
   */
  typedef std::unordered_map<numeric_index_type,numeric_index_type> GlobalToLocalMap;

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

  /**
   * Whether or not the data array has been manually retrieved using
   * get_array() or get_array_read()
   */
  mutable bool _values_manually_retrieved;

  /**
   * Whether or not the data array is for read only access
   */
  mutable bool _values_read_only;
};


/*----------------------- Inline functions ----------------------------------*/



template <typename T>
inline
PetscVector<T>::PetscVector (const Parallel::Communicator & comm_in, const ParallelType ptype) :
  NumericVector<T>(comm_in, ptype),
  _array_is_present(false),
  _first(0),
  _last(0),
  _local_form(nullptr),
  _values(nullptr),
  _global_to_local_map(),
  _destroy_vec_on_exit(true),
  _values_manually_retrieved(false),
  _values_read_only(false)
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
  _local_form(nullptr),
  _values(nullptr),
  _global_to_local_map(),
  _destroy_vec_on_exit(true),
  _values_manually_retrieved(false),
  _values_read_only(false)
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
  _local_form(nullptr),
  _values(nullptr),
  _global_to_local_map(),
  _destroy_vec_on_exit(true),
  _values_manually_retrieved(false),
  _values_read_only(false)
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
  _local_form(nullptr),
  _values(nullptr),
  _global_to_local_map(),
  _destroy_vec_on_exit(true),
  _values_manually_retrieved(false),
  _values_read_only(false)
{
  this->init(n, n_local, ghost, false, ptype);
}





template <typename T>
inline
PetscVector<T>::PetscVector (Vec v,
                             const Parallel::Communicator & comm_in) :
  NumericVector<T>(comm_in, AUTOMATIC),
  _array_is_present(false),
  _local_form(nullptr),
  _values(nullptr),
  _global_to_local_map(),
  _destroy_vec_on_exit(false),
  _values_manually_retrieved(false),
  _values_read_only(false)
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
  VecType ptype;
  ierr = VecGetType(_vec, &ptype);
  LIBMESH_CHKERR(ierr);

  if ((std::strcmp(ptype,VECSHARED) == 0) || (std::strcmp(ptype,VECMPI) == 0))
    {
      ISLocalToGlobalMapping mapping;
      ierr = VecGetLocalToGlobalMapping(_vec, &mapping);
      LIBMESH_CHKERR(ierr);

      Vec localrep;
      ierr = VecGhostGetLocalForm(_vec,&localrep);
      LIBMESH_CHKERR(ierr);
      // If is a sparsely stored vector, set up our new mapping
      // Only checking mapping is not enough to determine if a vec is ghosted
      // We need to check if vec has a local representation
      if (mapping && localrep)
        {
          const numeric_index_type my_local_size = static_cast<numeric_index_type>(petsc_local_size);
          const numeric_index_type ghost_begin = static_cast<numeric_index_type>(petsc_local_size);
          PetscInt n;
          ierr = ISLocalToGlobalMappingGetSize(mapping, &n);
          LIBMESH_CHKERR(ierr);

          const numeric_index_type ghost_end = static_cast<numeric_index_type>(n);
          const PetscInt * indices;
          ierr = ISLocalToGlobalMappingGetIndices(mapping,&indices);
          LIBMESH_CHKERR(ierr);

          for (numeric_index_type i=ghost_begin; i<ghost_end; i++)
            _global_to_local_map[indices[i]] = i-my_local_size;
          this->_type = GHOSTED;
          ierr = ISLocalToGlobalMappingRestoreIndices(mapping, &indices);
          LIBMESH_CHKERR(ierr);
        }
      else
        this->_type = PARALLEL;

      ierr = VecGhostRestoreLocalForm(_vec,&localrep);
      LIBMESH_CHKERR(ierr);
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
      ierr = VecCreate(PETSC_COMM_SELF, &_vec);CHKERRABORT(PETSC_COMM_SELF,ierr);
      ierr = VecSetSizes(_vec, petsc_n, petsc_n); CHKERRABORT(PETSC_COMM_SELF,ierr);
      ierr = VecSetFromOptions (_vec);
      CHKERRABORT(PETSC_COMM_SELF,ierr);
    }
  // otherwise create an MPI-enabled vector
  else if (this->_type == PARALLEL)
    {
#ifdef LIBMESH_HAVE_MPI
      PetscInt petsc_n_local=cast_int<PetscInt>(n_local);
      libmesh_assert_less_equal (n_local, n);
      // Use more generic function instead of VecCreateSeq/MPI
      ierr = VecCreate(this->comm().get(), &_vec);LIBMESH_CHKERR(ierr);
      ierr = VecSetSizes(_vec, petsc_n_local, petsc_n); LIBMESH_CHKERR(ierr);
#else
      libmesh_assert_equal_to (n_local, n);
      ierr = VecCreate(PETSC_COMM_SELF, &_vec);CHKERRABORT(PETSC_COMM_SELF,ierr);
      ierr = VecSetSizes(_vec, petsc_n, petsc_n); CHKERRABORT(PETSC_COMM_SELF,ierr);
#endif
      ierr = VecSetFromOptions (_vec);
      LIBMESH_CHKERR(ierr);
    }
  else
    libmesh_error_msg("Unsupported type " << Utility::enum_to_string(this->_type));

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
    const_cast<PetscInt *>(reinterpret_cast<const PetscInt *>(ghost.data()));

  // Clear initialized vectors
  if (this->initialized())
    this->clear();

  libmesh_assert(ptype == AUTOMATIC || ptype == GHOSTED);
  this->_type = GHOSTED;

  /* Make the global-to-local ghost cell map.  */
  for (auto i : index_range(ghost))
    {
      _global_to_local_map[ghost[i]] = i;
    }

  /* Create vector.  */
  ierr = VecCreateGhost (this->comm().get(), petsc_n_local, petsc_n,
                         petsc_n_ghost, petsc_ghost,
                         &_vec);
  LIBMESH_CHKERR(ierr);

  // Add a prefix so that we can at least distinguish a ghosted vector from a
  // nonghosted vector when using a petsc option.
  // PETSc does not fully support VecGhost on GPU yet. This change allows us to
  // trigger a nonghosted vector to use GPU without bothering the ghosted vectors.
  ierr = PetscObjectAppendOptionsPrefix((PetscObject)_vec,"ghost_");
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
  if (v.initialized())
    {
      v._restore_array();
    }

  this->_global_to_local_map = v._global_to_local_map;

  // Even if we're initializing sizes based on an uninitialized or
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

  VecAssemblyBeginEnd(this->comm(), _vec);

  if (this->type() == GHOSTED)
    VecGhostUpdateBeginEnd(this->comm(), _vec, INSERT_VALUES, SCATTER_FORWARD);

  this->_is_closed = true;
}



template <typename T>
inline
void PetscVector<T>::clear () noexcept
{
  exceptionless_parallel_object_only();

  if (this->initialized())
    this->_restore_array();

  if ((this->initialized()) && (this->_destroy_vec_on_exit))
    {
      // If we encounter an error here, print a warning but otherwise
      // keep going since we may be recovering from an exception.
      PetscErrorCode ierr = VecDestroy(&_vec);
      if (ierr)
        libmesh_warning("Warning: VecDestroy returned a non-zero error code which we ignored.");
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

  if (this->type() != GHOSTED)
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
std::unique_ptr<NumericVector<T>> PetscVector<T>::zero_clone () const
{
  NumericVector<T> * cloned_vector = new PetscVector<T>(this->comm(), this->type());
  cloned_vector->init(*this);
  return std::unique_ptr<NumericVector<T>>(cloned_vector);
}



template <typename T>
inline
std::unique_ptr<NumericVector<T>> PetscVector<T>::clone () const
{
  NumericVector<T> * cloned_vector = new PetscVector<T>(this->comm(), this->type());
  cloned_vector->init(*this, true);
  *cloned_vector = *this;
  return std::unique_ptr<NumericVector<T>>(cloned_vector);
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

  if (_array_is_present) // Can we use cached values?
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

  if (_array_is_present) // Can we use cached values?
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

  if (_array_is_present) // Can we use cached values?
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


  if ((i>=first) && (i<last))
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
  this->_get_array(true);

  const numeric_index_type local_index = this->map_global_to_local_index(i);

#ifndef NDEBUG
  if (this->type() == GHOSTED)
    {
      libmesh_assert_less (local_index, _local_size);
    }
#endif

  return static_cast<T>(_read_only_values[local_index]);
}



template <typename T>
inline
void PetscVector<T>::get(const std::vector<numeric_index_type> & index,
                         T * values) const
{
  this->_get_array(true);

  const std::size_t num = index.size();

  for (std::size_t i=0; i<num; i++)
    {
      const numeric_index_type local_index = this->map_global_to_local_index(index[i]);
#ifndef NDEBUG
      if (this->type() == GHOSTED)
        {
          libmesh_assert_less (local_index, _local_size);
        }
#endif
      values[i] = static_cast<T>(_read_only_values[local_index]);
    }
}


template <typename T>
inline
PetscScalar * PetscVector<T>::get_array()
{
  _values_manually_retrieved = true;
  _get_array(false);

  return _values;
}


template <typename T>
inline
const PetscScalar * PetscVector<T>::get_array_read() const
{
  _values_manually_retrieved = true;
  _get_array(true);

  return _read_only_values;
}

template <typename T>
inline
void PetscVector<T>::restore_array()
{
  // \note \p _values_manually_retrieved needs to be set to \p false
  // \e before calling \p _restore_array()!
  _values_manually_retrieved = false;
  _restore_array();
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



template <typename T>
inline
std::size_t PetscVector<T>::max_allowed_id () const
{
  // The PetscInt type is used for indexing, it may be either a signed
  // 4-byte or 8-byte integer depending on how PETSc is configured.
  return std::numeric_limits<PetscInt>::max();
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
