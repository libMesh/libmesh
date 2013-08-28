// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

/**
 * Petsc include files.
 */
EXTERN_C_FOR_PETSC_BEGIN
# include <petscvec.h>
EXTERN_C_FOR_PETSC_END

// C++ includes
#include <cstddef>
#include <cstring>
#include <vector>

namespace libMesh
{



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
  explicit
  PetscVector (const Parallel::Communicator &comm,
	       const ParallelType type = AUTOMATIC);

  /**
   * Constructor. Set dimension to \p n and initialize all elements with zero.
   */
  explicit
  PetscVector (const Parallel::Communicator &comm,
	       const numeric_index_type n,
               const ParallelType type = AUTOMATIC);

  /**
   * Constructor. Set local dimension to \p n_local, the global dimension
   * to \p n, and initialize all elements with zero.
   */
  PetscVector (const Parallel::Communicator &comm,
	       const numeric_index_type n,
	       const numeric_index_type n_local,
               const ParallelType type = AUTOMATIC);

  /**
   * Constructor. Set local dimension to \p n_local, the global
   * dimension to \p n, but additionally reserve memory for the
   * indices specified by the \p ghost argument.
   */
  PetscVector (const Parallel::Communicator &comm,
	       const numeric_index_type N,
	       const numeric_index_type n_local,
	       const std::vector<numeric_index_type>& ghost,
               const ParallelType type = AUTOMATIC);

  /**
   * Constructor.  Creates a PetscVector assuming you already have a
   * valid PETSc Vec object.  In this case, \p v is NOT destroyed by the
   * PetscVector constructor when this object goes out of scope.
   * This allows ownership of \p v to remain with the original creator,
   * and to simply provide additional functionality with the PetscVector.
   */
  PetscVector(Vec v,
	      const Parallel::Communicator &comm
	      LIBMESH_CAN_DEFAULT_TO_COMMWORLD);

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
   * Creates a vector which has the same type, size and partitioning
   * as this vector, but whose data is all zero.  Returns it in an \p
   * AutoPtr.
   */
  virtual AutoPtr<NumericVector<T> > zero_clone () const;

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

  void init (const numeric_index_type N,
	     const numeric_index_type n_local,
	     const bool         fast=false,
	     const ParallelType type=AUTOMATIC);

  /**
   * call init with n_local = N,
   */
  void init (const numeric_index_type N,
	     const bool         fast=false,
	     const ParallelType type=AUTOMATIC);

  /**
   * Create a vector that holds tha local indices plus those specified
   * in the \p ghost argument.
   */
  virtual void init (const numeric_index_type /*N*/,
		     const numeric_index_type /*n_local*/,
		     const std::vector<numeric_index_type>& /*ghost*/,
		     const bool /*fast*/ = false,
		     const ParallelType = AUTOMATIC);

  /**
   * Creates a vector that has the same dimension and storage type as
   * \p other, including ghost dofs.
   */
  virtual void init (const NumericVector<T>& other,
                     const bool fast = false);

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
  numeric_index_type size () const;

  /**
   * @returns the local size of the vector
   * (index_stop-index_start)
   */
  numeric_index_type local_size() const;

  /**
   * @returns the index of the first vector element
   * actually stored on this processor
   */
  numeric_index_type first_local_index() const;

  /**
   * @returns the index of the last vector element
   * actually stored on this processor
   */
  numeric_index_type last_local_index() const;

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
  T operator() (const numeric_index_type i) const;

  /**
   * Access multiple components at once.  Overloaded method that
   * should be faster (probably much faster) than calling \p
   * operator() individually for each index.
   */
  virtual void get(const std::vector<numeric_index_type>& index, std::vector<T>& values) const;

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
   * Replace each entry v_i of this vector by its reciprocal, 1/v_i.
   */
  virtual void reciprocal();

  /**
   * Replace each entry v_i = real(v_i) + imag(v_i)
   * of this vector by its complex conjugate, real(v_i) - imag(v_i)
   */
  virtual void conjugate();

  /**
   * v(i) = value
   */
  void set (const numeric_index_type i, const T value);

  /**
   * v(i) += value
   */
  void add (const numeric_index_type i, const T value);

  /**
   * \f$U(0-LIBMESH_DIM)+=s\f$.
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
		   const std::vector<numeric_index_type>& dof_indices);

  /**
   * \f$ U+=V \f$ where U and V are type
   * \p NumericVector<T> and you
   * want to specify WHERE to add
   * the \p NumericVector<T> V
   */
  void add_vector (const NumericVector<T>& V,
		   const std::vector<numeric_index_type>& dof_indices);


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
		   const std::vector<numeric_index_type>& dof_indices);

  /**
   * \f$U+=A^T*V\f$, add the product of the transpose
   * of \p SparseMatrix \p A_trans and a \p NumericVector \p V to
   * \p this, where \p this=U.
   */
  void add_vector_transpose (const NumericVector<T> &V,
		             const SparseMatrix<T> &A_trans);

  /**
   * \f$U+=A^H*V\f$, add the product of the conjugate-transpose
   * of \p SparseMatrix \p A_trans and a \p NumericVector \p V to
   * \p this, where \p this=U.
   */
  void add_vector_conjugate_transpose (const NumericVector<T> &V,
		                       const SparseMatrix<T> &A_trans);

  /**
   * \f$ U=v \f$ where v is a std::vector<T>
   * and you want to specify WHERE to insert it
   */
  virtual void insert (const std::vector<T>& v,
		       const std::vector<numeric_index_type>& dof_indices);

  /**
   * \f$U=V\f$, where U and V are type
   * NumericVector<T> and you
   * want to specify WHERE to insert
   * the NumericVector<T> V
   */
  virtual void insert (const NumericVector<T>& V,
		       const std::vector<numeric_index_type>& dof_indices);

  /**
   * \f$ U=V \f$ where V is type
   * DenseVector<T> and you
   * want to specify WHERE to insert it
   */
  virtual void insert (const DenseVector<T>& V,
		       const std::vector<numeric_index_type>& dof_indices);

  /**
   * \f$ U=V \f$ where V is type
   * DenseSubVector<T> and you
   * want to specify WHERE to insert it
   */
  virtual void insert (const DenseSubVector<T>& V,
		       const std::vector<numeric_index_type>& dof_indices);


  /**
   * Scale each element of the
   * vector by the given factor.
   */
  void scale (const T factor);

  /**
   * v = abs(v)... that is, each entry in v is replaced
   * by its absolute value.
   */
  virtual void abs();

  /**
   * Computes the dot product, p = U.V. Use complex-conjugate of V
   * in the complex-valued case.
   */
  virtual T dot(const NumericVector<T>&) const;

  /**
   * Computes the dot product, p = U.V. Do not use complex-conjugate of V
   * in the complex-valued case.
   */
  virtual T indefinite_dot(const NumericVector<T>&) const;

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
		 const std::vector<numeric_index_type>& send_list) const;

  /**
   * Updates a local vector with selected values from neighboring
   * processors, as defined by \p send_list.
   */
  void localize (const numeric_index_type first_local_idx,
		 const numeric_index_type last_local_idx,
		 const std::vector<numeric_index_type>& send_list);

  /**
   * Creates a local copy of the global vector in
   * \p v_local only on processor \p proc_id.  By
   * default the data is sent to processor 0.  This method
   * is useful for outputting data from one processor.
   */
  void localize_to_one (std::vector<T>& v_local,
			const processor_id_type proc_id=0) const;

  /**
   * Computes the pointwise (i.e. component-wise) product of \p vec1
   * and \p vec2 and stores the result in \p *this.
   */
  virtual void pointwise_mult (const NumericVector<T>& vec1,
			       const NumericVector<T>& vec2);

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
				const std::vector<numeric_index_type>& rows) const;

  /**
   * Swaps the raw PETSc vector context pointers.
   */
  virtual void swap (NumericVector<T> &v);

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
  mutable bool _array_is_present;

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
   */
  mutable PetscScalar* _values;

  /**
   * Queries the array (and the local form if the vector is ghosted)
   * from Petsc.
   */
  void _get_array(void) const;

  /**
   * Restores the array (and the local form if the vector is ghosted)
   * to Petsc.
   */
  void _restore_array(void) const;

  /**
   * Type for map that maps global to local ghost cells.
   */
  typedef std::map<numeric_index_type,numeric_index_type> GlobalToLocalMap;

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
PetscVector<T>::PetscVector (const Parallel::Communicator &comm, const ParallelType ptype)
    : NumericVector<T>(comm, ptype),
      _array_is_present(false),
      _local_form(NULL),
      _values(NULL),
      _global_to_local_map(),
      _destroy_vec_on_exit(true)
{
  this->_type = ptype;
}



template <typename T>
inline
PetscVector<T>::PetscVector (const Parallel::Communicator &comm,
			     const numeric_index_type n,
                             const ParallelType ptype)
  : NumericVector<T>(comm, ptype),
    _array_is_present(false),
    _local_form(NULL),
    _values(NULL),
    _global_to_local_map(),
    _destroy_vec_on_exit(true)
{
  this->init(n, n, false, ptype);
}



template <typename T>
inline
PetscVector<T>::PetscVector (const Parallel::Communicator &comm,
			     const numeric_index_type n,
			     const numeric_index_type n_local,
                             const ParallelType ptype)
  : NumericVector<T>(comm, ptype),
    _array_is_present(false),
    _local_form(NULL),
    _values(NULL),
    _global_to_local_map(),
    _destroy_vec_on_exit(true)
{
  this->init(n, n_local, false, ptype);
}



template <typename T>
inline
PetscVector<T>::PetscVector (const Parallel::Communicator &comm,
			     const numeric_index_type n,
			     const numeric_index_type n_local,
			     const std::vector<numeric_index_type>& ghost,
                             const ParallelType ptype)
  : NumericVector<T>(comm, ptype),
    _array_is_present(false),
    _local_form(NULL),
    _values(NULL),
    _global_to_local_map(),
    _destroy_vec_on_exit(true)
{
  this->init(n, n_local, ghost, false, ptype);
}





template <typename T>
inline
PetscVector<T>::PetscVector (Vec v,
			     const Parallel::Communicator &comm)
  : NumericVector<T>(comm, AUTOMATIC),
    _array_is_present(false),
    _local_form(NULL),
    _values(NULL),
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
  LIBMESH_CHKERRABORT(ierr);

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
  LIBMESH_CHKERRABORT(ierr);

  if((std::strcmp(ptype,VECSHARED) == 0) || (std::strcmp(ptype,VECMPI) == 0))
  {
#if PETSC_VERSION_RELEASE && PETSC_VERSION_LESS_THAN(3,1,1)
    ISLocalToGlobalMapping mapping = _vec->mapping;
#else
    ISLocalToGlobalMapping mapping;
    ierr = VecGetLocalToGlobalMapping(_vec, &mapping);
    LIBMESH_CHKERRABORT(ierr);
#endif

    // If is a sparsely stored vector, set up our new mapping
    if (mapping)
    {
      const numeric_index_type my_local_size = static_cast<numeric_index_type>(petsc_local_size);
      const numeric_index_type ghost_begin = static_cast<numeric_index_type>(petsc_local_size);
#if PETSC_VERSION_RELEASE && PETSC_VERSION_LESS_THAN(3,4,0)
      const numeric_index_type ghost_end = static_cast<numeric_index_type>(mapping->n);
#else
      PetscInt n;
      ierr = ISLocalToGlobalMappingGetSize(mapping, &n);
      LIBMESH_CHKERRABORT(ierr);
      const numeric_index_type ghost_end = static_cast<numeric_index_type>(n);
#endif
#if PETSC_VERSION_RELEASE && PETSC_VERSION_LESS_THAN(3,1,1)
      const PetscInt *indices = mapping->indices;
#else
      const PetscInt *indices;
      ierr = ISLocalToGlobalMappingGetIndices(mapping,&indices);
      LIBMESH_CHKERRABORT(ierr);
#endif
      for(numeric_index_type i=ghost_begin; i<ghost_end; i++)
        _global_to_local_map[indices[i]] = i-my_local_size;
      this->_type = GHOSTED;
#if !PETSC_VERSION_RELEASE || !PETSC_VERSION_LESS_THAN(3,1,1)
      ierr = ISLocalToGlobalMappingRestoreIndices(mapping, &indices);
      LIBMESH_CHKERRABORT(ierr);
#endif
    }
    else
      this->_type = PARALLEL;
  }
  else
    this->_type = SERIAL;

  this->close();
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
             LIBMESH_CHKERRABORT(ierr);
#else
      libmesh_assert_equal_to (n_local, n);
      ierr = VecCreateSeq (PETSC_COMM_SELF, petsc_n, &_vec);
             CHKERRABORT(PETSC_COMM_SELF,ierr);
#endif

      ierr = VecSetFromOptions (_vec);
             LIBMESH_CHKERRABORT(ierr);
    }
  else
    libmesh_error();

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
			   const std::vector<numeric_index_type>& ghost,
			   const bool fast,
                           const ParallelType libmesh_dbg_var(ptype))
{
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

  PetscInt* petsc_ghost = ghost.empty() ? PETSC_NULL :
    const_cast<PetscInt*>(reinterpret_cast<const PetscInt*>(&ghost[0]));

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
  LIBMESH_CHKERRABORT(ierr);

  ierr = VecSetFromOptions (_vec);
  LIBMESH_CHKERRABORT(ierr);

  this->_is_initialized = true;
  this->_is_closed = true;

  if (fast == false)
    this->zero ();
}



template <typename T>
inline
void PetscVector<T>::init (const NumericVector<T>& other,
                           const bool fast)
{
  // Clear initialized vectors
  if (this->initialized())
    this->clear();

  const PetscVector<T>& v = libmesh_cast_ref<const PetscVector<T>&>(other);

  // Other vector should restore array.
  if(v.initialized())
    {
      v._restore_array();
    }

  this->_global_to_local_map = v._global_to_local_map;
  this->_is_closed      = v._is_closed;
  this->_is_initialized = v._is_initialized;
  this->_type = v._type;

  if (v.size() != 0)
    {
      PetscErrorCode ierr = 0;

      ierr = VecDuplicate (v._vec, &this->_vec);
      LIBMESH_CHKERRABORT(ierr);
    }

  if (fast == false)
    this->zero ();
}



template <typename T>
inline
void PetscVector<T>::close ()
{
  this->_restore_array();

  PetscErrorCode ierr=0;

  ierr = VecAssemblyBegin(_vec);
  LIBMESH_CHKERRABORT(ierr);
  ierr = VecAssemblyEnd(_vec);
  LIBMESH_CHKERRABORT(ierr);

  if(this->type() == GHOSTED)
    {
      ierr = VecGhostUpdateBegin(_vec,INSERT_VALUES,SCATTER_FORWARD);
      LIBMESH_CHKERRABORT(ierr);
      ierr = VecGhostUpdateEnd(_vec,INSERT_VALUES,SCATTER_FORWARD);
      LIBMESH_CHKERRABORT(ierr);
    }

  this->_is_closed = true;
}



template <typename T>
inline
void PetscVector<T>::clear ()
{
  if (this->initialized())
    this->_restore_array();

  if ((this->initialized()) && (this->_destroy_vec_on_exit))
    {
      PetscErrorCode ierr=0;

      ierr = LibMeshVecDestroy(&_vec);
             LIBMESH_CHKERRABORT(ierr);
    }

  this->_is_closed = this->_is_initialized = false;

  _global_to_local_map.clear();
}



template <typename T>
inline
void PetscVector<T>::zero ()
{
  libmesh_assert(this->closed());

  this->_restore_array();

  PetscErrorCode ierr=0;

  PetscScalar z=0.;

  if(this->type() != GHOSTED)
    {
#if PETSC_VERSION_LESS_THAN(2,3,0)
      // 2.2.x & earlier style
      ierr = VecSet (&z, _vec);
      LIBMESH_CHKERRABORT(ierr);
#else
      // 2.3.x & newer
      ierr = VecSet (_vec, z);
      LIBMESH_CHKERRABORT(ierr);
#endif
    }
  else
    {
      /* Vectors that include ghost values require a special
	 handling.  */
      Vec loc_vec;
      ierr = VecGhostGetLocalForm (_vec,&loc_vec);
      LIBMESH_CHKERRABORT(ierr);
#if PETSC_VERSION_LESS_THAN(2,3,0)
      // 2.2.x & earlier style
      ierr = VecSet (&z, loc_vec);
      LIBMESH_CHKERRABORT(ierr);
#else
      // 2.3.x & newer
      ierr = VecSet (loc_vec, z);
      LIBMESH_CHKERRABORT(ierr);
#endif
      ierr = VecGhostRestoreLocalForm (_vec,&loc_vec);
      LIBMESH_CHKERRABORT(ierr);
    }
}



template <typename T>
inline
AutoPtr<NumericVector<T> > PetscVector<T>::zero_clone () const
{
  AutoPtr<NumericVector<T> > cloned_vector
    (new PetscVector<T>(this->comm(), this->type()));

  cloned_vector->init(*this);

  return cloned_vector;
}



template <typename T>
inline
AutoPtr<NumericVector<T> > PetscVector<T>::clone () const
{
  AutoPtr<NumericVector<T> > cloned_vector
    (new PetscVector<T>(this->comm(), this->type()));

  cloned_vector->init(*this, true);

  *cloned_vector = *this;

  return cloned_vector;
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
         LIBMESH_CHKERRABORT(ierr);

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
         LIBMESH_CHKERRABORT(ierr);

  return static_cast<numeric_index_type>(petsc_size);
}



template <typename T>
inline
numeric_index_type PetscVector<T>::first_local_index () const
{
  libmesh_assert (this->initialized());

  PetscErrorCode ierr=0;
  PetscInt petsc_first=0, petsc_last=0;

  ierr = VecGetOwnershipRange (_vec, &petsc_first, &petsc_last);
         LIBMESH_CHKERRABORT(ierr);

  return static_cast<numeric_index_type>(petsc_first);
}



template <typename T>
inline
numeric_index_type PetscVector<T>::last_local_index () const
{
  libmesh_assert (this->initialized());

  PetscErrorCode ierr=0;
  PetscInt petsc_first=0, petsc_last=0;

  ierr = VecGetOwnershipRange (_vec, &petsc_first, &petsc_last);
         LIBMESH_CHKERRABORT(ierr);

  return static_cast<numeric_index_type>(petsc_last);
}



template <typename T>
inline
numeric_index_type PetscVector<T>::map_global_to_local_index (const numeric_index_type i) const
{
  libmesh_assert (this->initialized());

  PetscErrorCode ierr=0;
  PetscInt petsc_first=0, petsc_last=0;
  ierr = VecGetOwnershipRange (_vec, &petsc_first, &petsc_last);
  LIBMESH_CHKERRABORT(ierr);
  const numeric_index_type first = static_cast<numeric_index_type>(petsc_first);
  const numeric_index_type last = static_cast<numeric_index_type>(petsc_last);

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
          for (b++; b != end; b++)
            error_message << ',' << b->first;
          error_message << "}\n";
        }

      libMesh::err << error_message.str();

      libmesh_error();
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
void PetscVector<T>::get(const std::vector<numeric_index_type>& index, std::vector<T>& values) const
{
  this->_get_array();

  const std::size_t num = index.size();
  values.resize(num);

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
  this->_restore_array();

  PetscErrorCode ierr=0;
  PetscInt index=0;
  PetscReal returnval=0.;

  ierr = VecMin (_vec, &index, &returnval);
         LIBMESH_CHKERRABORT(ierr);

  // this return value is correct: VecMin returns a PetscReal
  return static_cast<Real>(returnval);
}



template <typename T>
inline
Real PetscVector<T>::max() const
{
  this->_restore_array();

  PetscErrorCode ierr=0;
  PetscInt index=0;
  PetscReal returnval=0.;

  ierr = VecMax (_vec, &index, &returnval);
         LIBMESH_CHKERRABORT(ierr);

  // this return value is correct: VecMax returns a PetscReal
  return static_cast<Real>(returnval);
}



template <typename T>
inline
void PetscVector<T>::swap (NumericVector<T> &other)
{
  NumericVector<T>::swap(other);

  PetscVector<T>& v = libmesh_cast_ref<PetscVector<T>&>(other);

  std::swap(_vec, v._vec);
  std::swap(_destroy_vec_on_exit, v._destroy_vec_on_exit);
  std::swap(_global_to_local_map, v._global_to_local_map);
  std::swap(_array_is_present, v._array_is_present);
  std::swap(_local_form, v._local_form);
  std::swap(_values, v._values);
}



template <typename T>
inline
void PetscVector<T>::_get_array(void) const
{
  libmesh_assert (this->initialized());
  if(!_array_is_present)
    {
      PetscErrorCode ierr=0;
      if(this->type() != GHOSTED)
	{
	  ierr = VecGetArray(_vec, &_values);
	  LIBMESH_CHKERRABORT(ierr);
	}
      else
	{
	  ierr = VecGhostGetLocalForm (_vec,&_local_form);
	  LIBMESH_CHKERRABORT(ierr);
	  ierr = VecGetArray(_local_form, &_values);
	  LIBMESH_CHKERRABORT(ierr);
#ifndef NDEBUG
	  PetscInt my_local_size = 0;
	  ierr = VecGetLocalSize(_local_form, &my_local_size);
	  LIBMESH_CHKERRABORT(ierr);
	  _local_size = static_cast<numeric_index_type>(my_local_size);
#endif
	}
      _array_is_present = true;
    }
}



template <typename T>
inline
void PetscVector<T>::_restore_array(void) const
{
  libmesh_assert (this->initialized());
  if(_array_is_present)
    {
      PetscErrorCode ierr=0;
      if(this->type() != GHOSTED)
	{
	  ierr = VecRestoreArray (_vec, &_values);
	  LIBMESH_CHKERRABORT(ierr);
	  _values = NULL;
	}
      else
	{
	  ierr = VecRestoreArray (_local_form, &_values);
	  LIBMESH_CHKERRABORT(ierr);
	  _values = NULL;
	  ierr = VecGhostRestoreLocalForm (_vec,&_local_form);
	  LIBMESH_CHKERRABORT(ierr);
	  _local_form = NULL;
#ifndef NDEBUG
	  _local_size = 0;
#endif
	}
      _array_is_present = false;
    }
}


} // namespace libMesh


#endif // #ifdef LIBMESH_HAVE_PETSC
#endif // LIBMESH_PETSC_VECTOR_H
