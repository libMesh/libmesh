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

#ifndef LIBMESH_TRILINOS_EPETRA_VECTOR_H
#define LIBMESH_TRILINOS_EPETRA_VECTOR_H


#include "libmesh/libmesh_common.h"


#ifdef LIBMESH_HAVE_TRILINOS

// Local includes
#include "libmesh/numeric_vector.h"
#include "libmesh/parallel.h"

// Trilinos includes
#include <Epetra_CombineMode.h>
#include <Epetra_Map.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>
#include <Epetra_MpiComm.h>

// C++ includes
#include <cstddef>
#include <vector>

// Forward declarations
class Epetra_IntSerialDenseVector;
class Epetra_SerialDenseVector;

namespace libMesh
{

// forward declarations
template <typename T> class SparseMatrix;

/**
 * Epetra vector. Provides a nice interface to the
 * Trilinos Epetra data structures for parallel vectors.
 *
 * @author Derek R. Gaston, 2008
 */

template <typename T>
class EpetraVector : public NumericVector<T>
{
public:

  /**
   *  Dummy-Constructor. Dimension=0
   */
  explicit
  EpetraVector (const Parallel::Communicator &comm,
		const ParallelType type = AUTOMATIC);

  /**
   * Constructor. Set dimension to \p n and initialize all elements with zero.
   */
  explicit
  EpetraVector (const Parallel::Communicator &comm,
		const numeric_index_type n,
                const ParallelType type = AUTOMATIC);

  /**
   * Constructor. Set local dimension to \p n_local, the global dimension
   * to \p n, and initialize all elements with zero.
   */
  EpetraVector (const Parallel::Communicator &comm,
		const numeric_index_type n,
	        const numeric_index_type n_local,
                const ParallelType type = AUTOMATIC);

  /**
   * Constructor. Set local dimension to \p n_local, the global
   * dimension to \p n, but additionally reserve memory for the
   * indices specified by the \p ghost argument.
   */
  EpetraVector (const Parallel::Communicator &comm,
		const numeric_index_type N,
		const numeric_index_type n_local,
		const std::vector<numeric_index_type>& ghost,
                const ParallelType type = AUTOMATIC);

  /**
   * Constructor.  Creates a EpetraVector assuming you already have a
   * valid Epetra Vec object.  In this case, v is NOT destroyed by the
   * EpetraVector constructor when this object goes out of scope.
   * This allows ownership of v to remain with the original creator,
   * and to simply provide additional functionality with the EpetraVector.
   */
  EpetraVector(Epetra_Vector & v,
	       const Parallel::Communicator &comm
	       LIBMESH_CAN_DEFAULT_TO_COMMWORLD);

  /**
   * Destructor, deallocates memory. Made virtual to allow
   * for derived classes to behave properly.
   */
  ~EpetraVector ();

  /**
   * Call the assemble functions
   */
  void close ();

  /**
   * @returns the \p EpetraVector<T> to a pristine state.
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
   * Create a vector that holds the local indices plus those specified
   * in the \p ghost argument.
   */
  void init (const numeric_index_type /*N*/,
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
  EpetraVector<T> & operator= (const EpetraVector<T> &V);

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
   * was renamed to get the \p EpetraVector<T> class
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
   * Access components, returns \p U(i).
   */
  T operator() (const numeric_index_type i) const;

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
   * of this vector by its complex conjugate, real(v_i) - imag(v_i).
   * Epetra is real-valued only, rendering this a no-op.
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
   * \f$ U+=v \f$ where \p v is a \p std::vector<T>
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
   * \f$U+=A*V\f$, add the product of the transpose of a \p SparseMatrix \p A_trans
   * and a \p NumericVector \p V to \p this, where \p this=U.
   */
  void add_vector_transpose (const NumericVector<T> &V,
		             const SparseMatrix<T> &A_trans);

  /**
   * \f$ U=v \f$ where v is a \p std::vector<T>
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
   * Computes the dot product, p = U.V
   */
  virtual T dot(const NumericVector<T>& V) const;

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
  void print_matlab (const std::string name="NULL") const;

  /**
   * Creates a "subvector" from this vector using the rows indices
   * of the "rows" array.
   */
  virtual void create_subvector (NumericVector<T>& subvector,
				 const std::vector<numeric_index_type>& rows) const;

  /**
   * Swaps the raw Epetra vector context pointers.
   */
  virtual void swap (NumericVector<T> &v);

  /**
   * Returns the raw PETSc vector context pointer.  Note this is generally
   * not required in user-level code. Just don't do anything crazy like
   * calling LibMeshVecDestroy()!
   */
  Epetra_Vector * vec () { libmesh_assert(_vec); return _vec; }

private:

  /**
   * Actual Epetra vector datatype
   * to hold vector entries
   */
  Epetra_Vector * _vec;

  /**
   * Holds the distributed Map
   */
  Epetra_Map * _map;

  /**
   * This boolean value should only be set to false
   * for the constructor which takes a Epetra Vec object.
   */
  bool _destroy_vec_on_exit;



  /*********************************************************************
   * The following were copied (and slightly modified) from
   * Epetra_FEVector.h in order to allow us to use a standard
   * Epetra_Vector... which is more compatible with other Trilinos
   * packages such as NOX.  All of this code is originally under LGPL
   *********************************************************************/

  /** Accumulate values into the vector, adding them to any values that
      already exist for the specified indices.
  */
  int SumIntoGlobalValues(int numIDs, const int* GIDs, const double* values);

  /** Accumulate values into the vector, adding them to any values that
      already exist for the specified GIDs.

      @param GIDs List of global ids. Must be the same length as the
      accompanying list of values.

      @param values List of coefficient values. Must be the same length as
      the accompanying list of GIDs.
  */
  int SumIntoGlobalValues(const Epetra_IntSerialDenseVector& GIDs,
                          const Epetra_SerialDenseVector& values);

  /** Copy values into the vector overwriting any values that already exist
      for the specified indices.
  */
  int ReplaceGlobalValues(int numIDs, const int* GIDs, const double* values);

  /** Copy values into the vector, replacing any values that
      already exist for the specified GIDs.

      @param GIDs List of global ids. Must be the same length as the
      accompanying list of values.

      @param values List of coefficient values. Must be the same length as
      the accompanying list of GIDs.
  */
  int ReplaceGlobalValues(const Epetra_IntSerialDenseVector& GIDs,
                          const Epetra_SerialDenseVector& values);

  int SumIntoGlobalValues(int numIDs, const int* GIDs,
                          const int* numValuesPerID,
                          const double* values);

  int ReplaceGlobalValues(int numIDs, const int* GIDs,
                          const int* numValuesPerID,
                          const double* values);

  /** Gather any overlapping/shared data into the non-overlapping partitioning
      defined by the Map that was passed to this vector at construction time.
      Data imported from other processors is stored on the owning processor
      with a "sumInto" or accumulate operation.
      This is a collective method -- every processor must enter it before any
      will complete it.
  */
  int GlobalAssemble(Epetra_CombineMode mode = Add);

  /** Set whether or not non-local data values should be ignored.
   */
  void setIgnoreNonLocalEntries(bool flag) {
    ignoreNonLocalEntries_ = flag;
  }

  void FEoperatorequals(const EpetraVector& source);

  int inputValues(int numIDs,
                  const int* GIDs, const double* values,
                  bool accumulate);

  int inputValues(int numIDs,
                  const int* GIDs, const int* numValuesPerID,
		  const double* values,
                  bool accumulate);

  int inputNonlocalValue(int GID, double value, bool accumulate);

  int inputNonlocalValues(int GID, int numValues, const double* values,
			  bool accumulate);

  void destroyNonlocalData();

  int myFirstID_;
  int myNumIDs_;
  double* myCoefs_;

  int* nonlocalIDs_;
  int* nonlocalElementSize_;
  int numNonlocalIDs_;
  int allocatedNonlocalLength_;
  double** nonlocalCoefs_;

  /**
   * Keep track of whether the last write operation on this vector was
   * nothing (0) or a sum (1) or an add (2), so we can decide how to
   * do the GlobalAssemble()
   */
  unsigned char last_edit;

  bool ignoreNonLocalEntries_;
};


/*----------------------- Inline functions ----------------------------------*/



template <typename T>
inline
EpetraVector<T>::EpetraVector (const Parallel::Communicator &comm,
			       const ParallelType type)
: NumericVector<T>(comm, type),
  _destroy_vec_on_exit(true),
  myFirstID_(0),
  myNumIDs_(0),
  myCoefs_(NULL),
  nonlocalIDs_(NULL),
  nonlocalElementSize_(NULL),
  numNonlocalIDs_(0),
  allocatedNonlocalLength_(0),
  nonlocalCoefs_(NULL),
  last_edit(0),
  ignoreNonLocalEntries_(false)
{
  this->_type = type;
}



template <typename T>
inline
EpetraVector<T>::EpetraVector (const Parallel::Communicator &comm,
			       const numeric_index_type n,
                               const ParallelType type)
: NumericVector<T>(comm, type),
  _destroy_vec_on_exit(true),
  myFirstID_(0),
  myNumIDs_(0),
  myCoefs_(NULL),
  nonlocalIDs_(NULL),
  nonlocalElementSize_(NULL),
  numNonlocalIDs_(0),
  allocatedNonlocalLength_(0),
  nonlocalCoefs_(NULL),
  last_edit(0),
  ignoreNonLocalEntries_(false)

{
  this->init(n, n, false, type);
}



template <typename T>
inline
EpetraVector<T>::EpetraVector (const Parallel::Communicator &comm,
			       const numeric_index_type n,
			       const numeric_index_type n_local,
                               const ParallelType type)
: NumericVector<T>(comm, type),
  _destroy_vec_on_exit(true),
  myFirstID_(0),
  myNumIDs_(0),
  myCoefs_(NULL),
  nonlocalIDs_(NULL),
  nonlocalElementSize_(NULL),
  numNonlocalIDs_(0),
  allocatedNonlocalLength_(0),
  nonlocalCoefs_(NULL),
  last_edit(0),
  ignoreNonLocalEntries_(false)
{
  this->init(n, n_local, false, type);
}




template <typename T>
inline
EpetraVector<T>::EpetraVector(Epetra_Vector & v,
			      const Parallel::Communicator &comm)
  : NumericVector<T>(comm, AUTOMATIC),
    _destroy_vec_on_exit(false),
    myFirstID_(0),
    myNumIDs_(0),
    myCoefs_(NULL),
    nonlocalIDs_(NULL),
    nonlocalElementSize_(NULL),
    numNonlocalIDs_(0),
    allocatedNonlocalLength_(0),
    nonlocalCoefs_(NULL),
    last_edit(0),
    ignoreNonLocalEntries_(false)
{
  _vec = &v;

  this->_type = PARALLEL; // FIXME - need to determine this from v!

  myFirstID_ = _vec->Map().MinMyGID();
  myNumIDs_ = _vec->Map().NumMyElements();

  _map = new Epetra_Map(_vec->GlobalLength(),
                        _vec->MyLength(),
                        0,
                        Epetra_MpiComm (this->comm().get()));

  //Currently we impose the restriction that NumVectors==1, so we won't
  //need the LDA argument when calling ExtractView. Hence the "dummy" arg.
  int dummy;
  _vec->ExtractView(&myCoefs_, &dummy);

  this->_is_closed = true;
  this->_is_initialized = true;
}



template <typename T>
inline
EpetraVector<T>::EpetraVector (const Parallel::Communicator &comm,
			       const numeric_index_type n,
			       const numeric_index_type n_local,
		               const std::vector<numeric_index_type>& ghost,
                               const ParallelType type)
: NumericVector<T>(comm, AUTOMATIC),
  _destroy_vec_on_exit(true),
  myFirstID_(0),
  myNumIDs_(0),
  myCoefs_(NULL),
  nonlocalIDs_(NULL),
  nonlocalElementSize_(NULL),
  numNonlocalIDs_(0),
  allocatedNonlocalLength_(0),
  nonlocalCoefs_(NULL),
  last_edit(0),
  ignoreNonLocalEntries_(false)
{
  this->init(n, n_local, ghost, false, type);
}



/* Default implementation for solver packages for which ghosted
   vectors are not yet implemented.  */
template <class T>
void EpetraVector<T>::init (const NumericVector<T>& other,
                            const bool fast)
{
  this->init(other.size(),other.local_size(),fast,other.type());
}



template <typename T>
inline
EpetraVector<T>::~EpetraVector ()
{
  this->clear ();
}



template <typename T>
inline
void EpetraVector<T>::init (const numeric_index_type n,
			    const numeric_index_type n_local,
			    const bool fast,
                            const ParallelType type)
{
  // We default to allocating n_local local storage
  numeric_index_type my_n_local = n_local;

  if (type == AUTOMATIC)
    {
      if (n == n_local)
        this->_type = SERIAL;
      else
        this->_type = PARALLEL;
    }
  else if (type == GHOSTED)
    {
      // We don't yet support GHOSTED Epetra vectors, so to get the
      // same functionality we need a SERIAL vector with local
      // storage allocated for every entry.
      this->_type = SERIAL;
      my_n_local = n;
    }
  else
    this->_type = type;

  libmesh_assert ((this->_type==SERIAL && n==my_n_local) ||
                  this->_type==PARALLEL);

  _map = new Epetra_Map(static_cast<int>(n),
                        my_n_local,
                        0,
                        Epetra_MpiComm (this->comm().get()));

  _vec = new Epetra_Vector(*_map);

  myFirstID_ = _vec->Map().MinMyGID();
  myNumIDs_ = _vec->Map().NumMyElements();

  //Currently we impose the restriction that NumVectors==1, so we won't
  //need the LDA argument when calling ExtractView. Hence the "dummy" arg.
  int dummy;
  _vec->ExtractView(&myCoefs_, &dummy);

  this->_is_initialized = true;
  this->_is_closed = true;
  this->last_edit = 0;

  if (fast == false)
    this->zero ();
}


template <typename T>
inline
void EpetraVector<T>::init (const numeric_index_type n,
			    const numeric_index_type n_local,
		            const std::vector<numeric_index_type>& /*ghost*/,
			    const bool fast,
                            const ParallelType type)
{
  // TODO: we shouldn't ignore the ghost sparsity pattern
  this->init(n, n_local, fast, type);
}



template <typename T>
inline
void EpetraVector<T>::init (const numeric_index_type n,
			    const bool fast,
                            const ParallelType type)
{
  this->init(n,n,fast,type);
}



template <typename T>
inline
void EpetraVector<T>::close ()
{
  libmesh_assert (this->initialized());

  // Are we adding or inserting?
  unsigned char global_last_edit = last_edit;
  this->comm().max(global_last_edit);
  libmesh_assert(!last_edit || last_edit == global_last_edit);

  if (global_last_edit == 1)
    this->GlobalAssemble(Insert);
  else if (global_last_edit == 2)
    this->GlobalAssemble(Add);
  else
    libmesh_assert(!global_last_edit);

  this->_is_closed = true;
  this->last_edit = 0;
}



template <typename T>
inline
void EpetraVector<T>::clear ()
{
  if (this->initialized())
  {
    // We might just be an interface to a user-provided _vec
    if (this->_destroy_vec_on_exit)
      {
        delete _vec;
        _vec = NULL;
      }

    // But we currently always own our own _map
    delete _map;
    _map = NULL;
  }

  this->_is_closed = this->_is_initialized = false;
}



template <typename T>
inline
void EpetraVector<T>::zero ()
{
  libmesh_assert (this->initialized());
  libmesh_assert (this->closed());

  _vec->PutScalar(0.0);
}



template <typename T>
inline
AutoPtr<NumericVector<T> > EpetraVector<T>::zero_clone () const
{
  AutoPtr<NumericVector<T> > cloned_vector
    (new EpetraVector<T>(this->comm(), AUTOMATIC));

  cloned_vector->init(*this);

  return cloned_vector;
}



template <typename T>
inline
AutoPtr<NumericVector<T> > EpetraVector<T>::clone () const
{
  AutoPtr<NumericVector<T> > cloned_vector
    (new EpetraVector<T>(this->comm(), AUTOMATIC));

  cloned_vector->init(*this, true);

  *cloned_vector = *this;

  return cloned_vector;
}



template <typename T>
inline
numeric_index_type EpetraVector<T>::size () const
{
  libmesh_assert (this->initialized());

  return _vec->GlobalLength();
}



template <typename T>
inline
numeric_index_type EpetraVector<T>::local_size () const
{
  libmesh_assert (this->initialized());

  return _vec->MyLength();
}

template <typename T>
inline
numeric_index_type EpetraVector<T>::first_local_index () const
{
  libmesh_assert (this->initialized());

  return _vec->Map().MinMyGID();
}



template <typename T>
inline
numeric_index_type EpetraVector<T>::last_local_index () const
{
  libmesh_assert (this->initialized());

  return _vec->Map().MaxMyGID()+1;
}


template <typename T>
inline
T EpetraVector<T>::operator() (const numeric_index_type i) const
{
  libmesh_assert (this->initialized());
  libmesh_assert ( ((i >= this->first_local_index()) &&
		    (i <  this->last_local_index())) );

  return (*_vec)[i-this->first_local_index()];
}



template <typename T>
inline
Real EpetraVector<T>::min () const
{
  libmesh_assert (this->initialized());

  T value;

  _vec->MinValue(&value);

  return value;
}



template <typename T>
inline
Real EpetraVector<T>::max() const
{
  libmesh_assert (this->initialized());

  T value;

  _vec->MaxValue(&value);

  return value;
}



template <typename T>
inline
void EpetraVector<T>::swap (NumericVector<T> &other)
{
  NumericVector<T>::swap(other);

  EpetraVector<T>& v = libmesh_cast_ref<EpetraVector<T>&>(other);

  std::swap(_vec, v._vec);
  std::swap(_map, v._map);
  std::swap(_destroy_vec_on_exit, v._destroy_vec_on_exit);
  std::swap(myFirstID_, v.myFirstID_);
  std::swap(myNumIDs_, v.myNumIDs_);
  std::swap(myCoefs_, v.myCoefs_);
  std::swap(nonlocalIDs_, v.nonlocalIDs_);
  std::swap(nonlocalElementSize_, v.nonlocalElementSize_);
  std::swap(numNonlocalIDs_, v.numNonlocalIDs_);
  std::swap(allocatedNonlocalLength_, v.allocatedNonlocalLength_);
  std::swap(nonlocalCoefs_, v.nonlocalCoefs_);
  std::swap(last_edit, v.last_edit);
  std::swap(ignoreNonLocalEntries_, v.ignoreNonLocalEntries_);
}

} // namespace libMesh


#endif // #ifdef HAVE_EPETRA
#endif // LIBMESH_TRILINOS_EPETRA_VECTOR_H
