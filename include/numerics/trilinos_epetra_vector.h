// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifdef LIBMESH_TRILINOS_HAVE_EPETRA

// Local includes
#include "libmesh/numeric_vector.h"
#include "libmesh/parallel.h"

// Trilinos includes
#include "libmesh/ignore_warnings.h"
#include <Epetra_CombineMode.h>
#include <Epetra_Map.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>
#include <Epetra_MpiComm.h>
#include "libmesh/restore_warnings.h"

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
 * This class provides a nice interface to the Trilinos Epetra_Vector
 * object. All overridden virtual functions are documented in
 * numeric_vector.h.
 *
 * \author Derek R. Gaston
 * \date 2008
 */
template <typename T>
class EpetraVector libmesh_final : public NumericVector<T>
{
public:

  /**
   *  Dummy-Constructor. Dimension=0
   */
  explicit
  EpetraVector (const Parallel::Communicator & comm,
                const ParallelType type = AUTOMATIC);

  /**
   * Constructor. Set dimension to \p n and initialize all elements with zero.
   */
  explicit
  EpetraVector (const Parallel::Communicator & comm,
                const numeric_index_type n,
                const ParallelType type = AUTOMATIC);

  /**
   * Constructor. Set local dimension to \p n_local, the global dimension
   * to \p n, and initialize all elements with zero.
   */
  EpetraVector (const Parallel::Communicator & comm,
                const numeric_index_type n,
                const numeric_index_type n_local,
                const ParallelType type = AUTOMATIC);

  /**
   * Constructor. Set local dimension to \p n_local, the global
   * dimension to \p n, but additionally reserve memory for the
   * indices specified by the \p ghost argument.
   */
  EpetraVector (const Parallel::Communicator & comm,
                const numeric_index_type N,
                const numeric_index_type n_local,
                const std::vector<numeric_index_type> & ghost,
                const ParallelType type = AUTOMATIC);

  /**
   * Constructor.  Creates a EpetraVector assuming you already have a
   * valid Epetra Vec object.  In this case, v is NOT destroyed by the
   * EpetraVector constructor when this object goes out of scope.
   * This allows ownership of v to remain with the original creator,
   * and to simply provide additional functionality with the EpetraVector.
   */
  EpetraVector(Epetra_Vector & v,
               const Parallel::Communicator & comm
               LIBMESH_CAN_DEFAULT_TO_COMMWORLD);

  /**
   * Destructor, deallocates memory. Made virtual to allow
   * for derived classes to behave properly.
   */
  ~EpetraVector ();

  virtual void close () libmesh_override;

  virtual void clear () libmesh_override;

  virtual void zero () libmesh_override;

  virtual std::unique_ptr<NumericVector<T>> zero_clone () const libmesh_override;

  virtual std::unique_ptr<NumericVector<T>> clone () const libmesh_override;

  virtual void init (const numeric_index_type N,
                     const numeric_index_type n_local,
                     const bool fast=false,
                     const ParallelType type=AUTOMATIC) libmesh_override;

  virtual void init (const numeric_index_type N,
                     const bool fast=false,
                     const ParallelType type=AUTOMATIC) libmesh_override;

  virtual void init (const numeric_index_type N,
                     const numeric_index_type n_local,
                     const std::vector<numeric_index_type> & ghost,
                     const bool fast = false,
                     const ParallelType = AUTOMATIC) libmesh_override;

  virtual void init (const NumericVector<T> & other,
                     const bool fast = false) libmesh_override;

  virtual NumericVector<T> & operator= (const T s) libmesh_override;

  virtual NumericVector<T> & operator= (const NumericVector<T> & v) libmesh_override;

  /**
   * Sets (*this)(i) = v(i) for each entry of the vector.
   *
   * \returns A reference to *this as the derived type.
   */
  EpetraVector<T> & operator= (const EpetraVector<T> & v);

  virtual NumericVector<T> & operator= (const std::vector<T> & v) libmesh_override;

  virtual Real min () const libmesh_override;

  virtual Real max () const libmesh_override;

  virtual T sum () const libmesh_override;

  virtual Real l1_norm () const libmesh_override;

  virtual Real l2_norm () const libmesh_override;

  virtual Real linfty_norm () const libmesh_override;

  virtual numeric_index_type size () const libmesh_override;

  virtual numeric_index_type local_size() const libmesh_override;

  virtual numeric_index_type first_local_index() const libmesh_override;

  virtual numeric_index_type last_local_index() const libmesh_override;

  virtual T operator() (const numeric_index_type i) const libmesh_override;

  virtual NumericVector<T> & operator += (const NumericVector<T> & v) libmesh_override;

  virtual NumericVector<T> & operator -= (const NumericVector<T> & v) libmesh_override;

  virtual NumericVector<T> & operator /= (NumericVector<T> & v) libmesh_override;

  virtual void reciprocal() libmesh_override;

  virtual void conjugate() libmesh_override;

  virtual void set (const numeric_index_type i, const T value) libmesh_override;

  virtual void add (const numeric_index_type i, const T value) libmesh_override;

  virtual void add (const T s) libmesh_override;

  virtual void add (const NumericVector<T> & v) libmesh_override;

  virtual void add (const T a, const NumericVector<T> & v) libmesh_override;

  /**
   * We override two NumericVector<T>::add_vector() methods but don't
   * want to hide the other defaults.
   */
  using NumericVector<T>::add_vector;

  virtual void add_vector (const T * v,
                           const std::vector<numeric_index_type> & dof_indices) libmesh_override;

  virtual void add_vector (const NumericVector<T> & v,
                           const SparseMatrix<T> & A) libmesh_override;

  virtual void add_vector_transpose (const NumericVector<T> & v,
                                     const SparseMatrix<T> & A) libmesh_override;

  /**
   * We override one NumericVector<T>::insert() method but don't want
   * to hide the other defaults
   */
  using NumericVector<T>::insert;

  virtual void insert (const T * v,
                       const std::vector<numeric_index_type> & dof_indices) libmesh_override;

  virtual void scale (const T factor) libmesh_override;

  virtual void abs() libmesh_override;

  virtual T dot(const NumericVector<T> & v) const libmesh_override;

  virtual void localize (std::vector<T> & v_local) const libmesh_override;

  virtual void localize (NumericVector<T> & v_local) const libmesh_override;

  virtual void localize (NumericVector<T> & v_local,
                         const std::vector<numeric_index_type> & send_list) const libmesh_override;

  virtual void localize (std::vector<T> & v_local,
                         const std::vector<numeric_index_type> & indices) const libmesh_override;

  virtual void localize (const numeric_index_type first_local_idx,
                         const numeric_index_type last_local_idx,
                         const std::vector<numeric_index_type> & send_list) libmesh_override;

  virtual void localize_to_one (std::vector<T> & v_local,
                                const processor_id_type proc_id=0) const libmesh_override;

  virtual void pointwise_mult (const NumericVector<T> & vec1,
                               const NumericVector<T> & vec2) libmesh_override;

  virtual void create_subvector (NumericVector<T> & subvector,
                                 const std::vector<numeric_index_type> & rows) const libmesh_override;

  virtual void swap (NumericVector<T> & v) libmesh_override;

  /**
   * \returns The raw Epetra_Vector pointer.
   *
   * \note This is generally not required in user-level code.
   *
   * \note Don't do anything crazy like deleting the pointer, or very
   * bad things will likely happen!
   */
  Epetra_Vector * vec () { libmesh_assert(_vec); return _vec; }

private:

  /**
   * Actual Epetra vector datatype to hold vector entries.
   */
  Epetra_Vector * _vec;

  /**
   * Holds the distributed Map.
   */
  std::unique_ptr<Epetra_Map> _map;

  /**
   * This boolean value should only be set to false
   * for the constructor which takes a Epetra Vec object.
   */
  bool _destroy_vec_on_exit;

  // The following were copied (and slightly modified) from
  // Epetra_FEVector.h in order to allow us to use a standard
  // Epetra_Vector... which is more compatible with other Trilinos
  // packages such as NOX.  All of this code is originally under LGPL

  /**
   * Accumulate values into the vector, adding them to any values that
   * already exist for the specified indices.
   */
  int SumIntoGlobalValues(int numIDs,
                          const int * GIDs,
                          const double * values);

  /**
   * Accumulate values into the vector, adding them to any values that
   * already exist for the specified GIDs.
   *
   * \param GIDs List of global ids. Must be the same length as the
   * accompanying list of values.
   *
   * \param values List of coefficient values. Must be the same length as
   * the accompanying list of GIDs.
   */
  int SumIntoGlobalValues(const Epetra_IntSerialDenseVector & GIDs,
                          const Epetra_SerialDenseVector & values);

  /**
   * Copy values into the vector overwriting any values that already
   * exist for the specified indices.
   */
  int ReplaceGlobalValues(int numIDs,
                          const int * GIDs,
                          const double * values);

  /**
   * Copy values into the vector, replacing any values that
   * already exist for the specified GIDs.
   *
   * \param GIDs List of global ids. Must be the same length as the
   * accompanying list of values.
   *
   * \param values List of coefficient values. Must be the same length as
   * the accompanying list of GIDs.
   */
  int ReplaceGlobalValues(const Epetra_IntSerialDenseVector & GIDs,
                          const Epetra_SerialDenseVector & values);

  int SumIntoGlobalValues(int numIDs,
                          const int * GIDs,
                          const int * numValuesPerID,
                          const double * values);

  int ReplaceGlobalValues(int numIDs,
                          const int * GIDs,
                          const int * numValuesPerID,
                          const double * values);

  /**
   * Gather any overlapping/shared data into the non-overlapping
   * partitioning defined by the Map that was passed to this vector at
   * construction time.  Data imported from other processors is stored
   * on the owning processor with a "sumInto" or accumulate operation.
   * This is a collective method -- every processor must enter it
   * before any will complete it.
   */
  int GlobalAssemble(Epetra_CombineMode mode = Add);

  /**
   * Set whether or not non-local data values should be ignored.
   */
  void setIgnoreNonLocalEntries(bool flag) {
    ignoreNonLocalEntries_ = flag;
  }

  void FEoperatorequals(const EpetraVector & source);

  int inputValues(int numIDs,
                  const int * GIDs,
                  const double * values,
                  bool accumulate);

  int inputValues(int numIDs,
                  const int * GIDs,
                  const int * numValuesPerID,
                  const double * values,
                  bool accumulate);

  int inputNonlocalValue(int GID,
                         double value,
                         bool accumulate);

  int inputNonlocalValues(int GID,
                          int numValues,
                          const double * values,
                          bool accumulate);

  void destroyNonlocalData();

  int myFirstID_;
  int myNumIDs_;
  double * myCoefs_;

  int * nonlocalIDs_;
  int * nonlocalElementSize_;
  int numNonlocalIDs_;
  int allocatedNonlocalLength_;
  double ** nonlocalCoefs_;

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
EpetraVector<T>::EpetraVector (const Parallel::Communicator & comm,
                               const ParallelType type) :
  NumericVector<T>(comm, type),
  _destroy_vec_on_exit(true),
  myFirstID_(0),
  myNumIDs_(0),
  myCoefs_(libmesh_nullptr),
  nonlocalIDs_(libmesh_nullptr),
  nonlocalElementSize_(libmesh_nullptr),
  numNonlocalIDs_(0),
  allocatedNonlocalLength_(0),
  nonlocalCoefs_(libmesh_nullptr),
  last_edit(0),
  ignoreNonLocalEntries_(false)
{
  this->_type = type;
}



template <typename T>
inline
EpetraVector<T>::EpetraVector (const Parallel::Communicator & comm,
                               const numeric_index_type n,
                               const ParallelType type) :
  NumericVector<T>(comm, type),
  _destroy_vec_on_exit(true),
  myFirstID_(0),
  myNumIDs_(0),
  myCoefs_(libmesh_nullptr),
  nonlocalIDs_(libmesh_nullptr),
  nonlocalElementSize_(libmesh_nullptr),
  numNonlocalIDs_(0),
  allocatedNonlocalLength_(0),
  nonlocalCoefs_(libmesh_nullptr),
  last_edit(0),
  ignoreNonLocalEntries_(false)

{
  this->init(n, n, false, type);
}



template <typename T>
inline
EpetraVector<T>::EpetraVector (const Parallel::Communicator & comm,
                               const numeric_index_type n,
                               const numeric_index_type n_local,
                               const ParallelType type) :
  NumericVector<T>(comm, type),
  _destroy_vec_on_exit(true),
  myFirstID_(0),
  myNumIDs_(0),
  myCoefs_(libmesh_nullptr),
  nonlocalIDs_(libmesh_nullptr),
  nonlocalElementSize_(libmesh_nullptr),
  numNonlocalIDs_(0),
  allocatedNonlocalLength_(0),
  nonlocalCoefs_(libmesh_nullptr),
  last_edit(0),
  ignoreNonLocalEntries_(false)
{
  this->init(n, n_local, false, type);
}




template <typename T>
inline
EpetraVector<T>::EpetraVector(Epetra_Vector & v,
                              const Parallel::Communicator & comm) :
  NumericVector<T>(comm, AUTOMATIC),
  _destroy_vec_on_exit(false),
  myFirstID_(0),
  myNumIDs_(0),
  myCoefs_(libmesh_nullptr),
  nonlocalIDs_(libmesh_nullptr),
  nonlocalElementSize_(libmesh_nullptr),
  numNonlocalIDs_(0),
  allocatedNonlocalLength_(0),
  nonlocalCoefs_(libmesh_nullptr),
  last_edit(0),
  ignoreNonLocalEntries_(false)
{
  _vec = &v;

  this->_type = PARALLEL; // FIXME - need to determine this from v!

  myFirstID_ = _vec->Map().MinMyGID();
  myNumIDs_ = _vec->Map().NumMyElements();

  _map.reset(new Epetra_Map(_vec->GlobalLength(),
                            _vec->MyLength(),
                            0, // IndexBase = 0 for C/C++, 1 for Fortran.
                            Epetra_MpiComm (this->comm().get())));

  //Currently we impose the restriction that NumVectors==1, so we won't
  //need the LDA argument when calling ExtractView. Hence the "dummy" arg.
  int dummy;
  _vec->ExtractView(&myCoefs_, &dummy);

  this->_is_closed = true;
  this->_is_initialized = true;
}



template <typename T>
inline
EpetraVector<T>::EpetraVector (const Parallel::Communicator & comm,
                               const numeric_index_type n,
                               const numeric_index_type n_local,
                               const std::vector<numeric_index_type> & ghost,
                               const ParallelType type) :
  NumericVector<T>(comm, AUTOMATIC),
  _destroy_vec_on_exit(true),
  myFirstID_(0),
  myNumIDs_(0),
  myCoefs_(libmesh_nullptr),
  nonlocalIDs_(libmesh_nullptr),
  nonlocalElementSize_(libmesh_nullptr),
  numNonlocalIDs_(0),
  allocatedNonlocalLength_(0),
  nonlocalCoefs_(libmesh_nullptr),
  last_edit(0),
  ignoreNonLocalEntries_(false)
{
  this->init(n, n_local, ghost, false, type);
}



// Default implementation for solver packages for which ghosted
// vectors are not yet implemented.
template <class T>
void EpetraVector<T>::init (const NumericVector<T> & other,
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

  _map.reset(new Epetra_Map(static_cast<int>(n),
                            my_n_local,
                            0,
                            Epetra_MpiComm (this->comm().get())));

  _vec = new Epetra_Vector(*_map);

  myFirstID_ = _vec->Map().MinMyGID();
  myNumIDs_ = _vec->Map().NumMyElements();

  // Currently we impose the restriction that NumVectors==1, so we won't
  // need the LDA argument when calling ExtractView. Hence the "dummy" arg.
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
                            const std::vector<numeric_index_type> & /*ghost*/,
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
          _vec = libmesh_nullptr;
        }

      // But we currently always own our own _map
      _map.reset();
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
std::unique_ptr<NumericVector<T>> EpetraVector<T>::zero_clone () const
{
  NumericVector<T> * cloned_vector = new EpetraVector<T>(this->comm(), AUTOMATIC);
  cloned_vector->init(*this);
  return std::unique_ptr<NumericVector<T>>(cloned_vector);
}



template <typename T>
inline
std::unique_ptr<NumericVector<T>> EpetraVector<T>::clone () const
{
  NumericVector<T> * cloned_vector = new EpetraVector<T>(this->comm(), AUTOMATIC);
  cloned_vector->init(*this, true);
  *cloned_vector = *this;
  return std::unique_ptr<NumericVector<T>>(cloned_vector);
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
void EpetraVector<T>::swap (NumericVector<T> & other)
{
  NumericVector<T>::swap(other);

  EpetraVector<T> & v = cast_ref<EpetraVector<T> &>(other);

  std::swap(_vec, v._vec);
  _map.swap(v._map);
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


#endif // #ifdef LIBMESH_TRILINOS_HAVE_EPETRA
#endif // LIBMESH_TRILINOS_EPETRA_VECTOR_H
