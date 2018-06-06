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




#ifndef LIBMESH_EIGEN_SPARSE_VECTOR_H
#define LIBMESH_EIGEN_SPARSE_VECTOR_H



#include "libmesh/libmesh_common.h"

#ifdef LIBMESH_HAVE_EIGEN

// Local includes
#include "libmesh/eigen_core_support.h"
#include "libmesh/numeric_vector.h"

namespace libMesh
{

// Forward declarations
template <typename T> class EigenSparseMatrix;
template <typename T> class EigenSparseLinearSolver;
template <typename T> class SparseMatrix;

/**
 * This class provides a nice interface to the Eigen C++-based data
 * structures for serial vectors. All overridden virtual functions are
 * documented in numeric_vector.h.
 *
 * \author Benjamin S. Kirk
 * \date 2002
 */
template <typename T>
class EigenSparseVector libmesh_final : public NumericVector<T>
{
public:

  /**
   *  Dummy-Constructor. Dimension=0
   */
  explicit
  EigenSparseVector (const Parallel::Communicator & comm_in,
                     const ParallelType = AUTOMATIC);

  /**
   * Constructor. Set dimension to \p n and initialize all elements with zero.
   */
  explicit
  EigenSparseVector (const Parallel::Communicator & comm_in,
                     const numeric_index_type n,
                     const ParallelType = AUTOMATIC);

  /**
   * Constructor. Set local dimension to \p n_local, the global dimension
   * to \p n, and initialize all elements with zero.
   */
  EigenSparseVector (const Parallel::Communicator & comm_in,
                     const numeric_index_type n,
                     const numeric_index_type n_local,
                     const ParallelType = AUTOMATIC);

  /**
   * Constructor. Set local dimension to \p n_local, the global
   * dimension to \p n, but additionally reserve memory for the
   * indices specified by the \p ghost argument.
   */
  EigenSparseVector (const Parallel::Communicator & comm_in,
                     const numeric_index_type N,
                     const numeric_index_type n_local,
                     const std::vector<numeric_index_type> & ghost,
                     const ParallelType = AUTOMATIC);

  /**
   * Destructor, deallocates memory. Made virtual to allow
   * for derived classes to behave properly.
   */
  ~EigenSparseVector ();

  /**
   * Convenient typedefs
   */
  typedef EigenSV DataType;

  virtual void close () override;

  virtual void clear () override;

  virtual void zero () override;

  virtual std::unique_ptr<NumericVector<T>> zero_clone () const override;

  virtual std::unique_ptr<NumericVector<T>> clone () const override;

  virtual void init (const numeric_index_type N,
                     const numeric_index_type n_local,
                     const bool fast=false,
                     const ParallelType ptype=AUTOMATIC) override;

  virtual void init (const numeric_index_type N,
                     const bool fast=false,
                     const ParallelType ptype=AUTOMATIC) override;

  virtual void init (const numeric_index_type N,
                     const numeric_index_type n_local,
                     const std::vector<numeric_index_type> & ghost,
                     const bool fast = false,
                     const ParallelType = AUTOMATIC) override;

  virtual void init (const NumericVector<T> & other,
                     const bool fast = false) override;

  virtual NumericVector<T> & operator= (const T s) override;

  virtual NumericVector<T> & operator= (const NumericVector<T> & v) override;

  /**
   * Sets (*this)(i) = v(i) for each entry of the vector.
   *
   * \returns A reference to *this as the derived type.
   */
  EigenSparseVector<T> & operator= (const EigenSparseVector<T> & v);

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

  virtual T operator() (const numeric_index_type i) const override;

  virtual NumericVector<T> & operator += (const NumericVector<T> & v) override;

  virtual NumericVector<T> & operator -= (const NumericVector<T> & v) override;

  virtual NumericVector<T> & operator /= (NumericVector<T> & v_in) override;

  virtual void reciprocal() override;

  virtual void conjugate() override;

  virtual void set (const numeric_index_type i, const T value) override;

  virtual void add (const numeric_index_type i, const T value) override;

  virtual void add (const T s) override;

  virtual void add (const NumericVector<T> & v) override;

  virtual void add (const T a, const NumericVector<T> & v) override;

  /**
   * We override one NumericVector<T>::add_vector() method but don't
   * want to hide the other defaults.
   */
  using NumericVector<T>::add_vector;

  virtual void add_vector (const NumericVector<T> & v,
                           const SparseMatrix<T> & A) override;

  virtual void add_vector_transpose (const NumericVector<T> & v,
                                     const SparseMatrix<T> & A) override;

  virtual void scale (const T factor) override;

  virtual void abs() override;

  virtual T dot(const NumericVector<T> & v) const override;

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

  virtual void swap (NumericVector<T> & v) override;

  /**
   * References to the underlying Eigen data types.
   *
   * \note This is generally not required in user-level code.
   */
  DataType &       vec ()        { return _vec; }
  const DataType & vec () const  { return _vec; }

private:

  /**
   * Actual Eigen::SparseVector<> we are wrapping.
   */
  DataType _vec;

  /**
   * Make other Eigen datatypes friends
   */
  friend class EigenSparseMatrix<T>;
  friend class EigenSparseLinearSolver<T>;
};



// ---------------------------------------------------------
// EigenSparseVector inline methods
template <typename T>
inline
EigenSparseVector<T>::EigenSparseVector (const Parallel::Communicator & comm_in,
                                         const ParallelType ptype)
  : NumericVector<T>(comm_in, ptype)
{
  this->_type = ptype;
}



template <typename T>
inline
EigenSparseVector<T>::EigenSparseVector (const Parallel::Communicator & comm_in,
                                         const numeric_index_type n,
                                         const ParallelType ptype)
  : NumericVector<T>(comm_in, ptype)
{
  this->init(n, n, false, ptype);
}



template <typename T>
inline
EigenSparseVector<T>::EigenSparseVector (const Parallel::Communicator & comm_in,
                                         const numeric_index_type n,
                                         const numeric_index_type n_local,
                                         const ParallelType ptype)
  : NumericVector<T>(comm_in, ptype)
{
  this->init(n, n_local, false, ptype);
}



template <typename T>
inline
EigenSparseVector<T>::EigenSparseVector (const Parallel::Communicator & comm_in,
                                         const numeric_index_type N,
                                         const numeric_index_type n_local,
                                         const std::vector<numeric_index_type> & ghost,
                                         const ParallelType ptype)
  : NumericVector<T>(comm_in, ptype)
{
  this->init(N, n_local, ghost, false, ptype);
}



template <typename T>
inline
EigenSparseVector<T>::~EigenSparseVector ()
{
  this->clear ();
}



template <typename T>
inline
void EigenSparseVector<T>::init (const numeric_index_type n,
                                 const numeric_index_type n_local,
                                 const bool fast,
                                 const ParallelType)
{
  // Eigen vectors only for serial cases,
  // but can provide a "parallel" vector on one processor.
  if (n != n_local)
    libmesh_error_msg("Error: EigenSparseVectors can only be used in serial!");

  this->_type = SERIAL;

  // Clear initialized vectors
  if (this->initialized())
    this->clear();

  _vec.resize(n);

  this->_is_initialized = true;
#ifndef NDEBUG
  this->_is_closed = true;
#endif

  // Optionally zero out all components
  if (fast == false)
    this->zero ();

  return;
}



template <typename T>
inline
void EigenSparseVector<T>::init (const numeric_index_type n,
                                 const bool fast,
                                 const ParallelType ptype)
{
  this->init(n,n,fast,ptype);
}


template <typename T>
inline
void EigenSparseVector<T>::init (const numeric_index_type n,
                                 const numeric_index_type n_local,
                                 const std::vector<numeric_index_type> & libmesh_dbg_var(ghost),
                                 const bool fast,
                                 const ParallelType ptype)
{
  libmesh_assert(ghost.empty());
  this->init(n,n_local,fast,ptype);
}



/* Default implementation for solver packages for which ghosted
   vectors are not yet implemented.  */
template <class T>
void EigenSparseVector<T>::init (const NumericVector<T> & other,
                                 const bool fast)
{
  this->init(other.size(),other.local_size(),fast,other.type());
}



template <typename T>
inline
void EigenSparseVector<T>::close ()
{
  libmesh_assert (this->initialized());

#ifndef NDEBUG
  this->_is_closed = true;
#endif
}



template <typename T>
inline
void EigenSparseVector<T>::clear ()
{
  _vec.resize(0);

  this->_is_initialized = false;
#ifndef NDEBUG
  this->_is_closed = false;
#endif
}



template <typename T> inline
void EigenSparseVector<T>::zero ()
{
  libmesh_assert (this->initialized());
  libmesh_assert (this->closed());

  _vec.setZero();
}



template <typename T>
inline
std::unique_ptr<NumericVector<T>> EigenSparseVector<T>::zero_clone () const
{
  NumericVector<T> * cloned_vector = new EigenSparseVector<T>(this->comm());
  cloned_vector->init(*this);
  return std::unique_ptr<NumericVector<T>>(cloned_vector);
}



template <typename T>
inline
std::unique_ptr<NumericVector<T>> EigenSparseVector<T>::clone () const
{
  NumericVector<T> * cloned_vector = new EigenSparseVector<T>(this->comm());
  cloned_vector->init(*this, true);
  *cloned_vector = *this;
  return std::unique_ptr<NumericVector<T>>(cloned_vector);
}



template <typename T>
inline
numeric_index_type EigenSparseVector<T>::size () const
{
  libmesh_assert (this->initialized());

  return static_cast<numeric_index_type>(_vec.size());
}



template <typename T>
inline
numeric_index_type EigenSparseVector<T>::local_size () const
{
  libmesh_assert (this->initialized());

  return this->size();
}



template <typename T>
inline
numeric_index_type EigenSparseVector<T>::first_local_index () const
{
  libmesh_assert (this->initialized());

  return 0;
}



template <typename T>
inline
numeric_index_type EigenSparseVector<T>::last_local_index () const
{
  libmesh_assert (this->initialized());

  return this->size();
}



template <typename T>
inline
void EigenSparseVector<T>::set (const numeric_index_type i, const T value)
{
  libmesh_assert (this->initialized());
  libmesh_assert_less (i, this->size());

  _vec[static_cast<eigen_idx_type>(i)] = value;

#ifndef NDEBUG
  this->_is_closed = false;
#endif
}



template <typename T>
inline
void EigenSparseVector<T>::add (const numeric_index_type i, const T value)
{
  libmesh_assert (this->initialized());
  libmesh_assert_less (i, this->size());

  _vec[static_cast<eigen_idx_type>(i)] += value;

#ifndef NDEBUG
  this->_is_closed = false;
#endif
}



template <typename T>
inline
T EigenSparseVector<T>::operator() (const numeric_index_type i) const
{
  libmesh_assert (this->initialized());
  libmesh_assert ( ((i >= this->first_local_index()) &&
                    (i <  this->last_local_index())) );

  return _vec[static_cast<eigen_idx_type>(i)];
}



template <typename T>
inline
void EigenSparseVector<T>::swap (NumericVector<T> & other)
{
  EigenSparseVector<T> & v = cast_ref<EigenSparseVector<T> &>(other);

  _vec.swap(v._vec);

  std::swap (this->_is_closed,      v._is_closed);
  std::swap (this->_is_initialized, v._is_initialized);
  std::swap (this->_type,           v._type);
}


} // namespace libMesh


#endif // #ifdef LIBMESH_HAVE_EIGEN
#endif // LIBMESH_EIGEN_SPARSE_VECTOR_H
