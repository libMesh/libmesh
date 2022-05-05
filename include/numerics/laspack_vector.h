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




#ifndef LIBMESH_LASPACK_VECTOR_H
#define LIBMESH_LASPACK_VECTOR_H



#include "libmesh/libmesh_common.h"

#ifdef LIBMESH_HAVE_LASPACK

// Local includes
#include "libmesh/numeric_vector.h"

// Laspack includes
#include <operats.h>
#include <qvector.h>

// C++ includes
#include <cstdio> // for std::sprintf
#include <limits>
#include <mutex>

namespace libMesh
{

// Forward declarations
template <typename T> class LaspackLinearSolver;
template <typename T> class SparseMatrix;

/**
 * This class provides a nice interface to the Laspack C-based data
 * structures for serial vectors.  All overridden virtual functions
 * are documented in numeric_vector.h.
 *
 * \author Benjamin S. Kirk
 * \date 2002
 */
template <typename T>
class LaspackVector final : public NumericVector<T>
{
public:

  /**
   *  Dummy-Constructor. Dimension=0
   */
  explicit
  LaspackVector (const Parallel::Communicator & comm,
                 const ParallelType = AUTOMATIC);

  /**
   * Constructor. Set dimension to \p n and initialize all elements with zero.
   */
  explicit
  LaspackVector (const Parallel::Communicator & comm,
                 const numeric_index_type n,
                 const ParallelType = AUTOMATIC);

  /**
   * Constructor. Set local dimension to \p n_local, the global dimension
   * to \p n, and initialize all elements with zero.
   */
  LaspackVector (const Parallel::Communicator & comm,
                 const numeric_index_type n,
                 const numeric_index_type n_local,
                 const ParallelType = AUTOMATIC);

  /**
   * Constructor. Set local dimension to \p n_local, the global
   * dimension to \p n, but additionally reserve memory for the
   * indices specified by the \p ghost argument.
   */
  LaspackVector (const Parallel::Communicator & comm,
                 const numeric_index_type N,
                 const numeric_index_type n_local,
                 const std::vector<numeric_index_type> & ghost,
                 const ParallelType = AUTOMATIC);

  /**
   * Copy assignment operator.
   * Calls Asgn_VV() to assign the contents of one vector to another.
   * \returns A reference to *this as the derived type.
   */
  LaspackVector<T> & operator= (const LaspackVector<T> & v);

  /**
   * This class manages a C-style struct (QVector) manually, so we
   * don't want to allow any automatic copy/move functions to be
   * generated, and we can't default the destructor.
   */
  LaspackVector (LaspackVector &&) = delete;
  LaspackVector (const LaspackVector &) = delete;
  LaspackVector & operator= (LaspackVector &&) = delete;
  virtual ~LaspackVector ();

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

  virtual NumericVector<T> & operator *= (const NumericVector<T> & v) override;

  virtual NumericVector<T> & operator /= (const NumericVector<T> & v) override;

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

  virtual std::size_t max_allowed_id() const override;

private:

  /**
   * Actual Laspack vector datatype
   * to hold vector entries
   */
  QVector _vec;

  /**
   * Make other Laspack datatypes friends
   */
  friend class LaspackLinearSolver<T>;
};



//----------------------------------------------------------
// LaspackVector inline methods
template <typename T>
inline
LaspackVector<T>::LaspackVector (const Parallel::Communicator & comm,
                                 const ParallelType ptype)
  : NumericVector<T>(comm, ptype)
{
  this->_type = ptype;
}



template <typename T>
inline
LaspackVector<T>::LaspackVector (const Parallel::Communicator & comm,
                                 const numeric_index_type n,
                                 const ParallelType ptype)
  : NumericVector<T>(comm, ptype)
{
  this->init(n, n, false, ptype);
}



template <typename T>
inline
LaspackVector<T>::LaspackVector (const Parallel::Communicator & comm,
                                 const numeric_index_type n,
                                 const numeric_index_type n_local,
                                 const ParallelType ptype)
  : NumericVector<T>(comm, ptype)
{
  this->init(n, n_local, false, ptype);
}



template <typename T>
inline
LaspackVector<T>::LaspackVector (const Parallel::Communicator & comm,
                                 const numeric_index_type N,
                                 const numeric_index_type n_local,
                                 const std::vector<numeric_index_type> & ghost,
                                 const ParallelType ptype)
  : NumericVector<T>(comm, ptype)
{
  this->init(N, n_local, ghost, false, ptype);
}



template <typename T>
inline
LaspackVector<T>::~LaspackVector ()
{
  this->clear ();
}



template <typename T>
inline
void LaspackVector<T>::init (const numeric_index_type n,
                             const numeric_index_type libmesh_dbg_var(n_local),
                             const bool fast,
                             const ParallelType)
{
  // Laspack vectors only for serial cases,
  // but can provide a "parallel" vector on one processor.
  libmesh_assert_equal_to (n, n_local);

  this->_type = SERIAL;

  // Clear initialized vectors
  if (this->initialized())
    this->clear();

  // create a sequential vector

  static int cnt = 0;
  char foo[80];
  std::sprintf(foo,  "Vec-%d", cnt++);

  V_Constr(&_vec, const_cast<char *>(foo), n, Normal, _LPTrue);

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
void LaspackVector<T>::init (const numeric_index_type n,
                             const bool fast,
                             const ParallelType ptype)
{
  this->init(n,n,fast,ptype);
}


template <typename T>
inline
void LaspackVector<T>::init (const numeric_index_type n,
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
void LaspackVector<T>::init (const NumericVector<T> & other,
                             const bool fast)
{
  this->init(other.size(),other.local_size(),fast,other.type());
}



template <typename T>
inline
void LaspackVector<T>::close ()
{
  libmesh_assert (this->initialized());

#ifndef NDEBUG
  this->_is_closed = true;
#endif
}



template <typename T>
inline
void LaspackVector<T>::clear ()
{
  if (this->initialized())
    {
      V_Destr (&_vec);
    }

  this->_is_initialized = false;
#ifndef NDEBUG
  this->_is_closed = false;
#endif
}



template <typename T> inline
void LaspackVector<T>::zero ()
{
  libmesh_assert (this->initialized());
  libmesh_assert (this->closed());

  V_SetAllCmp (&_vec, 0.);
}



template <typename T>
inline
std::unique_ptr<NumericVector<T>> LaspackVector<T>::zero_clone () const
{
  NumericVector<T> * cloned_vector = new LaspackVector<T>(this->comm());

  cloned_vector->init(*this);

  return std::unique_ptr<NumericVector<T>>(cloned_vector);
}



template <typename T>
inline
std::unique_ptr<NumericVector<T>> LaspackVector<T>::clone () const
{
  NumericVector<T> * cloned_vector = new LaspackVector<T>(this->comm());

  cloned_vector->init(*this, true);

  *cloned_vector = *this;

  return std::unique_ptr<NumericVector<T>>(cloned_vector);
}



template <typename T>
inline
numeric_index_type LaspackVector<T>::size () const
{
  libmesh_assert (this->initialized());

  return static_cast<numeric_index_type>(V_GetDim(const_cast<QVector*>(&_vec)));
}



template <typename T>
inline
numeric_index_type LaspackVector<T>::local_size () const
{
  libmesh_assert (this->initialized());

  return this->size();
}



template <typename T>
inline
numeric_index_type LaspackVector<T>::first_local_index () const
{
  libmesh_assert (this->initialized());

  return 0;
}



template <typename T>
inline
numeric_index_type LaspackVector<T>::last_local_index () const
{
  libmesh_assert (this->initialized());

  return this->size();
}



template <typename T>
inline
void LaspackVector<T>::set (const numeric_index_type i, const T value)
{
  libmesh_assert (this->initialized());
  libmesh_assert_less (i, this->size());

  std::scoped_lock lock(this->_numeric_vector_mutex);
  V_SetCmp (&_vec, i+1, value);

#ifndef NDEBUG
  this->_is_closed = false;
#endif
}



template <typename T>
inline
void LaspackVector<T>::add (const numeric_index_type i, const T value)
{
  libmesh_assert (this->initialized());
  libmesh_assert_less (i, this->size());

  std::scoped_lock lock(this->_numeric_vector_mutex);
  V_AddCmp (&_vec, i+1, value);

#ifndef NDEBUG
  this->_is_closed = false;
#endif
}



template <typename T>
inline
T LaspackVector<T>::operator() (const numeric_index_type i) const
{
  libmesh_assert (this->initialized());
  libmesh_assert ( ((i >= this->first_local_index()) &&
                    (i <  this->last_local_index())) );


  return static_cast<T>(V_GetCmp(const_cast<QVector*>(&_vec), i+1));
}



template <typename T>
inline
void LaspackVector<T>::swap (NumericVector<T> & other)
{
  LaspackVector<T> & v = cast_ref<LaspackVector<T> &>(other);

  // This is all grossly dependent on Laspack version...

  std::swap(_vec.Name, v._vec.Name);
  std::swap(_vec.Dim, v._vec.Dim);
  std::swap(_vec.Instance, v._vec.Instance);
  std::swap(_vec.LockLevel, v._vec.LockLevel);
  std::swap(_vec.Multipl, v._vec.Multipl);
  std::swap(_vec.OwnData, v._vec.OwnData);

  // This should still be O(1), since _vec.Cmp is just a pointer to
  // data on the heap

  std::swap(_vec.Cmp, v._vec.Cmp);
}



template <typename T>
inline
std::size_t LaspackVector<T>::max_allowed_id () const
{
  // The QVector type declares a "size_t Dim;"
  return std::numeric_limits<std::size_t>::max();
}


} // namespace libMesh


#endif // #ifdef LIBMESH_HAVE_LASPACK
#endif // LIBMESH_LASPACK_VECTOR_H
