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




#ifndef LIBMESH_EIGEN_SPARSE_VECTOR_H
#define LIBMESH_EIGEN_SPARSE_VECTOR_H



#include "libmesh/libmesh_common.h"

#ifdef LIBMESH_HAVE_EIGEN

// Local includes
#include "libmesh/eigen_core_support.h"
#include "libmesh/numeric_vector.h"

// C++ includes

// Eigen includes

namespace libMesh
{


// Forward declarations
template <typename T> class EigenSparseMatrix;
template <typename T> class EigenSparseLinearSolver;
template <typename T> class SparseMatrix;


/**
 * Eigen vector.  Provides a nice interface to the
 * Eigen C-based data structures for serial vectors.
 *
 * @author Benjamin S. Kirk, 2002
 */

template <typename T>
class EigenSparseVector : public NumericVector<T>
{
 public:

  /**
   *  Dummy-Constructor. Dimension=0
   */
  explicit
  EigenSparseVector (const Parallel::Communicator &comm,
		     const ParallelType = AUTOMATIC);

  /**
   * Constructor. Set dimension to \p n and initialize all elements with zero.
   */
  explicit
  EigenSparseVector (const Parallel::Communicator &comm,
		     const numeric_index_type n,
                     const ParallelType = AUTOMATIC);

  /**
   * Constructor. Set local dimension to \p n_local, the global dimension
   * to \p n, and initialize all elements with zero.
   */
  EigenSparseVector (const Parallel::Communicator &comm,
		     const numeric_index_type n,
		     const numeric_index_type n_local,
		     const ParallelType = AUTOMATIC);

  /**
   * Constructor. Set local dimension to \p n_local, the global
   * dimension to \p n, but additionally reserve memory for the
   * indices specified by the \p ghost argument.
   */
  EigenSparseVector (const Parallel::Communicator &comm,
		     const numeric_index_type N,
		     const numeric_index_type n_local,
		     const std::vector<numeric_index_type>& ghost,
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

  /**
   * Call the assemble functions
   */
  void close ();

  /**
   * @returns the \p EigenSparseVector to a pristine state.
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
             const ParallelType ptype=AUTOMATIC);

  /**
   * call init with n_local = N,
   */
  void init (const numeric_index_type N,
	     const bool         fast=false,
             const ParallelType ptype=AUTOMATIC);

  /**
   * Create a vector that holds tha local indices plus those specified
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
  EigenSparseVector<T> & operator= (const EigenSparseVector<T> &V);

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
   * was renamed to get the \p EigenSparseVector class
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
   * \f$ U+=V \f$.
   * Simple vector addition, equal to the
   * \p operator +=.
   */
  void add (const NumericVector<T>& V);

  /**
   * \f$ U+=a*V \f$.
   * Simple vector addition, equal to the
   * \p operator +=.
   */
  void add (const T a, const NumericVector<T>& v);

  /**
   * \f$ U+=v \f$ where v is a \p std::vector<T>
   * and you
   * want to specify WHERE to add it
   */
  void add_vector (const std::vector<T>& v,
		   const std::vector<numeric_index_type>& dof_indices);

  /**
   * \f$ U+=V \f$ where U and V are type
   * NumericVector<T> and you
   * want to specify WHERE to add
   * the NumericVector<T> V
   */
  void add_vector (const NumericVector<T>& V,
		   const std::vector<numeric_index_type>& dof_indices);

  /**
   * \f$U+=A*V\f$, add the product of a \p SparseMatrix \p A
   * and a \p NumericVector \p V to \p this, where \p this=U.
   */
  void add_vector (const NumericVector<T> &,
		   const SparseMatrix<T> &);

  /**
   * \f$U+=V \f$ where U and V are type
   * DenseVector<T> and you
   * want to specify WHERE to add
   * the DenseVector<T> V
   */
  void add_vector (const DenseVector<T>& V,
		   const std::vector<numeric_index_type>& dof_indices);

  /**
   * \f$U+=A^T*V\f$, add the product of the transpose of a \p SparseMatrix \p A_trans
   * and a \p NumericVector \p V to \p this, where \p this=U.
   */
  void add_vector_transpose (const NumericVector<T> &,
		             const SparseMatrix<T> &);

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
   * Swaps the contents.
   */
  virtual void swap (NumericVector<T> &v);

  /**
   * References to the underlying Eigen data types. Note this is generally
   * not required in user-level code.
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



//----------------------- ----------------------------------
// EigenSparseVector inline methods
template <typename T>
inline
EigenSparseVector<T>::EigenSparseVector (const Parallel::Communicator &comm,
					 const ParallelType ptype)
  : NumericVector<T>(comm, ptype)
{
  this->_type = ptype;
}



template <typename T>
inline
EigenSparseVector<T>::EigenSparseVector (const Parallel::Communicator &comm,
					 const numeric_index_type n,
					 const ParallelType ptype)
  : NumericVector<T>(comm, ptype)
{
  this->init(n, n, false, ptype);
}



template <typename T>
inline
EigenSparseVector<T>::EigenSparseVector (const Parallel::Communicator &comm,
					 const numeric_index_type n,
					 const numeric_index_type n_local,
					 const ParallelType ptype)
  : NumericVector<T>(comm, ptype)
{
  this->init(n, n_local, false, ptype);
}



template <typename T>
inline
EigenSparseVector<T>::EigenSparseVector (const Parallel::Communicator &comm,
					 const numeric_index_type N,
					 const numeric_index_type n_local,
					 const std::vector<numeric_index_type>& ghost,
					 const ParallelType ptype)
  : NumericVector<T>(comm, ptype)
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
				 const numeric_index_type libmesh_dbg_var(n_local),
				 const bool fast,
				 const ParallelType)
{
  // Eigen vectors only for serial cases,
  // but can provide a "parallel" vector on one processor.
  libmesh_assert_equal_to (n, n_local);

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
				 const std::vector<numeric_index_type>& libmesh_dbg_var(ghost),
				 const bool fast,
				 const ParallelType ptype)
{
  libmesh_assert(ghost.empty());
  this->init(n,n_local,fast,ptype);
}



/* Default implementation for solver packages for which ghosted
   vectors are not yet implemented.  */
template <class T>
void EigenSparseVector<T>::init (const NumericVector<T>& other,
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
AutoPtr<NumericVector<T> > EigenSparseVector<T>::zero_clone () const
{
  AutoPtr<NumericVector<T> > cloned_vector
    (new EigenSparseVector<T>(this->comm()));

  cloned_vector->init(*this);

  return cloned_vector;
}



template <typename T>
inline
AutoPtr<NumericVector<T> > EigenSparseVector<T>::clone () const
{
  AutoPtr<NumericVector<T> > cloned_vector
    (new EigenSparseVector<T>(this->comm()));

  cloned_vector->init(*this, true);

  *cloned_vector = *this;

  return cloned_vector;
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
void EigenSparseVector<T>::swap (NumericVector<T> &other)
{
  EigenSparseVector<T>& v = libmesh_cast_ref<EigenSparseVector<T>&>(other);

  _vec.swap(v._vec);

  std::swap (this->_is_closed,      v._is_closed);
  std::swap (this->_is_initialized, v._is_initialized);
  std::swap (this->_type,           v._type);
}


} // namespace libMesh


#endif // #ifdef LIBMESH_HAVE_EIGEN
#endif // LIBMESH_EIGEN_SPARSE_VECTOR_H
