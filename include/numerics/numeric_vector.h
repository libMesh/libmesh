// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_NUMERIC_VECTOR_H
#define LIBMESH_NUMERIC_VECTOR_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/auto_ptr.h" // deprecated
#include "libmesh/enum_parallel_type.h"
#include "libmesh/id_types.h"
#include "libmesh/reference_counted_object.h"
#include "libmesh/libmesh.h"
#include "libmesh/parallel_object.h"
#include "libmesh/dense_subvector.h"
#include "libmesh/dense_vector.h"

#ifdef LIBMESH_FORWARD_DECLARE_ENUMS
namespace libMesh
{
enum SolverPackage : int;
}
#else
#include "libmesh/enum_solver_package.h"
#endif

// C++ includes
#include <cstddef>
#include <set>
#include <vector>
#include <memory>

namespace libMesh
{


// forward declarations
template <typename T> class NumericVector;
template <typename T> class DenseVector;
template <typename T> class DenseSubVector;
template <typename T> class SparseMatrix;
template <typename T> class ShellMatrix;

/**
 * \brief Provides a uniform interface to vector storage schemes for different
 * linear algebra libraries.
 *
 * \note This class is the abstract base class for different implementations
 * of numeric vectors.  Most of the time you should use a System object to
 * create numeric vectors.  If this is not desired, you can instantiate one of
 * the derived classes (PetscVector, EigenSparseVector, etc.) or use the
 * NumericVector::build method. When creating the vector yourself, make sure
 * that you initialize the vector properly (NumericVector::init).
 *
 * \author Benjamin S. Kirk
 * \date 2003
 */
template <typename T>
class NumericVector : public ReferenceCountedObject<NumericVector<T>>,
                      public ParallelObject
{
public:

  /**
   * Dummy-Constructor. Dimension=0
   */
  explicit
  NumericVector (const Parallel::Communicator & comm_in,
                 const ParallelType ptype = AUTOMATIC);

  /**
   * Constructor. Set dimension to \p n and initialize all elements with zero.
   */
  explicit
  NumericVector (const Parallel::Communicator & comm_in,
                 const numeric_index_type n,
                 const ParallelType ptype = AUTOMATIC);

  /**
   * Constructor. Set local dimension to \p n_local, the global dimension
   * to \p n, and initialize all elements with zero.
   */
  NumericVector (const Parallel::Communicator & comm_in,
                 const numeric_index_type n,
                 const numeric_index_type n_local,
                 const ParallelType ptype = AUTOMATIC);

  /**
   * Constructor. Set local dimension to \p n_local, the global
   * dimension to \p n, but additionally reserve memory for the
   * indices specified by the \p ghost argument.
   */
  NumericVector (const Parallel::Communicator & comm_in,
                 const numeric_index_type N,
                 const numeric_index_type n_local,
                 const std::vector<numeric_index_type> & ghost,
                 const ParallelType ptype = AUTOMATIC);

  /**
   * This _looks_ like a copy assignment operator, but note that,
   * unlike normal copy assignment operators, it is pure virtual. This
   * function should be overridden in derived classes so that they can
   * be copied correctly via references to the base class. This design
   * usually isn't a good idea in general, but in this context it
   * works because we usually don't have a mix of different kinds of
   * NumericVectors active in the library at a single time.
   *
   * \returns A reference to *this as the base type.
   */
  virtual NumericVector<T> & operator= (const NumericVector<T> & v) = 0;

  /**
   * The 5 special functions can be defaulted for this class, as it
   * does not manage any memory itself.
   */
  NumericVector (NumericVector &&) = default;
  NumericVector (const NumericVector &) = default;
  NumericVector & operator= (NumericVector &&) = default;

  /**
   * Builds a \p NumericVector on the processors in communicator
   * \p comm using the linear solver package specified by
   * \p solver_package
   */
  static std::unique_ptr<NumericVector<T>>
  build(const Parallel::Communicator & comm,
        const SolverPackage solver_package = libMesh::default_solver_package());

  /**
   * \returns \p true if the vector has been initialized,
   * false otherwise.
   */
  virtual bool initialized() const { return _is_initialized; }

  /**
   * \returns The type (SERIAL, PARALLEL, GHOSTED) of the vector.
   */
  ParallelType type() const { return _type; }

  /**
   * \returns The type (SERIAL, PARALLEL, GHOSTED) of the vector.
   */
  ParallelType & type() { return _type; }

  /**
   * \returns \p true if the vector is closed and ready for
   * computation, false otherwise.
   */
  virtual bool closed() const { return _is_closed; }

  /**
   * Calls the NumericVector's internal assembly routines, ensuring
   * that the values are consistent across processors.
   */
  virtual void close () = 0;

  /**
   * Restores the \p NumericVector<T> to a pristine state.
   */
  virtual void clear ();

  /**
   * Set all entries to zero. Equivalent to \p v = 0, but more obvious and
   * faster.
   */
  virtual void zero () = 0;

  /**
   * \returns A smart pointer to a copy of this vector with the same
   * type, size, and partitioning, but with all zero entries.
   *
   * \note This must be overridden in the derived classes.
   */
  virtual std::unique_ptr<NumericVector<T>> zero_clone () const = 0;

  /**
   * \returns A copy of this vector wrapped in a smart pointer.
   *
   * \note This must be overridden in the derived classes.
   */
  virtual std::unique_ptr<NumericVector<T>> clone () const = 0;

  /**
   * Change the dimension of the vector to \p n. The reserved memory
   * for this vector remains unchanged if possible.  If \p n==0, all
   * memory is freed. Therefore, if you want to resize the vector and
   * release the memory not needed, you have to first call \p init(0)
   * and then \p init(n). This behaviour is analogous to that of the
   * STL containers.
   *
   * On \p fast==false, the vector is filled by zeros.
   */
  virtual void init (const numeric_index_type n,
                     const numeric_index_type n_local,
                     const bool fast = false,
                     const ParallelType ptype = AUTOMATIC) = 0;

  /**
   * Call \p init() with n_local = N.
   */
  virtual void init (const numeric_index_type n,
                     const bool fast = false,
                     const ParallelType ptype = AUTOMATIC) = 0;

  /**
   * Create a vector that holds tha local indices plus those specified
   * in the \p ghost argument.
   */
  virtual void init (const numeric_index_type n,
                     const numeric_index_type n_local,
                     const std::vector<numeric_index_type> & ghost,
                     const bool fast = false,
                     const ParallelType ptype = AUTOMATIC) = 0;

  /**
   * Creates a vector that has the same dimension and storage type as
   * \p other, including ghost dofs.
   */
  virtual void init (const NumericVector<T> & other,
                     const bool fast = false) = 0;

  /**
   * Sets all entries of the vector to the value \p s.
   *
   * \returns A reference to *this.
   */
  virtual NumericVector<T> & operator= (const T s) = 0;

  /**
   * Sets (*this)(i) = v(i) for each entry of the vector.
   *
   * \returns A reference to *this as the base type.
   */
  virtual NumericVector<T> & operator= (const std::vector<T> & v) = 0;

  /**
   * \returns The minimum entry in the vector, or the minimum real
   * part in the case of complex numbers.
   */
  virtual Real min () const = 0;

  /**
   * \returns The maximum entry in the vector, or the maximum real
   * part in the case of complex numbers.
   */
  virtual Real max () const = 0;

  /**
   * \returns The sum of all values in the vector.
   */
  virtual T sum() const = 0;

  /**
   * \returns The \f$ \ell_1 \f$-norm of the vector, i.e. the sum of the
   * absolute values of the entries.
   */
  virtual Real l1_norm () const = 0;

  /**
   * \returns The \f$ \ell_2 \f$-norm of the vector, i.e. the square root
   * of the sum of the squares of the entries.
   */
  virtual Real l2_norm () const = 0;

  /**
   * \returns The \f$ \ell_{\infty} \f$-norm of the vector, i.e. the maximum
   * absolute value of the entries of the vector.
   */
  virtual Real linfty_norm () const = 0;

  /**
   * \returns The \f$ \ell_1 \f$-norm of the vector, i.e. the sum of the
   * absolute values for the specified entries in the vector.
   *
   * \note The indices must necessarily live on this processor.
   */
  virtual Real subset_l1_norm (const std::set<numeric_index_type> & indices) const;

  /**
   * \returns The \f$ \ell_2 \f$-norm of the vector, i.e. the square root
   * of the sum of the squares of the elements for the specified
   * entries in the vector.
   *
   * \note The indices must necessarily live on this processor.
   */
  virtual Real subset_l2_norm (const std::set<numeric_index_type> & indices) const;

  /**
   * \returns The maximum absolute value of the specified entries of
   * this vector, which is the \f$ \ell_{\infty} \f$-norm of a vector.
   *
   * \note The indices must necessarily live on this processor.
   */
  virtual Real subset_linfty_norm (const std::set<numeric_index_type> & indices) const;

  /**
   * \returns The size of the vector.
   */
  virtual numeric_index_type size () const = 0;

  /**
   * \returns The local size of the vector, i.e. \p index_stop - \p index_start.
   */
  virtual numeric_index_type local_size() const = 0;

  /**
   * \returns The index of the first vector element actually stored on
   * this processor.
   *
   * \note The minimum for this index is \p 0.
   */
  virtual numeric_index_type first_local_index() const = 0;

  /**
   * \returns The index+1 of the last vector element actually stored
   * on this processor.
   *
   * \note The maximum for this index is \p size().
   */
  virtual numeric_index_type last_local_index() const = 0;

  /**
   * \returns A copy of the ith entry of the vector.
   */
  virtual T operator() (const numeric_index_type i) const = 0;

  /**
   * \returns (*this)(i).
   */
  virtual T el(const numeric_index_type i) const { return (*this)(i); }

  /**
   * Access multiple components at once.  \p values will *not* be
   * reallocated; it should already have enough space.  The default
   * implementation calls \p operator() for each index, but some
   * implementations may supply faster methods here.
   */
  virtual void get(const std::vector<numeric_index_type> & index,
                   T * values) const;

  /**
   * Access multiple components at once.  \p values will be resized,
   * if necessary, and filled.  The default implementation calls \p
   * operator() for each index, but some implementations may supply
   * faster methods here.
   */
  void get(const std::vector<numeric_index_type> & index,
           std::vector<T> & values) const;

  /**
   * Adds \p v to *this,
   * \f$ \vec{u} \leftarrow \vec{u} + \vec{v} \f$.
   * Equivalent to \p u.add(1, v).
   *
   * \returns A reference to *this.
   */
  virtual NumericVector<T> & operator += (const NumericVector<T> & v) = 0;

  /**
   * Subtracts \p v from *this,
   * \f$ \vec{u} \leftarrow \vec{u} - \vec{v} \f$.
   * Equivalent to \p u.add(-1, v).
   *
   * \returns A reference to *this.
   */
  virtual NumericVector<T> & operator -= (const NumericVector<T> & v) = 0;

  /**
   * Scales the vector by \p a,
   * \f$ \vec{u} \leftarrow a\vec{u} \f$.
   * Equivalent to \p u.scale(a)
   *
   * \returns A reference to *this.
   */
  NumericVector<T> & operator *= (const T a) { this->scale(a); return *this; }

  /**
   * Scales the vector by \p 1/a,
   * \f$ \vec{u} \leftarrow \frac{1}{a}\vec{u} \f$.
   * Equivalent to \p u.scale(1./a)
   *
   * \returns A reference to *this.
   */
  NumericVector<T> & operator /= (const T a) { this->scale(1./a); return *this; }

  /**
   * Computes the pointwise division of this vector's entries by another's,
   * \f$ u_i \leftarrow \frac{u_i}{v_i} \, \forall i\f$
   *
   * \returns A reference to *this.
   */
  virtual NumericVector<T> & operator /= (const NumericVector<T> & /*v*/) = 0;

  /**
   * Computes the pointwise reciprocal,
   * \f$ u_i \leftarrow \frac{1}{u_i} \, \forall i\f$
   */
  virtual void reciprocal() = 0;

  /**
   * Negates the imaginary component of each entry in the vector.
   */
  virtual void conjugate() = 0;

  /**
   * Sets v(i) = \p value.
   */
  virtual void set (const numeric_index_type i, const T value) = 0;

  /**
   * Adds \p value to each entry of the vector.
   */
  virtual void add (const numeric_index_type i, const T value) = 0;

  /**
   * Adds \p s to each entry of the vector,
   * \f$ u_i \leftarrow u_i + s \f$
   */
  virtual void add (const T s) = 0;

  /**
   * Adds \p v to \p this,
   * \f$ \vec{u} \leftarrow \vec{u} + \vec{v} \f$.
   * Equivalent to calling \p operator+=().
   */
  virtual void add (const NumericVector<T> & v) = 0;

  /**
   * Vector addition with a scalar multiple,
   * \f$ \vec{u} \leftarrow \vec{u} + a\vec{v} \f$.
   * Equivalent to calling \p operator+=().
   */
  virtual void add (const T a, const NumericVector<T> & v) = 0;

  /**
   * Computes \f$ \vec{u} \leftarrow \vec{u} + \vec{v} \f$,
   * where \p v is a pointer and each \p dof_indices[i] specifies where
   * to add value \p v[i].  This should be overridden in subclasses
   * for efficiency.
   */
  virtual void add_vector (const T * v,
                           const std::vector<numeric_index_type> & dof_indices);

  /**
   * Computes \f$ \vec{u} \leftarrow \vec{u} + \vec{v} \f$,
   * where \p v is a std::vector and each \p dof_indices[i] specifies
   * where to add value \p v[i].
   */
  void add_vector (const std::vector<T> & v,
                   const std::vector<numeric_index_type> & dof_indices);

  /**
   * Computes \f$ \vec{u} \leftarrow \vec{u} + \vec{v} \f$,
   * where \p v is a NumericVector and each \p dof_indices[i]
   * specifies where to add value \p v(i).
   */
  virtual void add_vector (const NumericVector<T> & v,
                           const std::vector<numeric_index_type> & dof_indices);

  /**
   * Computes \f$ \vec{u} \leftarrow \vec{u} + \vec{v} \f$,
   * where \p v is a DenseVector and each \p dof_indices[i] specifies
   * where to add value \p v(i).
   */
  void add_vector (const DenseVector<T> & v,
                   const std::vector<numeric_index_type> & dof_indices);

  /**
   * Computes \f$ \vec{u} \leftarrow \vec{u} + A \vec{v} \f$,
   * i.e. adds the product of a \p SparseMatrix \p A and a \p
   * NumericVector \p v to \p this.
   */
  virtual void add_vector (const NumericVector<T> & v,
                           const SparseMatrix<T> & A) = 0;

  /**
   * Computes \f$ \vec{u} \leftarrow \vec{u} + A \vec{v} \f$,
   * i.e. adds the product of a \p ShellMatrix \p A and a \p
   * NumericVector \p v to \p this.
   */
  void add_vector (const NumericVector<T> & v,
                   const ShellMatrix<T> & A);

  /**
   * Computes \f$ \vec{u} \leftarrow \vec{u} + A^T \vec{v} \f$,
   * i.e. adds the product of the transpose of a \p SparseMatrix \p A
   * and a \p NumericVector \p v to \p this.
   */
  virtual void add_vector_transpose (const NumericVector<T> & v,
                                     const SparseMatrix<T> & A) = 0;

  /**
   * Inserts the entries of \p v in *this at the locations specified by \p v.
   */
  virtual void insert (const T * v,
                       const std::vector<numeric_index_type> & dof_indices);

  /**
   * Inserts the entries of \p v in *this at the locations specified by \p v.
   */
  void insert (const std::vector<T> & v,
               const std::vector<numeric_index_type> & dof_indices);

  /**
   * Inserts the entries of \p v in *this at the locations specified by \p v.
   */
  virtual void insert (const NumericVector<T> & v,
                       const std::vector<numeric_index_type> & dof_indices);

  /**
   * Inserts the entries of \p v in *this at the locations specified by \p v.
   */
  void insert (const DenseVector<T> & v,
               const std::vector<numeric_index_type> & dof_indices);

  /**
   * Inserts the entries of \p v in *this at the locations specified by \p v.
   */
  void insert (const DenseSubVector<T> & v,
               const std::vector<numeric_index_type> & dof_indices);

  /**
   * Scale each element of the vector by the given \p factor.
   */
  virtual void scale (const T factor) = 0;

  /**
   * Sets \f$ u_i \leftarrow |u_i| \f$ for each entry in the vector.
   */
  virtual void abs() = 0;

  /**
   * \returns \f$ \vec{u} \cdot \vec{v} \f$, the dot product of
   * (*this) with the vector \p v.
   *
   * Uses the complex-conjugate of \p v in the complex-valued case.
   */
  virtual T dot(const NumericVector<T> & v) const = 0;

  /**
   * Creates a copy of the global vector in the local vector \p
   * v_local.
   */
  virtual void localize (std::vector<T> & v_local) const = 0;

  /**
   * Same, but fills a \p NumericVector<T> instead of a \p
   * std::vector.
   */
  virtual void localize (NumericVector<T> & v_local) const = 0;

  /**
   * Creates a local vector \p v_local containing only information
   * relevant to this processor, as defined by the \p send_list.
   */
  virtual void localize (NumericVector<T> & v_local,
                         const std::vector<numeric_index_type> & send_list) const = 0;

  /**
   * Fill in the local std::vector "v_local" with the global indices
   * given in "indices".
   *
   * \note The indices can be different on every processor, and the
   * same index can be localized to more than one processor.  The
   * resulting v_local can be shorter than the original, and the
   * entries will be in the order specified by indices.
   *
   * Example:
   * \verbatim
   *   On 4 procs *this = {a, b, c, d, e, f, g, h, i} is a parallel vector.
   *   On each proc, the indices arrays are set up as:
   *   proc0, indices = {1,2,4,5}
   *   proc1, indices = {2,5,6,8}
   *   proc2, indices = {2,3,6,7}
   *   proc3, indices = {0,1,2,3}
   *
   *   After calling this version of localize, the v_local vectors are:
   *   proc0, v_local = {b,c,e,f}
   *   proc1, v_local = {c,f,g,i}
   *   proc2, v_local = {c,d,g,h}
   *   proc3, v_local = {a,b,c,d}
   * \endverbatim
   *
   * This function is useful in parallel I/O routines, when you have a
   * parallel vector of solution values which you want to write a
   * subset of.
   */
  virtual void localize (std::vector<T> & v_local,
                         const std::vector<numeric_index_type> & indices) const = 0;

  /**
   * Updates a local vector with selected values from neighboring
   * processors, as defined by \p send_list.
   */
  virtual void localize (const numeric_index_type first_local_idx,
                         const numeric_index_type last_local_idx,
                         const std::vector<numeric_index_type> & send_list) = 0;

  /**
   * Creates a local copy of the global vector in \p v_local only on
   * processor \p proc_id.  By default the data is sent to processor
   * 0.  This method is useful for outputting data from one processor.
   */
  virtual void localize_to_one (std::vector<T> & v_local,
                                const processor_id_type proc_id=0) const = 0;

  /**
   * \returns \p -1 when \p this is equivalent to \p other_vector
   * (up to the given \p threshold), or the first index where
   * \p abs(a[i]-b[i]) exceeds the threshold.
   */
  virtual int compare (const NumericVector<T> & other_vector,
                       const Real threshold = TOLERANCE) const;

  /**
   * \returns \p -1 when \p this is equivalent to \p other_vector, (up
   * to the given local relative \p threshold), or the first index
   * where \p abs(a[i]-b[i])/max(a[i],b[i]) exceeds the threshold.
   */
  virtual int local_relative_compare (const NumericVector<T> & other_vector,
                                      const Real threshold = TOLERANCE) const;

  /**
   * \returns \p -1 when \p this is equivalent to \p other_vector (up
   * to the given global relative \p threshold), or the first index
   * where \p abs(a[i]-b[i])/max_j(a[j],b[j]) exceeds the threshold.
   */
  virtual int global_relative_compare (const NumericVector<T> & other_vector,
                                       const Real threshold = TOLERANCE) const;

  /**
   * Computes \f$ u_i \leftarrow u_i v_i \f$ (summation not implied)
   * i.e. the pointwise (component-wise) product of \p vec1 and
   * \p vec2, and stores the result in \p *this.
   */
  virtual void pointwise_mult (const NumericVector<T> & vec1,
                               const NumericVector<T> & vec2) = 0;

  /**
   * Prints the local contents of the vector, by default to
   * libMesh::out
   */
  virtual void print(std::ostream & os=libMesh::out) const;

  /**
   * Prints the global contents of the vector, by default to
   * libMesh::out
   */
  virtual void print_global(std::ostream & os=libMesh::out) const;

  /**
   * Same as above but allows you to use stream syntax.
   */
  friend std::ostream & operator << (std::ostream & os, const NumericVector<T> & v)
  {
    v.print_global(os);
    return os;
  }

  /**
   * Print the contents of the vector in Matlab's sparse matrix
   * format. Optionally prints the vector to the file named \p name.
   * If \p name is not specified it is dumped to the screen.
   */
  virtual void print_matlab(const std::string & /*name*/ = "") const
  {
    libmesh_not_implemented();
  }

  /**
   * Fills in \p subvector from this vector using the indices in \p
   * rows.  Similar to the \p create_submatrix() routine for the
   * SparseMatrix class, it is currently only implemented for
   * PetscVectors.
   */
  virtual void create_subvector(NumericVector<T> & ,
                                const std::vector<numeric_index_type> &) const
  {
    libmesh_not_implemented();
  }

  /**
   * Swaps the contents of this with \p v.  There should be enough
   * indirection in subclasses to make this an O(1) header-swap
   * operation.
   */
  virtual void swap (NumericVector<T> & v);

protected:

  /**
   * Flag which tracks whether the vector's values are consistent on
   * all processors after insertion or addition of values has occurred
   * on some or all processors.
   */
  bool _is_closed;

  /**
   * \p true once init() has been called.
   */
  bool _is_initialized;

  /**
   * Type of vector.
   */
  ParallelType _type;
};


/*----------------------- Inline functions ----------------------------------*/



template <typename T>
inline
NumericVector<T>::NumericVector (const Parallel::Communicator & comm_in,
                                 const ParallelType ptype) :
  ParallelObject(comm_in),
  _is_closed(false),
  _is_initialized(false),
  _type(ptype)
{
}



template <typename T>
inline
NumericVector<T>::NumericVector (const Parallel::Communicator & comm_in,
                                 const numeric_index_type /*n*/,
                                 const ParallelType ptype) :
  ParallelObject(comm_in),
  _is_closed(false),
  _is_initialized(false),
  _type(ptype)
{
  libmesh_not_implemented(); // Abstract base class!
  // init(n, n, false, ptype);
}



template <typename T>
inline
NumericVector<T>::NumericVector (const Parallel::Communicator & comm_in,
                                 const numeric_index_type /*n*/,
                                 const numeric_index_type /*n_local*/,
                                 const ParallelType ptype) :
  ParallelObject(comm_in),
  _is_closed(false),
  _is_initialized(false),
  _type(ptype)
{
  libmesh_not_implemented(); // Abstract base class!
  // init(n, n_local, false, ptype);
}



template <typename T>
inline
NumericVector<T>::NumericVector (const Parallel::Communicator & comm_in,
                                 const numeric_index_type /*n*/,
                                 const numeric_index_type /*n_local*/,
                                 const std::vector<numeric_index_type> & /*ghost*/,
                                 const ParallelType ptype) :
  ParallelObject(comm_in),
  _is_closed(false),
  _is_initialized(false),
  _type(ptype)
{
  libmesh_not_implemented(); // Abstract base class!
  // init(n, n_local, ghost, false, ptype);
}



template <typename T>
inline
void NumericVector<T>::clear ()
{
  _is_closed      = false;
  _is_initialized = false;
}



template <typename T>
inline
void NumericVector<T>::get(const std::vector<numeric_index_type> & index,
                           T * values) const
{
  const std::size_t num = index.size();
  for (std::size_t i=0; i<num; i++)
    {
      values[i] = (*this)(index[i]);
    }
}



template <typename T>
inline
void NumericVector<T>::get(const std::vector<numeric_index_type> & index,
                           std::vector<T> & values) const
{
  const std::size_t num = index.size();
  values.resize(num);
  if (!num)
    return;

  this->get(index, values.data());
}



template <typename T>
inline
void NumericVector<T>::add_vector(const std::vector<T> & v,
                                  const std::vector<numeric_index_type> & dof_indices)
{
  libmesh_assert(v.size() == dof_indices.size());
  if (!v.empty())
    this->add_vector(v.data(), dof_indices);
}



template <typename T>
inline
void NumericVector<T>::add_vector(const DenseVector<T> & v,
                                  const std::vector<numeric_index_type> & dof_indices)
{
  libmesh_assert(v.size() == dof_indices.size());
  if (!v.empty())
    this->add_vector(&v(0), dof_indices);
}



template <typename T>
inline
void NumericVector<T>::insert(const std::vector<T> & v,
                              const std::vector<numeric_index_type> & dof_indices)
{
  libmesh_assert(v.size() == dof_indices.size());
  if (!v.empty())
    this->insert(v.data(), dof_indices);
}



template <typename T>
inline
void NumericVector<T>::insert(const DenseVector<T> & v,
                              const std::vector<numeric_index_type> & dof_indices)
{
  libmesh_assert(v.size() == dof_indices.size());
  if (!v.empty())
    this->insert(&v(0), dof_indices);
}



template <typename T>
inline
void NumericVector<T>::insert(const DenseSubVector<T> & v,
                              const std::vector<numeric_index_type> & dof_indices)
{
  libmesh_assert(v.size() == dof_indices.size());
  if (!v.empty())
    this->insert(&v(0), dof_indices);
}



// Full specialization of the print() member for complex
// variables.  This must precede the non-specialized
// version, at least according to icc v7.1
template <>
inline
void NumericVector<Complex>::print(std::ostream & os) const
{
  libmesh_assert (this->initialized());
  os << "Size\tglobal =  " << this->size()
     << "\t\tlocal =  " << this->local_size() << std::endl;

  // std::complex<>::operator<<() is defined, but use this form
  os << "#\tReal part\t\tImaginary part" << std::endl;
  for (numeric_index_type i=this->first_local_index(); i<this->last_local_index(); i++)
    os << i << "\t"
       << (*this)(i).real() << "\t\t"
       << (*this)(i).imag() << std::endl;
}



template <typename T>
inline
void NumericVector<T>::print(std::ostream & os) const
{
  libmesh_assert (this->initialized());
  os << "Size\tglobal =  " << this->size()
     << "\t\tlocal =  " << this->local_size() << std::endl;

  os << "#\tValue" << std::endl;
  for (numeric_index_type i=this->first_local_index(); i<this->last_local_index(); i++)
    os << i << "\t" << (*this)(i) << std::endl;
}



template <>
inline
void NumericVector<Complex>::print_global(std::ostream & os) const
{
  libmesh_assert (this->initialized());

  std::vector<Complex> v(this->size());
  this->localize(v);

  // Right now we only want one copy of the output
  if (this->processor_id())
    return;

  os << "Size\tglobal =  " << this->size() << std::endl;
  os << "#\tReal part\t\tImaginary part" << std::endl;
  for (numeric_index_type i=0; i!=v.size(); i++)
    os << i << "\t"
       << v[i].real() << "\t\t"
       << v[i].imag() << std::endl;
}


template <typename T>
inline
void NumericVector<T>::print_global(std::ostream & os) const
{
  libmesh_assert (this->initialized());

  std::vector<T> v(this->size());
  this->localize(v);

  // Right now we only want one copy of the output
  if (this->processor_id())
    return;

  os << "Size\tglobal =  " << this->size() << std::endl;
  os << "#\tValue" << std::endl;
  for (numeric_index_type i=0; i!=v.size(); i++)
    os << i << "\t" << v[i] << std::endl;
}



template <typename T>
inline
void  NumericVector<T>::swap (NumericVector<T> & v)
{
  std::swap(_is_closed, v._is_closed);
  std::swap(_is_initialized, v._is_initialized);
  std::swap(_type, v._type);
}


} // namespace libMesh


#endif  // LIBMESH_NUMERIC_VECTOR_H
