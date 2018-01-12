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



// C++ includes

// Local includes
#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_HAVE_EIGEN

#include "libmesh/eigen_sparse_vector.h"
#include "libmesh/eigen_sparse_matrix.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparsity_pattern.h"

namespace libMesh
{


//-----------------------------------------------------------------------
// EigenSparseMatrix members
template <typename T>
void EigenSparseMatrix<T>::init (const numeric_index_type m_in,
                                 const numeric_index_type n_in,
                                 const numeric_index_type libmesh_dbg_var(m_l),
                                 const numeric_index_type libmesh_dbg_var(n_l),
                                 const numeric_index_type nnz,
                                 const numeric_index_type,
                                 const numeric_index_type)
{
  // noz ignored...  only used for multiple processors!
  libmesh_assert_equal_to (m_in, m_l);
  libmesh_assert_equal_to (n_in, n_l);
  libmesh_assert_equal_to (m_in, n_in);
  libmesh_assert_greater  (nnz, 0);

  _mat.resize(m_in, n_in);
  _mat.reserve(Eigen::Matrix<numeric_index_type, Eigen::Dynamic, 1>::Constant(m_in,nnz));

  this->_is_initialized = true;
}



template <typename T>
void EigenSparseMatrix<T>::init ()
{
  // Ignore calls on initialized objects
  if (this->initialized())
    return;

  // We need the DofMap for this!
  libmesh_assert(this->_dof_map);

  // Clear initialized matrices
  if (this->initialized())
    this->clear();

  const numeric_index_type n_rows   = this->_dof_map->n_dofs();
  const numeric_index_type n_cols   = n_rows;

#ifndef NDEBUG
  // The following variables are only used for assertions,
  // so avoid declaring them when asserts are inactive.
  const numeric_index_type n_l = this->_dof_map->n_dofs_on_processor(0);
  const numeric_index_type m_l = n_l;
#endif

  // Laspack Matrices only work for uniprocessor cases
  libmesh_assert_equal_to (n_rows, n_cols);
  libmesh_assert_equal_to (m_l, n_rows);
  libmesh_assert_equal_to (n_l, n_cols);

  const std::vector<numeric_index_type> & n_nz = this->_dof_map->get_n_nz();

#ifndef NDEBUG
  // The following variables are only used for assertions,
  // so avoid declaring them when asserts are inactive.
  const std::vector<numeric_index_type> & n_oz = this->_dof_map->get_n_oz();
#endif

  // Make sure the sparsity pattern isn't empty
  libmesh_assert_equal_to (n_nz.size(), n_l);
  libmesh_assert_equal_to (n_oz.size(), n_l);

  if (n_rows==0)
    {
      _mat.resize(0,0);
      return;
    }

  _mat.resize(n_rows,n_cols);
  _mat.reserve(n_nz);

  this->_is_initialized = true;

  libmesh_assert_equal_to (n_rows, this->m());
  libmesh_assert_equal_to (n_cols, this->n());
}



template <typename T>
void EigenSparseMatrix<T>::add_matrix(const DenseMatrix<T> & dm,
                                      const std::vector<numeric_index_type> & rows,
                                      const std::vector<numeric_index_type> & cols)

{
  libmesh_assert (this->initialized());
  unsigned int n_rows = cast_int<unsigned int>(rows.size());
  unsigned int n_cols = cast_int<unsigned int>(cols.size());
  libmesh_assert_equal_to (dm.m(), n_rows);
  libmesh_assert_equal_to (dm.n(), n_cols);


  for (unsigned int i=0; i<n_rows; i++)
    for (unsigned int j=0; j<n_cols; j++)
      this->add(rows[i],cols[j],dm(i,j));
}



template <typename T>
void EigenSparseMatrix<T>::get_diagonal (NumericVector<T> & dest_in) const
{
  EigenSparseVector<T> & dest = cast_ref<EigenSparseVector<T> &>(dest_in);

  dest._vec = _mat.diagonal();
}



template <typename T>
void EigenSparseMatrix<T>::get_transpose (SparseMatrix<T> & dest_in) const
{
  EigenSparseMatrix<T> & dest = cast_ref<EigenSparseMatrix<T> &>(dest_in);

  dest._mat = _mat.transpose();
}



template <typename T>
EigenSparseMatrix<T>::EigenSparseMatrix (const Parallel::Communicator & comm_in) :
  SparseMatrix<T>(comm_in),
  _closed (false)
{
}



template <typename T>
EigenSparseMatrix<T>::~EigenSparseMatrix ()
{
  this->clear ();
}



template <typename T>
void EigenSparseMatrix<T>::clear ()
{
  _mat.resize(0,0);

  _closed = false;
  this->_is_initialized = false;
}



template <typename T>
void EigenSparseMatrix<T>::zero ()
{
  _mat.setZero();
}



template <typename T>
numeric_index_type EigenSparseMatrix<T>::m () const
{
  libmesh_assert (this->initialized());

  return _mat.rows();
}



template <typename T>
numeric_index_type EigenSparseMatrix<T>::n () const
{
  libmesh_assert (this->initialized());

  return _mat.cols();
}



template <typename T>
numeric_index_type EigenSparseMatrix<T>::row_start () const
{
  return 0;
}



template <typename T>
numeric_index_type EigenSparseMatrix<T>::row_stop () const
{
  return this->m();
}



template <typename T>
void EigenSparseMatrix<T>::set (const numeric_index_type i,
                                const numeric_index_type j,
                                const T value)
{
  libmesh_assert (this->initialized());
  libmesh_assert_less (i, this->m());
  libmesh_assert_less (j, this->n());

  _mat.coeffRef(i,j) = value;
}



template <typename T>
void EigenSparseMatrix<T>::add (const numeric_index_type i,
                                const numeric_index_type j,
                                const T value)
{
  libmesh_assert (this->initialized());
  libmesh_assert_less (i, this->m());
  libmesh_assert_less (j, this->n());

  _mat.coeffRef(i,j) += value;
}



template <typename T>
void EigenSparseMatrix<T>::add_matrix(const DenseMatrix<T> & dm,
                                      const std::vector<numeric_index_type> & dof_indices)
{
  this->add_matrix (dm, dof_indices, dof_indices);
}



template <typename T>
void EigenSparseMatrix<T>::add (const T a_in, SparseMatrix<T> & X_in)
{
  libmesh_assert (this->initialized());
  libmesh_assert_equal_to (this->m(), X_in.m());
  libmesh_assert_equal_to (this->n(), X_in.n());

  EigenSparseMatrix<T> & X = cast_ref<EigenSparseMatrix<T> &> (X_in);

  _mat += X._mat*a_in;
}




template <typename T>
T EigenSparseMatrix<T>::operator () (const numeric_index_type i,
                                     const numeric_index_type j) const
{
  libmesh_assert (this->initialized());
  libmesh_assert_less (i, this->m());
  libmesh_assert_less (j, this->n());

  return _mat.coeff(i,j);
}



template <typename T>
Real EigenSparseMatrix<T>::l1_norm () const
{
  // There does not seem to be a straightforward way to iterate over
  // the columns of an EigenSparseMatrix.  So we use some extra
  // storage and keep track of the column sums while going over the
  // row entries...
  std::vector<Real> abs_col_sums(this->n());

  // For a row-major Eigen SparseMatrix like we're using, the
  // InnerIterator iterates over the non-zero entries of rows.
  for (unsigned row=0; row<this->m(); ++row)
    {
      EigenSM::InnerIterator it(_mat, row);
      for (; it; ++it)
        abs_col_sums[it.col()] += std::abs(it.value());
    }

  return *(std::max_element(abs_col_sums.begin(), abs_col_sums.end()));
}



template <typename T>
Real EigenSparseMatrix<T>::linfty_norm () const
{
  Real max_abs_row_sum = 0.;

  // For a row-major Eigen SparseMatrix like we're using, the
  // InnerIterator iterates over the non-zero entries of rows.
  for (unsigned row=0; row<this->m(); ++row)
    {
      Real current_abs_row_sum = 0.;
      EigenSM::InnerIterator it(_mat, row);
      for (; it; ++it)
        current_abs_row_sum += std::abs(it.value());

      max_abs_row_sum = std::max(max_abs_row_sum, current_abs_row_sum);
    }

  return max_abs_row_sum;
}


//------------------------------------------------------------------
// Explicit instantiations
template class EigenSparseMatrix<Number>;

} // namespace libMesh


#endif // #ifdef LIBMESH_HAVE_EIGEN
