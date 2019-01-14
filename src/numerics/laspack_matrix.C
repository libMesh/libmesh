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



// C++ includes

// Local includes
#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_HAVE_LASPACK

#include "libmesh/laspack_matrix.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparsity_pattern.h"

namespace libMesh
{


//-----------------------------------------------------------------------
// LaspackMatrix members
template <typename T>
void LaspackMatrix<T>::update_sparsity_pattern (const SparsityPattern::Graph & sparsity_pattern)
{
  // clear data, start over
  this->clear ();

  // big trouble if this fails!
  libmesh_assert(this->_dof_map);

  const numeric_index_type n_rows = sparsity_pattern.size();

  // Initialize the _row_start data structure,
  // allocate storage for the _csr array
  {
    std::size_t size = 0;

    for (numeric_index_type row=0; row<n_rows; row++)
      size += sparsity_pattern[row].size();

    _csr.resize       (size);
    _row_start.reserve(n_rows + 1);
  }


  // Initialize the _csr data structure.
  {
    std::vector<numeric_index_type>::iterator position = _csr.begin();

    _row_start.push_back (position);

    for (numeric_index_type row=0; row<n_rows; row++)
      {
        // insert the row indices
        for (const auto & col : sparsity_pattern[row])
          {
            libmesh_assert (position != _csr.end());
            *position = col;
            ++position;
          }

        _row_start.push_back (position);
      }
  }


  // Initialize the matrix
  libmesh_assert (!this->initialized());
  this->init ();
  libmesh_assert (this->initialized());
  //libMesh::out << "n_rows=" << n_rows << std::endl;
  //libMesh::out << "m()=" << m() << std::endl;
  libmesh_assert_equal_to (n_rows, this->m());

  // Tell the matrix about its structure.  Initialize it
  // to zero.
  for (numeric_index_type i=0; i<n_rows; i++)
    {
      auto rs = _row_start[i];

      const numeric_index_type length = _row_start[i+1] - rs;

      Q_SetLen (&_QMat, i+1, length);

      for (numeric_index_type l=0; l<length; l++)
        {
          const numeric_index_type j = *(rs+l);

          // sanity check
          //libMesh::out << "m()=" << m() << std::endl;
          //libMesh::out << "(i,j,l) = (" << i
          //          << "," << j
          //          << "," << l
          //           << ")" << std::endl;
          //libMesh::out << "pos(i,j)=" << pos(i,j)
          //              << std::endl;
          libmesh_assert_equal_to (this->pos(i,j), l);
          Q_SetEntry (&_QMat, i+1, l, j+1, 0.);
        }
    }

  // That's it!
  //libmesh_here();
}



template <typename T>
void LaspackMatrix<T>::init (const numeric_index_type libmesh_dbg_var(m_in),
                             const numeric_index_type libmesh_dbg_var(n_in),
                             const numeric_index_type libmesh_dbg_var(m_l),
                             const numeric_index_type libmesh_dbg_var(n_l),
                             const numeric_index_type libmesh_dbg_var(nnz),
                             const numeric_index_type,
                             const numeric_index_type)
{
  // noz ignored...  only used for multiple processors!
  libmesh_assert_equal_to (m_in, m_l);
  libmesh_assert_equal_to (n_in, n_l);
  libmesh_assert_equal_to (m_in, n_in);
  libmesh_assert_greater (nnz, 0);

  libmesh_error_msg("ERROR: Only the init() member that uses the DofMap is implemented for Laspack matrices!");

  this->_is_initialized = true;
}



template <typename T>
void LaspackMatrix<T>::init ()
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
#ifndef NDEBUG
  // The following variables are only used for assertions,
  // so avoid declaring them when asserts are inactive.
  const numeric_index_type n_cols   = n_rows;
  const numeric_index_type n_l = this->_dof_map->n_dofs_on_processor(0);
  const numeric_index_type m_l = n_l;
#endif

  // Laspack Matrices only work for uniprocessor cases
  libmesh_assert_equal_to (n_rows, n_cols);
  libmesh_assert_equal_to (m_l, n_rows);
  libmesh_assert_equal_to (n_l, n_cols);

#ifndef NDEBUG
  // The following variables are only used for assertions,
  // so avoid declaring them when asserts are inactive.
  const std::vector<numeric_index_type> & n_nz = this->_dof_map->get_n_nz();
  const std::vector<numeric_index_type> & n_oz = this->_dof_map->get_n_oz();
#endif

  // Make sure the sparsity pattern isn't empty
  libmesh_assert_equal_to (n_nz.size(), n_l);
  libmesh_assert_equal_to (n_oz.size(), n_l);

  if (n_rows==0)
    return;

  Q_Constr(&_QMat, const_cast<char *>("Mat"), n_rows, _LPFalse, Rowws, Normal, _LPTrue);

  this->_is_initialized = true;

  libmesh_assert_equal_to (n_rows, this->m());
}



template <typename T>
void LaspackMatrix<T>::add_matrix(const DenseMatrix<T> & dm,
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
void LaspackMatrix<T>::get_diagonal (NumericVector<T> & /*dest*/) const
{
  libmesh_not_implemented();
}



template <typename T>
void LaspackMatrix<T>::get_transpose (SparseMatrix<T> & /*dest*/) const
{
  libmesh_not_implemented();
}



template <typename T>
LaspackMatrix<T>::LaspackMatrix (const Parallel::Communicator & comm) :
  SparseMatrix<T>(comm),
  _closed (false)
{
}



template <typename T>
LaspackMatrix<T>::~LaspackMatrix ()
{
  this->clear ();
}



template <typename T>
void LaspackMatrix<T>::clear ()
{
  if (this->initialized())
    {
      Q_Destr(&_QMat);
    }

  _csr.clear();
  _row_start.clear();
  _closed = false;
  this->_is_initialized = false;
}



template <typename T>
void LaspackMatrix<T>::zero ()
{
  const numeric_index_type n_rows = this->m();

  for (numeric_index_type row=0; row<n_rows; row++)
    {
      auto r_start = _row_start[row];

      const numeric_index_type len = (_row_start[row+1] - _row_start[row]);

      // Make sure we agree on the row length
      libmesh_assert_equal_to (len, Q_GetLen(&_QMat, row+1));

      for (numeric_index_type l=0; l<len; l++)
        {
          const numeric_index_type j = *(r_start + l);

          // Make sure the data structures are working
          libmesh_assert_equal_to ((j+1), Q_GetPos (&_QMat, row+1, l));

          Q_SetEntry (&_QMat, row+1, l, j+1, 0.);
        }
    }

  this->close();
}



template <typename T>
numeric_index_type LaspackMatrix<T>::m () const
{
  libmesh_assert (this->initialized());

  return static_cast<numeric_index_type>(Q_GetDim(const_cast<QMatrix*>(&_QMat)));
}



template <typename T>
numeric_index_type LaspackMatrix<T>::n () const
{
  libmesh_assert (this->initialized());

  return static_cast<numeric_index_type>(Q_GetDim(const_cast<QMatrix*>(&_QMat)));
}



template <typename T>
numeric_index_type LaspackMatrix<T>::row_start () const
{
  return 0;
}



template <typename T>
numeric_index_type LaspackMatrix<T>::row_stop () const
{
  return this->m();
}



template <typename T>
void LaspackMatrix<T>::set (const numeric_index_type i,
                            const numeric_index_type j,
                            const T value)
{
  libmesh_assert (this->initialized());
  libmesh_assert_less (i, this->m());
  libmesh_assert_less (j, this->n());

  const numeric_index_type position = this->pos(i,j);

  // Sanity check
  libmesh_assert_equal_to (*(_row_start[i]+position), j);
  libmesh_assert_equal_to ((j+1), Q_GetPos (&_QMat, i+1, position));

  Q_SetEntry (&_QMat, i+1, position, j+1, value);
}



template <typename T>
void LaspackMatrix<T>::add (const numeric_index_type i,
                            const numeric_index_type j,
                            const T value)
{
  libmesh_assert (this->initialized());
  libmesh_assert_less (i, this->m());
  libmesh_assert_less (j, this->n());

  const numeric_index_type position = this->pos(i,j);

  // Sanity check
  libmesh_assert_equal_to (*(_row_start[i]+position), j);

  Q_AddVal (&_QMat, i+1, position, value);
}



template <typename T>
void LaspackMatrix<T>::add_matrix(const DenseMatrix<T> & dm,
                                  const std::vector<numeric_index_type> & dof_indices)
{
  this->add_matrix (dm, dof_indices, dof_indices);
}



template <typename T>
void LaspackMatrix<T>::add (const T a_in, const SparseMatrix<T> & X_in)
{
  libmesh_assert (this->initialized());
  libmesh_assert_equal_to (this->m(), X_in.m());
  libmesh_assert_equal_to (this->n(), X_in.n());

  const LaspackMatrix<T> * X =
    cast_ptr<const LaspackMatrix<T> *> (&X_in);

  _LPNumber a = static_cast<_LPNumber> (a_in);

  libmesh_assert(X);

  // loops taken from LaspackMatrix<T>::zero ()

  const numeric_index_type n_rows = this->m();

  for (numeric_index_type row=0; row<n_rows; row++)
    {
      auto r_start = _row_start[row];

      const numeric_index_type len = (_row_start[row+1] - _row_start[row]);

      // Make sure we agree on the row length
      libmesh_assert_equal_to (len, Q_GetLen(&_QMat, row+1));
      // compare matrix sparsity structures
      libmesh_assert_equal_to (len, Q_GetLen(&(X->_QMat), row+1));


      for (numeric_index_type l=0; l<len; l++)
        {
          const numeric_index_type j = *(r_start + l);

          // Make sure the data structures are working
          libmesh_assert_equal_to ((j+1), Q_GetPos (&_QMat, row+1, l));

          const _LPNumber value = a * Q_GetEl(const_cast<QMatrix*>(&(X->_QMat)), row+1, j+1);
          Q_AddVal   (&_QMat, row+1, l, value);
        }
    }
}




template <typename T>
T LaspackMatrix<T>::operator () (const numeric_index_type i,
                                 const numeric_index_type j) const
{
  libmesh_assert (this->initialized());
  libmesh_assert_less (i, this->m());
  libmesh_assert_less (j, this->n());

  return Q_GetEl (const_cast<QMatrix*>(&_QMat), i+1, j+1);
}



template <typename T>
numeric_index_type LaspackMatrix<T>::pos (const numeric_index_type i,
                                          const numeric_index_type j) const
{
  libmesh_assert_less (i, this->m());
  libmesh_assert_less (j, this->n());
  libmesh_assert_less (i+1, _row_start.size());
  libmesh_assert (_row_start.back() == _csr.end());

  // note this requires the _csr to be sorted
  auto p = std::equal_range (_row_start[i], _row_start[i+1], j);

  // Make sure the row contains the element j
  libmesh_assert (p.first != p.second);

  // Make sure the values match
  libmesh_assert (*p.first == j);

  // Return the position in the compressed row
  return std::distance (_row_start[i], p.first);
}



template <typename T>
void LaspackMatrix<T>::close()
{
  libmesh_assert(this->initialized());

  this->_closed = true;

  // We've probably changed some entries so we need to tell LASPACK
  // that cached data is now invalid.
  *_QMat.DiagElAlloc = _LPFalse;
  *_QMat.ElSorted = _LPFalse;
  if (*_QMat.ILUExists)
    {
      *_QMat.ILUExists = _LPFalse;
      Q_Destr(_QMat.ILU);
    }
}



//------------------------------------------------------------------
// Explicit instantiations
template class LaspackMatrix<Number>;

} // namespace libMesh


#endif // #ifdef LIBMESH_HAVE_LASPACK
