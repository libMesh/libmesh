// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_TRILINOS_HAVE_EPETRA

// Local includes
#include "libmesh/trilinos_epetra_matrix.h"
#include "libmesh/trilinos_epetra_vector.h"
#include "libmesh/dof_map.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/parallel.h"
#include "libmesh/sparsity_pattern.h"
#include "libmesh/int_range.h"

namespace libMesh
{



//-----------------------------------------------------------------------
//EpetraMatrix members

template <typename T>
void EpetraMatrix<T>::update_sparsity_pattern (const SparsityPattern::Graph & sparsity_pattern)
{
  // clear data, start over
  this->clear ();

  // big trouble if this fails!
  libmesh_assert(this->_dof_map);

  const numeric_index_type n_rows = cast_int<numeric_index_type>
    (sparsity_pattern.size());

  const numeric_index_type m   = this->_dof_map->n_dofs();
  const numeric_index_type n   = m;
  const numeric_index_type n_l = this->_dof_map->n_local_dofs();
  const numeric_index_type m_l = n_l;

  // error checking
#ifndef NDEBUG
  {
    libmesh_assert_equal_to (n, m);
    libmesh_assert_equal_to (n_l, m_l);

    numeric_index_type
      summed_m_l = m_l,
      summed_n_l = n_l;

    this->comm().sum (summed_m_l);
    this->comm().sum (summed_n_l);

    libmesh_assert_equal_to (m, summed_m_l);
    libmesh_assert_equal_to (n, summed_n_l);
  }
#endif

  // build a map defining the data distribution
  _map = new Epetra_Map (static_cast<int>(m),
                         m_l,
                         0,
                         Epetra_MpiComm (this->comm().get()));

  libmesh_assert_equal_to (static_cast<numeric_index_type>(_map->NumGlobalPoints()), m);
  libmesh_assert_equal_to (static_cast<numeric_index_type>(_map->MaxAllGID()+1), m);

  const std::vector<numeric_index_type> & n_nz = this->_sp->get_n_nz();
  const std::vector<numeric_index_type> & n_oz = this->_sp->get_n_oz();

  // Make sure the sparsity pattern isn't empty
  libmesh_assert_equal_to (n_nz.size(), n_l);
  libmesh_assert_equal_to (n_oz.size(), n_l);

  // Epetra wants the total number of nonzeros, both local and remote.
  std::vector<int> n_nz_tot; /**/ n_nz_tot.reserve(n_nz.size());

  for (auto i : index_range(n_nz))
    n_nz_tot.push_back(std::min(n_nz[i] + n_oz[i], n));

  if (m==0)
    return;

  _graph = new Epetra_CrsGraph(Copy, *_map, n_nz_tot.data());

  // Tell the matrix about its structure.  Initialize it
  // to zero.
  for (numeric_index_type i=0; i<n_rows; i++)
    _graph->InsertGlobalIndices(_graph->GRID(i),
                                cast_int<numeric_index_type>(sparsity_pattern[i].size()),
                                const_cast<int *>(reinterpret_cast<const int *>(sparsity_pattern[i].data())));

  _graph->FillComplete();

  //Initialize the matrix
  libmesh_assert (!this->initialized());
  this->init ();
  libmesh_assert (this->initialized());
}



template <typename T>
void EpetraMatrix<T>::init (const numeric_index_type m,
                            const numeric_index_type n,
                            const numeric_index_type m_l,
                            const numeric_index_type libmesh_dbg_var(n_l),
                            const numeric_index_type nnz,
                            const numeric_index_type noz,
                            const numeric_index_type /* blocksize */)
{
  if ((m==0) || (n==0))
    return;

  {
    // Clear initialized matrices
    if (this->initialized())
      this->clear();

    libmesh_assert (!this->_mat);
    libmesh_assert (!this->_map);

    this->_is_initialized = true;
  }

  // error checking
#ifndef NDEBUG
  {
    libmesh_assert_equal_to (n, m);
    libmesh_assert_equal_to (n_l, m_l);

    numeric_index_type
      summed_m_l = m_l,
      summed_n_l = n_l;

    this->comm().sum (summed_m_l);
    this->comm().sum (summed_n_l);

    libmesh_assert_equal_to (m, summed_m_l);
    libmesh_assert_equal_to (n, summed_n_l);
  }
#endif

  // build a map defining the data distribution
  _map = new Epetra_Map (static_cast<int>(m),
                         m_l,
                         0,
                         Epetra_MpiComm (this->comm().get()));

  libmesh_assert_equal_to (static_cast<numeric_index_type>(_map->NumGlobalPoints()), m);
  libmesh_assert_equal_to (static_cast<numeric_index_type>(_map->MaxAllGID()+1), m);

  _mat = new Epetra_FECrsMatrix (Copy, *_map, nnz + noz);
}




template <typename T>
void EpetraMatrix<T>::init (const ParallelType)
{
  libmesh_assert(this->_dof_map);

  {
    // Clear initialized matrices
    if (this->initialized())
      this->clear();

    this->_is_initialized = true;
  }


  _mat = new Epetra_FECrsMatrix (Copy, *_graph);
}



template <typename T>
void EpetraMatrix<T>::zero ()
{
  libmesh_assert (this->initialized());

  _mat->Scale(0.0);
}



template <typename T>
std::unique_ptr<SparseMatrix<T>> EpetraMatrix<T>::zero_clone () const
{
  // This function is marked as "not implemented" since it hasn't been
  // tested, the code below might serve as a possible implementation.
  libmesh_not_implemented();

  // Make empty copy with matching comm, initialize, and return.
  auto mat_copy = std::make_unique<EpetraMatrix<T>>(this->comm());
  mat_copy->init();
  mat_copy->zero();

  return mat_copy;
}



template <typename T>
std::unique_ptr<SparseMatrix<T>> EpetraMatrix<T>::clone () const
{
  // We don't currently have a faster implementation than making a
  // zero clone and then filling in the values.
  auto mat_copy = this->zero_clone();
  mat_copy->add(1., *this);

  return mat_copy;
}



template <typename T>
void EpetraMatrix<T>::clear () noexcept
{
  // FIXME: clear() doesn't actually free the memory managed by this
  // class, so it probably leaks memory.
  // delete _mat;
  // delete _map;

  this->_is_initialized = false;
}



template <typename T>
Real EpetraMatrix<T>::l1_norm () const
{
  libmesh_assert (this->initialized());

  libmesh_assert(_mat);

  return static_cast<Real>(_mat->NormOne());
}



template <typename T>
Real EpetraMatrix<T>::linfty_norm () const
{
  libmesh_assert (this->initialized());


  libmesh_assert(_mat);

  return static_cast<Real>(_mat->NormInf());
}



template <typename T>
void EpetraMatrix<T>::add_matrix(const DenseMatrix<T> & dm,
                                 const std::vector<numeric_index_type> & rows,
                                 const std::vector<numeric_index_type> & cols)
{
  libmesh_assert (this->initialized());

  const numeric_index_type m = dm.m();
  const numeric_index_type n = dm.n();

  libmesh_assert_equal_to (rows.size(), m);
  libmesh_assert_equal_to (cols.size(), n);

  _mat->SumIntoGlobalValues(m, numeric_trilinos_cast(rows.data()),
                            n, numeric_trilinos_cast(cols.data()),
                            dm.get_values().data());
}






template <typename T>
void EpetraMatrix<T>::get_diagonal (NumericVector<T> & dest) const
{
  // Convert vector to EpetraVector.
  EpetraVector<T> * epetra_dest = cast_ptr<EpetraVector<T> *>(&dest);

  // Call Epetra function.
  _mat->ExtractDiagonalCopy(*(epetra_dest->vec()));
}



template <typename T>
void EpetraMatrix<T>::get_transpose (SparseMatrix<T> & dest) const
{
  // Make sure the SparseMatrix passed in is really a EpetraMatrix
  EpetraMatrix<T> & epetra_dest = cast_ref<EpetraMatrix<T> &>(dest);

  // We currently only support calling get_transpose() with ourself
  // as the destination. Previously, this called the default copy
  // constructor which was not safe because this class manually
  // manages memory.
  if (&epetra_dest != this)
    libmesh_not_implemented();

  epetra_dest._use_transpose = !epetra_dest._use_transpose;
  epetra_dest._mat->SetUseTranspose(epetra_dest._use_transpose);
}



template <typename T>
void EpetraMatrix<T>::get_row(numeric_index_type i,
                              std::vector<numeric_index_type> & indices,
                              std::vector<T> & values) const
{
  libmesh_assert (this->initialized());
  libmesh_assert(this->_mat);
  libmesh_assert (this->_mat->MyGlobalRow(static_cast<int>(i)));
  libmesh_assert_greater_equal (i, this->row_start());
  libmesh_assert_less (i, this->row_stop());

  int row_length;
  int * row_indices;
  double * row_values;

  _mat->ExtractMyRowView (i-this->row_start(),
                          row_length,
                          row_values,
                          row_indices);

  indices.resize(row_length);
  values.resize(row_length);

  for (auto i : make_range(row_length))
    {
      indices[i] = row_indices[i];
      values[i] = row_values[i];
    }
}



template <typename T>
EpetraMatrix<T>::EpetraMatrix(const Parallel::Communicator & comm) :
  SparseMatrix<T>(comm),
  _destroy_mat_on_exit(true),
  _use_transpose(false)
{}




template <typename T>
EpetraMatrix<T>::EpetraMatrix(Epetra_FECrsMatrix * m,
                              const Parallel::Communicator & comm) :
  SparseMatrix<T>(comm),
  _destroy_mat_on_exit(false),
  _use_transpose(false) // dumb guess is the best we can do...
{
  this->_mat = m;
  this->_is_initialized = true;
}




template <typename T>
EpetraMatrix<T>::~EpetraMatrix()
{
  this->clear();
}



template <typename T>
void EpetraMatrix<T>::close ()
{
  libmesh_assert(_mat);

  _mat->GlobalAssemble();
}



template <typename T>
numeric_index_type EpetraMatrix<T>::m () const
{
  libmesh_assert (this->initialized());

  return static_cast<numeric_index_type>(_mat->NumGlobalRows());
}



template <typename T>
numeric_index_type EpetraMatrix<T>::n () const
{
  libmesh_assert (this->initialized());

  return static_cast<numeric_index_type>(_mat->NumGlobalCols());
}



template <typename T>
numeric_index_type EpetraMatrix<T>::row_start () const
{
  libmesh_assert (this->initialized());
  libmesh_assert(_map);

  return static_cast<numeric_index_type>(_map->MinMyGID());
}



template <typename T>
numeric_index_type EpetraMatrix<T>::row_stop () const
{
  libmesh_assert (this->initialized());
  libmesh_assert(_map);

  return static_cast<numeric_index_type>(_map->MaxMyGID())+1;
}



template <typename T>
numeric_index_type EpetraMatrix<T>::col_start () const
{
  libmesh_assert (this->initialized());
  libmesh_assert(_map);

  return static_cast<numeric_index_type>(_map->MinMyGID());
}



template <typename T>
numeric_index_type EpetraMatrix<T>::col_stop () const
{
  libmesh_assert (this->initialized());
  libmesh_assert(_map);

  return static_cast<numeric_index_type>(_map->MaxMyGID())+1;
}



template <typename T>
void EpetraMatrix<T>::set (const numeric_index_type i,
                           const numeric_index_type j,
                           const T value)
{
  libmesh_assert (this->initialized());

  int
    epetra_i = static_cast<int>(i),
    epetra_j = static_cast<int>(j);

  T epetra_value = value;

  if (_mat->Filled())
    _mat->ReplaceGlobalValues (epetra_i, 1, &epetra_value, &epetra_j);
  else
    _mat->InsertGlobalValues (epetra_i, 1, &epetra_value, &epetra_j);
}



template <typename T>
void EpetraMatrix<T>::add (const numeric_index_type i,
                           const numeric_index_type j,
                           const T value)
{
  libmesh_assert (this->initialized());

  int
    epetra_i = static_cast<int>(i),
    epetra_j = static_cast<int>(j);

  T epetra_value = value;

  _mat->SumIntoGlobalValues (epetra_i, 1, &epetra_value, &epetra_j);
}



template <typename T>
void EpetraMatrix<T>::add_matrix(const DenseMatrix<T> & dm,
                                 const std::vector<numeric_index_type> & dof_indices)
{
  this->add_matrix (dm, dof_indices, dof_indices);
}



template <typename T>
void EpetraMatrix<T>::add (const T a_in, const SparseMatrix<T> & X_in)
{
#ifdef LIBMESH_TRILINOS_HAVE_EPETRAEXT
  libmesh_assert (this->initialized());

  // sanity check. but this cannot avoid
  // crash due to incompatible sparsity structure...
  libmesh_assert_equal_to (this->m(), X_in.m());
  libmesh_assert_equal_to (this->n(), X_in.n());

  const EpetraMatrix<T> * X =
    cast_ptr<const EpetraMatrix<T> *> (&X_in);

  EpetraExt::MatrixMatrix::Add (*X->_mat, false, a_in, *_mat, 1.);
#else
  libmesh_error_msg("ERROR: EpetraExt is required for EpetraMatrix::add()!");
#endif
}




template <typename T>
T EpetraMatrix<T>::operator () (const numeric_index_type i,
                                const numeric_index_type j) const
{
  libmesh_assert (this->initialized());
  libmesh_assert(this->_mat);
  libmesh_assert (this->_mat->MyGlobalRow(static_cast<int>(i)));
  libmesh_assert_greater_equal (i, this->row_start());
  libmesh_assert_less (i, this->row_stop());


  int row_length;
  int * row_indices;
  double * values;

  _mat->ExtractMyRowView (i-this->row_start(),
                          row_length,
                          values,
                          row_indices);

  //libMesh::out << "row_length=" << row_length << std::endl;

  int * index = std::lower_bound (row_indices, row_indices+row_length, j);

  libmesh_assert_less (*index, row_length);
  libmesh_assert_equal_to (static_cast<numeric_index_type>(row_indices[*index]), j);

  //libMesh::out << "val=" << values[*index] << std::endl;

  return values[*index];
}




template <typename T>
bool EpetraMatrix<T>::closed() const
{
  libmesh_assert (this->initialized());
  libmesh_assert(this->_mat);

  return this->_mat->Filled();
}


template <typename T>
void EpetraMatrix<T>::swap(EpetraMatrix<T> & m)
{
  std::swap(_mat, m._mat);
  std::swap(_destroy_mat_on_exit, m._destroy_mat_on_exit);
}





template <typename T>
void EpetraMatrix<T>::print_personal(std::ostream & os) const
{
  libmesh_assert (this->initialized());
  libmesh_assert(_mat);

  os << *_mat;
}



//------------------------------------------------------------------
// Explicit instantiations
template class LIBMESH_EXPORT EpetraMatrix<Number>;

} // namespace libMesh


#endif // LIBMESH_TRILINOS_HAVE_EPETRA
