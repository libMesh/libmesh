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
void LaspackMatrix<T>::update_sparsity_pattern (const SparsityPattern::Graph &sparsity_pattern)
{
  // clear data, start over
  this->clear ();

  // big trouble if this fails!
  libmesh_assert(this->_dof_map);

  const unsigned int n_rows = sparsity_pattern.size();

  // Initialize the _row_start data structure,
  // allocate storage for the _csr array
  {
    unsigned int size = 0;

    for (unsigned int row=0; row<n_rows; row++)
      size += sparsity_pattern[row].size();

    _csr.resize       (size);
    _row_start.reserve(n_rows + 1);
  }


  // Initize the _csr data structure.
  {
    std::vector<unsigned int>::iterator pos = _csr.begin();

    _row_start.push_back (pos);

    for (unsigned int row=0; row<n_rows; row++)
      {
	// insert the row indices
	for (SparsityPattern::Row::const_iterator col = sparsity_pattern[row].begin();
	     col != sparsity_pattern[row].end(); ++col)
	  {
	    libmesh_assert (pos != _csr.end());
	    *pos = *col;
	    ++pos;
	  }

	_row_start.push_back (pos);
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
  for (unsigned int i=0; i<n_rows; i++)
    {
      const std::vector<unsigned int>::const_iterator
	rs = _row_start[i];

      const unsigned int length = _row_start[i+1] - rs;

      Q_SetLen (&_QMat, i+1, length);

      for (unsigned int l=0; l<length; l++)
	{
	  const unsigned int j = *(rs+l);

	  // sanity check
	  //libMesh::out << "m()=" << m() << std::endl;
	  //libMesh::out << "(i,j,l) = (" << i
	  //	          << "," << j
	  //	          << "," << l
	  // 	          << ")" << std::endl;
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
void LaspackMatrix<T>::init (const unsigned int libmesh_dbg_var(m),
			     const unsigned int libmesh_dbg_var(n),
			     const unsigned int libmesh_dbg_var(m_l),
			     const unsigned int libmesh_dbg_var(n_l),
			     const unsigned int libmesh_dbg_var(nnz),
			     const unsigned int)
{
  // noz ignored...  only used for multiple processors!
  libmesh_assert_equal_to (m, m_l);
  libmesh_assert_equal_to (n, n_l);
  libmesh_assert_equal_to (m, n);
  libmesh_assert_greater (nnz, 0);


  libMesh::err << "ERROR: Only the init() member that uses the" << std::endl
	        << "DofMap is implemented for Laspack matrices!" << std::endl;
  libmesh_error();

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

  // Clear intialized matrices
  if (this->initialized())
    this->clear();

  const unsigned int m   = this->_dof_map->n_dofs();
#ifndef NDEBUG
  // The following variables are only used for assertions,
  // so avoid declaring them when asserts are inactive.
  const unsigned int n   = m;
  const unsigned int n_l = this->_dof_map->n_dofs_on_processor(0);
  const unsigned int m_l = n_l;
#endif

  // Laspack Matrices only work for uniprocessor cases
  libmesh_assert_equal_to (m, n);
  libmesh_assert_equal_to (m_l, m);
  libmesh_assert_equal_to (n_l, n);

#ifndef NDEBUG
  // The following variables are only used for assertions,
  // so avoid declaring them when asserts are inactive.
  const std::vector<unsigned int>& n_nz = this->_dof_map->get_n_nz();
  const std::vector<unsigned int>& n_oz = this->_dof_map->get_n_oz();
#endif

  // Make sure the sparsity pattern isn't empty
  libmesh_assert_equal_to (n_nz.size(), n_l);
  libmesh_assert_equal_to (n_oz.size(), n_l);

  if (m==0)
    return;

  Q_Constr(&_QMat, const_cast<char*>("Mat"), m, _LPFalse, Rowws, Normal, _LPTrue);

  this->_is_initialized = true;

  libmesh_assert_equal_to (m, this->m());
}



template <typename T>
void LaspackMatrix<T>::add_matrix(const DenseMatrix<T>& dm,
				  const std::vector<unsigned int>& rows,
				  const std::vector<unsigned int>& cols)

{
  libmesh_assert (this->initialized());
  libmesh_assert_equal_to (dm.m(), rows.size());
  libmesh_assert_equal_to (dm.n(), cols.size());


  for (unsigned int i=0; i<rows.size(); i++)
    for (unsigned int j=0; j<cols.size(); j++)
      this->add(rows[i],cols[j],dm(i,j));
}



template <typename T>
void LaspackMatrix<T>::get_diagonal (NumericVector<T>& /*dest*/) const
{
  libmesh_not_implemented();
}



template <typename T>
void LaspackMatrix<T>::get_transpose (SparseMatrix<T>& /*dest*/) const
{
  libmesh_not_implemented();
}



template <typename T>
LaspackMatrix<T>::LaspackMatrix () :
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
  const unsigned int n_rows = this->m();

  for (unsigned int row=0; row<n_rows; row++)
    {
      const std::vector<unsigned int>::const_iterator
	r_start = _row_start[row];

      const unsigned int len = (_row_start[row+1] - _row_start[row]);

      // Make sure we agree on the row length
      libmesh_assert_equal_to (len, Q_GetLen(&_QMat, row+1));

      for (unsigned int l=0; l<len; l++)
	{
	  const unsigned int j = *(r_start + l);

	  // Make sure the data structures are working
	  libmesh_assert_equal_to ((j+1), Q_GetPos (&_QMat, row+1, l));

	  Q_SetEntry (&_QMat, row+1, l, j+1, 0.);
	}
    }
}



template <typename T>
unsigned int LaspackMatrix<T>::m () const
{
  libmesh_assert (this->initialized());

  return static_cast<unsigned int>(Q_GetDim(const_cast<QMatrix*>(&_QMat)));
}



template <typename T>
unsigned int LaspackMatrix<T>::n () const
{
  libmesh_assert (this->initialized());

  return static_cast<unsigned int>(Q_GetDim(const_cast<QMatrix*>(&_QMat)));
}



template <typename T>
unsigned int LaspackMatrix<T>::row_start () const
{
  return 0;
}



template <typename T>
unsigned int LaspackMatrix<T>::row_stop () const
{
  return this->m();
}



template <typename T>
void LaspackMatrix<T>::set (const unsigned int i,
			    const unsigned int j,
			    const T value)
{
  libmesh_assert (this->initialized());
  libmesh_assert_less (i, this->m());
  libmesh_assert_less (j, this->n());

  const unsigned int position = this->pos(i,j);

  // Sanity check
  libmesh_assert_equal_to (*(_row_start[i]+position), j);
  libmesh_assert_equal_to ((j+1), Q_GetPos (&_QMat, i+1, position));

  Q_SetEntry (&_QMat, i+1, position, j+1, value);
}



template <typename T>
void LaspackMatrix<T>::add (const unsigned int i,
			    const unsigned int j,
			    const T value)
{
  libmesh_assert (this->initialized());
  libmesh_assert_less (i, this->m());
  libmesh_assert_less (j, this->n());

  const unsigned int position = this->pos(i,j);

  // Sanity check
  libmesh_assert_equal_to (*(_row_start[i]+position), j);

  Q_AddVal (&_QMat, i+1, position, value);
}



template <typename T>
void LaspackMatrix<T>::add_matrix(const DenseMatrix<T>& dm,
				  const std::vector<unsigned int>& dof_indices)
{
  this->add_matrix (dm, dof_indices, dof_indices);
}



template <typename T>
void LaspackMatrix<T>::add (const T a_in, SparseMatrix<T> &X_in)
{
  libmesh_assert (this->initialized());
  libmesh_assert_equal_to (this->m(), X_in.m());
  libmesh_assert_equal_to (this->n(), X_in.n());

  LaspackMatrix<T>* X = libmesh_cast_ptr<LaspackMatrix<T>*> (&X_in);
  _LPNumber         a = static_cast<_LPNumber>          (a_in);

  libmesh_assert(X);

  // loops taken from LaspackMatrix<T>::zero ()

  const unsigned int n_rows = this->m();

  for (unsigned int row=0; row<n_rows; row++)
    {
      const std::vector<unsigned int>::const_iterator
	r_start = _row_start[row];

      const unsigned int len = (_row_start[row+1] - _row_start[row]);

      // Make sure we agree on the row length
      libmesh_assert_equal_to (len, Q_GetLen(&_QMat, row+1));
      // compare matrix sparsity structures
      libmesh_assert_equal_to (len, Q_GetLen(&(X->_QMat), row+1));


      for (unsigned int l=0; l<len; l++)
	{
	  const unsigned int j = *(r_start + l);

	  // Make sure the data structures are working
	  libmesh_assert_equal_to ((j+1), Q_GetPos (&_QMat, row+1, l));

	  const _LPNumber value = a * Q_GetEl(const_cast<QMatrix*>(&(X->_QMat)), row+1, j+1);
	  Q_AddVal   (&_QMat, row+1, l, value);
	}
    }
}




template <typename T>
T LaspackMatrix<T>::operator () (const unsigned int i,
				 const unsigned int j) const
{
  libmesh_assert (this->initialized());
  libmesh_assert_less (i, this->m());
  libmesh_assert_less (j, this->n());

  return Q_GetEl (const_cast<QMatrix*>(&_QMat), i+1, j+1);
}



template <typename T>
unsigned int LaspackMatrix<T>::pos (const unsigned int i,
				    const unsigned int j) const
{
  libmesh_assert_less (i, this->m());
  libmesh_assert_less (j, this->n());
  libmesh_assert_less (i+1, _row_start.size());
  libmesh_assert (_row_start.back() == _csr.end());

  // note this requires the _csr to be
  std::pair<std::vector<unsigned int>::const_iterator,
	    std::vector<unsigned int>::const_iterator> p =
    std::equal_range (_row_start[i],
		      _row_start[i+1],
		      j);

  // Make sure the row contains the element j
  libmesh_assert (p.first != p.second);

  // Make sure the values match
  libmesh_assert (*p.first == j);

  // Return the position in the compressed row
  return std::distance (_row_start[i], p.first);
}



//------------------------------------------------------------------
// Explicit instantiations
template class LaspackMatrix<Number>;

} // namespace libMesh


#endif // #ifdef LIBMESH_HAVE_LASPACK
