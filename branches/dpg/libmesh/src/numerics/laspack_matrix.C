// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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
#include "libmesh_config.h"

#ifdef LIBMESH_HAVE_LASPACK

#include "laspack_matrix.h"
#include "dense_matrix.h"
#include "dof_map.h"
#include "sparsity_pattern.h"

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
  libmesh_assert (this->_dof_map != NULL);
  
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
  libmesh_assert (n_rows == this->m());

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
	  libmesh_assert (this->pos(i,j) == l);
	  Q_SetEntry (&_QMat, i+1, l, j+1, 0.);
	}
    }
  
  // That's it!
  //libmesh_here();
}



template <typename T> 
void LaspackMatrix<T>::init (const unsigned int m,
			     const unsigned int n,
			     const unsigned int m_l,
			     const unsigned int n_l,
			     const unsigned int nnz,
			     const unsigned int)
{
  // noz ignored...  only used for multiple processors!
  libmesh_assert (m == m_l);
  libmesh_assert (n == n_l);
  libmesh_assert (m == n);
  libmesh_assert (nnz > 0);


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
  libmesh_assert (this->_dof_map != NULL);

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
  libmesh_assert (m   == n);
  libmesh_assert (m_l == m);
  libmesh_assert (n_l == n);

#ifndef NDEBUG
  // The following variables are only used for assertions,
  // so avoid declaring them when asserts are inactive.
  const std::vector<unsigned int>& n_nz = this->_dof_map->get_n_nz();
  const std::vector<unsigned int>& n_oz = this->_dof_map->get_n_oz();
#endif
  
  // Make sure the sparsity pattern isn't empty
  libmesh_assert (n_nz.size() == n_l);
  libmesh_assert (n_oz.size() == n_l);
  
  if (m==0)
    return;

  Q_Constr(&_QMat, const_cast<char*>("Mat"), m, _LPFalse, Rowws, Normal, _LPTrue);

  this->_is_initialized = true;
  
  libmesh_assert (m == this->m());
}



template <typename T> 
void LaspackMatrix<T>::add_matrix(const DenseMatrix<T>& dm,
				  const std::vector<unsigned int>& rows,
				  const std::vector<unsigned int>& cols)
		    
{
  libmesh_assert (this->initialized());
  libmesh_assert (dm.m() == rows.size());
  libmesh_assert (dm.n() == cols.size());

  
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



//------------------------------------------------------------------
// Explicit instantiations
template class LaspackMatrix<Number>;

} // namespace libMesh
 

#endif // #ifdef LIBMESH_HAVE_LASPACK
