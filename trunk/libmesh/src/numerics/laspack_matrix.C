// $Id: laspack_matrix.C,v 1.7 2003-03-14 09:56:41 ddreyer Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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
#include <math.h>

// Local includes
#include "laspack_matrix.h"


#ifdef HAVE_LASPACK




//-----------------------------------------------------------------------
// LaspackMatrix members
template <typename T> 
void LaspackMatrix<T>::update_sparsity_pattern (std::vector<std::set<unsigned int> >&
					       sparsity_pattern)
{
  // clear data, start over
  clear ();    

  // big trouble if this fails!
  assert (_dof_map != NULL);
  
  const unsigned int n_rows = sparsity_pattern.size();

  //assert (n_rows == n());
  
  _row_start.resize(n_rows+1);

  // Initialize the _row_start data structure,
  // allocate storage for the _csr array
  {
    unsigned int size = 0;
 
    for (unsigned int row=0; row<n_rows; row++)
      {
	_row_start[row] = size;  
	size += sparsity_pattern[row].size();
	_row_start[row+1] = size;
      }
    
    _csr.resize(size);
  }


  // Initize the _csr data structure.
  {
    unsigned int pos = 0;
    
    for (unsigned int row=0; row<n_rows; row++)
      {
	// make sure that _row_start is set up properly
	assert (pos == _row_start[row]);
	 
	const std::set<unsigned int>& sparsity_row =
	  sparsity_pattern[row];
	
	// make sure that _row_start is set up properly
	assert (pos+sparsity_row.size() == _row_start[row+1]);

	std::set<unsigned int>::const_iterator it_end =
	  sparsity_row.end();

	for (std::set<unsigned int>::const_iterator it =
	       sparsity_row.begin(); it != it_end; ++it)
	  _csr[pos++] = *it;
      }
    
    // make sure that _row_start is set up properly
    assert (pos == _row_start[n_rows]);
  }


  // Initialize the matrix
  assert (!initialized());
  init ();
  assert (initialized());
  //std::cout << "n_rows=" << n_rows << std::endl;
  //std::cout << "m()=" << m() << std::endl;
  assert (n_rows == m());

  // Tell the matrix about its structure.  Initialize it
  // to zero.
  for (unsigned int i=0; i<n_rows; i++)
    {
      const unsigned int rs     = _row_start[i];
      const unsigned int length = _row_start[i+1] - rs;
      
      //std::cout << "m()=" << m() << std::endl;
      Q_SetLen (&_QMat, i+1, length);
      //std::cout << "m()=" << m() << std::endl;

      for (unsigned int l=0; l<length; l++)
	{
	  const unsigned int j = _csr[rs+l];

	  // sanity check
	  //std::cout << "m()=" << m() << std::endl;
	  //std::cout << "(i,j,l) = (" << i
	  //	    << "," << j
	  //	    << "," << l
	  // 	    << ")" << std::endl;
	  //std::cout << "pos(i,j)=" << pos(i,j)
	  //          << std::endl;	  
	  assert (pos(i,j) == l);
	  Q_SetEntry (&_QMat, i+1, l, j+1, 0.);
	}
    }
  
  // That's it!
  //here();
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
  assert (m == m_l);
  assert (n == n_l);
  assert (m == n);
  assert (nnz > 0);


  std::cerr << "ERROR: Only the init() member that uses the" << std::endl
	    << "DofMap is implemented for Laspack matrices!" << std::endl;
  error();

  _is_initialized = true;
}



template <typename T> 
void LaspackMatrix<T>::init ()
{
  // Ignore calls on initialized objects
  if (initialized())
    return;
  
  
  // We need the DofMap for this!
  assert (_dof_map != NULL);
  assert (!initialized());  

  const unsigned int m   = _dof_map->n_dofs();
  const unsigned int n   = m;
  const unsigned int n_l = _dof_map->n_dofs_on_processor(0); 
  const unsigned int m_l = n_l;

  // Laspack Matrices only work for uniprocessor cases
  assert (m   == n);
  assert (m_l == m);
  assert (n_l == n);

  const std::vector<unsigned int>& n_nz = _dof_map->get_n_nz();
  const std::vector<unsigned int>& n_oz = _dof_map->get_n_oz();

  // Make sure the sparsity pattern isn't empty
  assert (n_nz.size() == n_l);
  assert (n_oz.size() == n_l);
  
  if (m==0)
    return;

  Q_Constr(&_QMat, const_cast<char*>("Mat"), m, _LPFalse, Rowws, Normal, _LPTrue);

  _is_initialized = true;
  
  assert (m == this->m());
}


//------------------------------------------------------------------
// Explicit instantiations
template class LaspackMatrix<Number>;
 

#endif // #ifdef HAVE_LASPACK
