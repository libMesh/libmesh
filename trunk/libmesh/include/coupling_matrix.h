// $Id: coupling_matrix.h,v 1.1 2003-01-20 16:31:21 jwpeterson Exp $

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



#ifndef __coupling_matrix_h__
#define __coupling_matrix_h__

// C++ includes
#include <vector>

// Local Includes
#include "mesh_common.h"




/**
 * This class defines a coupling matrix.  A coupling
 * matrix is simply a matrix of ones and zeros describing
 * how different components in a system couple with each
 * other.  A coupling matrix is necessarily square but not
 * necessarily symmetric.
 */
class CouplingMatrix
{
public:
  
  /**
   * Constructor.
   */
  CouplingMatrix (const unsigned int n=0);

  /**
   * @returns the (i,j) entry of the matrix.
   */
  unsigned char operator() (const unsigned int i,
			    const unsigned int j) const;
  
  /**
   * @returns the (i,j) entry of the matrix as
   * a writeable reference.
   */
  unsigned char & operator() (const unsigned int i,
			      const unsigned int j);
  
  /**
   * @returns the size of the matrix, i.e. N for an
   * NxN matrix.
   */ 
  unsigned int size() const;

  /**
   * Resizes the matrix and initializes
   * all entries to be 0.
   */
  void resize(const unsigned int n);

  /**
   * Clears the matrix.
   */
  void clear();

  /**
   * @returns true if the matrix is empty.
   */
  bool empty() const;
  
private:
  
  /**
   * The actual matrix values.  These
   * are stored as unsigned chars because
   * a vector of bools is not what you
   * think.
   */
  std::vector<unsigned char> values;  

  /**
   * The size of the matrix.
   */
  unsigned int _size;
};


  



//--------------------------------------------------
// CouplingMatrix inline methods
inline
CouplingMatrix::CouplingMatrix (const unsigned int n) :
  _size(n)
{
  resize(n);
};



inline
unsigned char CouplingMatrix::operator() (const unsigned int i,
					  const unsigned int j) const
{
  assert (i < _size);
  assert (j < _size);

  return values[i*_size + j];
};



inline
unsigned char & CouplingMatrix::operator() (const unsigned int i,
					    const unsigned int j)
{
  assert (i < _size);
  assert (j < _size);

  return values[i*_size + j];
};



inline
unsigned int CouplingMatrix::size() const
{
  return _size;
};



inline
void CouplingMatrix::resize(const unsigned int n)
{
  _size = n;

  values.resize(_size*_size);

  for (unsigned int i=0; i<values.size(); i++)
    values[i] = 0;
};



inline
void CouplingMatrix::clear()
{
  _size = 0;

  values.clear();
};



inline
bool CouplingMatrix::empty() const
{
  return (_size == 0);
};



#endif
