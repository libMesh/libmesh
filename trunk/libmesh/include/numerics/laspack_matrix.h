// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __laspack_matrix_h__
#define __laspack_matrix_h__

#include "libmesh_config.h"

#ifdef HAVE_LASPACK 

// C++ includes
#include <algorithm>

// Local includes
#include "sparse_matrix.h"
#include "dense_matrix.h"


#include <qmatrix.h>



// Forward declarations
template <typename T> class LaspackVector;
template <typename T> class LaspackLinearSolver;



/**
 * Generic laspack matrix. This class contains
 * pure virtual members that must be overloaded
 * in derived classes.  Using a derived class
 * allows for uniform access to laspack matrices
 * from various different solver packages in
 * different formats.
 * Currently Laspack only supports real datatypes, so
 * this class is a full specialization of \p SparseMatrix<>
 * with \p T = \p Real

 *
 * @author Benjamin S. Kirk, 2003
 */

template <typename T>
class LaspackMatrix : public SparseMatrix<T>
{

public:
  /**
   * Constructor; initializes the matrix to
   * be empty, without any structure, i.e.
   * the matrix is not usable at all. This
   * constructor is therefore only useful
   * for matrices which are members of a
   * class. All other matrices should be
   * created at a point in the data flow
   * where all necessary information is
   * available.
   *
   * You have to initialize
   * the matrix before usage with
   * \p init(...).
   */
  LaspackMatrix ();

  /**
   * Destructor. Free all memory, but do not
   * release the memory of the sparsity
   * structure.
   */
  ~LaspackMatrix ();

  /**
   * The \p LaspackMatrix needs the full sparsity pattern.
   */ 
  bool need_full_sparsity_pattern() const 
  { return true; }

  /**
   * Updates the matrix sparsity pattern.  This will
   * tell the underlying matrix storage scheme how
   * to map the \f$ (i,j) \f$ elements.
   */
  void update_sparsity_pattern (const SparsityPattern::Graph &);
  
  /**
   * Initialize a Laspack matrix that is of global
   * dimension \f$ m \times  n \f$ with local dimensions
   * \f$ m_l \times n_l \f$.  \p nnz is the number of on-processor
   * nonzeros per row (defaults to 30).
   * \p noz is the number of on-processor
   * nonzeros per row (defaults to 30).
   */
  void init (const unsigned int m,
	     const unsigned int n,
	     const unsigned int m_l,
	     const unsigned int n_l,
	     const unsigned int nnz=30,
	     const unsigned int noz=10);

  /**
   * Initialize using sparsity structure computed by \p dof_map.
   */   
  void init ();
  
  /**
   * Release all memory and return
   * to a state just like after
   * having called the default
   * constructor. 
   */
  void clear ();

  /**
   * Set all entries to 0.
   */
  void zero ();
  
  /**
   * Close the matrix.  Dummy routine.  After calling
   * this method \p closed() is true and the matrix can
   * be used in computations.
   */
  void close () const { const_cast<LaspackMatrix<T>*>(this)->_closed = true; }
  
  /**
   * @returns \p m, the row-dimension of
   * the matrix where the marix is \f$ M \times N \f$.
   */  
  unsigned int m () const;

  /**
   * @returns \p n, the column-dimension of
   * the matrix where the marix is \f$ M \times N \f$.
   */  
  unsigned int n () const;

  /**
   * return row_start, the index of the first
   * matrix row stored on this processor
   */
  unsigned int row_start () const;

  /**
   * return row_stop, the index of the last
   * matrix row (+1) stored on this processor
   */
  unsigned int row_stop () const;

  /**
   * Set the element \p (i,j) to \p value.
   * Throws an error if the entry does
   * not exist. Still, it is allowed to store
   * zero values in non-existent fields.
   */
  void set (const unsigned int i,
	    const unsigned int j,
	    const T value);
    
  /**
   * Add \p value to the element
   * \p (i,j).  Throws an error if
   * the entry does not
   * exist. Still, it is allowed to
   * store zero values in
   * non-existent fields.
   */
  void add (const unsigned int i,
	    const unsigned int j,
	    const T value);

  /**
   * Add the full matrix to the
   * Laspack matrix.  This is useful
   * for adding an element matrix
   * at assembly time
   */
    
  void add_matrix (const DenseMatrix<T> &dm,
		   const std::vector<unsigned int> &rows,
		   const std::vector<unsigned int> &cols);
  
  /**
   * Same, but assumes the row and column maps are the same.
   * Thus the matrix \p dm must be square.
   */
  void add_matrix (const DenseMatrix<T> &dm,
		   const std::vector<unsigned int> &dof_indices);
      
  /**
   * Add a Sparse matrix \p X, scaled with \p a, to \p this,
   * stores the result in \p this: \f$\texttt{this} += a*X \f$.
   * \p LASPACK does not provide a true \p axpy for matrices,
   * so a hand-coded version with hopefully acceptable performance
   * is provided.
   */
  void add (const T a, SparseMatrix<T> &X);
      
  /**
   * Return the value of the entry
   * \p (i,j).  This may be an
   * expensive operation and you
   * should always take care where
   * to call this function.  In
   * order to avoid abuse, this
   * function throws an exception
   * if the required element does
   * not exist in the matrix.
   *
   * In case you want a function
   * that returns zero instead (for
   * entries that are not in the
   * sparsity pattern of the
   * matrix), use the \p el
   * function.
   */
  T operator () (const unsigned int i,
		 const unsigned int j) const;

  /**
   * Return the l1-norm of the matrix, that is
   * \f$|M|_1=max_{all columns j}\sum_{all 
   * rows i} |M_ij|\f$,
   * (max. sum of columns).
   * This is the
   * natural matrix norm that is compatible
   * to the l1-norm for vectors, i.e.
   * \f$|Mv|_1\leq |M|_1 |v|_1\f$.
   */
  Real l1_norm () const { error(); return 0.; }

  /**
   * Return the linfty-norm of the
   * matrix, that is
   * \f$|M|_\infty=max_{all rows i}\sum_{all 
   * columns j} |M_ij|\f$,
   * (max. sum of rows).
   * This is the
   * natural matrix norm that is compatible
   * to the linfty-norm of vectors, i.e.
   * \f$|Mv|_\infty \leq |M|_\infty |v|_\infty\f$.
   */
  Real linfty_norm () const { error(); return 0.; }

  /**
   * see if Laspack matrix has been closed
   * and fully assembled yet
   */
  bool closed() const { return _closed; }
  
  /**
   * Print the contents of the matrix to the screen,
   * currently identical to \p print().
   */
  void print_personal(std::ostream& os=std::cout) const { this->print(os); }

  
private:

  /**
   * @returns the position in the compressed row
   * storage scheme of the \f$ (i,j) \f$ element.
   */
  unsigned int pos (const unsigned int i,
		    const unsigned int j) const;
  
  /**
   *  The Laspack sparse matrix pointer.
   */
  QMatrix _QMat;

  /**
   * The compressed row indices.
   */
  std::vector<unsigned int> _csr;

  /**
   * The start of each row in the compressed
   * row index data structure.
   */
  std::vector<std::vector<unsigned int>::const_iterator> _row_start;

  /**
   * Flag indicating if the matrix has been closed yet.
   */
  bool _closed;

  /**
   * Make other Laspack datatypes friends
   */
  friend class LaspackVector<T>;
  friend class LaspackLinearSolver<T>;
};



//-----------------------------------------------------------------------
// LaspackMatrix class inline members
template <typename T>
inline
LaspackMatrix<T>::LaspackMatrix () :
  _closed (false)
{
}



template <typename T>
inline
LaspackMatrix<T>::~LaspackMatrix ()
{
  this->clear ();
}



template <typename T>
inline
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
inline
void LaspackMatrix<T>::zero ()
{
  const unsigned int n_rows = this->m();

  for (unsigned int row=0; row<n_rows; row++)
    {
      const std::vector<unsigned int>::const_iterator
	r_start = _row_start[row];
      
      const unsigned int len = (_row_start[row+1] - _row_start[row]);
	
      // Make sure we agree on the row length
      assert (len == Q_GetLen(&_QMat, row+1));
      
      for (unsigned int l=0; l<len; l++)
	{
	  const unsigned int j = *(r_start + l);

	  // Make sure the data structures are working
	  assert ((j+1) == Q_GetPos (&_QMat, row+1, l));
	  
	  Q_SetEntry (&_QMat, row+1, l, j+1, 0.);
	}
    }    
}



template <typename T> 
inline
unsigned int LaspackMatrix<T>::m () const
{
  assert (this->initialized());

  return static_cast<unsigned int>(Q_GetDim(const_cast<QMatrix*>(&_QMat)));
}



template <typename T> 
inline
unsigned int LaspackMatrix<T>::n () const
{
  assert (this->initialized());
  
  return static_cast<unsigned int>(Q_GetDim(const_cast<QMatrix*>(&_QMat)));
}



template <typename T> 
inline
unsigned int LaspackMatrix<T>::row_start () const
{
  return 0;
}



template <typename T> 
inline
unsigned int LaspackMatrix<T>::row_stop () const
{
  return this->m();
}



template <typename T> 
inline
void LaspackMatrix<T>::set (const unsigned int i,
			    const unsigned int j,
			    const T value)
{
  assert (this->initialized());
  assert (i < this->m());
  assert (j < this->n());
  
  const unsigned int position = this->pos(i,j);

  // Sanity check
  assert (*(_row_start[i]+position) == j);
  assert ((j+1) == Q_GetPos (&_QMat, i+1, position));

  Q_SetEntry (&_QMat, i+1, position, j+1, value);
}



template <typename T> 
inline
void LaspackMatrix<T>::add (const unsigned int i,
			    const unsigned int j,
			    const T value)
{
  assert (this->initialized());
  assert (i < this->m());
  assert (j < this->n());
  
  const unsigned int position = this->pos(i,j);

  // Sanity check
  assert (*(_row_start[i]+position) == j);

  Q_AddVal (&_QMat, i+1, position, value);
}



template <typename T> 
inline
void LaspackMatrix<T>::add_matrix(const DenseMatrix<T>& dm,
				  const std::vector<unsigned int>& dof_indices)
{
  this->add_matrix (dm, dof_indices, dof_indices);
}



template <typename T> 
inline
void LaspackMatrix<T>::add_matrix(const DenseMatrix<T>& dm,
				  const std::vector<unsigned int>& rows,
				  const std::vector<unsigned int>& cols)
		    
{
  assert (this->initialized());
  assert (dm.m() == rows.size());
  assert (dm.n() == cols.size());

  
  for (unsigned int i=0; i<rows.size(); i++)
    for (unsigned int j=0; j<cols.size(); j++)
      this->add(rows[i],cols[j],dm(i,j));
}



template <typename T>
void LaspackMatrix<T>::add (const T a_in, SparseMatrix<T> &X_in)
{
  assert (this->initialized());
  assert (this->m() == X_in.m());
  assert (this->n() == X_in.n());

  LaspackMatrix<T>* X = dynamic_cast<LaspackMatrix<T>*> (&X_in);
  _LPNumber         a = static_cast<_LPNumber>          (a_in);
  
  assert(X != NULL);
  
  // loops taken from LaspackMatrix<T>::zero ()

  const unsigned int n_rows = this->m();

  for (unsigned int row=0; row<n_rows; row++)
    {
      const std::vector<unsigned int>::const_iterator
	r_start = _row_start[row];
      
      const unsigned int len = (_row_start[row+1] - _row_start[row]);

      // Make sure we agree on the row length
      assert (len == Q_GetLen(&_QMat, row+1));
      // compare matrix sparsity structures
      assert (len == Q_GetLen(&(X->_QMat), row+1));
	
      
      for (unsigned int l=0; l<len; l++)
	{
	  const unsigned int j = *(r_start + l);

	  // Make sure the data structures are working
	  assert ((j+1) == Q_GetPos (&_QMat, row+1, l));

	  const _LPNumber value = a * Q_GetEl(const_cast<QMatrix*>(&(X->_QMat)), row+1, j+1);
	  Q_AddVal   (&_QMat, row+1, l, value);
	}
    }    
}




template <typename T> 
inline
T LaspackMatrix<T>::operator () (const unsigned int i,
				 const unsigned int j) const
{
  assert (this->initialized());
  assert (i < this->m());
  assert (j < this->n());
  
  return Q_GetEl (const_cast<QMatrix*>(&_QMat), i+1, j+1);
}



template <typename T> 
inline
unsigned int LaspackMatrix<T>::pos (const unsigned int i,
				    const unsigned int j) const
{
  //std::cout << "m()=" << m() << std::endl;
  assert (i < this->m());
  assert (j < this->n());
  assert (i+1 < _row_start.size());
  assert (_row_start.back() == _csr.end());

  // note this requires the _csr to be 
  std::pair<std::vector<unsigned int>::const_iterator,
	    std::vector<unsigned int>::const_iterator> p =
    std::equal_range (_row_start[i],
		      _row_start[i+1],
		      j);

  // Make sure the row contains the element j
  assert (p.first != p.second);

  // Make sure the values match
  assert (*p.first == j);

  // Return the position in the compressed row
  return std::distance (_row_start[i], p.first);
}


#endif // #ifdef HAVE_LASPACK
#endif // #ifdef __laspack_matrix_h__
