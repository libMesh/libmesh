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



#ifndef __sparse_matrix_h__
#define __sparse_matrix_h__


// C++ includes
#include <iomanip>
#include <vector>

// Local includes
#include "libmesh_common.h"
#include "auto_ptr.h"
#include "reference_counted_object.h"
#include "libmesh.h"

// forward declarations
template <typename T> class SparseMatrix;
template <typename T> class DenseMatrix;
template <typename T> inline std::ostream& operator << (std::ostream& os, const SparseMatrix<T>& m);
class DofMap;
namespace SparsityPattern { class Graph; }
template <typename T> class NumericVector;


/**
 * Generic sparse matrix. This class contains
 * pure virtual members that must be overloaded
 * in derived classes.  Using a derived class
 * allows for uniform access to sparse matrices
 * from various different solver packages in
 * different formats.
 *
 * @author Benjamin S. Kirk, 2003
 */

template <typename T>
class SparseMatrix : public ReferenceCountedObject<SparseMatrix<T> >
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
  SparseMatrix ();

  /**
   * Destructor. Free all memory, but do not
   * release the memory of the sparsity
   * structure.
   */
  virtual ~SparseMatrix ();

  /**
   * Builds a \p SparseMatrix<T> using the linear solver package specified by
   * \p solver_package
   */
  static AutoPtr<SparseMatrix<T> >
  build(const SolverPackage solver_package = libMesh::default_solver_package());
  
  /**
   * @returns true if the matrix has been initialized,
   * false otherwise.
   */
  virtual bool initialized() const { return _is_initialized; }

  /**
   * Get a pointer to the \p DofMap to use.
   */
  void attach_dof_map (const DofMap& dof_map)
  { _dof_map = &dof_map; }

  /**
   * \p returns true if this sparse matrix format needs to be fed the 
   * graph of the sparse matrix.  This is true in the case of the \p LaspackMatrix, 
   * but not for the \p PetscMatrix.  In the case where the full graph is not
   * required we can efficiently approximate it to provide a good estimate of the 
   * required size of the sparse matrix.
   */ 
  virtual bool need_full_sparsity_pattern() const 
  { return false; }  

  /**
   * Updates the matrix sparsity pattern. When your \p SparseMatrix<T>
   * implementation does not need this data simply do
   * not overload this method.
   */
  virtual void update_sparsity_pattern (const SparsityPattern::Graph &) {}
  
  /**
   * Initialize a Sparse matrix that is of global
   * dimension \f$ m \times  n \f$ with local dimensions
   * \f$ m_l \times n_l \f$.  \p nnz is the number of on-processor
   * nonzeros per row (defaults to 30).
   * \p noz is the number of on-processor
   * nonzeros per row (defaults to 10).
   */
  virtual void init (const unsigned int m,
		     const unsigned int n,
		     const unsigned int m_l,
		     const unsigned int n_l,
		     const unsigned int nnz=30,
		     const unsigned int noz=10) = 0;

  /**
   * Initialize using sparsity structure computed by \p dof_map.
   */   
  virtual void init () = 0;
  
  /**
   * Release all memory and return
   * to a state just like after
   * having called the default
   * constructor. 
   */
  virtual void clear () = 0;

  /**
   * Set all entries to 0.
   */
  virtual void zero () = 0;

  /**
   * Set all row entries to 0 then puts diag_value in the diagonal entry
   */
  virtual void zero_rows (std::vector<int> & rows, T diag_value = 0.0);
  
  /**
   * Call the Sparse assemble routines.
   * sends necessary messages to other
   * processors
   */
  virtual void close () const = 0;
  
  /**
   * @returns \p m, the row-dimension of
   * the matrix where the marix is \f$ M \times N \f$.
   */  
  virtual unsigned int m () const = 0;

  /**
   * @returns \p n, the column-dimension of
   * the matrix where the marix is \f$ M \times N \f$.
   */  
  virtual unsigned int n () const = 0;

  /**
   * return row_start, the index of the first
   * matrix row stored on this processor
   */
  virtual unsigned int row_start () const = 0;

  /**
   * return row_stop, the index of the last
   * matrix row (+1) stored on this processor
   */
  virtual unsigned int row_stop () const = 0;

  /**
   * Set the element \p (i,j) to \p value.
   * Throws an error if the entry does
   * not exist. Still, it is allowed to store
   * zero values in non-existent fields.
   */
  virtual void set (const unsigned int i,
		    const unsigned int j,
		    const T value) = 0;
    
  /**
   * Add \p value to the element
   * \p (i,j).  Throws an error if
   * the entry does not
   * exist. Still, it is allowed to
   * store zero values in
   * non-existent fields.
   */
  virtual void add (const unsigned int i,
		    const unsigned int j,
		    const T value) = 0;

  /**
   * Add the full matrix to the
   * Sparse matrix.  This is useful
   * for adding an element matrix
   * at assembly time
   */
  virtual void add_matrix (const DenseMatrix<T> &dm,
			   const std::vector<unsigned int> &rows,
			   const std::vector<unsigned int> &cols) = 0;
  
  /**
   * Same, but assumes the row and column maps are the same.
   * Thus the matrix \p dm must be square.
   */
  virtual void add_matrix (const DenseMatrix<T> &dm,
			   const std::vector<unsigned int> &dof_indices) = 0;
      
  /**
   * Add a Sparse matrix \p _X, scaled with \p _a, to \p this,
   * stores the result in \p this: 
   * \f$\texttt{this} = \_a*\_X + \texttt{this} \f$.
   */
  virtual void add (const T, SparseMatrix<T> &) = 0;

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
  virtual T operator () (const unsigned int i,
			 const unsigned int j) const = 0;

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
  virtual Real l1_norm () const = 0;

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
  virtual Real linfty_norm () const = 0;

  /**
   * see if Sparse matrix has been closed
   * and fully assembled yet
   */
  virtual bool closed() const = 0;

  /**
   * Print the contents of the matrix to the screen
   * in a uniform style, regardless of matrix/solver
   * package being used.
   */
  void print(std::ostream& os=libMesh::out) const;

  /**
   * Same as the print method above, but allows you
   * to print to a stream in the standard syntax.
   */
  template <typename U>
  friend std::ostream& operator << (std::ostream& os, const SparseMatrix<U>& m);
  
  /**
   * Print the contents of the matrix to the screen
   * in a package-personalized style, if available.
   */
  virtual void print_personal(std::ostream& os=libMesh::out) const = 0;
  
  /**
   * Print the contents of the matrix in Matlab's
   * sparse matrix format. Optionally prints the
   * matrix to the file named \p name.  If \p name
   * is not specified it is dumped to the screen.
x   */
  virtual void print_matlab(const std::string name="NULL") const
  {
    libMesh::err << "ERROR: Not Implemented in base class yet!" << std::endl;
    libMesh::err << "ERROR writing MATLAB file " << name << std::endl;
    libmesh_error();
  }

  /**
   * This function creates a matrix called "submatrix" which is defined
   * by the row and column indices given in the "rows" and "cols" entries.
   * Currently this operation is only defined for the PetscMatrix type.
   */
  virtual void create_submatrix(SparseMatrix<T>& submatrix,
				const std::vector<unsigned int>& rows,
				const std::vector<unsigned int>& cols) const
  {
    this->_get_submatrix(submatrix,
			 rows,
			 cols,
			 false); // false means DO NOT REUSE submatrix
  }

  /**
   * This function is similar to the one above, but it allows you to reuse
   * the existing sparsity pattern of "submatrix" instead of reallocating
   * it again.  This should hopefully be more efficient if you are frequently
   * extracting submatrices of the same size.
   */
  virtual void reinit_submatrix(SparseMatrix<T>& submatrix,
				const std::vector<unsigned int>& rows,
				const std::vector<unsigned int>& cols) const
  {
    this->_get_submatrix(submatrix,
			 rows,
			 cols,
			 true); // true means REUSE submatrix
  }
  
  /**
   * Multiplies the matrix with \p arg and stores the result in \p
   * dest.
   */
  void vector_mult (NumericVector<T>& dest,
		    const NumericVector<T>& arg) const;

  /**
   * Multiplies the matrix with \p arg and adds the result to \p dest.
   */
  void vector_mult_add (NumericVector<T>& dest,
			const NumericVector<T>& arg) const;

  /**
   * Copies the diagonal part of the matrix into \p dest.
   */
  virtual void get_diagonal (NumericVector<T>& dest) const = 0;

  /**
   * Copies the transpose of the matrix into \p dest, which may be
   * *this.
   */
  virtual void get_transpose (SparseMatrix<T>& dest) const = 0;

protected:

  /**
   * Protected implementation of the create_submatrix and reinit_submatrix
   * routines.  Note that this function must be redefined in derived classes
   * for it to work properly!
   */
  virtual void _get_submatrix(SparseMatrix<T>& ,
			      const std::vector<unsigned int>& ,
			      const std::vector<unsigned int>& ,
			      const bool) const
  {
    libMesh::err << "Error! This function is not yet implemented in the base class!"
	          << std::endl;
    libmesh_error();
  }
  
  /**
   * The \p DofMap object associated with this object.
   */
  DofMap const *_dof_map;
  
  /**
   * Flag indicating whether or not the matrix
   * has been initialized.
   */
  bool _is_initialized;
};



//-----------------------------------------------------------------------
// SparseMatrix inline members
template <typename T>
inline
SparseMatrix<T>::SparseMatrix () :
  _dof_map(NULL),
  _is_initialized(false)
{}



template <typename T>
inline
SparseMatrix<T>::~SparseMatrix ()
{}





// Full specialization for Complex datatypes
template <>
inline
void SparseMatrix<Complex>::print(std::ostream& os) const
{
  // std::complex<>::operator<<() is defined, but use this form

  os << "Real part:" << std::endl;
  for (unsigned int i=0; i<this->m(); i++)
    {
      for (unsigned int j=0; j<this->n(); j++)
	os << std::setw(8) << (*this)(i,j).real() << " ";
      os << std::endl;
    }

  os << std::endl << "Imaginary part:" << std::endl;
  for (unsigned int i=0; i<this->m(); i++)
    {
      for (unsigned int j=0; j<this->n(); j++)
	os << std::setw(8) << (*this)(i,j).imag() << " ";
      os << std::endl;
    }
}



// For SGI MIPSpro this implementation must occur after
// the partial specialization of the print() member.
template <typename T>
inline
std::ostream& operator << (std::ostream& os, const SparseMatrix<T>& m)
{
  m.print(os);
  return os;
}



#endif // #ifndef __sparse_matrix_h__
