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



#ifndef __numeric_vector_h__
#define __numeric_vector_h__


// C++ includes
#include <vector>
#include <set>

// Local includes
#include "libmesh_common.h"
#include "auto_ptr.h"
#include "enum_parallel_type.h"
#include "enum_solver_package.h"
#include "reference_counted_object.h"
#include "libmesh.h"


// forward declarations
template <typename T> class NumericVector;
template <typename T> class DenseVector;
template <typename T> class DenseSubVector;
template <typename T> class SparseMatrix;
template <typename T> class ShellMatrix;


/**
 * Numeric vector. Provides a uniform interface
 * to vector storage schemes for different linear
 * algebra libraries.
 *
 * @author Benjamin S. Kirk, 2003
 */
template <typename T>
class NumericVector : public ReferenceCountedObject<NumericVector<T> >
{
public:

  /**
   *  Dummy-Constructor. Dimension=0
   */
  explicit
  NumericVector (const ParallelType type = AUTOMATIC);
    
  /**
   * Constructor. Set dimension to \p n and initialize all elements with zero.
   */
  explicit
  NumericVector (const unsigned int n,
                 const ParallelType type = AUTOMATIC);
    
  /**
   * Constructor. Set local dimension to \p n_local, the global dimension
   * to \p n, and initialize all elements with zero.
   */
  NumericVector (const unsigned n,
		 const unsigned int n_local,
                 const ParallelType type = AUTOMATIC);

  /**
   * Constructor. Set local dimension to \p n_local, the global
   * dimension to \p n, but additionally reserve memory for the
   * indices specified by the \p ghost argument.
   */
  NumericVector (const unsigned int N,
		 const unsigned int n_local,
		 const std::vector<unsigned int>& ghost,
                 const ParallelType type = AUTOMATIC);
    
public:

  /**
   * Destructor, deallocates memory. Made virtual to allow
   * for derived classes to behave properly.
   */
  virtual ~NumericVector ();

  /**
   * Builds a \p NumericVector using the linear solver package specified by
   * \p solver_package
   */
  static AutoPtr<NumericVector<T> >
  build(const SolverPackage solver_package = libMesh::default_solver_package());
  
  /**
   * @returns true if the vector has been initialized,
   * false otherwise.
   */
  virtual bool initialized() const { return _is_initialized; }

  /**
   * @returns the type (SERIAL, PARALLEL, GHOSTED) of the vector.
   */
  ParallelType type() const { return _type; }

  /**
   * @returns true if the vector is closed and ready for
   * computation, false otherwise.
   */
  virtual bool closed() const { return _is_closed; }
  
  /**
   * Call the assemble functions
   */
  virtual void close () = 0; 

  /**
   * @returns the \p NumericVector<T> to a pristine state.
   */
  virtual void clear ();

  /**
   * Set all entries to zero. Equivalent to \p v = 0, but more obvious and
   * faster. 
   */
  virtual void zero () = 0;    

  /**
   * Creates a copy of this vector and returns it in an \p AutoPtr.
   * This must be overloaded in the derived classes.
   */
  virtual AutoPtr<NumericVector<T> > clone () const = 0;
  
  /**
   * Change the dimension of the vector to \p N. The reserved memory for
   * this vector remains unchanged if possible, to make things faster, but
   * this may waste some memory, so take this in the back of your head.
   * However, if \p N==0 all memory is freed, i.e. if you want to resize
   * the vector and release the memory not needed, you have to first call
   * \p init(0) and then \p init(N). This cited behaviour is analogous
   * to that of the STL containers.
   *
   * On \p fast==false, the vector is filled by
   * zeros.
   */
    
  virtual void init (const unsigned int,
		     const unsigned int,
		     const bool = false,
                     const ParallelType = AUTOMATIC) = 0;
  
  /**
   * call init with n_local = N,
   */
  virtual void init (const unsigned int,
		     const bool = false,
                     const ParallelType = AUTOMATIC) = 0;
    
  /**
   * Create a vector that holds tha local indices plus those specified
   * in the \p ghost argument.
   */
  virtual void init (const unsigned int /*N*/,
		     const unsigned int /*n_local*/,
		     const std::vector<unsigned int>& /*ghost*/,
		     const bool /*fast*/ = false,
                     const ParallelType = AUTOMATIC) = 0;

  /**
   * Creates a vector that has the same dimension and storage type as
   * \p other, including ghost dofs.
   */
  virtual void init (const NumericVector<T>& other,
                     const bool fast = false);
    
  /**
   * \f$U(0-N) = s\f$: fill all components.
   */
  virtual NumericVector<T> & operator= (const T s) = 0;
  
  /**
   *  \f$U = V\f$: copy all components.
   */
  virtual NumericVector<T> & operator= (const NumericVector<T> &V) = 0;

  /**
   *  \f$U = V\f$: copy all components.
   */
  virtual NumericVector<T> & operator= (const std::vector<T> &v) = 0;

  /**
   * @returns the minimum element in the vector.
   * In case of complex numbers, this returns the minimum
   * Real part.
   */
  virtual Real min () const = 0;
  
  /**
   * @returns the maximum element in the vector.
   * In case of complex numbers, this returns the maximum
   * Real part.
   */
  virtual Real max () const = 0;
 
  /**
   * returns the sum of the elements in a vector
   */
  virtual T sum() const = 0;

  /**
   * @returns the \f$l_1\f$-norm of the vector, i.e.
   * the sum of the absolute values.
   */
  virtual Real l1_norm () const = 0;

  /**
   * @returns the \f$l_2\f$-norm of the vector, i.e.
   * the square root of the sum of the
   * squares of the elements.
   */
  virtual Real l2_norm () const = 0;

  /**
   * @returns the maximum absolute value of the
   * elements of this vector, which is the
   * \f$l_\infty\f$-norm of a vector.
   */
  virtual Real linfty_norm () const = 0;

  /**
   * @returns the \f$l_1\f$-norm of the vector, i.e.
   * the sum of the absolute values for the specified
   * entries in the vector.
   *
   * Note that the indices must necessary live on this
   * processor.
   */
  virtual Real subset_l1_norm (const std::set<unsigned int> & indices);

  /**
   * @returns the \f$l_2\f$-norm of the vector, i.e.
   * the square root of the sum of the
   * squares of the elements for the specified entries
   * in the vector.
   *
   * Note that the indices must necessary live on this
   * processor.
   */
  virtual Real subset_l2_norm (const std::set<unsigned int> & indices);

  /**
   * @returns the maximum absolute value of the
   * specified entries of this vector, which is the
   * \f$l_\infty\f$-norm of a vector.
   *
   * Note that the indices must necessary live on this
   * processor.
   */
  virtual Real subset_linfty_norm (const std::set<unsigned int> & indices);

  /**
   * @returns dimension of the vector. This
   * function was formerly called \p n(), but
   * was renamed to get the \p NumericVector<T> class
   * closer to the C++ standard library's
   * \p std::vector container.
   */
  virtual unsigned int size () const = 0;

  /**
   * @returns the local size of the vector
   * (index_stop-index_start).
   * In ghost cell mode, this does *not* include the ghost cells.
   */
  virtual unsigned int local_size() const = 0;

  /**
   * @returns the index of the first vector element
   * actually stored on this processor.  Hint: the
   * minimum for this index is \p 0.
   */
  virtual unsigned int first_local_index() const = 0;

  /**
   * @returns the index+1 of the last vector element
   * actually stored on this processor.  Hint: the
   * maximum for this index is \p size().
   */
  virtual unsigned int last_local_index() const = 0;
    
  /**
   * Access components, returns \p U(i).
   */
  virtual T operator() (const unsigned int i) const = 0;
    
  /**
   * Addition operator.
   * Fast equivalent to \p U.add(1, V).
   */
  virtual NumericVector<T> & operator += (const NumericVector<T> &V) = 0;

  /**
   * Subtraction operator.
   * Fast equivalent to \p U.add(-1, V).
   */
  virtual NumericVector<T> & operator -= (const NumericVector<T> &V) = 0;
    
  /**
   * v(i) = value
   */
  virtual void set (const unsigned int i, const T value) = 0;
    
  /**
   * v(i) += value
   */
  virtual void add (const unsigned int i, const T value) = 0;
    
  /**
   * \f$U(0-LIBMESH_DIM)+=s\f$.
   * Addition of \p s to all components. Note
   * that \p s is a scalar and not a vector.
   */
  virtual void add (const T s) = 0;
    
  /**
   * \f$U+=V\f$:
   * Simple vector addition, equal to the
   * \p operator +=.
   */
  virtual void add (const NumericVector<T>& V) = 0;

  /**
   * \f$U+=a*V\f$.
   * Simple vector addition, equal to the
   * \p operator +=.
   */
  virtual void add (const T a, const NumericVector<T>& v) = 0;

  /**
   * \f$ U+=v \f$ where v is a DenseVector<T> 
   * and you
   * want to specify WHERE to add it
   */
  virtual void add_vector (const std::vector<T>& v,
			   const std::vector<unsigned int>& dof_indices) = 0;

  /**
   * \f$U+=V\f$, where U and V are type 
   * NumericVector<T> and you
   * want to specify WHERE to add
   * the NumericVector<T> V 
   */
  virtual void add_vector (const NumericVector<T>& V,
			   const std::vector<unsigned int>& dof_indices) = 0;

  /**
   * \f$U+=A*V\f$, add the product of a \p SparseMatrix \p A
   * and a \p NumericVector \p V to \p this, where \p this=U.
   */
  virtual void add_vector (const NumericVector<T>&,
			   const SparseMatrix<T>&) = 0;
      
  /**
   * \f$U+=A*V\f$, add the product of a \p ShellMatrix \p A
   * and a \p NumericVector \p V to \p this, where \p this=U.
   */
  void add_vector (const NumericVector<T>& v,
		   const ShellMatrix<T>& a);
      
  /**
   * \f$ U+=V \f$ where U and V are type 
   * DenseVector<T> and you
   * want to specify WHERE to add
   * the DenseVector<T> V 
   */
  virtual void add_vector (const DenseVector<T>& V,
			   const std::vector<unsigned int>& dof_indices) = 0;

  /**
   * \f$ U=v \f$ where v is a std::vector<T> 
   * and you want to specify WHERE to insert it
   */
  virtual void insert (const std::vector<T>& v,
		       const std::vector<unsigned int>& dof_indices) = 0;

  /**
   * \f$U=V\f$, where U and V are type 
   * NumericVector<T> and you
   * want to specify WHERE to insert
   * the NumericVector<T> V 
   */
  virtual void insert (const NumericVector<T>& V,
		       const std::vector<unsigned int>& dof_indices) = 0;
      
  /**
   * \f$ U=V \f$ where U and V are type 
   * DenseVector<T> and you
   * want to specify WHERE to insert
   * the DenseVector<T> V 
   */
  virtual void insert (const DenseVector<T>& V,
		       const std::vector<unsigned int>& dof_indices) = 0;

  /**
   * \f$ U=V \f$ where V is a
   * DenseSubVector<T> and you
   * want to specify WHERE to insert it
   */
  virtual void insert (const DenseSubVector<T>& V,
		       const std::vector<unsigned int>& dof_indices) = 0;
    
  /**
   * Scale each element of the
   * vector by the given factor.
   */
  virtual void scale (const T factor) = 0;

  /**
   * Computes the dot product, p = U.V
   */
  virtual T dot(const NumericVector<T>&) const = 0;
  
  /**
   * Creates a copy of the global vector in the
   * local vector \p v_local.
   */
  virtual void localize (std::vector<T>& v_local) const = 0;

  /**
   * Same, but fills a \p NumericVector<T> instead of
   * a \p std::vector.
   */
  virtual void localize (NumericVector<T>& v_local) const = 0;

  /**
   * Creates a local vector \p v_local containing
   * only information relevant to this processor, as
   * defined by the \p send_list.
   */
  virtual void localize (NumericVector<T>& v_local,
			 const std::vector<unsigned int>& send_list) const = 0;
  
  /**
   * Updates a local vector with selected values from neighboring
   * processors, as defined by \p send_list.
   */
  virtual void localize (const unsigned int first_local_idx,
			 const unsigned int last_local_idx,
			 const std::vector<unsigned int>& send_list) = 0;

  /**
   * Creates a local copy of the global vector in
   * \p v_local only on processor \p proc_id.  By
   * default the data is sent to processor 0.  This method
   * is useful for outputting data from one processor.
   */
  virtual void localize_to_one (std::vector<T>& v_local,
				const unsigned int proc_id=0) const = 0;
    
  /**
   * @returns \p -1 when \p this is equivalent to \p other_vector,
   * up to the given \p threshold.  When differences occur,
   * the return value contains the first index where
   * the difference exceeded the threshold.  When
   * no threshold is given, the \p libMesh \p TOLERANCE
   * is used.
   */
  virtual int compare (const NumericVector<T> &other_vector,
		       const Real threshold = TOLERANCE) const;

  /**
   * Computes the pointwise (i.e. component-wise) product of \p vec1
   * and \p vec2 and stores the result in \p *this.
   */
  virtual void pointwise_mult (const NumericVector<T>& vec1,
			       const NumericVector<T>& vec2) = 0;

  /**
   * Prints the local contents of the vector to the screen.
   */
  virtual void print(std::ostream& os=std::cout) const;

  /**
   * Prints the global contents of the vector to the screen.
   */
  virtual void print_global(std::ostream& os=std::cout) const;

  /**
   * Same as above but allows you to use stream syntax.
   */
  friend std::ostream& operator << (std::ostream& os, const NumericVector<T>& v)
  {
    v.print_global(os);
    return os;
  }
  
  /**
   * Print the contents of the matrix in Matlab's
   * sparse matrix format. Optionally prints the
   * matrix to the file named \p name.  If \p name
   * is not specified it is dumped to the screen.
   */
  virtual void print_matlab(const std::string name="NULL") const
  {
    std::cerr << "ERROR: Not Implemented in base class yet!" << std::endl;
    std::cerr << "ERROR writing MATLAB file " << name << std::endl;
    libmesh_error();
  }

  /**
   * Creates the subvector "subvector" from the indices in the
   * "rows" array.  Similar to the create_submatrix routine for
   * the SparseMatrix class, it is currently only implemented for
   * PetscVectors.
   */
  virtual void create_subvector(NumericVector<T>& ,
				const std::vector<unsigned int>& ) const
  {
    std::cerr << "ERROR: Not Implemented in base class yet!" << std::endl;
    libmesh_error();
  }
    
protected:
  
  /**
   * Flag to see if the Numeric
   * assemble routines have been called yet
   */
  bool _is_closed;
  
  /**
   * Flag to tell if init 
   * has been called yet
   */
  bool _is_initialized;

  /**
   * Type of vector
   */
  ParallelType _type;
};


/*----------------------- Inline functions ----------------------------------*/



template <typename T>
inline
NumericVector<T>::NumericVector (const ParallelType type) :
  _is_closed(false),
  _is_initialized(false),
  _type(type)
{
}



template <typename T>
inline
NumericVector<T>::NumericVector (const unsigned int /*n*/,
                                 const ParallelType type) :
  _is_closed(false),
  _is_initialized(false),
  _type(type)
{
  libmesh_error(); // Abstract base class!
  // init(n, n, false, type);
}



template <typename T>
inline
NumericVector<T>::NumericVector (const unsigned int /*n*/,
				 const unsigned int /*n_local*/,
                                 const ParallelType type) :
  _is_closed(false),
  _is_initialized(false),
  _type(type)
{
  libmesh_error(); // Abstract base class!
  // init(n, n_local, false, type);
}



template <typename T>
inline
NumericVector<T>::NumericVector (const unsigned int /*n*/,
				 const unsigned int /*n_local*/,
				 const std::vector<unsigned int>& /*ghost*/,
                                 const ParallelType type) :
  _is_closed(false),
  _is_initialized(false),
  _type(type)
{
  libmesh_error(); // Abstract base class!
  // init(n, n_local, ghost, false, type);
}



template <typename T>
inline
NumericVector<T>::~NumericVector ()
{
  clear ();
}



// These should be pure virtual, not bugs waiting to happen - RHS
/*
template <typename T>
inline
NumericVector<T> & NumericVector<T>::operator= (const T) 
{
  //  libmesh_error();

  return *this;
}



template <typename T>
inline
NumericVector<T> & NumericVector<T>::operator= (const NumericVector<T>&) 
{
  //  libmesh_error();

  return *this;
}



template <typename T>
inline
NumericVector<T> & NumericVector<T>::operator= (const std::vector<T>&) 
{
  //  libmesh_error();

  return *this;
}
*/



template <typename T>
inline
void NumericVector<T>::clear ()
{
  _is_closed      = false;
  _is_initialized = false;
}



// Full specialization of the print() member for complex
// variables.  This must precede the non-specialized
// version, at least according to icc v7.1
template <>
inline
void NumericVector<Complex>::print(std::ostream& os) const
{
  libmesh_assert (this->initialized());
  os << "Size\tglobal =  " << this->size()
     << "\t\tlocal =  " << this->local_size() << std::endl;
  
  // std::complex<>::operator<<() is defined, but use this form
  os << "#\tReal part\t\tImaginary part" << std::endl;
  for (unsigned int i=this->first_local_index(); i<this->last_local_index(); i++)
    os << i << "\t" 
       << (*this)(i).real() << "\t\t" 
       << (*this)(i).imag() << std::endl;
}



template <typename T>
inline
void NumericVector<T>::print(std::ostream& os) const
{
  libmesh_assert (this->initialized());
  os << "Size\tglobal =  " << this->size()
     << "\t\tlocal =  " << this->local_size() << std::endl;

  os << "#\tValue" << std::endl;
  for (unsigned int i=this->first_local_index(); i<this->last_local_index(); i++)
    os << i << "\t" << (*this)(i) << std::endl;
}



template <>
inline
void NumericVector<Complex>::print_global(std::ostream& os) const
{
  libmesh_assert (this->initialized());

  std::vector<Complex> v(this->size());
  this->localize(v);

  // Right now we only want one copy of the output
  if (libMesh::processor_id())
    return;
  
  os << "Size\tglobal =  " << this->size() << std::endl;
  os << "#\tReal part\t\tImaginary part" << std::endl;
  for (unsigned int i=0; i!=v.size(); i++)
    os << i << "\t" 
       << v[i].real() << "\t\t" 
       << v[i].imag() << std::endl;
}


template <typename T>
inline
void NumericVector<T>::print_global(std::ostream& os) const
{
  libmesh_assert (this->initialized());

  std::vector<T> v(this->size());
  this->localize(v);

  // Right now we only want one copy of the output
  if (libMesh::processor_id())
    return;

  os << "Size\tglobal =  " << this->size() << std::endl;
  os << "#\tValue" << std::endl;
  for (unsigned int i=0; i!=v.size(); i++)
    os << i << "\t" << v[i] << std::endl;
}



#endif  // #ifdef __numeric_vector_h__
