// $Id: distributed_vector.h,v 1.11 2005-12-28 13:47:10 spetersen Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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



#include "libmesh_common.h"



#ifndef __distributed_vector_h__
#define __distributed_vector_h__




// C++ includes
#include <vector>
#include <algorithm>

#ifdef HAVE_MPI
# include <mpi.h>
#endif

// Local includes
#include "numeric_vector.h"



/**
 * Distributed vector. Provides an interface for simple
 * parallel, distributed vectors. Offers some collective
 * communication capabilities.  Note that the class will
 * sill function without MPI, but only on one processor.
 * This lets us keep the parallel details behind the scenes.
 *
 * @author Benjamin S. Kirk, 2003
 */

template <typename T>
class DistributedVector : public NumericVector<T>
{
public:

  /**
   *  Dummy-Constructor. Dimension=0
   */
  DistributedVector ();
  
  /**
   * Constructor. Set dimension to \p n and initialize all elements with zero.
   */
  DistributedVector (const unsigned int n);
    
  /**
   * Constructor. Set local dimension to \p n_local, the global dimension
   * to \p n, and initialize all elements with zero.
   */
  DistributedVector (const unsigned int n,
		     const unsigned int n_local);
    
  /**
   * Destructor, deallocates memory. Made virtual to allow
   * for derived classes to behave properly.
   */
  ~DistributedVector ();

  /**
   * Call the assemble functions
   */
  void close (); 

  /**
   * @returns the \p DistributedVector to a pristine state.
   */
  void clear ();
  
  /**
   * Set all entries to zero. Equivalent to \p v = 0, but more obvious and
   * faster. 
   */
  void zero ();    

  /**
   * Creates a copy of this vector and returns it in an \p AutoPtr.
   */
  AutoPtr<NumericVector<T> > clone () const;

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
  void init (const unsigned int N,
	     const unsigned int n_local,
	     const bool         fast=false);
    
  /**
   * call init with n_local = N,
   */
  void init (const unsigned int N,
	     const bool         fast=false);
    
  /**
   * \f$U(0-N) = s\f$: fill all components.
   */
  NumericVector<T> & operator= (const T s);
    
  /**
   *  \f$U = V\f$: copy all components.
   */
  NumericVector<T> & operator= (const NumericVector<T> &V);

  /**
   *  \f$U = V\f$: copy all components.
   */
  DistributedVector<T> & operator= (const DistributedVector<T> &V);

  /**
   *  \f$U = V\f$: copy all components.
   */
  NumericVector<T> & operator= (const std::vector<T> &v);
  
  /**
   * @returns the minimum element in the vector.
   * In case of complex numbers, this returns the minimum
   * Real part.
   */
  Real min () const;
  
  /**
   * @returns the maximum element in the vector.
   * In case of complex numbers, this returns the maximum
   * Real part.
   */
  Real max () const;
  
  /**
   * @returns the \f$l_1\f$-norm of the vector, i.e.
   * the sum of the absolute values.
   */
  Real l1_norm () const;

  /**
   * @returns the \f$l_2\f$-norm of the vector, i.e.
   * the square root of the sum of the
   * squares of the elements.
   */
  Real l2_norm () const;

  /**
   * @returns the maximum absolute value of the
   * elements of this vector, which is the
   * \f$l_\infty\f$-norm of a vector.
   */
  Real linfty_norm () const;

  /**
   * @returns dimension of the vector. This
   * function was formerly called \p n(), but
   * was renamed to get the \p DistributedVector class
   * closer to the C++ standard library's
   * \p std::vector container.
   */
  unsigned int size () const;

  /**
   * @returns the local size of the vector
   * (index_stop-index_start)
   */
  unsigned int local_size() const;

  /**
   * @returns the index of the first vector element
   * actually stored on this processor
   */
  unsigned int first_local_index() const;

  /**
   * @returns the index of the last vector element
   * actually stored on this processor
   */
  unsigned int last_local_index() const;
    
  /**
   * Access components, returns \p U(i).
   */
  T operator() (const unsigned int i) const;
    
  /**
   * Addition operator.
   * Fast equivalent to \p U.add(1, V).
   */
  NumericVector<T> & operator += (const NumericVector<T> &V);

  /**
   * Subtraction operator.
   * Fast equivalent to \p U.add(-1, V).
   */
  NumericVector<T> & operator -= (const NumericVector<T> &V);
    
  /**
   * v(i) = value
   */
  void set (const unsigned int i, const T value);
    
  /**
   * v(i) += value
   */
  void add (const unsigned int i, const T value);
    
  /**
   * \f$U(0-DIM)+=s\f$.
   * Addition of \p s to all components. Note
   * that \p s is a scalar and not a vector.
   */
  void add (const T s);
    
  /**
   * \f$U+=V\f$.
   * Simple vector addition, equal to the
   * \p operator +=.
   */
  void add (const NumericVector<T>& V);

  /**
   * \f$U+=a*V\f$.
   * Simple vector addition, equal to the
   * \p operator +=.
   */
  void add (const T a, const NumericVector<T>& v);
  
  /**
   * \f$U+=v\f$ where v is a \p std::vector<T> 
   * and you
   * want to specify WHERE to add it
   */
  void add_vector (const std::vector<T>& v,
		   const std::vector<unsigned int>& dof_indices);

  /**
   * \f$U+=V\f$ where U and V are type 
   * \p NumericVector<T> and you
   * want to specify WHERE to add
   * the \p NumericVector<T> V 
   */
  void add_vector (const NumericVector<T>& V,
		   const std::vector<unsigned int>& dof_indices);
  
  /**
   * \f$U+=A*V\f$.
   * Add the product of a Sparse matrix \p A
   * and a Numeric vector \p V to this Numeric vector.
   * @e Not @e implemented.
   */
  void add_vector (const NumericVector<T>&,
		   const SparseMatrix<T>&)
  { error(); }
  
  /**
   * \f$U+=V\f$ where U and V are type 
   * \p DenseVector<T> and you
   * want to specify WHERE to add
   * the \p DenseVector<T> V 
   */
  void add_vector (const DenseVector<T>& V,
		   const std::vector<unsigned int>& dof_indices);
  
  /**
   * \f$ U=v \f$ where v is a DenseVector<T> 
   * and you want to specify WHERE to insert it
   */
  virtual void insert (const std::vector<T>& v,
		       const std::vector<unsigned int>& dof_indices);

  /**
   * \f$U=V\f$, where U and V are type 
   * NumericVector<T> and you
   * want to specify WHERE to insert
   * the NumericVector<T> V 
   */
  virtual void insert (const NumericVector<T>& V,
		       const std::vector<unsigned int>& dof_indices);
      
  /**
   * \f$ U+=V \f$ where U and V are type 
   * DenseVector<T> and you
   * want to specify WHERE to insert
   * the DenseVector<T> V 
   */
  virtual void insert (const DenseVector<T>& V,
		       const std::vector<unsigned int>& dof_indices);
    
  /**
   * Scale each element of the
   * vector by the given factor.
   */
  void scale (const T factor);

  /**
   * Computes the dot product, p = U.V
   */
  virtual Number dot(const NumericVector<T>& V) const;
  
  /**
   * Creates a copy of the global vector in the
   * local vector \p v_local.
   */
  void localize (std::vector<T>& v_local) const;

  /**
   * Same, but fills a \p NumericVector<T> instead of
   * a \p std::vector.
   */
  void localize (NumericVector<T>& v_local) const;

  /**
   * Creates a local vector \p v_local containing
   * only information relevant to this processor, as
   * defined by the \p send_list.
   */
  void localize (NumericVector<T>& v_local,
		 const std::vector<unsigned int>& send_list) const;
  
  /**
   * Updates a local vector with selected values from neighboring
   * processors, as defined by \p send_list.
   */
  void localize (const unsigned int first_local_idx,
		 const unsigned int last_local_idx,
		 const std::vector<unsigned int>& send_list);

  /**
   * Creates a local copy of the global vector in
   * \p v_local only on processor \p proc_id.  By
   * default the data is sent to processor 0.  This method
   * is useful for outputting data from one processor.
   */
  void localize_to_one (std::vector<T>& v_local,
			const unsigned int proc_id=0) const;
    
private:

  /**
   * Actual vector datatype
   * to hold vector entries
   */
  std::vector<T> _values;

  /**
   * The global vector size
   */
  unsigned int _global_size;

  /**
   * The local vector size
   */
  unsigned int _local_size;

  /**
   * The first component stored locally
   */
  unsigned int _first_local_index;

  /**
   * The last component (+1) stored locally
   */
  unsigned int _last_local_index;
};



//--------------------------------------------------------------------------
// DistributedVector inline methods
template <typename T>
inline
DistributedVector<T>::DistributedVector () :
  _global_size      (0),
  _local_size       (0),
  _first_local_index(0),
  _last_local_index (0)
{
}



template <typename T>
inline
DistributedVector<T>::DistributedVector (const unsigned int n)
{
  this->init(n, n, false);
}



template <typename T>
inline
DistributedVector<T>::DistributedVector (const unsigned int n,
					 const unsigned int n_local)
{
  this->init(n, n_local, false);
}



template <typename T>
inline
DistributedVector<T>::~DistributedVector ()
{
  this->clear ();
}



template <typename T>
inline
void DistributedVector<T>::init (const unsigned int n,
				 const unsigned int n_local,
				 const bool fast)
{
  assert (n_local <= n);

  // Clear the data structures if already initialized
  if (this->initialized())
    this->clear();
    
  // Initialize data structures
  _values.resize(n_local);
  _local_size  = n_local;
  _global_size = n;

  _first_local_index = 0;
  
#ifdef HAVE_MPI

  int n_proc=0, proc_id=0;
  
  MPI_Comm_rank (libMesh::COMM_WORLD, &proc_id);
  MPI_Comm_size (libMesh::COMM_WORLD, &n_proc);
  
  std::vector<int> local_sizes     (n_proc, 0);
  std::vector<int> local_sizes_send(n_proc, 0);
  
  local_sizes_send[proc_id] = n_local;

  MPI_Allreduce (&local_sizes_send[0], &local_sizes[0], local_sizes.size(),
		 MPI_INT, MPI_SUM, libMesh::COMM_WORLD);
  
  // _first_local_index is the sum of _local_size
  // for all processor ids less than ours
  for (int p=0; p<proc_id; p++)
    _first_local_index += local_sizes[p];


#  ifdef DEBUG
  // Make sure all the local sizes sum up to the global
  // size, otherwise there is big trouble!
  int sum=0;

  for (int p=0; p<n_proc; p++)
    sum += local_sizes[p];

  assert (sum == static_cast<int>(n));
  
#  endif
  
#else
  
  // No other options without MPI!
  if (n != n_local)
    {
      std::cerr << "ERROR:  MPI is required for n != n_local!"
		<< std::endl;
      error();
    }
  
#endif

  _last_local_index = _first_local_index + n_local;
  
  // Set the initialized flag
  this->_is_initialized = true;

  // Zero the components unless directed otherwise
  if (!fast)
    this->zero();
}



template <typename T>
inline
void DistributedVector<T>::init (const unsigned int n,
				 const bool fast)
{
  this->init(n,n,fast);
}



template <typename T>
inline
void DistributedVector<T>::close ()
{
  assert (this->initialized());
  
  this->_is_closed = true;
}



template <typename T>
inline
void DistributedVector<T>::clear ()
{
  _values.clear();
  
  _global_size =
    _local_size =
    _first_local_index =
    _last_local_index = 0;
  
  
  this->_is_closed = this->_is_initialized = false;
}



template <typename T>
inline
void DistributedVector<T>::zero ()
{
  assert (this->initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);
  
  std::fill (_values.begin(),
	     _values.end(),
	     0.);
}



template <typename T>
inline
AutoPtr<NumericVector<T> > DistributedVector<T>::clone () const
{
  AutoPtr<NumericVector<T> > cloned_vector (new DistributedVector<T>);

  *cloned_vector = *this;

  return cloned_vector;
}



template <typename T>
inline
unsigned int DistributedVector<T>::size () const
{
  assert (this->initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  return _global_size;
}



template <typename T>
inline
unsigned int DistributedVector<T>::local_size () const
{
  assert (this->initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  return _local_size;
}



template <typename T>
inline
unsigned int DistributedVector<T>::first_local_index () const
{
  assert (this->initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  return _first_local_index;
}



template <typename T>
inline
unsigned int DistributedVector<T>::last_local_index () const
{
  assert (this->initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  return _last_local_index;
}



template <typename T>
inline
T DistributedVector<T>::operator() (const unsigned int i) const
{
  assert (this->initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);
  assert ( ((i >= first_local_index()) &&
	    (i <  last_local_index())) );

  return _values[i - _first_local_index];
}



template <typename T>
inline
void DistributedVector<T>::set (const unsigned int i, const T value)
{
  assert (this->initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);
  assert (i<size());
  assert (i-first_local_index() < local_size());
  
  _values[i - _first_local_index] = value;
}



template <typename T>
inline
void DistributedVector<T>::add (const unsigned int i, const T value)
{
  assert (this->initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);
  assert (i<size());
  assert (i-first_local_index() < local_size());
  
  _values[i - _first_local_index] += value;
}



template <typename T>
inline
Real DistributedVector<T>::min () const
{
  assert (this->initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  double local_min = static_cast<double>(*std::min(_values.begin(),
						   _values.end()));
  
  double global_min = local_min;


#ifdef HAVE_MPI

  MPI_Allreduce (&local_min, &global_min, 1,
		 MPI_DOUBLE, MPI_MIN, libMesh::COMM_WORLD);

#endif

  return global_min;
}



template <typename T>
inline
Real DistributedVector<T>::max() const
{
  assert (this->initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  double local_max = static_cast<double>(*std::max(_values.begin(),
						   _values.end()));

  double global_max = local_max;


#ifdef HAVE_MPI

  MPI_Allreduce (&local_max, &global_max, 1,
		 MPI_DOUBLE, MPI_MAX, libMesh::COMM_WORLD);
  
#endif

  return global_max;
}


#endif  // #ifdef __distributed_vector_h__
