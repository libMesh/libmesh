// $Id: distributed_vector.C,v 1.1 2003-02-18 19:43:38 benkirk Exp $

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



#include "mesh_common.h"

// C++ includes
#include <math.h>

// Local Includes
#include "distributed_vector.h"




//--------------------------------------------------------------------------
// DistributedVector methods
template <typename Tp>
Real DistributedVector<Tp>::l1_norm () const
{
  assert (initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  double local_l1 = 0.;

  for (unsigned int i=0; i<local_size(); i++)
    local_l1 += fabs(_values[i]);
  
  double global_l1 = local_l1;


#ifdef HAVE_MPI

  using namespace Mpi;

  MPI_Allreduce (&local_l1, &global_l1, 1,
		 MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
#endif

  return global_l1;
}



template <typename Tp>
Real DistributedVector<Tp>::l2_norm () const
{
  assert (initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  double local_l2 = 0.;
  
  for (unsigned int i=0; i<local_size(); i++)
    local_l2 += _values[i]*_values[i];
  
  double global_l2 = local_l2;


#ifdef HAVE_MPI

  using namespace Mpi;

  MPI_Allreduce (&local_l2, &global_l2, 1,
		 MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
#endif

  return global_l2;
}



template <typename Tp>
Real DistributedVector<Tp>::linfty_norm () const
{
  assert (initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  double local_linfty = 0.;
  
  for (unsigned int i=0; i<local_size(); i++)
    local_linfty  = std::max(local_linfty, fabs(_values[i]));
  
  double global_linfty = local_linfty;


#ifdef HAVE_MPI

  using namespace Mpi;

  MPI_Allreduce (&local_linfty, &global_linfty, 1,
		 MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
#endif

  return global_linfty;
}



template <typename Tp>
NumericVector& DistributedVector<Tp>::operator += (const NumericVector& v)
{
  assert(closed());
  assert (initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);
  
  add(1., v);
  
  return *this;
}



template <typename Tp>
NumericVector& DistributedVector<Tp>::operator -= (const NumericVector& v)
{
  assert(closed());
  assert (initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);
  
  add(-1., v);
  
  return *this;
}



template <typename Tp>
void DistributedVector<Tp>::add_vector (const std::vector<Tp>& v,
					const std::vector<unsigned int>& dof_indices)
{
  assert (!v.empty());
  assert (v.size() == dof_indices.size());
  assert (initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);
  
  for (unsigned int i=0; i<v.size(); i++)
    add (dof_indices[i], v[i]);
}



template <typename Tp>
void DistributedVector<Tp>::add_vector (const NumericVector& V,
					const std::vector<unsigned int>& dof_indices)
{
  assert (V.size() == dof_indices.size());
  assert (initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  for (unsigned int i=0; i<V.size(); i++)
    add (dof_indices[i], V(i));
}



template <typename Tp>
void DistributedVector<Tp>::add (const Tp v)
{
  assert (initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  for (unsigned int i=0; i<local_size(); i++)
    _values[i] += v;
}



template <typename Tp>
void DistributedVector<Tp>::add (const NumericVector& v)
{
  assert (initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  add (1., v);
}



template <typename Tp>
void DistributedVector<Tp>::add (const Tp a, const NumericVector& v)
{
  assert (initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  add(a, v);
}



template <typename Tp>
void DistributedVector<Tp>::scale (const Tp factor)
{
  assert (initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  for (unsigned int i=0; i<local_size(); i++)
    _values[i] *= factor;  
}



template <typename Tp>
NumericVector& 
DistributedVector<Tp>::operator = (const Tp s)
{
  assert (initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  for (unsigned int i=0; i<local_size(); i++)
    _values[i] = s;
  
  return *this;
}



template <typename Tp>
NumericVector&
DistributedVector<Tp>::operator = (const NumericVector& v_in)
{
  assert (initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  const DistributedVector<Tp>& v = reinterpret_cast<const DistributedVector<Tp>&>(v_in);
  
  *this = v;
  
  return *this;
}



template <typename Tp>
DistributedVector<Tp>&
DistributedVector<Tp>::operator = (const DistributedVector<Tp>& v)
{
  assert (initialized());
  assert (v.local_size() == local_size());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  if (v.local_size() == local_size())
    {
      for (unsigned int i=0; i<local_size(); i++)
	_values[i] = v._values[i];
    }
  else
    {
      error();
    }
  
  return *this;
}



template <typename Tp>
NumericVector&
DistributedVector<Tp>::operator = (const std::vector<Tp>& v)
{
  assert (initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  if (v.size() == local_size())
    _values = v;

  else if (v.size() == size())
    for (unsigned int i=first_local_index(); i<last_local_index(); i++)
      _values[i-first_local_index()] = v[i];

  else
    {
      error();
    }

  
  return *this;
}



template <typename Tp>
void DistributedVector<Tp>::localize (NumericVector& v_local_in) const

{
  assert (initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  DistributedVector<Tp>& v_local = reinterpret_cast<DistributedVector<Tp>&>(v_local_in);

  v_local._first_local_index = 0;
  
  v_local._global_size =
    v_local._local_size =
    v_local._last_local_index = size();

  v_local._is_initialized =
    v_local._is_closed = true;
  
  // Call localize on the vector's values.  This will help
  // prevent code duplication
  localize (v_local._values);    
  
#ifndef HAVE_MPI

  assert (local_size() == size());
  
#endif
}



template <typename Tp>
void DistributedVector<Tp>::localize (NumericVector& v_local_in,
				      const std::vector<unsigned int>&) const
{
  assert (initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  // We don't support the send list.  Call the less efficient localize(v_local_in)
  localize (v_local_in);
}



template <typename Tp>
void DistributedVector<Tp>::localize (std::vector<Tp>& v_local) const

{
  assert (initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  v_local.resize(size());
  
  std::fill (v_local.begin(),
	     v_local.end(),
	     0.);

  for (unsigned int i=0; i<local_size(); i++)
    v_local[i+first_local_index()] = _values[i];

#ifdef HAVE_MPI

  using namespace Mpi;

  if (sizeof(Tp) == sizeof(double))
    MPI_Allreduce (&v_local[0], &v_local[0], v_local.size(),
		   MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  else if (sizeof(Tp) == sizeof(float))
    MPI_Allreduce (&v_local[0], &v_local[0], v_local.size(),
		   MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

  else
    {
      error();
    }
    
  
#else

  assert (local_size() == size());
  
#endif  
}



template <typename Tp>
void DistributedVector<Tp>::localize_to_one (std::vector<Tp>& v_local,
					     const unsigned int pid) const
{
  assert (initialized());
  assert (_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  v_local.resize(size());
  
  std::fill (v_local.begin(),
	     v_local.end(),
	     0.);

  for (unsigned int i=0; i<local_size(); i++)
    v_local[i+first_local_index()] = _values[i];

#ifdef HAVE_MPI

  using namespace Mpi;

  if (sizeof(Tp) == sizeof(double))
    MPI_Reduce (&v_local[0], &v_local[0], v_local.size(),
		MPI_DOUBLE, MPI_SUM, pid, MPI_COMM_WORLD);

  else if (sizeof(Tp) == sizeof(float))
    MPI_Reduce (&v_local[0], &v_local[0], v_local.size(),
		   MPI_FLOAT, MPI_SUM, pid, MPI_COMM_WORLD);

  else
    {
      error();
    }
    
  
#else

  assert (local_size() == size());
  
#endif  
}


//--------------------------------------------------------------
// Explicit instantiations
template class DistributedVector<float>;
template class DistributedVector<double>;
