
// $Id: parallel.h,v 1.1 2007-06-15 20:33:39 roystgnr Exp $

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


#ifndef __parallel_h__
#define __parallel_h__

// System includes
#include <vector>

// Local includes
#include "libmesh_common.h" // for Real


// ------------------------------------------------------------
// The Parallel namespace is for wrapper functions
// for common general parallel synchronization tasks

namespace Parallel
{
  //-------------------------------------------------------------------
  /**
   * Take a local variable and replace it with the sum of it's value
   * on all processors
   */
  template <typename T>
  void sum(T &r);

  //-------------------------------------------------------------------
  /**
   * Take a vector of local variables and replace each entry with the sum of
   * it's value on all processors
   */
  template <typename T>
  void sum(std::vector<T> &r);



//-----------------------------------------------------------------------
// Parallel members


template <>
void sum(Real &r)
{
#ifdef HAVE_MPI
  if (libMesh::n_processors() > 1)
    {
      Real temp;
      MPI_Allreduce (&r,
                     &temp,
                     1,
                     MPI_REAL,
                     MPI_SUM,
                     libMesh::COMM_WORLD);
      r = temp;
    }
#endif // #ifdef HAVE_MPI
}


template <>
void sum(std::vector<Real> &r)
{
#ifdef HAVE_MPI
  if (libMesh::n_processors() > 1)
    {
      std::vector<Real> temp(r.size());
      MPI_Allreduce (&r[0],
                     &temp[0],
                     r.size(),
                     MPI_REAL,
                     MPI_SUM,
                     libMesh::COMM_WORLD);
      r = temp;
    }
#endif // #ifdef HAVE_MPI
}

}

#endif // #define __parallel_h__
