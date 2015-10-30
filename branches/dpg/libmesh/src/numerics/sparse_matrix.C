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

// Local Includes
#include "dof_map.h"
#include "laspack_matrix.h"
#include "parallel.h"
#include "petsc_matrix.h"
#include "sparse_matrix.h"
#include "trilinos_epetra_matrix.h"
#include "numeric_vector.h"

namespace libMesh
{


//------------------------------------------------------------------
// SparseMatrix Methods

// Full specialization for Real datatypes
template <typename T>
AutoPtr<SparseMatrix<T> >
SparseMatrix<T>::build(const SolverPackage solver_package)
{
  // Build the appropriate vector
  switch (solver_package)
    {


#ifdef LIBMESH_HAVE_LASPACK
    case LASPACK_SOLVERS:
      {
	AutoPtr<SparseMatrix<T> > ap(new LaspackMatrix<T>);
	return ap;
      }
#endif


#ifdef LIBMESH_HAVE_PETSC
    case PETSC_SOLVERS:
      {
	AutoPtr<SparseMatrix<T> > ap(new PetscMatrix<T>);
	return ap;
      }
#endif


#ifdef LIBMESH_HAVE_TRILINOS
    case TRILINOS_SOLVERS:
      {
	AutoPtr<SparseMatrix<T> > ap(new EpetraMatrix<T>);
	return ap;
      }
#endif


    default:
      libMesh::err << "ERROR:  Unrecognized solver package: "
		    << solver_package
		    << std::endl;
      libmesh_error();
    }

  AutoPtr<SparseMatrix<T> > ap(NULL);
  return ap;    
}


template <typename T>
void SparseMatrix<T>::vector_mult (NumericVector<T>& dest,
				   const NumericVector<T>& arg) const
{
  dest.zero();
  this->vector_mult_add(dest,arg);
}



template <typename T>
void SparseMatrix<T>::vector_mult_add (NumericVector<T>& dest,
				       const NumericVector<T>& arg) const
{
  /* This functionality is actually implemented in the \p
     NumericVector class.  */
  dest.add_vector(arg,*this);
}



template <typename T>
void SparseMatrix<T>::zero_rows (std::vector<int> &, T)
{
  /* This functionality isn't implemented or stubbed in every subclass yet */
  libmesh_not_implemented();
}

  

template <typename T>
inline
void SparseMatrix<T>::print(std::ostream& os, const bool sparse) const
{
  parallel_only();

  libmesh_assert (this->initialized());

  if(!this->_dof_map)
  {
    os << std::endl << "Error!  Trying to print a matrix with no dof_map set!" << std::endl << std::endl;
    libmesh_error();
  }  

  // We'll print the matrix from processor 0 to make sure
  // it's serialized properly
  if (libMesh::processor_id() == 0)
    {
      libmesh_assert(this->_dof_map->first_dof() == 0);
      for (unsigned int i=this->_dof_map->first_dof();
           i!=this->_dof_map->end_dof(); ++i)
        {
	  if(sparse)
	    {
	      for (unsigned int j=0; j<this->n(); j++)
		{
		  T c = (*this)(i,j);
		  if (c != static_cast<T>(0.0))
		    {
		      os << i << " " << j << " " << c << std::endl;
		    }
		}
	    }
	  else
	    {
	      for (unsigned int j=0; j<this->n(); j++)
		os << (*this)(i,j) << " ";
	      os << std::endl;
	    }
        }

      std::vector<unsigned int> ibuf, jbuf;
      std::vector<T> cbuf;
      unsigned int currenti = this->_dof_map->end_dof();
      for (unsigned int p=1; p < libMesh::n_processors(); ++p)
        {
          Parallel::receive(p, ibuf);
          Parallel::receive(p, jbuf);
          Parallel::receive(p, cbuf);
          libmesh_assert(ibuf.size() == jbuf.size());
          libmesh_assert(ibuf.size() == cbuf.size());

          if (ibuf.empty())
            continue;
          libmesh_assert(ibuf.front() >= currenti);
          libmesh_assert(ibuf.back() >= ibuf.front());

          unsigned int currentb = 0;
          for (;currenti <= ibuf.back(); ++currenti)
            {
	      if(sparse)
		{
		  for (unsigned int j=0; j<this->n(); j++)
		    {
		      if (currentb < ibuf.size() &&
			  ibuf[currentb] == currenti &&
			  jbuf[currentb] == j)
			{
			  os << currenti << " " << j << " " << cbuf[currentb] << std::endl;
			  currentb++;
			}
		    }
		}
	      else
		{
		  for (unsigned int j=0; j<this->n(); j++)
		    {
		      if (currentb < ibuf.size() &&
			  ibuf[currentb] == currenti &&
			  jbuf[currentb] == j)
			{
			  os << cbuf[currentb] << " ";
			  currentb++;
			}
		      else
			os << static_cast<T>(0.0) << " ";
		    }
		  os << std::endl;
		}
            }
        }
      if(!sparse)
	{
	  for (; currenti != this->m(); ++currenti)
	    {
	      for (unsigned int j=0; j<this->n(); j++)
		os << static_cast<T>(0.0) << " ";
	      os << std::endl;
	    }
	}
    }
  else
    {
      std::vector<unsigned int> ibuf, jbuf;
      std::vector<T> cbuf;

      // We'll assume each processor has access to entire
      // matrix rows, so (*this)(i,j) is valid if i is a local index.
      for (unsigned int i=this->_dof_map->first_dof();
           i!=this->_dof_map->end_dof(); ++i)
        {
          for (unsigned int j=0; j<this->n(); j++)
	    {
              T c = (*this)(i,j);
              if (c != static_cast<T>(0.0))
                {
                  ibuf.push_back(i);
                  jbuf.push_back(j);
                  cbuf.push_back(c);
                }
	    }
        }
      Parallel::send(0,ibuf);
      Parallel::send(0,jbuf);
      Parallel::send(0,cbuf);
    }
}



//------------------------------------------------------------------
// Explicit instantiations
template class SparseMatrix<Number>;

} // namespace libMesh
