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


#include "libmesh_config.h"
#ifdef LIBMESH_HAVE_SLEPC

// C++ includes

// Local Includes
#include "eigen_solver.h"
#include "slepc_eigen_solver.h"
#include "auto_ptr.h"


//------------------------------------------------------------------
// EigenSolver members
template <typename T>
AutoPtr<EigenSolver<T> >
EigenSolver<T>::build(const SolverPackage solver_package)
{
  // Build the appropriate solver
  switch (solver_package)
    {



#ifdef LIBMESH_HAVE_SLEPC
	case SLEPC_SOLVERS:
      {
	AutoPtr<EigenSolver<T> > ap(new SlepcEigenSolver<T>);
	return ap;
      }
#endif


    default:
      std::cerr << "ERROR:  Unrecognized eigen solver package: "
		<< solver_package
		<< std::endl;
      libmesh_error();
    }
    
  AutoPtr<EigenSolver<T> > ap(NULL);
  return ap;    
}



//------------------------------------------------------------------
// Explicit instantiations
template class EigenSolver<Number>;


#endif // LIBMESH_HAVE_SLEPC
