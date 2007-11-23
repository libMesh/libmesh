// $Id: quadrature_build.C,v 1.11 2007-10-21 20:48:53 benkirk Exp $

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


// C++ includes



// Local includes
#include "quadrature_clough.h"
#include "quadrature_gauss.h"
#include "quadrature_jacobi.h"
#include "quadrature_simpson.h"
#include "quadrature_trap.h"



AutoPtr<QBase> QBase::build(const QuadratureType _qt,
			    const unsigned int _dim,
			    const Order _order)
{
  switch (_qt)
    {
      
    case QCLOUGH:
      {
#ifdef DEBUG
	if (_order > TWENTYTHIRD)
	  {
	    std::cout << "WARNING: Clough quadrature implemented" << std::endl
		      << " up to TWENTYTHIRD order." << std::endl;
	  }
#endif

	AutoPtr<QBase> ap(new QClough(_dim, _order));
	return ap;
      }

    case QGAUSS:
      {

#ifdef DEBUG
	if (_order > FORTYTHIRD)
	  {
	    std::cout << "WARNING: Gauss quadrature implemented" << std::endl
		      << " up to FORTYTHIRD order." << std::endl;
	  }
#endif

	AutoPtr<QBase> ap(new QGauss(_dim, _order));
	return ap;
      }

    case QJACOBI_1_0:
      {

#ifdef DEBUG
	if (_order > TWENTYTHIRD)
	  {
	    std::cout << "WARNING: Jacobi(1,0) quadrature implemented" << std::endl
		      << " up to TWENTYTHIRD order." << std::endl;
	  }

	if (_dim > 1)
	  {
	    std::cout << "WARNING: Jacobi(1,0) quadrature implemented" << std::endl
		      << " in 1D only." << std::endl;
	  }
#endif

	AutoPtr<QBase> ap(new QJacobi(_dim, _order, 1, 0));
	return ap;
      }

    case QJACOBI_2_0:
      {

#ifdef DEBUG
	if (_order > TWENTYTHIRD)
	  {
	    std::cout << "WARNING: Jacobi(2,0) quadrature implemented" << std::endl
		      << " up to TWENTYTHIRD order." << std::endl;
	  }

	if (_dim > 1)
	  {
	    std::cout << "WARNING: Jacobi(2,0) quadrature implemented" << std::endl
		      << " in 1D only." << std::endl;
	  }
#endif

	AutoPtr<QBase> ap(new QJacobi(_dim, _order, 2, 0));
	return ap;
      }

    case QSIMPSON:
      {

#ifdef DEBUG
	if (_order > THIRD)
	  {
	    std::cout << "WARNING: Simpson rule provides only" << std::endl
		      << " THIRD order!" << std::endl;
	  }
#endif

	AutoPtr<QBase> ap(new QSimpson(_dim));
	return ap;
      }

    case QTRAP:
      {

#ifdef DEBUG
	if (_order > FIRST)
	  {
	    std::cout << "WARNING: Trapezoidal rule provides only" << std::endl
		      << " FIRST order!" << std::endl;
	  }
#endif

	AutoPtr<QBase> ap(new QTrap(_dim));
	return ap;
      }


    default:
      { 
	std::cerr << "ERROR: Bad qt=" << _qt << std::endl;
	error();
      }
    }


  error();
  AutoPtr<QBase> ap(NULL);
  return ap;
}

