// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


// Local includes
#include "libmesh/quadrature_clough.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/quadrature_gm.h"
#include "libmesh/quadrature_grid.h"
#include "libmesh/quadrature_jacobi.h"
#include "libmesh/quadrature_monomial.h"
#include "libmesh/quadrature_simpson.h"
#include "libmesh/quadrature_trap.h"
#include "libmesh/quadrature_gauss_lobatto.h"
#include "libmesh/quadrature_conical.h"
#include "libmesh/quadrature_nodal.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/auto_ptr.h" // libmesh_make_unique
#include "libmesh/enum_quadrature_type.h"

namespace libMesh
{



//---------------------------------------------------------------
std::unique_ptr<QBase> QBase::build (const std::string & type,
                                     const unsigned int _dim,
                                     const Order _order)
{
  return QBase::build (Utility::string_to_enum<QuadratureType> (type),
                       _dim,
                       _order);
}



std::unique_ptr<QBase> QBase::build(const QuadratureType _qt,
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
            libMesh::out << "WARNING: Clough quadrature implemented" << std::endl
                         << " up to TWENTYTHIRD order." << std::endl;
          }
#endif

        return libmesh_make_unique<QClough>(_dim, _order);
      }

    case QGAUSS:
      {

#ifdef DEBUG
        if (_order > FORTYTHIRD)
          {
            libMesh::out << "WARNING: Gauss quadrature implemented" << std::endl
                         << " up to FORTYTHIRD order." << std::endl;
          }
#endif

        return libmesh_make_unique<QGauss>(_dim, _order);
      }

    case QJACOBI_1_0:
      {

#ifdef DEBUG
        if (_order > FORTYTHIRD)
          {
            libMesh::out << "WARNING: Jacobi(1,0) quadrature implemented" << std::endl
                         << " up to FORTYTHIRD order." << std::endl;
          }

        if (_dim > 1)
          {
            libMesh::out << "WARNING: Jacobi(1,0) quadrature implemented" << std::endl
                         << " in 1D only." << std::endl;
          }
#endif

        return libmesh_make_unique<QJacobi>(_dim, _order, 1, 0);
      }

    case QJACOBI_2_0:
      {

#ifdef DEBUG
        if (_order > FORTYTHIRD)
          {
            libMesh::out << "WARNING: Jacobi(2,0) quadrature implemented" << std::endl
                         << " up to FORTYTHIRD order." << std::endl;
          }

        if (_dim > 1)
          {
            libMesh::out << "WARNING: Jacobi(2,0) quadrature implemented" << std::endl
                         << " in 1D only." << std::endl;
          }
#endif

        return libmesh_make_unique<QJacobi>(_dim, _order, 2, 0);
      }

    case QSIMPSON:
      {

#ifdef DEBUG
        if (_order > THIRD)
          {
            libMesh::out << "WARNING: Simpson rule provides only" << std::endl
                         << " THIRD order!" << std::endl;
          }
#endif

        return libmesh_make_unique<QSimpson>(_dim);
      }

    case QTRAP:
      {

#ifdef DEBUG
        if (_order > FIRST)
          {
            libMesh::out << "WARNING: Trapezoidal rule provides only" << std::endl
                         << " FIRST order!" << std::endl;
          }
#endif

        return libmesh_make_unique<QTrap>(_dim);
      }

    case QGRID:
      return libmesh_make_unique<QGrid>(_dim, _order);

    case QGRUNDMANN_MOLLER:
      return libmesh_make_unique<QGrundmann_Moller>(_dim, _order);

    case QMONOMIAL:
      return libmesh_make_unique<QMonomial>(_dim, _order);

    case QGAUSS_LOBATTO:
      return libmesh_make_unique<QGaussLobatto>(_dim, _order);

    case QCONICAL:
      return libmesh_make_unique<QConical>(_dim, _order);

    case QNODAL:
      return libmesh_make_unique<QNodal>(_dim, _order);

    default:
      libmesh_error_msg("ERROR: Bad qt=" << _qt);
    }
}

} // namespace libMesh
