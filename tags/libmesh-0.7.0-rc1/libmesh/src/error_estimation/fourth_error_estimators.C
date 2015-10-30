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
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES


// C++ includes
#include <algorithm> // for std::fill
#include <cmath>    // for sqrt


// Local Includes
#include "libmesh_common.h"
#include "fourth_error_estimators.h"
#include "error_vector.h"
#include "fe_base.h"
#include "libmesh_logging.h"
#include "elem.h"
#include "system.h"

#include "dense_vector.h"


namespace libMesh
{


void
LaplacianErrorEstimator::initialize(const System&,
                                    ErrorVector&,
                                    bool)
{
  // We'll need second derivatives for Laplacian jump computation
  fe_fine->get_d2phi();
  fe_coarse->get_d2phi();
}



void
LaplacianErrorEstimator::internal_side_integration ()
{
  Real error = 1.e-30;
  unsigned int n_qp = fe_fine->n_quadrature_points();
  unsigned int n_fine_dofs = Ufine.size();
  unsigned int n_coarse_dofs = Ucoarse.size();

  unsigned int dim = fine_elem->dim();

  std::vector<std::vector<RealTensor> > d2phi_coarse = fe_coarse->get_d2phi();
  std::vector<std::vector<RealTensor> > d2phi_fine = fe_fine->get_d2phi();
  std::vector<Real> JxW_face = fe_fine->get_JxW();

  for (unsigned int qp=0; qp != n_qp; ++qp)
    {
      // Calculate solution gradients on fine and coarse elements
      // at this quadrature point
      Number laplacian_fine = 0., laplacian_coarse = 0.;

      for (unsigned int i=0; i != n_coarse_dofs; ++i)
        {
          laplacian_coarse += d2phi_coarse[i][qp](0,0) * Ucoarse(i);
          if (dim > 1)
            laplacian_coarse += d2phi_coarse[i][qp](1,1) * Ucoarse(i);
          if (dim > 2)
            laplacian_coarse += d2phi_coarse[i][qp](2,2) * Ucoarse(i);
        }

      for (unsigned int i=0; i != n_fine_dofs; ++i)
        {
          laplacian_fine += d2phi_fine[i][qp](0,0) * Ufine(i);
          if (dim > 1)
            laplacian_fine += d2phi_fine[i][qp](1,1) * Ufine(i);
          if (dim > 2)
            laplacian_fine += d2phi_fine[i][qp](2,2) * Ufine(i);
        }

                                
      // Find the jump in the Laplacian
      // at this quadrature point
      const Number jump = laplacian_fine - laplacian_coarse;
      const Real jump2 = libmesh_norm(jump);

      // Accumulate the jump integral
      error += JxW_face[qp] * jump2;
    }

  // Add the h-weighted jump integral to each error term
  fine_error =
    error * fine_elem->hmax() * error_norm.weight(var);
  coarse_error =
    error * coarse_elem->hmax() * error_norm.weight(var);
}

} // namespace libMesh

#else // defined (LIBMESH_ENABLE_SECOND_DERIVATIVES)

#include "fourth_error_estimators.h"

namespace libMesh
{

void
LaplacianErrorEstimator::initialize(const System&,
                                    ErrorVector&,
                                    bool)
{
  libMesh::err << "Error: LaplacianErrorEstimator requires second "
                << "derivative support; try configuring libmesh with "
                << "--enable-second" << std::endl;
  libmesh_error();
}


void
LaplacianErrorEstimator::internal_side_integration ()
{
  libMesh::err << "Error: LaplacianErrorEstimator requires second "
                << "derivative support; try configuring libmesh with "
                << "--enable-second" << std::endl;
  libmesh_error();
}

} // namespace libMesh

#endif // defined (LIBMESH_ENABLE_SECOND_DERIVATIVES)


