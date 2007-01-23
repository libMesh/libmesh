// $Id: discontinuity_measure.C,v 1.6 2007-01-23 21:03:15 roystgnr Exp $

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


// C++ includes
#include <algorithm> // for std::fill
#include <cmath>    // for sqrt


// Local Includes
#include "libmesh_common.h"
#include "discontinuity_measure.h"
#include "dof_map.h"
#include "error_vector.h"
#include "fe.h"
#include "fe_interface.h"
#include "libmesh_logging.h"
#include "elem.h"
#include "mesh.h"
#include "system.h"

#include "dense_vector.h"



void
DiscontinuityMeasure::initialize(const System& system,
                                 ErrorVector&,
                                 bool)
{
  // Hang onto the system - we may need it for variable names later.
  my_system = &system;

  // We'll need values for jump computation
  fe_fine->get_phi();
  fe_coarse->get_phi();
}



void
DiscontinuityMeasure::internal_side_integration ()
{
  Real error = 1.e-30;
  unsigned int n_qp = fe_fine->n_quadrature_points();
  unsigned int n_fine_dofs = Ufine.size();
  unsigned int n_coarse_dofs = Ucoarse.size();

  std::vector<std::vector<Real> > phi_coarse = fe_coarse->get_phi();
  std::vector<std::vector<Real> > phi_fine = fe_fine->get_phi();
  std::vector<Real> JxW_face = fe_fine->get_JxW();

  for (unsigned int qp=0; qp != n_qp; ++qp)
    {
      // Calculate solution values on fine and coarse elements
      // at this quadrature point
      Number u_fine, u_coarse;
      for (unsigned int i=0; i != n_coarse_dofs; ++i)
        u_coarse += phi_coarse[i][qp] * Ucoarse(i);

      for (unsigned int i=0; i != n_fine_dofs; ++i)
        u_fine += phi_fine[i][qp] * Ufine(i);
                                
      // Find the jump in the value
      // at this quadrature point
      const Number jump = u_fine - u_coarse;
#ifndef USE_COMPLEX_NUMBERS
      const Real jump2 = jump*jump;
#else
      const Real jump2 = std::norm(jump);
#endif
      // Accumulate the jump integral
      error += JxW_face[qp] * jump2;
    }

  // Add the h-weighted jump integral to each error term
  fine_error =
    error * fine_elem->hmax() * component_scale[var];
  coarse_error =
    error * coarse_elem->hmax() * component_scale[var];
}


bool
DiscontinuityMeasure::boundary_side_integration ()
{
  const std::string &var_name = my_system->variable_name(var);

  std::vector<std::vector<Real> > phi_fine = fe_fine->get_phi();
  std::vector<Real> JxW_face = fe_fine->get_JxW();
  std::vector<Point> qface_point = fe_fine->get_xyz();

  // The reinitialization also recomputes the locations of
  // the quadrature points on the side.  By checking if the
  // first quadrature point on the side is on an essential boundary
  // for a particular variable, we will determine if the whole
  // element is on an essential boundary (assuming quadrature points
  // are strictly contained in the side).
  if (this->_bc_function(*my_system, qface_point[0], var_name).first)
    {
      const Real h = fine_elem->hmax();
		    
      // The number of quadrature points
      const unsigned int n_qp = fe_fine->n_quadrature_points();

      // The error contribution from this face
      Real error = 1.e-30;
		    
      // loop over the integration points on the face.
      for (unsigned int qp=0; qp<n_qp; qp++)
        {
          // Value of the imposed essential BC at this quadrature point.
          const std::pair<bool,Real> essential_bc =
            this->_bc_function(*my_system, qface_point[qp], var_name);

          // Be sure the BC function still thinks we're on the 
          // essential boundary.
          assert (essential_bc.first == true);

          // The solution gradient from each element
          Number u_fine;

          // Compute the solution gradient on element e
          for (unsigned int i=0; i != Ufine.size(); i++)
            u_fine += phi_fine[i][qp] * Ufine(i);

          // The difference between the desired BC and the approximate solution. 
          const Number jump = essential_bc.second - u_fine;

          // The flux jump squared.  If using complex numbers,
          // std::norm(z) returns |z|^2, where |z| is the modulus of z.
#ifndef USE_COMPLEX_NUMBERS
          const Real jump2 = jump*jump;
#else
          const Real jump2 = std::norm(jump);
#endif

          // Integrate the error on the face.  The error is
          // scaled by an additional power of h, where h is
          // the maximum side length for the element.  This
          // arises in the definition of the indicator.
          error += JxW_face[qp]*jump2;			
			
        } // End quadrature point loop

      fine_error = error*h*component_scale[var];

      return true;
    } // end if side on flux boundary
  return false;
}



void
DiscontinuityMeasure::attach_essential_bc_function
  (std::pair<bool,Real> fptr(const System& system,
   const Point& p,
   const std::string& var_name))
{
  _bc_function = fptr;

// We may be turning boundary side integration on or off
  if (fptr)
    integrate_boundary_sides = true;
  else
    integrate_boundary_sides = false;
}
