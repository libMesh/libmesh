// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include <cstdlib> // *must* precede <cmath> for proper std:abs() on PGI, Sun Studio CC
#include <cmath>    // for sqrt


// Local Includes
#include "libmesh/libmesh_common.h"
#include "libmesh/discontinuity_measure.h"
#include "libmesh/error_vector.h"
#include "libmesh/fe_base.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/elem.h"
#include "libmesh/system.h"

#include "libmesh/dense_vector.h"
#include "libmesh/tensor_tools.h"


namespace libMesh
{

void
DiscontinuityMeasure::init_context(FEMContext & c)
{
  const unsigned int n_vars = c.n_vars();
  for (unsigned int v=0; v<n_vars; v++)
    {
      // Possibly skip this variable
      if (error_norm.weight(v) == 0.0) continue;

      // FIXME: Need to generalize this to vector-valued elements. [PB]
      FEBase * side_fe = libmesh_nullptr;

      const std::set<unsigned char> & elem_dims =
        c.elem_dimensions();

      for (std::set<unsigned char>::const_iterator dim_it =
             elem_dims.begin(); dim_it != elem_dims.end(); ++dim_it)
        {
          const unsigned char dim = *dim_it;

          fine_context->get_side_fe( v, side_fe, dim );

          // We'll need values on both sides for discontinuity computation
          side_fe->get_phi();
        }
    }
}



void
DiscontinuityMeasure::internal_side_integration ()
{
  const Elem & coarse_elem = coarse_context->get_elem();
  const Elem & fine_elem = fine_context->get_elem();

  FEBase * fe_fine = libmesh_nullptr;
  fine_context->get_side_fe( var, fe_fine, fine_elem.dim() );

  FEBase * fe_coarse = libmesh_nullptr;
  coarse_context->get_side_fe( var, fe_coarse, fine_elem.dim() );

  Real error = 1.e-30;
  unsigned int n_qp = fe_fine->n_quadrature_points();

  std::vector<std::vector<Real>> phi_coarse = fe_coarse->get_phi();
  std::vector<std::vector<Real>> phi_fine = fe_fine->get_phi();
  std::vector<Real> JxW_face = fe_fine->get_JxW();

  for (unsigned int qp=0; qp != n_qp; ++qp)
    {
      // Calculate solution values on fine and coarse elements
      // at this quadrature point
      Number
        u_fine   = fine_context->side_value(var, qp),
        u_coarse = coarse_context->side_value(var, qp);

      // Find the jump in the value
      // at this quadrature point
      const Number jump = u_fine - u_coarse;
      const Real jump2 = TensorTools::norm_sq(jump);
      // Accumulate the jump integral
      error += JxW_face[qp] * jump2;
    }

  // Add the h-weighted jump integral to each error term
  fine_error =
    error * fine_elem.hmax() * error_norm.weight(var);
  coarse_error =
    error * coarse_elem.hmax() * error_norm.weight(var);
}


bool
DiscontinuityMeasure::boundary_side_integration ()
{
  const Elem & fine_elem = fine_context->get_elem();

  FEBase * fe_fine = libmesh_nullptr;
  fine_context->get_side_fe( var, fe_fine, fine_elem.dim() );

  const std::string & var_name =
    fine_context->get_system().variable_name(var);

  std::vector<std::vector<Real>> phi_fine = fe_fine->get_phi();
  std::vector<Real> JxW_face = fe_fine->get_JxW();
  std::vector<Point> qface_point = fe_fine->get_xyz();

  // The reinitialization also recomputes the locations of
  // the quadrature points on the side.  By checking if the
  // first quadrature point on the side is on an essential boundary
  // for a particular variable, we will determine if the whole
  // element is on an essential boundary (assuming quadrature points
  // are strictly contained in the side).
  if (this->_bc_function(fine_context->get_system(),
                         qface_point[0], var_name).first)
    {
      const Real h = fine_elem.hmax();

      // The number of quadrature points
      const unsigned int n_qp = fe_fine->n_quadrature_points();

      // The error contribution from this face
      Real error = 1.e-30;

      // loop over the integration points on the face.
      for (unsigned int qp=0; qp<n_qp; qp++)
        {
          // Value of the imposed essential BC at this quadrature point.
          const std::pair<bool,Real> essential_bc =
            this->_bc_function(fine_context->get_system(), qface_point[qp],
                               var_name);

          // Be sure the BC function still thinks we're on the
          // essential boundary.
          libmesh_assert_equal_to (essential_bc.first, true);

          // The solution value on each point
          Number u_fine = fine_context->side_value(var, qp);

          // The difference between the desired BC and the approximate solution.
          const Number jump = essential_bc.second - u_fine;

          // The flux jump squared.  If using complex numbers,
          // norm_sq(z) returns |z|^2, where |z| is the modulus of z.
          const Real jump2 = TensorTools::norm_sq(jump);

          // Integrate the error on the face.  The error is
          // scaled by an additional power of h, where h is
          // the maximum side length for the element.  This
          // arises in the definition of the indicator.
          error += JxW_face[qp]*jump2;

        } // End quadrature point loop

      fine_error = error*h*error_norm.weight(var);

      return true;
    } // end if side on flux boundary
  return false;
}



void
DiscontinuityMeasure::attach_essential_bc_function (std::pair<bool,Real> fptr(const System & system,
                                                                              const Point & p,
                                                                              const std::string & var_name))
{
  _bc_function = fptr;

  // We may be turning boundary side integration on or off
  if (fptr)
    integrate_boundary_sides = true;
  else
    integrate_boundary_sides = false;
}

} // namespace libMesh
