// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include <cstdlib> // *must* precede <cmath> for proper std:abs() on PGI, Sun Studio CC
#include <cmath> // for std::sqrt


// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature.h"
#include "libmesh/elem.h"
#include "libmesh/libmesh_logging.h"

namespace libMesh
{




//-------------------------------------------------------
// Full specialization for 0D when this is a useless method
template <>
void FEXYZ<0>::reinit(const Elem *,
                      const unsigned int,
                      const Real,
                      const std::vector<Point> * const,
                      const std::vector<Real> * const)
{
  libmesh_error_msg("ERROR: This method only makes sense for 2D/3D elements!");
}


//-------------------------------------------------------
// Full specialization for 1D when this is a useless method
template <>
void FEXYZ<1>::reinit(const Elem *,
                      const unsigned int,
                      const Real,
                      const std::vector<Point> * const,
                      const std::vector<Real> * const)
{
  libmesh_error_msg("ERROR: This method only makes sense for 2D/3D elements!");
}


//-------------------------------------------------------
// Method for 2D, 3D
template <unsigned int Dim>
void FEXYZ<Dim>::reinit(const Elem * elem,
                        const unsigned int s,
                        const Real,
                        const std::vector<Point> * const pts,
                        const std::vector<Real> * const weights)
{
  libmesh_assert(elem);
  libmesh_assert (this->qrule != libmesh_nullptr || pts != libmesh_nullptr);
  // We don't do this for 1D elements!
  libmesh_assert_not_equal_to (Dim, 1);

  // Build the side of interest
  const std::unique_ptr<const Elem> side(elem->build_side_ptr(s));

  // Initialize the shape functions at the user-specified
  // points
  if (pts != libmesh_nullptr)
    {
      // We can't get away without recomputing shape functions next
      // time
      this->shapes_on_quadrature = false;

      // Set the element type
      this->elem_type = elem->type();

      // Initialize the face shape functions
      this->_fe_map->template init_face_shape_functions<Dim>(*pts,  side.get());
      if (weights != libmesh_nullptr)
        {
          this->compute_face_values (elem, side.get(), *weights);
        }
      else
        {
          std::vector<Real> dummy_weights (pts->size(), 1.);
          // Compute data on the face for integration
          this->compute_face_values (elem, side.get(), dummy_weights);
        }
    }
  else
    {
      // initialize quadrature rule
      this->qrule->init(side->type(), elem->p_level());

      {
        // Set the element type
        this->elem_type = elem->type();

        // Initialize the face shape functions
        this->_fe_map->template init_face_shape_functions<Dim>(this->qrule->get_points(),  side.get());
      }
      // We can't get away without recomputing shape functions next
      // time
      this->shapes_on_quadrature = false;
      // Compute data on the face for integration
      this->compute_face_values (elem, side.get(), this->qrule->get_weights());
    }
}



template <unsigned int Dim>
void FEXYZ<Dim>::compute_face_values(const Elem * elem,
                                     const Elem * side,
                                     const std::vector<Real> & qw)
{
  libmesh_assert(elem);
  libmesh_assert(side);

  LOG_SCOPE("compute_face_values()", "FEXYZ");

  // The number of quadrature points.
  const std::size_t n_qp = qw.size();

  // Number of shape functions in the finite element approximation
  // space.
  const unsigned int n_approx_shape_functions =
    this->n_shape_functions(this->get_type(),
                            this->get_order());

  // Resize the shape functions and their gradients
  this->phi.resize    (n_approx_shape_functions);
  this->dphi.resize   (n_approx_shape_functions);
  this->dphidx.resize (n_approx_shape_functions);
  this->dphidy.resize (n_approx_shape_functions);
  this->dphidz.resize (n_approx_shape_functions);

  for (unsigned int i=0; i<n_approx_shape_functions; i++)
    {
      this->phi[i].resize    (n_qp);
      this->dphi[i].resize   (n_qp);
      this->dphidx[i].resize (n_qp);
      this->dphidy[i].resize (n_qp);
      this->dphidz[i].resize (n_qp);
    }

  this->_fe_map->compute_face_map(this->dim, qw, side);

  const std::vector<libMesh::Point> & xyz = this->_fe_map->get_xyz();

  switch (this->dim)
    {
      // A 2D finite element living in either 2D or 3D space.
      // This means the boundary is a 1D finite element, i.e.
      // and EDGE2 or EDGE3.
    case 2:
      {
        // compute the shape function values & gradients
        for (unsigned int i=0; i<n_approx_shape_functions; i++)
          for (std::size_t p=0; p<n_qp; p++)
            {
              this->phi[i][p] = FE<Dim,XYZ>::shape (elem, this->fe_type.order, i, xyz[p]);

              this->dphi[i][p](0) =
                this->dphidx[i][p] = FE<Dim,XYZ>::shape_deriv (elem, this->fe_type.order, i, 0, xyz[p]);

              this->dphi[i][p](1) =
                this->dphidy[i][p] = FE<Dim,XYZ>::shape_deriv (elem, this->fe_type.order, i, 1, xyz[p]);

#if LIBMESH_DIM == 3
              this->dphi[i][p](2) = // can only assign to the Z component if LIBMESH_DIM==3
#endif
                this->dphidz[i][p] = 0.;
            }

        // done computing face values
        break;
      }

      // A 3D finite element living in 3D space.
    case 3:
      {
        // compute the shape function values & gradients
        for (unsigned int i=0; i<n_approx_shape_functions; i++)
          for (std::size_t p=0; p<n_qp; p++)
            {
              this->phi[i][p] = FE<Dim,XYZ>::shape (elem, this->fe_type.order, i, xyz[p]);

              this->dphi[i][p](0) =
                this->dphidx[i][p] = FE<Dim,XYZ>::shape_deriv (elem, this->fe_type.order, i, 0, xyz[p]);

              this->dphi[i][p](1) =
                this->dphidy[i][p] = FE<Dim,XYZ>::shape_deriv (elem, this->fe_type.order, i, 1, xyz[p]);

              this->dphi[i][p](2) =
                this->dphidz[i][p] = FE<Dim,XYZ>::shape_deriv (elem, this->fe_type.order, i, 2, xyz[p]);
            }

        // done computing face values
        break;
      }

    default:
      libmesh_error_msg("Invalid dim " << this->dim);
    }
}




//--------------------------------------------------------------
// Explicit instantiations (doesn't make sense in 1D!) using fe_macro.h's macro

template class FEXYZ<2>;
template class FEXYZ<3>;


} // namespace libMesh
