// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include <cmath> // for std::sqrt


// Local includes
#include "libmesh_common.h"
#include "fe.h"
#include "quadrature.h"
#include "elem.h"
#include "libmesh_logging.h"

namespace libMesh
{




//-------------------------------------------------------
// Full specialization for 0D when this is a useless method
template <>
void FEXYZ<0>::reinit(const Elem*,
		      const unsigned int,
		      const Real,
                      const std::vector<Point>* const,
                      const std::vector<Real>* const)
{
  libMesh::err << "ERROR: This method only makes sense for 2D/3D elements!"
	        << std::endl;
  libmesh_error();
}


//-------------------------------------------------------
// Full specialization for 1D when this is a useless method
template <>
void FEXYZ<1>::reinit(const Elem*,
		      const unsigned int,
		      const Real,
                      const std::vector<Point>* const,
                      const std::vector<Real>* const)
{
  libMesh::err << "ERROR: This method only makes sense for 2D/3D elements!"
	        << std::endl;
  libmesh_error();
}


//-------------------------------------------------------
// Method for 2D, 3D
template <unsigned int Dim>
void FEXYZ<Dim>::reinit(const Elem* elem,
			const unsigned int s,
			const Real,
                        const std::vector<Point>* const pts,
                        const std::vector<Real>* const weights)
{
  libmesh_assert (elem  != NULL);
  libmesh_assert (this->qrule != NULL || pts != NULL);
  // We don't do this for 1D elements!
  libmesh_assert (Dim != 1);

  // Build the side of interest
  const AutoPtr<Elem> side(elem->build_side(s));

  // Initialize the shape functions at the user-specified
  // points
  if (pts != NULL)
    {
      // We can't get away without recomputing shape functions next
      // time
      this->shapes_on_quadrature = false;

      // Initialize the face shape functions
      this->init_face_shape_functions (*pts,  side.get());

    }
  else
    {
      // initialize quadrature rule
      this->qrule->init(side->type(), elem->p_level());

        {
          // Set the element type
          this->elem_type = elem->type();

          // Initialize the face shape functions
          this->init_face_shape_functions (this->qrule->get_points(),  side.get());
        }
    }

  // Compute data on the face for integration
  this->compute_face_values (elem, side.get(),
    weights ? *weights : this->qrule->get_weights());
}



template <unsigned int Dim>
void FEXYZ<Dim>::compute_face_values(const Elem* elem,
				     const Elem* side,
                                     const std::vector<Real>& qw)
{
  libmesh_assert (elem != NULL);
  libmesh_assert (side != NULL);

  START_LOG("compute_face_values()", "FEXYZ");

  // The number of quadrature points.
  const unsigned int n_qp = qw.size();

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

  switch (this->dim)
    {

    case 2:
      {
	// A 2D finite element living in either 2D or 3D space.
	// This means the boundary is a 1D finite element, i.e.
	// and EDGE2 or EDGE3.
	// Resize the vectors to hold data at the quadrature points
	{
	  this->xyz.resize(n_qp);
	  this->dxyzdxi_map.resize(n_qp);
	  this->d2xyzdxi2_map.resize(n_qp);
	  this->tangents.resize(n_qp);
	  this->normals.resize(n_qp);
	  this->curvatures.resize(n_qp);

	  this->JxW.resize(n_qp);
	}

	// Clear the entities that will be summed
	// Compute the tangent & normal at the quadrature point
	for (unsigned int p=0; p<n_qp; p++)
	  {
	    this->tangents[p].resize(LIBMESH_DIM-1); // 1 Tangent in 2D, 2 in 3D
	    this->xyz[p].zero();
	    this->dxyzdxi_map[p].zero();
	    this->d2xyzdxi2_map[p].zero();
	  }

	// compute x, dxdxi at the quadrature points
	for (unsigned int i=0; i<this->psi_map.size(); i++) // sum over the nodes
	  {
	    const Point& side_point = side->point(i);

	    for (unsigned int p=0; p<n_qp; p++) // for each quadrature point...
	      {
		this->xyz[p].add_scaled          (side_point, this->psi_map[i][p]);
		this->dxyzdxi_map[p].add_scaled  (side_point, this->dpsidxi_map[i][p]);
		this->d2xyzdxi2_map[p].add_scaled(side_point, this->d2psidxi2_map[i][p]);
	      }
	  }

	// Compute the tangent & normal at the quadrature point
	for (unsigned int p=0; p<n_qp; p++)
	  {
	    const Point n(this->dxyzdxi_map[p](1), -this->dxyzdxi_map[p](0), 0.);

	    this->normals[p]     = n.unit();
	    this->tangents[p][0] = this->dxyzdxi_map[p].unit();
#if LIBMESH_DIM == 3  // Only good in 3D space
	    this->tangents[p][1] = this->dxyzdxi_map[p].cross(n).unit();
#endif
	    // The curvature is computed via the familiar Frenet formula:
	    // curvature = [d^2(x) / d (xi)^2] dot [normal]
	    // For a reference, see:
	    // F.S. Merritt, Mathematics Manual, 1962, McGraw-Hill, p. 310
	    //
	    // Note: The sign convention here is different from the
	    // 3D case.  Concave-upward curves (smiles) have a positive
	    // curvature.  Concave-downward curves (frowns) have a
	    // negative curvature.  Be sure to take that into account!
	    const Real numerator   = this->d2xyzdxi2_map[p] * this->normals[p];
	    const Real denominator = this->dxyzdxi_map[p].size_sq();
	    libmesh_assert (denominator != 0);
	    this->curvatures[p] = numerator / denominator;
	  }

	// compute the jacobian at the quadrature points
	for (unsigned int p=0; p<n_qp; p++)
	  {
	    const Real jac = std::sqrt(this->dxdxi_map(p)*this->dxdxi_map(p) +
				       this->dydxi_map(p)*this->dydxi_map(p));

	    libmesh_assert (jac > 0.);

	    this->JxW[p] = jac*qw[p];
	  }

	// compute the shape function values & gradients
	for (unsigned int i=0; i<n_approx_shape_functions; i++)
	  for (unsigned int p=0; p<n_qp; p++)
	    {
	      this->phi[i][p] = FE<Dim,XYZ>::shape (elem, this->fe_type.order, i, this->xyz[p]);

	      this->dphi[i][p](0) =
		this->dphidx[i][p] = FE<Dim,XYZ>::shape_deriv (elem, this->fe_type.order, i, 0, this->xyz[p]);

	      this->dphi[i][p](1) =
		this->dphidy[i][p] = FE<Dim,XYZ>::shape_deriv (elem, this->fe_type.order, i, 1, this->xyz[p]);

#if LIBMESH_DIM == 3
	      this->dphi[i][p](2) = // can only assign to the Z component if LIBMESH_DIM==3
#endif
		this->dphidz[i][p] = 0.;
	    }

	// done computing face values
	break;
      }



    case 3:
      {
	// A 3D finite element living in 3D space.
	// Resize the vectors to hold data at the quadrature points
	{
	  this->xyz.resize(n_qp);
	  this->dxyzdxi_map.resize(n_qp);
	  this->dxyzdeta_map.resize(n_qp);
	  this->d2xyzdxi2_map.resize(n_qp);
	  this->d2xyzdxideta_map.resize(n_qp);
	  this->d2xyzdeta2_map.resize(n_qp);
	  this->tangents.resize(n_qp);
	  this->normals.resize(n_qp);
	  this->curvatures.resize(n_qp);

	  this->JxW.resize(n_qp);
	}

	// Clear the entities that will be summed
	for (unsigned int p=0; p<n_qp; p++)
	  {
	    this->tangents[p].resize(LIBMESH_DIM-1); // 1 Tangent in 2D, 2 in 3D
	    this->xyz[p].zero();
	    this->dxyzdxi_map[p].zero();
	    this->dxyzdeta_map[p].zero();
	    this->d2xyzdxi2_map[p].zero();
	    this->d2xyzdxideta_map[p].zero();
	    this->d2xyzdeta2_map[p].zero();
	  }

	// compute x, dxdxi at the quadrature points
	for (unsigned int i=0; i<this->psi_map.size(); i++) // sum over the nodes
	  {
	    const Point& side_point = side->point(i);

	    for (unsigned int p=0; p<n_qp; p++) // for each quadrature point...
	      {
		this->xyz[p].add_scaled             (side_point, this->psi_map[i][p]);
		this->dxyzdxi_map[p].add_scaled     (side_point, this->dpsidxi_map[i][p]);
		this->dxyzdeta_map[p].add_scaled    (side_point, this->dpsideta_map[i][p]);
		this->d2xyzdxi2_map[p].add_scaled   (side_point, this->d2psidxi2_map[i][p]);
		this->d2xyzdxideta_map[p].add_scaled(side_point, this->d2psidxideta_map[i][p]);
		this->d2xyzdeta2_map[p].add_scaled  (side_point, this->d2psideta2_map[i][p]);
	      }
	  }

	// Compute the tangents, normal, and curvature at the quadrature point
	for (unsigned int p=0; p<n_qp; p++)
	  {
	    const Point n  = this->dxyzdxi_map[p].cross(this->dxyzdeta_map[p]);
	    this->normals[p]     = n.unit();
	    this->tangents[p][0] = this->dxyzdxi_map[p].unit();
	    this->tangents[p][1] = n.cross(this->dxyzdxi_map[p]).unit();

	    // Compute curvature using the typical nomenclature
	    // of the first and second fundamental forms.
	    // For reference, see:
	    // 1) http://mathworld.wolfram.com/MeanCurvature.html
	    //    (note -- they are using inward normal)
	    // 2) F.S. Merritt, Mathematics Manual, 1962, McGraw-Hill
	    const Real L  = -this->d2xyzdxi2_map[p]    * this->normals[p];
	    const Real M  = -this->d2xyzdxideta_map[p] * this->normals[p];
	    const Real N  = -this->d2xyzdeta2_map[p]   * this->normals[p];
	    const Real E  =  this->dxyzdxi_map[p].size_sq();
	    const Real F  =  this->dxyzdxi_map[p]      * this->dxyzdeta_map[p];
	    const Real G  =  this->dxyzdeta_map[p].size_sq();

	    const Real numerator   = E*N -2.*F*M + G*L;
	    const Real denominator = E*G - F*F;
	    libmesh_assert (denominator != 0.);
	    this->curvatures[p] = 0.5*numerator/denominator;
	  }

	// compute the jacobian at the quadrature points, see
	// http://sp81.msi.umn.edu:999/fluent/fidap/help/theory/thtoc.htm
	for (unsigned int p=0; p<n_qp; p++)
	  {
	    const Real g11 = (this->dxdxi_map(p)*this->dxdxi_map(p) +
			      this->dydxi_map(p)*this->dydxi_map(p) +
			      this->dzdxi_map(p)*this->dzdxi_map(p));

	    const Real g12 = (this->dxdxi_map(p)*this->dxdeta_map(p) +
			      this->dydxi_map(p)*this->dydeta_map(p) +
			      this->dzdxi_map(p)*this->dzdeta_map(p));

	    const Real g21 = g12;

	    const Real g22 = (this->dxdeta_map(p)*this->dxdeta_map(p) +
			      this->dydeta_map(p)*this->dydeta_map(p) +
			      this->dzdeta_map(p)*this->dzdeta_map(p));


	    const Real jac = std::sqrt(g11*g22 - g12*g21);

	    libmesh_assert (jac > 0.);

	    this->JxW[p] = jac*qw[p];
	  }

	// compute the shape function values & gradients
	for (unsigned int i=0; i<n_approx_shape_functions; i++)
	  for (unsigned int p=0; p<n_qp; p++)
	    {
	      this->phi[i][p] = FE<Dim,XYZ>::shape (elem, this->fe_type.order, i, this->xyz[p]);

	      this->dphi[i][p](0) =
		this->dphidx[i][p] = FE<Dim,XYZ>::shape_deriv (elem, this->fe_type.order, i, 0, this->xyz[p]);

	      this->dphi[i][p](1) =
		this->dphidy[i][p] = FE<Dim,XYZ>::shape_deriv (elem, this->fe_type.order, i, 1, this->xyz[p]);

	      this->dphi[i][p](2) =
		this->dphidz[i][p] = FE<Dim,XYZ>::shape_deriv (elem, this->fe_type.order, i, 2, this->xyz[p]);
	    }

	// done computing face values
	break;
      }


    default:
      libmesh_error();

    }
  STOP_LOG("compute_face_values()", "FEXYZ");
}




//--------------------------------------------------------------
// Explicit instantiations (doesn't make sense in 1D!) using fe_macro.h's macro

template class FEXYZ<2>;
template class FEXYZ<3>;


} // namespace libMesh
