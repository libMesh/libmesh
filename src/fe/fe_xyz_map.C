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

#include "libmesh/fe_xyz_map.h"
#include "libmesh/elem.h"

namespace libMesh
{

void FEXYZMap::compute_face_map(int dim, const std::vector<Real> & qw, const Elem * side)
{
  libmesh_assert(side);

  LOG_SCOPE("compute_face_map()", "FEXYZMap");

  // The number of quadrature points.
  const unsigned int n_qp = cast_int<unsigned int>(qw.size());

  switch(dim)
    {
    case 2:
      {

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
        for (std::size_t i=0; i<this->psi_map.size(); i++) // sum over the nodes
          {
            const Point & side_point = side->point(i);

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
            const Real denominator = this->dxyzdxi_map[p].norm_sq();
            libmesh_assert_not_equal_to (denominator, 0);
            this->curvatures[p] = numerator / denominator;
          }

        // compute the jacobian at the quadrature points
        for (unsigned int p=0; p<n_qp; p++)
          {
            const Real the_jac = std::sqrt(this->dxdxi_map(p)*this->dxdxi_map(p) +
                                           this->dydxi_map(p)*this->dydxi_map(p));

            libmesh_assert_greater (the_jac, 0.);

            this->JxW[p] = the_jac*qw[p];
          }

        break;
      }

    case 3:
      {
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
        for (std::size_t i=0; i<this->psi_map.size(); i++) // sum over the nodes
          {
            const Point & side_point = side->point(i);

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
            const Real E  =  this->dxyzdxi_map[p].norm_sq();
            const Real F  =  this->dxyzdxi_map[p]      * this->dxyzdeta_map[p];
            const Real G  =  this->dxyzdeta_map[p].norm_sq();

            const Real numerator   = E*N -2.*F*M + G*L;
            const Real denominator = E*G - F*F;
            libmesh_assert_not_equal_to (denominator, 0.);
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


            const Real the_jac = std::sqrt(g11*g22 - g12*g21);

            libmesh_assert_greater (the_jac, 0.);

            this->JxW[p] = the_jac*qw[p];
          }

        break;
      }
    default:
      libmesh_error_msg("Invalid dim = " << dim);

    } // switch(dim)
}

} // namespace libMesh
