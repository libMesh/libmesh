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

#ifndef LIBMESH_FE_MAP_IMPL_H
#define LIBMESH_FE_MAP_IMPL_H

// C++ includes
#include <cstdlib> // *must* precede <cmath> for proper std:abs() on PGI, Sun Studio CC
#include <cmath> // for std::sqrt, std::abs


// Local includes
#include "libmesh/fe.h"
#include "libmesh/elem.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fe_macro.h"
#include "libmesh/fe_map.h"
#include "libmesh/fe_xyz_map.h"
#include "libmesh/inf_fe_map.h"
#include "libmesh/mesh_subdivision_support.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/tensor_value.h"
#include "libmesh/auto_ptr.h" // libmesh_make_unique
#include "libmesh/enum_elem_type.h"
#include "libmesh/int_range.h"

namespace libMesh
{

template <typename RealType>
inline
FEFamily
FEMapTempl<RealType>::map_fe_type(const Elem & elem)
{
  switch (elem.mapping_type())
  {
  case RATIONAL_BERNSTEIN_MAP:
    return RATIONAL_BERNSTEIN;
  case LAGRANGE_MAP:
    return LAGRANGE;
  default:
    libmesh_error_msg("Unknown mapping type " << elem.mapping_type());
  }
  return LAGRANGE;
}



// Constructor
template <typename RealType>
FEMapTempl<RealType>::FEMapTempl(Real jtol) :
  calculations_started(false),
  calculate_xyz(false),
  calculate_dxyz(false),
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  calculate_d2xyz(false),
#endif
  jacobian_tolerance(jtol)
{}



template <typename RealType>
std::unique_ptr<FEMapTempl<RealType>> FEMapTempl<RealType>::build( FEType fe_type )
{
  switch( fe_type.family )
    {
    case XYZ:
      return libmesh_make_unique<FEXYZMap>();

    default:
      return libmesh_make_unique<FEMap>();
    }
}



template <typename RealType>
template<unsigned int Dim>
void FEMapTempl<RealType>::init_reference_to_physical_map(const std::vector<Point> & qp,
                                           const Elem * elem)
{
  // Start logging the reference->physical map initialization
  LOG_SCOPE("init_reference_to_physical_map()", "FEMap");

  // We're calculating now!
  this->determine_calculations();

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
  if (elem->infinite())
    {
      //This mainly requires to change the FE<>-calls
      // to FEInterface in this function.
      libmesh_not_implemented();
    }
#endif

  // The number of quadrature points.
  const std::size_t n_qp = qp.size();

  // The element type and order to use in
  // the map
  const FEFamily mapping_family = FEMap::map_fe_type(*elem);
  const Order    mapping_order     (elem->default_order());
  const ElemType mapping_elem_type (elem->type());

  const FEType map_fe_type(mapping_order, mapping_family);

  // Number of shape functions used to construct the map
  // (Lagrange shape functions are used for mapping)
  const unsigned int n_mapping_shape_functions =
    FEInterface::n_shape_functions (Dim, map_fe_type,
                                    mapping_elem_type);

  if (calculate_xyz)
    this->phi_map.resize         (n_mapping_shape_functions);
  if (Dim > 0)
    {
      if (calculate_dxyz)
        this->dphidxi_map.resize     (n_mapping_shape_functions);
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
      if (calculate_d2xyz)
        this->d2phidxi2_map.resize   (n_mapping_shape_functions);
#endif // ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
    }

  if (Dim > 1)
    {
      if (calculate_dxyz)
        this->dphideta_map.resize  (n_mapping_shape_functions);
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
      if (calculate_d2xyz)
        {
          this->d2phidxideta_map.resize   (n_mapping_shape_functions);
          this->d2phideta2_map.resize     (n_mapping_shape_functions);
        }
#endif // ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
    }

  if (Dim > 2)
    {
      if (calculate_dxyz)
        this->dphidzeta_map.resize (n_mapping_shape_functions);
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
      if (calculate_d2xyz)
        {
          this->d2phidxidzeta_map.resize  (n_mapping_shape_functions);
          this->d2phidetadzeta_map.resize (n_mapping_shape_functions);
          this->d2phidzeta2_map.resize    (n_mapping_shape_functions);
        }
#endif // ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
    }


  for (unsigned int i=0; i<n_mapping_shape_functions; i++)
    {
      if (calculate_xyz)
        this->phi_map[i].resize         (n_qp);
      if (Dim > 0)
        {
          if (calculate_dxyz)
            this->dphidxi_map[i].resize     (n_qp);
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
          if (calculate_d2xyz)
            {
              this->d2phidxi2_map[i].resize     (n_qp);
              if (Dim > 1)
                {
                  this->d2phidxideta_map[i].resize (n_qp);
                  this->d2phideta2_map[i].resize (n_qp);
                }
              if (Dim > 2)
                {
                  this->d2phidxidzeta_map[i].resize  (n_qp);
                  this->d2phidetadzeta_map[i].resize (n_qp);
                  this->d2phidzeta2_map[i].resize    (n_qp);
                }
            }
#endif // ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

          if (Dim > 1 && calculate_dxyz)
            this->dphideta_map[i].resize  (n_qp);

          if (Dim > 2 && calculate_dxyz)
            this->dphidzeta_map[i].resize (n_qp);
        }
    }

  // Optimize for the *linear* geometric elements case:
  bool is_linear = elem->is_linear();

  FEInterface::shape_ptr<RealType> shape_ptr =
    FEInterface::shape_function<RealType>(Dim, map_fe_type);

  FEInterface::shape_deriv_ptr<RealType> shape_deriv_ptr =
    FEInterface::shape_deriv_function<RealType>(Dim, map_fe_type);

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  FEInterface::shape_second_deriv_ptr<RealType> shape_second_deriv_ptr =
    FEInterface::shape_second_deriv_function<RealType>(Dim, map_fe_type);
#endif // ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

  switch (Dim)
    {
      //------------------------------------------------------------
      // 0D
    case 0:
      {
        if (calculate_xyz)
          for (unsigned int i=0; i<n_mapping_shape_functions; i++)
            for (std::size_t p=0; p<n_qp; p++)
              this->phi_map[i][p] =
                shape_ptr(elem, mapping_order, i, qp[p], false);

        break;
      }

      //------------------------------------------------------------
      // 1D
    case 1:
      {
        // Compute the value of the mapping shape function i at quadrature point p
        // (Lagrange shape functions are used for mapping)
        if (is_linear)
          {
            for (unsigned int i=0; i<n_mapping_shape_functions; i++)
              {
                if (calculate_xyz)
                  this->phi_map[i][0] =
                    shape_ptr(elem, mapping_order, i, qp[0], false);

                if (calculate_dxyz)
                  this->dphidxi_map[i][0] =
                    shape_deriv_ptr(elem, mapping_order, i, 0, qp[0], false);

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
                if (calculate_d2xyz)
                  this->d2phidxi2_map[i][0] =
                    shape_second_deriv_ptr(elem, mapping_order, i, 0, qp[0], false);
#endif // ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
                for (std::size_t p=1; p<n_qp; p++)
                  {
                    if (calculate_xyz)
                      this->phi_map[i][p] =
                        shape_ptr(elem, mapping_order, i, qp[p], false);
                    if (calculate_dxyz)
                      this->dphidxi_map[i][p]  = this->dphidxi_map[i][0];
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
                    if (calculate_d2xyz)
                      this->d2phidxi2_map[i][p] = this->d2phidxi2_map[i][0];
#endif // ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
                  }
              }
          }
        else
          for (unsigned int i=0; i<n_mapping_shape_functions; i++)
            for (std::size_t p=0; p<n_qp; p++)
              {
                if (calculate_xyz)
                  this->phi_map[i][p] =
                    shape_ptr (elem, mapping_order, i, qp[p], false);
                if (calculate_dxyz)
                  this->dphidxi_map[i][p]  = shape_deriv_ptr (elem, mapping_order, i, 0, qp[p], false);
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
                if (calculate_d2xyz)
                  this->d2phidxi2_map[i][p] = shape_second_deriv_ptr (elem, mapping_order, i, 0, qp[p], false);
#endif // ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
              }

        break;
      }
      //------------------------------------------------------------
      // 2D
    case 2:
      {
        // Compute the value of the mapping shape function i at quadrature point p
        // (Lagrange shape functions are used for mapping)
        if (is_linear)
          {
            for (unsigned int i=0; i<n_mapping_shape_functions; i++)
              {
                if (calculate_xyz)
                  this->phi_map[i][0] =
                    shape_ptr (elem, mapping_order, i, qp[0], false);
                if (calculate_dxyz)
                  {
                    this->dphidxi_map[i][0]  = shape_deriv_ptr (elem, mapping_order, i, 0, qp[0], false);
                    this->dphideta_map[i][0] = shape_deriv_ptr (elem, mapping_order, i, 1, qp[0], false);
                  }
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
                if (calculate_d2xyz)
                  {
                    this->d2phidxi2_map[i][0]    = shape_second_deriv_ptr (elem, mapping_order, i, 0, qp[0], false);
                    this->d2phidxideta_map[i][0] = shape_second_deriv_ptr (elem, mapping_order, i, 1, qp[0], false);
                    this->d2phideta2_map[i][0]   = shape_second_deriv_ptr (elem, mapping_order, i, 2, qp[0], false);
                  }
#endif // ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
                for (std::size_t p=1; p<n_qp; p++)
                  {
                    if (calculate_xyz)
                      this->phi_map[i][p] =
                        shape_ptr (elem, mapping_order, i, qp[p], false);
                    if (calculate_dxyz)
                      {
                        this->dphidxi_map[i][p]  = this->dphidxi_map[i][0];
                        this->dphideta_map[i][p] = this->dphideta_map[i][0];
                      }
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
                    if (calculate_d2xyz)
                      {
                        this->d2phidxi2_map[i][p] = this->d2phidxi2_map[i][0];
                        this->d2phidxideta_map[i][p] = this->d2phidxideta_map[i][0];
                        this->d2phideta2_map[i][p] = this->d2phideta2_map[i][0];
                      }
#endif // ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
                  }
              }
          }
        else
          for (unsigned int i=0; i<n_mapping_shape_functions; i++)
            for (std::size_t p=0; p<n_qp; p++)
              {
                if (calculate_xyz)
                  this->phi_map[i][p] =
                    shape_ptr(elem, mapping_order, i, qp[p], false);
                if (calculate_dxyz)
                  {
                    this->dphidxi_map[i][p]  = shape_deriv_ptr (elem, mapping_order, i, 0, qp[p], false);
                    this->dphideta_map[i][p] = shape_deriv_ptr (elem, mapping_order, i, 1, qp[p], false);
                  }
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
                if (calculate_d2xyz)
                  {
                    this->d2phidxi2_map[i][p] = shape_second_deriv_ptr (elem, mapping_order, i, 0, qp[p], false);
                    this->d2phidxideta_map[i][p] = shape_second_deriv_ptr (elem, mapping_order, i, 1, qp[p], false);
                    this->d2phideta2_map[i][p] = shape_second_deriv_ptr (elem, mapping_order, i, 2, qp[p], false);
                  }
#endif // ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
              }

        break;
      }

      //------------------------------------------------------------
      // 3D
    case 3:
      {
        // Compute the value of the mapping shape function i at quadrature point p
        // (Lagrange shape functions are used for mapping)
        if (is_linear)
          {
            for (unsigned int i=0; i<n_mapping_shape_functions; i++)
              {
                if (calculate_xyz)
                  this->phi_map[i][0] =
                    shape_ptr (elem, mapping_order, i, qp[0], false);
                if (calculate_dxyz)
                  {
                    this->dphidxi_map[i][0]  = shape_deriv_ptr (elem, mapping_order, i, 0, qp[0], false);
                    this->dphideta_map[i][0] = shape_deriv_ptr (elem, mapping_order, i, 1, qp[0], false);
                    this->dphidzeta_map[i][0] = shape_deriv_ptr (elem, mapping_order, i, 2, qp[0], false);
                  }
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
                if (calculate_d2xyz)
                  {
                    this->d2phidxi2_map[i][0]      = shape_second_deriv_ptr (elem, mapping_order, i, 0, qp[0], false);
                    this->d2phidxideta_map[i][0]   = shape_second_deriv_ptr (elem, mapping_order, i, 1, qp[0], false);
                    this->d2phideta2_map[i][0]     = shape_second_deriv_ptr (elem, mapping_order, i, 2, qp[0], false);
                    this->d2phidxidzeta_map[i][0]  = shape_second_deriv_ptr (elem, mapping_order, i, 3, qp[0], false);
                    this->d2phidetadzeta_map[i][0] = shape_second_deriv_ptr (elem, mapping_order, i, 4, qp[0], false);
                    this->d2phidzeta2_map[i][0]    = shape_second_deriv_ptr (elem, mapping_order, i, 5, qp[0], false);
                  }
#endif // ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
                for (std::size_t p=1; p<n_qp; p++)
                  {
                    if (calculate_xyz)
                      this->phi_map[i][p] =
                        shape_ptr (elem, mapping_order, i, qp[p], false);
                    if (calculate_dxyz)
                      {
                        this->dphidxi_map[i][p]  = this->dphidxi_map[i][0];
                        this->dphideta_map[i][p] = this->dphideta_map[i][0];
                        this->dphidzeta_map[i][p] = this->dphidzeta_map[i][0];
                      }
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
                    if (calculate_d2xyz)
                      {
                        this->d2phidxi2_map[i][p] = this->d2phidxi2_map[i][0];
                        this->d2phidxideta_map[i][p] = this->d2phidxideta_map[i][0];
                        this->d2phideta2_map[i][p] = this->d2phideta2_map[i][0];
                        this->d2phidxidzeta_map[i][p] = this->d2phidxidzeta_map[i][0];
                        this->d2phidetadzeta_map[i][p] = this->d2phidetadzeta_map[i][0];
                        this->d2phidzeta2_map[i][p] = this->d2phidzeta2_map[i][0];
                      }
#endif // ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
                  }
              }
          }
        else
          for (unsigned int i=0; i<n_mapping_shape_functions; i++)
            for (std::size_t p=0; p<n_qp; p++)
              {
                if (calculate_xyz)
                  this->phi_map[i][p] =
                    shape_ptr(elem, mapping_order, i, qp[p], false);
                if (calculate_dxyz)
                  {
                    this->dphidxi_map[i][p]   = shape_deriv_ptr (elem, mapping_order, i, 0, qp[p], false);
                    this->dphideta_map[i][p]  = shape_deriv_ptr (elem, mapping_order, i, 1, qp[p], false);
                    this->dphidzeta_map[i][p] = shape_deriv_ptr (elem, mapping_order, i, 2, qp[p], false);
                  }
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
                if (calculate_d2xyz)
                  {
                    this->d2phidxi2_map[i][p]      = shape_second_deriv_ptr (elem, mapping_order, i, 0, qp[p], false);
                    this->d2phidxideta_map[i][p]   = shape_second_deriv_ptr (elem, mapping_order, i, 1, qp[p], false);
                    this->d2phideta2_map[i][p]     = shape_second_deriv_ptr (elem, mapping_order, i, 2, qp[p], false);
                    this->d2phidxidzeta_map[i][p]  = shape_second_deriv_ptr (elem, mapping_order, i, 3, qp[p], false);
                    this->d2phidetadzeta_map[i][p] = shape_second_deriv_ptr (elem, mapping_order, i, 4, qp[p], false);
                    this->d2phidzeta2_map[i][p]    = shape_second_deriv_ptr (elem, mapping_order, i, 5, qp[p], false);
                  }
#endif // ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
              }

        break;
      }

    default:
      libmesh_error_msg("Invalid Dim = " << Dim);
    }
}



template <typename RealType>
void FEMapTempl<RealType>::compute_single_point_map(const unsigned int dim,
                                     const std::vector<Real> & qw,
                                     const Elem * elem,
                                     unsigned int p,
                                     const std::vector<const Node *> & elem_nodes,
                                     bool compute_second_derivatives)
{
  libmesh_assert(elem);
  libmesh_assert(calculations_started);
#ifndef LIBMESH_ENABLE_SECOND_DERIVATIVES
  libmesh_assert(!compute_second_derivatives);
#endif

  if (calculate_xyz)
    libmesh_assert_equal_to(phi_map.size(), elem_nodes.size());

  switch (dim)
    {
      //--------------------------------------------------------------------
      // 0D
    case 0:
      {
        libmesh_assert(elem_nodes[0]);
        if (calculate_xyz)
          xyz[p] = *elem_nodes[0];
        if (calculate_dxyz)
          {
            jac[p] = 1.0;
            JxW[p] = qw[p];
          }
        break;
      }

      //--------------------------------------------------------------------
      // 1D
    case 1:
      {
        // Clear the entities that will be summed
        if (calculate_xyz)
          xyz[p].zero();
        if (calculate_dxyz)
          dxyzdxi_map[p].zero();
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
        if (calculate_d2xyz)
          {
            d2xyzdxi2_map[p].zero();
            // Inverse map second derivatives
            d2xidxyz2_map[p].assign(6, 0.);
          }
#endif

        // compute x, dx, d2x at the quadrature point
        for (auto i : index_range(elem_nodes)) // sum over the nodes
          {
            // Reference to the point, helps eliminate
            // excessive temporaries in the inner loop
            libmesh_assert(elem_nodes[i]);
            const Point & elem_point = *elem_nodes[i];

            if (calculate_xyz)
              xyz[p].add_scaled          (elem_point, phi_map[i][p]    );
            if (calculate_dxyz)
              dxyzdxi_map[p].add_scaled  (elem_point, dphidxi_map[i][p]);
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
            if (calculate_d2xyz)
              d2xyzdxi2_map[p].add_scaled(elem_point, d2phidxi2_map[i][p]);
#endif
          }

        // Compute the jacobian
        //
        // 1D elements can live in 2D or 3D space.
        // The transformation matrix from local->global
        // coordinates is
        //
        // T = | dx/dxi |
        //     | dy/dxi |
        //     | dz/dxi |
        //
        // The generalized determinant of T (from the
        // so-called "normal" eqns.) is
        // jac = "det(T)" = sqrt(det(T'T))
        //
        // where T'= transpose of T, so
        //
        // jac = sqrt( (dx/dxi)^2 + (dy/dxi)^2 + (dz/dxi)^2 )

        if (calculate_dxyz)
          {
            jac[p] = dxyzdxi_map[p].norm();

            if (jac[p] <= jacobian_tolerance)
              {
                // Don't call print_info() recursively if we're already
                // failing.  print_info() calls Elem::volume() which may
                // call FE::reinit() and trigger the same failure again.
                static bool failing = false;
                if (!failing)
                  {
                    failing = true;
                    elem->print_info(libMesh::err);
                    failing = false;
                    if (calculate_xyz)
                      {
                        libmesh_error_msg("ERROR: negative Jacobian " \
                                          << jac[p] \
                                          << " at point " \
                                          << xyz[p] \
                                          << " in element " \
                                          << elem->id());
                      }
                    else
                      {
                        // In this case xyz[p] is not defined, so don't
                        // try to print it out.
                        libmesh_error_msg("ERROR: negative Jacobian " \
                                          << jac[p] \
                                          << " at point index " \
                                          << p \
                                          << " in element " \
                                          << elem->id());
                      }
                  }
                else
                  {
                    // We were already failing when we called this, so just
                    // stop the current computation and return with
                    // incomplete results.
                    return;
                  }
              }

            // The inverse Jacobian entries also come from the
            // generalized inverse of T (see also the 2D element
            // living in 3D code).
            const auto jacm2 = 1./jac[p]/jac[p];
            dxidx_map[p] = jacm2*dxdxi_map(p);
#if LIBMESH_DIM > 1
            dxidy_map[p] = jacm2*dydxi_map(p);
#endif
#if LIBMESH_DIM > 2
            dxidz_map[p] = jacm2*dzdxi_map(p);
#endif

            JxW[p] = jac[p]*qw[p];
          }

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

        if (calculate_d2xyz)
          {
#if LIBMESH_DIM == 1
            // Compute inverse map second derivatives for 1D-element-living-in-1D case
            this->compute_inverse_map_second_derivs(p);
#elif LIBMESH_DIM == 2
            // Compute inverse map second derivatives for 1D-element-living-in-2D case
            // See JWP notes for details

            // numer = x_xi*x_{xi xi} + y_xi*y_{xi xi}
            auto numer =
              dxyzdxi_map[p](0)*d2xyzdxi2_map[p](0) +
              dxyzdxi_map[p](1)*d2xyzdxi2_map[p](1);

            // denom = (x_xi)^2 + (y_xi)^2 must be >= 0.0
            auto denom =
              dxyzdxi_map[p](0)*dxyzdxi_map[p](0) +
              dxyzdxi_map[p](1)*dxyzdxi_map[p](1);

            if (denom <= 0.0)
              {
                // Don't call print_info() recursively if we're already
                // failing.  print_info() calls Elem::volume() which may
                // call FE::reinit() and trigger the same failure again.
                static bool failing = false;
                if (!failing)
                  {
                    failing = true;
                    elem->print_info(libMesh::err);
                    failing = false;
                    libmesh_error_msg("Encountered invalid 1D element!");
                  }
                else
                  {
                    // We were already failing when we called this, so just
                    // stop the current computation and return with
                    // incomplete results.
                    return;
                  }
              }

            // xi_{x x}
            d2xidxyz2_map[p][0] = -numer * dxidx_map[p]*dxidx_map[p] / denom;

            // xi_{x y}
            d2xidxyz2_map[p][1] = -numer * dxidx_map[p]*dxidy_map[p] / denom;

            // xi_{y y}
            d2xidxyz2_map[p][3] = -numer * dxidy_map[p]*dxidy_map[p] / denom;

#elif LIBMESH_DIM == 3
            // Compute inverse map second derivatives for 1D-element-living-in-3D case
            // See JWP notes for details

            // numer = x_xi*x_{xi xi} + y_xi*y_{xi xi} + z_xi*z_{xi xi}
            auto numer =
              dxyzdxi_map[p](0)*d2xyzdxi2_map[p](0) +
              dxyzdxi_map[p](1)*d2xyzdxi2_map[p](1) +
              dxyzdxi_map[p](2)*d2xyzdxi2_map[p](2);

            // denom = (x_xi)^2 + (y_xi)^2 + (z_xi)^2 must be >= 0.0
            auto denom =
              dxyzdxi_map[p](0)*dxyzdxi_map[p](0) +
              dxyzdxi_map[p](1)*dxyzdxi_map[p](1) +
              dxyzdxi_map[p](2)*dxyzdxi_map[p](2);

            if (denom <= 0.0)
              {
                // Don't call print_info() recursively if we're already
                // failing.  print_info() calls Elem::volume() which may
                // call FE::reinit() and trigger the same failure again.
                static bool failing = false;
                if (!failing)
                  {
                    failing = true;
                    elem->print_info(libMesh::err);
                    failing = false;
                    libmesh_error_msg("Encountered invalid 1D element!");
                  }
                else
                  {
                    // We were already failing when we called this, so just
                    // stop the current computation and return with
                    // incomplete results.
                    return;
                  }
              }

            // xi_{x x}
            d2xidxyz2_map[p][0] = -numer * dxidx_map[p]*dxidx_map[p] / denom;

            // xi_{x y}
            d2xidxyz2_map[p][1] = -numer * dxidx_map[p]*dxidy_map[p] / denom;

            // xi_{x z}
            d2xidxyz2_map[p][2] = -numer * dxidx_map[p]*dxidz_map[p] / denom;

            // xi_{y y}
            d2xidxyz2_map[p][3] = -numer * dxidy_map[p]*dxidy_map[p] / denom;

            // xi_{y z}
            d2xidxyz2_map[p][4] = -numer * dxidy_map[p]*dxidz_map[p] / denom;

            // xi_{z z}
            d2xidxyz2_map[p][5] = -numer * dxidz_map[p]*dxidz_map[p] / denom;
#endif //LIBMESH_DIM == 3
          } // calculate_d2xyz

#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES

        // done computing the map
        break;
      }


      //--------------------------------------------------------------------
      // 2D
    case 2:
      {
        //------------------------------------------------------------------
        // Compute the (x,y) values at the quadrature points,
        // the Jacobian at the quadrature points

        if (calculate_xyz)
          xyz[p].zero();

        if (calculate_dxyz)
          {
            dxyzdxi_map[p].zero();
            dxyzdeta_map[p].zero();
          }
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
        if (calculate_d2xyz)
          {
            d2xyzdxi2_map[p].zero();
            d2xyzdxideta_map[p].zero();
            d2xyzdeta2_map[p].zero();
            // Inverse map second derivatives
            d2xidxyz2_map[p].assign(6, 0.);
            d2etadxyz2_map[p].assign(6, 0.);
          }
#endif


        // compute (x,y) at the quadrature points, derivatives once
        for (auto i : index_range(elem_nodes)) // sum over the nodes
          {
            // Reference to the point, helps eliminate
            // excessive temporaries in the inner loop
            libmesh_assert(elem_nodes[i]);
            const Point & elem_point = *elem_nodes[i];

            if (calculate_xyz)
              xyz[p].add_scaled          (elem_point, phi_map[i][p]     );

            if (calculate_dxyz)
              {
                dxyzdxi_map[p].add_scaled      (elem_point, dphidxi_map[i][p] );
                dxyzdeta_map[p].add_scaled     (elem_point, dphideta_map[i][p]);
              }
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
            if (calculate_d2xyz)
              {
                d2xyzdxi2_map[p].add_scaled    (elem_point, d2phidxi2_map[i][p]);
                d2xyzdxideta_map[p].add_scaled (elem_point, d2phidxideta_map[i][p]);
                d2xyzdeta2_map[p].add_scaled   (elem_point, d2phideta2_map[i][p]);
              }
#endif
          }

        if (calculate_dxyz)
          {
            // compute the jacobian once
            const auto dx_dxi = dxdxi_map(p),
              dx_deta = dxdeta_map(p),
              dy_dxi = dydxi_map(p),
              dy_deta = dydeta_map(p);

#if LIBMESH_DIM == 2
            // Compute the Jacobian.  This assumes the 2D face
            // lives in 2D space
            //
            // Symbolically, the matrix determinant is
            //
            //         | dx/dxi  dx/deta |
            // jac =   | dy/dxi  dy/deta |
            //
            // jac = dx/dxi*dy/deta - dx/deta*dy/dxi
            jac[p] = (dx_dxi*dy_deta - dx_deta*dy_dxi);

            if (jac[p] <= jacobian_tolerance)
              {
                // Don't call print_info() recursively if we're already
                // failing.  print_info() calls Elem::volume() which may
                // call FE::reinit() and trigger the same failure again.
                static bool failing = false;
                if (!failing)
                  {
                    failing = true;
                    elem->print_info(libMesh::err);
                    failing = false;
                    if (calculate_xyz)
                      {
                        libmesh_error_msg("ERROR: negative Jacobian " \
                                          << jac[p] \
                                          << " at point " \
                                          << xyz[p] \
                                          << " in element " \
                                          << elem->id());
                      }
                    else
                      {
                        // In this case xyz[p] is not defined, so don't
                        // try to print it out.
                        libmesh_error_msg("ERROR: negative Jacobian " \
                                          << jac[p] \
                                          << " at point index " \
                                          << p \
                                          << " in element " \
                                          << elem->id());
                      }
                  }
                else
                  {
                    // We were already failing when we called this, so just
                    // stop the current computation and return with
                    // incomplete results.
                    return;
                  }
              }

            JxW[p] = jac[p]*qw[p];

            // Compute the shape function derivatives wrt x,y at the
            // quadrature points
            const auto inv_jac = 1./jac[p];

            dxidx_map[p]  =  dy_deta*inv_jac; //dxi/dx  =  (1/J)*dy/deta
            dxidy_map[p]  = -dx_deta*inv_jac; //dxi/dy  = -(1/J)*dx/deta
            detadx_map[p] = -dy_dxi* inv_jac; //deta/dx = -(1/J)*dy/dxi
            detady_map[p] =  dx_dxi* inv_jac; //deta/dy =  (1/J)*dx/dxi

            dxidz_map[p] = detadz_map[p] = 0.;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
            if (compute_second_derivatives)
              this->compute_inverse_map_second_derivs(p);
#endif
#else // LIBMESH_DIM == 3

            const auto dz_dxi = dzdxi_map(p),
              dz_deta = dzdeta_map(p);

            // Compute the Jacobian.  This assumes a 2D face in
            // 3D space.
            //
            // The transformation matrix T from local to global
            // coordinates is
            //
            //         | dx/dxi  dx/deta |
            //     T = | dy/dxi  dy/deta |
            //         | dz/dxi  dz/deta |
            // note det(T' T) = det(T')det(T) = det(T)det(T)
            // so det(T) = std::sqrt(det(T' T))
            //
            //----------------------------------------------
            // Notes:
            //
            //       dX = R dXi -> R'dX = R'R dXi
            // (R^-1)dX =   dXi    [(R'R)^-1 R']dX = dXi
            //
            // so R^-1 = (R'R)^-1 R'
            //
            // and R^-1 R = (R'R)^-1 R'R = I.
            //
            const auto g11 = (dx_dxi*dx_dxi +
                              dy_dxi*dy_dxi +
                              dz_dxi*dz_dxi);

            const auto g12 = (dx_dxi*dx_deta +
                              dy_dxi*dy_deta +
                              dz_dxi*dz_deta);

            const auto g21 = g12;

            const auto g22 = (dx_deta*dx_deta +
                              dy_deta*dy_deta +
                              dz_deta*dz_deta);

            const auto det = (g11*g22 - g12*g21);

            if (det <= 0.)
              {
                // Don't call print_info() recursively if we're already
                // failing.  print_info() calls Elem::volume() which may
                // call FE::reinit() and trigger the same failure again.
                static bool failing = false;
                if (!failing)
                  {
                    failing = true;
                    elem->print_info(libMesh::err);
                    failing = false;
                    if (calculate_xyz)
                      {
                        libmesh_error_msg("ERROR: negative Jacobian " \
                                          << det \
                                          << " at point " \
                                          << xyz[p] \
                                          << " in element " \
                                          << elem->id());
                      }
                    else
                      {
                        // In this case xyz[p] is not defined, so don't
                        // try to print it out.
                        libmesh_error_msg("ERROR: negative Jacobian " \
                                          << det \
                                          << " at point index " \
                                          << p \
                                          << " in element " \
                                          << elem->id());
                      }
                  }
                else
                  {
                    // We were already failing when we called this, so just
                    // stop the current computation and return with
                    // incomplete results.
                    return;
                  }
              }

            const auto inv_det = 1./det;
            jac[p] = std::sqrt(det);

            JxW[p] = jac[p]*qw[p];

            const auto g11inv =  g22*inv_det;
            const auto g12inv = -g12*inv_det;
            const auto g21inv = -g21*inv_det;
            const auto g22inv =  g11*inv_det;

            dxidx_map[p]  = g11inv*dx_dxi + g12inv*dx_deta;
            dxidy_map[p]  = g11inv*dy_dxi + g12inv*dy_deta;
            dxidz_map[p]  = g11inv*dz_dxi + g12inv*dz_deta;

            detadx_map[p] = g21inv*dx_dxi + g22inv*dx_deta;
            detady_map[p] = g21inv*dy_dxi + g22inv*dy_deta;
            detadz_map[p] = g21inv*dz_dxi + g22inv*dz_deta;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

            if (calculate_d2xyz)
              {
                // Compute inverse map second derivative values for
                // 2D-element-living-in-3D case.  We pursue a least-squares
                // solution approach for this "non-square" case, see JWP notes
                // for details.

                // A = [ x_{xi xi} x_{eta eta} ]
                //     [ y_{xi xi} y_{eta eta} ]
                //     [ z_{xi xi} z_{eta eta} ]
                DenseMatrix<RealType> A(3,2);
                A(0,0) = d2xyzdxi2_map[p](0);  A(0,1) = d2xyzdeta2_map[p](0);
                A(1,0) = d2xyzdxi2_map[p](1);  A(1,1) = d2xyzdeta2_map[p](1);
                A(2,0) = d2xyzdxi2_map[p](2);  A(2,1) = d2xyzdeta2_map[p](2);

                // J^T, the transpose of the Jacobian matrix
                DenseMatrix<RealType> JT(2,3);
                JT(0,0) = dx_dxi;   JT(0,1) = dy_dxi;   JT(0,2) = dz_dxi;
                JT(1,0) = dx_deta;  JT(1,1) = dy_deta;  JT(1,2) = dz_deta;

                // (J^T J)^(-1), this has already been computed for us above...
                DenseMatrix<RealType> JTJinv(2,2);
                JTJinv(0,0) = g11inv;  JTJinv(0,1) = g12inv;
                JTJinv(1,0) = g21inv;  JTJinv(1,1) = g22inv;

                // Some helper variables
                RealVectorValue
                  dxi  (dxidx_map[p],   dxidy_map[p],   dxidz_map[p]),
                  deta (detadx_map[p],  detady_map[p],  detadz_map[p]);

                // To be filled in below
                DenseVector<RealType> tmp1(2);
                DenseVector<RealType> tmp2(3);
                DenseVector<RealType> tmp3(2);

                // For (s,t) in {(x,x), (x,y), (x,z), (y,y), (y,z), (z,z)}, compute the
                // vector of inverse map second derivatives [xi_{s t}, eta_{s t}]
                unsigned ctr=0;
                for (unsigned s=0; s<3; ++s)
                  for (unsigned t=s; t<3; ++t)
                    {
                      // Construct tmp1 = [xi_s*xi_t, eta_s*eta_t]
                      tmp1(0) = dxi(s)*dxi(t);
                      tmp1(1) = deta(s)*deta(t);

                      // Compute tmp2 = A * tmp1
                      A.vector_mult(tmp2, tmp1);

                      // Compute scalar value "alpha"
                      auto alpha = dxi(s)*deta(t) + deta(s)*dxi(t);

                      // Compute tmp2 <- tmp2 + alpha * x_{xi eta}
                      for (unsigned i=0; i<3; ++i)
                        tmp2(i) += alpha*d2xyzdxideta_map[p](i);

                      // Compute tmp3 = J^T * tmp2
                      JT.vector_mult(tmp3, tmp2);

                      // Compute tmp1 = (J^T J)^(-1) * tmp3.  tmp1 is available for us to reuse.
                      JTJinv.vector_mult(tmp1, tmp3);

                      // Fill in appropriate entries, don't forget to multiply by -1!
                      d2xidxyz2_map[p][ctr]  = -tmp1(0);
                      d2etadxyz2_map[p][ctr] = -tmp1(1);

                      // Increment the counter
                      ctr++;
                    }
              }

#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES

#endif // LIBMESH_DIM == 3
          }
        // done computing the map
        break;
      }



      //--------------------------------------------------------------------
      // 3D
    case 3:
      {
        //------------------------------------------------------------------
        // Compute the (x,y,z) values at the quadrature points,
        // the Jacobian at the quadrature point

        // Clear the entities that will be summed
        if (calculate_xyz)
          xyz[p].zero           ();
        if (calculate_dxyz)
          {
            dxyzdxi_map[p].zero   ();
            dxyzdeta_map[p].zero  ();
            dxyzdzeta_map[p].zero ();
          }
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
        if (calculate_d2xyz)
          {
            d2xyzdxi2_map[p].zero();
            d2xyzdxideta_map[p].zero();
            d2xyzdxidzeta_map[p].zero();
            d2xyzdeta2_map[p].zero();
            d2xyzdetadzeta_map[p].zero();
            d2xyzdzeta2_map[p].zero();
            // Inverse map second derivatives
            d2xidxyz2_map[p].assign(6, 0.);
            d2etadxyz2_map[p].assign(6, 0.);
            d2zetadxyz2_map[p].assign(6, 0.);
          }
#endif


        // compute (x,y,z) at the quadrature points,
        // dxdxi,   dydxi,   dzdxi,
        // dxdeta,  dydeta,  dzdeta,
        // dxdzeta, dydzeta, dzdzeta  all once
        for (auto i : index_range(elem_nodes)) // sum over the nodes
          {
            // Reference to the point, helps eliminate
            // excessive temporaries in the inner loop
            libmesh_assert(elem_nodes[i]);
            const Point & elem_point = *elem_nodes[i];

            if (calculate_xyz)
              xyz[p].add_scaled           (elem_point, phi_map[i][p]      );
            if (calculate_dxyz)
              {
                dxyzdxi_map[p].add_scaled   (elem_point, dphidxi_map[i][p]  );
                dxyzdeta_map[p].add_scaled  (elem_point, dphideta_map[i][p] );
                dxyzdzeta_map[p].add_scaled (elem_point, dphidzeta_map[i][p]);
              }
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
            if (calculate_d2xyz)
              {
                d2xyzdxi2_map[p].add_scaled      (elem_point,
                                                  d2phidxi2_map[i][p]);
                d2xyzdxideta_map[p].add_scaled   (elem_point,
                                                  d2phidxideta_map[i][p]);
                d2xyzdxidzeta_map[p].add_scaled  (elem_point,
                                                  d2phidxidzeta_map[i][p]);
                d2xyzdeta2_map[p].add_scaled     (elem_point,
                                                  d2phideta2_map[i][p]);
                d2xyzdetadzeta_map[p].add_scaled (elem_point,
                                                  d2phidetadzeta_map[i][p]);
                d2xyzdzeta2_map[p].add_scaled    (elem_point,
                                                  d2phidzeta2_map[i][p]);
              }
#endif
          }

        if (calculate_dxyz)
          {
            // compute the jacobian
            const auto
              dx_dxi   = dxdxi_map(p),   dy_dxi   = dydxi_map(p),   dz_dxi   = dzdxi_map(p),
              dx_deta  = dxdeta_map(p),  dy_deta  = dydeta_map(p),  dz_deta  = dzdeta_map(p),
              dx_dzeta = dxdzeta_map(p), dy_dzeta = dydzeta_map(p), dz_dzeta = dzdzeta_map(p);

            // Symbolically, the matrix determinant is
            //
            //         | dx/dxi   dy/dxi   dz/dxi   |
            // jac =   | dx/deta  dy/deta  dz/deta  |
            //         | dx/dzeta dy/dzeta dz/dzeta |
            //
            // jac = dx/dxi*(dy/deta*dz/dzeta - dz/deta*dy/dzeta) +
            //       dy/dxi*(dz/deta*dx/dzeta - dx/deta*dz/dzeta) +
            //       dz/dxi*(dx/deta*dy/dzeta - dy/deta*dx/dzeta)

            jac[p] = (dx_dxi*(dy_deta*dz_dzeta - dz_deta*dy_dzeta)  +
                      dy_dxi*(dz_deta*dx_dzeta - dx_deta*dz_dzeta)  +
                      dz_dxi*(dx_deta*dy_dzeta - dy_deta*dx_dzeta));

            if (jac[p] <= jacobian_tolerance)
              {
                // Don't call print_info() recursively if we're already
                // failing.  print_info() calls Elem::volume() which may
                // call FE::reinit() and trigger the same failure again.
                static bool failing = false;
                if (!failing)
                  {
                    failing = true;
                    elem->print_info(libMesh::err);
                    failing = false;
                    if (calculate_xyz)
                      {
                        libmesh_error_msg("ERROR: negative Jacobian " \
                                          << jac[p] \
                                          << " at point " \
                                          << xyz[p] \
                                          << " in element " \
                                          << elem->id());
                      }
                    else
                      {
                        // In this case xyz[p] is not defined, so don't
                        // try to print it out.
                        libmesh_error_msg("ERROR: negative Jacobian " \
                                          << jac[p] \
                                          << " at point index " \
                                          << p \
                                          << " in element " \
                                          << elem->id());
                      }
                  }
                else
                  {
                    // We were already failing when we called this, so just
                    // stop the current computation and return with
                    // incomplete results.
                    return;
                  }
              }

            JxW[p] = jac[p]*qw[p];

            // Compute the shape function derivatives wrt x,y at the
            // quadrature points
            const auto inv_jac  = 1./jac[p];

            dxidx_map[p]   = (dy_deta*dz_dzeta - dz_deta*dy_dzeta)*inv_jac;
            dxidy_map[p]   = (dz_deta*dx_dzeta - dx_deta*dz_dzeta)*inv_jac;
            dxidz_map[p]   = (dx_deta*dy_dzeta - dy_deta*dx_dzeta)*inv_jac;

            detadx_map[p]  = (dz_dxi*dy_dzeta  - dy_dxi*dz_dzeta )*inv_jac;
            detady_map[p]  = (dx_dxi*dz_dzeta  - dz_dxi*dx_dzeta )*inv_jac;
            detadz_map[p]  = (dy_dxi*dx_dzeta  - dx_dxi*dy_dzeta )*inv_jac;

            dzetadx_map[p] = (dy_dxi*dz_deta   - dz_dxi*dy_deta  )*inv_jac;
            dzetady_map[p] = (dz_dxi*dx_deta   - dx_dxi*dz_deta  )*inv_jac;
            dzetadz_map[p] = (dx_dxi*dy_deta   - dy_dxi*dx_deta  )*inv_jac;
          }

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
        if (compute_second_derivatives)
          this->compute_inverse_map_second_derivs(p);
#endif
        // done computing the map
        break;
      }

    default:
      libmesh_error_msg("Invalid dim = " << dim);
    }
}



template <typename RealType>
void FEMapTempl<RealType>::resize_quadrature_map_vectors(const unsigned int dim, unsigned int n_qp)
{
  // We're calculating now!
  this->determine_calculations();

  // Resize the vectors to hold data at the quadrature points
  if (calculate_xyz)
    xyz.resize(n_qp);
  if (calculate_dxyz)
    {
      dxyzdxi_map.resize(n_qp);
      dxidx_map.resize(n_qp);
      dxidy_map.resize(n_qp); // 1D element may live in 2D ...
      dxidz_map.resize(n_qp); // ... or 3D
    }
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  if (calculate_d2xyz)
    {
      d2xyzdxi2_map.resize(n_qp);

      // Inverse map second derivatives
      d2xidxyz2_map.resize(n_qp);
      for (auto i : index_range(d2xidxyz2_map))
        d2xidxyz2_map[i].assign(6, 0.);
    }
#endif
  if (dim > 1)
    {
      if (calculate_dxyz)
        {
          dxyzdeta_map.resize(n_qp);
          detadx_map.resize(n_qp);
          detady_map.resize(n_qp);
          detadz_map.resize(n_qp);
        }
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
      if (calculate_d2xyz)
        {
          d2xyzdxideta_map.resize(n_qp);
          d2xyzdeta2_map.resize(n_qp);

          // Inverse map second derivatives
          d2etadxyz2_map.resize(n_qp);
          for (auto i : index_range(d2etadxyz2_map))
            d2etadxyz2_map[i].assign(6, 0.);
        }
#endif
      if (dim > 2)
        {
          if (calculate_dxyz)
            {
              dxyzdzeta_map.resize (n_qp);
              dzetadx_map.resize   (n_qp);
              dzetady_map.resize   (n_qp);
              dzetadz_map.resize   (n_qp);
            }
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
          if (calculate_d2xyz)
            {
              d2xyzdxidzeta_map.resize(n_qp);
              d2xyzdetadzeta_map.resize(n_qp);
              d2xyzdzeta2_map.resize(n_qp);

              // Inverse map second derivatives
              d2zetadxyz2_map.resize(n_qp);
              for (auto i : index_range(d2zetadxyz2_map))
                d2zetadxyz2_map[i].assign(6, 0.);
            }
#endif
        }
    }

  if (calculate_dxyz)
    {
      jac.resize(n_qp);
      JxW.resize(n_qp);
    }
}



template <typename RealType>
void FEMapTempl<RealType>::compute_affine_map(const unsigned int dim,
                               const std::vector<Real> & qw,
                               const Elem * elem)
{
  // Start logging the map computation.
  LOG_SCOPE("compute_affine_map()", "FEMap");

  libmesh_assert(elem);

  const unsigned int n_qp = cast_int<unsigned int>(qw.size());

  // Resize the vectors to hold data at the quadrature points
  this->resize_quadrature_map_vectors(dim, n_qp);

  // Determine the nodes contributing to element elem
  unsigned int n_nodes = elem->n_nodes();
  _elem_nodes.resize(elem->n_nodes());
  for (unsigned int i=0; i<n_nodes; i++)
    _elem_nodes[i] = elem->node_ptr(i);

  // Compute map at quadrature point 0
  this->compute_single_point_map(dim, qw, elem, 0, _elem_nodes, /*compute_second_derivatives=*/false);

  // Compute xyz at all other quadrature points
  if (calculate_xyz)
    for (unsigned int p=1; p<n_qp; p++)
      {
        xyz[p].zero();
        for (auto i : index_range(phi_map)) // sum over the nodes
          xyz[p].add_scaled (*_elem_nodes[i], phi_map[i][p]);
      }

  // Copy other map data from quadrature point 0
  if (calculate_dxyz)
    for (unsigned int p=1; p<n_qp; p++) // for each extra quadrature point
      {
        dxyzdxi_map[p] = dxyzdxi_map[0];
        dxidx_map[p] = dxidx_map[0];
        dxidy_map[p] = dxidy_map[0];
        dxidz_map[p] = dxidz_map[0];
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
        // The map should be affine, so second derivatives are zero
        if (calculate_d2xyz)
          d2xyzdxi2_map[p] = 0.;
#endif
        if (dim > 1)
          {
            dxyzdeta_map[p] = dxyzdeta_map[0];
            detadx_map[p] = detadx_map[0];
            detady_map[p] = detady_map[0];
            detadz_map[p] = detadz_map[0];
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
            if (calculate_d2xyz)
              {
                d2xyzdxideta_map[p] = 0.;
                d2xyzdeta2_map[p] = 0.;
              }
#endif
            if (dim > 2)
              {
                dxyzdzeta_map[p] = dxyzdzeta_map[0];
                dzetadx_map[p] = dzetadx_map[0];
                dzetady_map[p] = dzetady_map[0];
                dzetadz_map[p] = dzetadz_map[0];
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
                if (calculate_d2xyz)
                  {
                    d2xyzdxidzeta_map[p] = 0.;
                    d2xyzdetadzeta_map[p] = 0.;
                    d2xyzdzeta2_map[p] = 0.;
                  }
#endif
              }
          }
        jac[p] = jac[0];
        JxW[p] = JxW[0] / qw[0] * qw[p];
      }
}



template <typename RealType>
void FEMapTempl<RealType>::compute_null_map(const unsigned int dim,
                             const std::vector<Real> & qw)
{
  // Start logging the map computation.
  LOG_SCOPE("compute_null_map()", "FEMap");

  const unsigned int n_qp = cast_int<unsigned int>(qw.size());

  // Resize the vectors to hold data at the quadrature points
  this->resize_quadrature_map_vectors(dim, n_qp);

  // Compute "fake" xyz
  for (unsigned int p=1; p<n_qp; p++)
    {
      if (calculate_xyz)
        xyz[p].zero();

      if (calculate_dxyz)
        {
          dxyzdxi_map[p] = 0;
          dxidx_map[p] = 0;
          dxidy_map[p] = 0;
          dxidz_map[p] = 0;
        }
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
      if (calculate_d2xyz)
        {
          d2xyzdxi2_map[p] = 0;
        }
#endif
      if (dim > 1)
        {
          if (calculate_dxyz)
            {
              dxyzdeta_map[p] = 0;
              detadx_map[p] = 0;
              detady_map[p] = 0;
              detadz_map[p] = 0;
            }
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
          if (calculate_d2xyz)
            {
              d2xyzdxideta_map[p] = 0.;
              d2xyzdeta2_map[p] = 0.;
            }
#endif
          if (dim > 2)
            {
              if (calculate_dxyz)
                {
                  dxyzdzeta_map[p] = 0;
                  dzetadx_map[p] = 0;
                  dzetady_map[p] = 0;
                  dzetadz_map[p] = 0;
                }
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
              if (calculate_d2xyz)
                {
                  d2xyzdxidzeta_map[p] = 0;
                  d2xyzdetadzeta_map[p] = 0;
                  d2xyzdzeta2_map[p] = 0;
                }
#endif
            }
        }
      if (calculate_dxyz)
        {
          jac[p] = 1;
          JxW[p] = qw[p];
        }
    }
}



template <typename RealType>
void FEMapTempl<RealType>::compute_map(const unsigned int dim,
                        const std::vector<Real> & qw,
                        const Elem * elem,
                        bool calculate_d2phi)
{
  if (!elem)
    {
      compute_null_map(dim, qw);
      return;
    }

  if (elem->has_affine_map())
    {
      compute_affine_map(dim, qw, elem);
      return;
    }
#ifndef LIBMESH_ENABLE_SECOND_DERIVATIVES
    libmesh_assert(!calculate_d2phi);
#endif

  // Start logging the map computation.
  LOG_SCOPE("compute_map()", "FEMap");

  libmesh_assert(elem);

  const unsigned int n_qp = cast_int<unsigned int>(qw.size());

  // Resize the vectors to hold data at the quadrature points
  this->resize_quadrature_map_vectors(dim, n_qp);

  // Determine the nodes contributing to element elem
  if (elem->type() == TRI3SUBDIVISION)
    {
      // Subdivision surface FE require the 1-ring around elem
      libmesh_assert_equal_to (dim, 2);
      const Tri3Subdivision * sd_elem = static_cast<const Tri3Subdivision *>(elem);
      MeshTools::Subdivision::find_one_ring(sd_elem, _elem_nodes);
    }
  else
    {
      // All other FE use only the nodes of elem itself
      _elem_nodes.resize(elem->n_nodes(), nullptr);
      for (auto i : elem->node_index_range())
        _elem_nodes[i] = elem->node_ptr(i);
    }

  // Compute map at all quadrature points
  for (unsigned int p=0; p!=n_qp; p++)
    this->compute_single_point_map(dim, qw, elem, p, _elem_nodes, calculate_d2phi);
}



template <typename RealType>
void FEMapTempl<RealType>::print_JxW(std::ostream & os) const
{
  for (auto i : index_range(JxW))
    os << " [" << i << "]: " << JxW[i] << std::endl;
}



template <typename RealType>
void FEMapTempl<RealType>::print_xyz(std::ostream & os) const
{
  for (auto i : index_range(xyz))
    os << " [" << i << "]: " << xyz[i];
}



template <typename RealType>
void FEMapTempl<RealType>::compute_inverse_map_second_derivs(unsigned p)
{
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  // Only certain second derivatives are valid depending on the
  // dimension...
  std::set<unsigned> valid_indices;

  // Construct J^{-1}, A, and B matrices (see JWP's notes for details)
  // for cases in which the element dimension matches LIBMESH_DIM.
#if LIBMESH_DIM==1
  RealTensor
    Jinv(dxidx_map[p],  0.,  0.,
         0.,            0.,  0.,
         0.,            0.,  0.),

    A(d2xyzdxi2_map[p](0), 0., 0.,
      0.,                  0., 0.,
      0.,                  0., 0.),

    B(0., 0., 0.,
      0., 0., 0.,
      0., 0., 0.);

  RealVectorValue
    dxi  (dxidx_map[p], 0., 0.),
    deta (0.,           0., 0.),
    dzeta(0.,           0., 0.);

  // In 1D, we have only the xx second derivative
  valid_indices.insert(0);

#elif LIBMESH_DIM==2
  RealTensor
    Jinv(dxidx_map[p],  dxidy_map[p],  0.,
         detadx_map[p], detady_map[p], 0.,
         0.,            0.,            0.),

    A(d2xyzdxi2_map[p](0), d2xyzdeta2_map[p](0), 0.,
      d2xyzdxi2_map[p](1), d2xyzdeta2_map[p](1), 0.,
      0.,                  0.,                   0.),

    B(d2xyzdxideta_map[p](0), 0., 0.,
      d2xyzdxideta_map[p](1), 0., 0.,
      0.,                     0., 0.);

  RealVectorValue
    dxi  (dxidx_map[p],  dxidy_map[p],  0.),
    deta (detadx_map[p], detady_map[p], 0.),
    dzeta(0.,            0.,            0.);

  // In 2D, we have xx, xy, and yy second derivatives
  const unsigned tmp[3] = {0,1,3};
  valid_indices.insert(tmp, tmp+3);

#elif LIBMESH_DIM==3
  RealTensor
    Jinv(dxidx_map[p],   dxidy_map[p],   dxidz_map[p],
         detadx_map[p],  detady_map[p],  detadz_map[p],
         dzetadx_map[p], dzetady_map[p], dzetadz_map[p]),

    A(d2xyzdxi2_map[p](0), d2xyzdeta2_map[p](0), d2xyzdzeta2_map[p](0),
      d2xyzdxi2_map[p](1), d2xyzdeta2_map[p](1), d2xyzdzeta2_map[p](1),
      d2xyzdxi2_map[p](2), d2xyzdeta2_map[p](2), d2xyzdzeta2_map[p](2)),

    B(d2xyzdxideta_map[p](0), d2xyzdxidzeta_map[p](0), d2xyzdetadzeta_map[p](0),
      d2xyzdxideta_map[p](1), d2xyzdxidzeta_map[p](1), d2xyzdetadzeta_map[p](1),
      d2xyzdxideta_map[p](2), d2xyzdxidzeta_map[p](2), d2xyzdetadzeta_map[p](2));

  RealVectorValue
    dxi  (dxidx_map[p],   dxidy_map[p],   dxidz_map[p]),
    deta (detadx_map[p],  detady_map[p],  detadz_map[p]),
    dzeta(dzetadx_map[p], dzetady_map[p], dzetadz_map[p]);

  // In 3D, we have xx, xy, xz, yy, yz, and zz second derivatives
  const unsigned tmp[6] = {0,1,2,3,4,5};
  valid_indices.insert(tmp, tmp+6);

#endif

  // For (s,t) in {(x,x), (x,y), (x,z), (y,y), (y,z), (z,z)}, compute the
  // vector of inverse map second derivatives [xi_{s t}, eta_{s t}, zeta_{s t}]
  unsigned ctr=0;
  for (unsigned s=0; s<3; ++s)
    for (unsigned t=s; t<3; ++t)
      {
        if (valid_indices.count(ctr))
          {
            RealVectorValue
              v1(dxi(s)*dxi(t),
                 deta(s)*deta(t),
                 dzeta(s)*dzeta(t)),

              v2(dxi(s)*deta(t) + deta(s)*dxi(t),
                 dxi(s)*dzeta(t) + dzeta(s)*dxi(t),
                 deta(s)*dzeta(t) + dzeta(s)*deta(t));

            // Compute the inverse map second derivatives
            RealVectorValue v3 = -Jinv*(A*v1 + B*v2);

            // Store them in the appropriate locations in the class data structures
            d2xidxyz2_map[p][ctr] = v3(0);

            if (LIBMESH_DIM > 1)
              d2etadxyz2_map[p][ctr] = v3(1);

            if (LIBMESH_DIM > 2)
              d2zetadxyz2_map[p][ctr] = v3(2);
          }

        // Increment the counter
        ctr++;
      }
#else
   // to avoid compiler warnings:
   libmesh_ignore(p);
#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES
}



template <typename RealType>
PointTempl<RealType> FEMapTempl<RealType>::inverse_map (const unsigned int dim,
                                         const Elem * elem,
                                         const Point & physical_point,
                                         const Real tolerance,
                                         const bool secure)
{
  libmesh_assert(elem);
  libmesh_assert_greater_equal (tolerance, 0.);

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
  if (elem->infinite())
    return InfFEMap::inverse_map(dim, elem, physical_point, tolerance,
                                 secure);
#endif

  // Start logging the map inversion.
  LOG_SCOPE("inverse_map()", "FEMap");

  // How much did the point on the reference
  // element change by in this Newton step?
  Real inverse_map_error = 0.;

  //  The point on the reference element.  This is
  //  the "initial guess" for Newton's method.  The
  //  centroid seems like a good idea, but computing
  //  it is a little more intensive than, say taking
  //  the zero point.
  //
  //  Convergence should be insensitive of this choice
  //  for "good" elements.
  Point p; // the zero point.  No computation required

  //  The number of iterations in the map inversion process.
  unsigned int cnt = 0;

  //  The number of iterations after which we give up and declare
  //  divergence
  const unsigned int max_cnt = 10;

  //  The distance (in master element space) beyond which we give up
  //  and declare divergence.  This is no longer used...
  // Real max_step_length = 4.;



  //  Newton iteration loop.
  do
    {
      //  Where our current iterate \p p maps to.
      const Point physical_guess = map(dim, elem, p);

      //  How far our current iterate is from the actual point.
      const Point delta = physical_point - physical_guess;

      //  Increment in current iterate \p p, will be computed.
      Point dp;


      //  The form of the map and how we invert it depends
      //  on the dimension that we are in.
      switch (dim)
        {
          // ------------------------------------------------------------------
          //  0D map inversion is trivial
        case 0:
          {
            break;
          }

          // ------------------------------------------------------------------
          //  1D map inversion
          //
          //  Here we find the point on a 1D reference element that maps to
          //  the point \p physical_point in the domain.  This is a bit tricky
          //  since we do not want to assume that the point \p physical_point
          //  is also in a 1D domain.  In particular, this method might get
          //  called on the edge of a 3D element, in which case
          //  \p physical_point actually lives in 3D.
        case 1:
          {
            const Point dxi = map_deriv (dim, elem, 0, p);

            //  Newton's method in this case looks like
            //
            //  {X} - {X_n} = [J]*dp
            //
            //  Where {X}, {X_n} are 3x1 vectors, [J] is a 3x1 matrix
            //  d(x,y,z)/dxi, and we seek dp, a scalar.  Since the above
            //  system is either overdetermined or rank-deficient, we will
            //  solve the normal equations for this system
            //
            //  [J]^T ({X} - {X_n}) = [J]^T [J] {dp}
            //
            //  which involves the trivial inversion of the scalar
            //  G = [J]^T [J]
            const auto G = dxi*dxi;

            if (secure)
              libmesh_assert_greater (G, 0.);

            const auto Ginv = 1./G;

            const auto  dxidelta = dxi*delta;

            dp(0) = Ginv*dxidelta;

            // No master elements have radius > 4, but sometimes we
            // can take a step that big while still converging
            // if (secure)
            // libmesh_assert_less (dp.size(), max_step_length);

            break;
          }



          // ------------------------------------------------------------------
          //  2D map inversion
          //
          //  Here we find the point on a 2D reference element that maps to
          //  the point \p physical_point in the domain.  This is a bit tricky
          //  since we do not want to assume that the point \p physical_point
          //  is also in a 2D domain.  In particular, this method might get
          //  called on the face of a 3D element, in which case
          //  \p physical_point actually lives in 3D.
        case 2:
          {
            const Point dxi  = map_deriv (dim, elem, 0, p);
            const Point deta = map_deriv (dim, elem, 1, p);

            //  Newton's method in this case looks like
            //
            //  {X} - {X_n} = [J]*{dp}
            //
            //  Where {X}, {X_n} are 3x1 vectors, [J] is a 3x2 matrix
            //  d(x,y,z)/d(xi,eta), and we seek {dp}, a 2x1 vector.  Since
            //  the above system is either over-determined or rank-deficient,
            //  we will solve the normal equations for this system
            //
            //  [J]^T ({X} - {X_n}) = [J]^T [J] {dp}
            //
            //  which involves the inversion of the 2x2 matrix
            //  [G] = [J]^T [J]
            const auto
              G11 = dxi*dxi,  G12 = dxi*deta,
              G21 = dxi*deta, G22 = deta*deta;


            const auto det = (G11*G22 - G12*G21);

            if (secure)
              libmesh_assert_not_equal_to (det, 0.);

            const auto inv_det = 1./det;

            const auto
              Ginv11 =  G22*inv_det,
              Ginv12 = -G12*inv_det,

              Ginv21 = -G21*inv_det,
              Ginv22 =  G11*inv_det;


            const auto  dxidelta  = dxi*delta;
            const auto  detadelta = deta*delta;

            dp(0) = (Ginv11*dxidelta + Ginv12*detadelta);
            dp(1) = (Ginv21*dxidelta + Ginv22*detadelta);

            // No master elements have radius > 4, but sometimes we
            // can take a step that big while still converging
            // if (secure)
            // libmesh_assert_less (dp.size(), max_step_length);

            break;
          }



          // ------------------------------------------------------------------
          //  3D map inversion
          //
          //  Here we find the point in a 3D reference element that maps to
          //  the point \p physical_point in a 3D domain. Nothing special
          //  has to happen here, since (unless the map is singular because
          //  you have a BAD element) the map will be invertible and we can
          //  apply Newton's method directly.
        case 3:
          {
            const Point dxi   = map_deriv (dim, elem, 0, p);
            const Point deta  = map_deriv (dim, elem, 1, p);
            const Point dzeta = map_deriv (dim, elem, 2, p);

            //  Newton's method in this case looks like
            //
            //  {X} = {X_n} + [J]*{dp}
            //
            //  Where {X}, {X_n} are 3x1 vectors, [J] is a 3x3 matrix
            //  d(x,y,z)/d(xi,eta,zeta), and we seek {dp}, a 3x1 vector.
            //  Since the above system is nonsingular for invertible maps
            //  we will solve
            //
            //  {dp} = [J]^-1 ({X} - {X_n})
            //
            //  which involves the inversion of the 3x3 matrix [J]
            libmesh_try
              {
                RealTensorValue(dxi(0), deta(0), dzeta(0),
                                dxi(1), deta(1), dzeta(1),
                                dxi(2), deta(2), dzeta(2)).solve(delta, dp);
              }
            libmesh_catch (ConvergenceFailure &)
              {
                // We encountered a singular Jacobian.  The value of
                // dp is zero, since it was never changed during the
                // call to RealTensorValue::solve().  We don't want to
                // continue iterating until max_cnt since there is no
                // update to the Newton iterate, and we don't want to
                // print the inverse_map_error value since it will
                // confusingly be 0.  Therefore, in the secure case we
                // need to throw an error message while in the !secure
                // case we can just return a far away point.
                if (secure)
                  {
                    libMesh::err << "ERROR: Newton scheme encountered a singular Jacobian in element: "
                                 << elem->id()
                                 << std::endl;

                    elem->print_info(libMesh::err);

                    libmesh_error_msg("Exiting...");
                  }
                else
                  {
                    for (unsigned int i=0; i != dim; ++i)
                      p(i) = 1e6;
                    return p;
                  }
              }

            // No master elements have radius > 4, but sometimes we
            // can take a step that big while still converging
            // if (secure)
            // libmesh_assert_less (dp.size(), max_step_length);

            break;
          }


          //  Some other dimension?
        default:
          libmesh_error_msg("Invalid dim = " << dim);
        } // end switch(Dim), dp now computed



      //  ||P_n+1 - P_n||
      inverse_map_error = dp.norm();

      //  P_n+1 = P_n + dp
      p.add (dp);

      //  Increment the iteration count.
      cnt++;

      //  Watch for divergence of Newton's
      //  method.  Here's how it goes:
      //  (1) For good elements, we expect convergence in 10
      //      iterations, with no too-large steps.
      //      - If called with (secure == true) and we have not yet converged
      //        print out a warning message.
      //      - If called with (secure == true) and we have not converged in
      //        20 iterations abort
      //  (2) This method may be called in cases when the target point is not
      //      inside the element and we have no business expecting convergence.
      //      For these cases if we have not converged in 10 iterations forget
      //      about it.
      if (cnt > max_cnt)
        {
          //  Warn about divergence when secure is true - this
          //  shouldn't happen
          if (secure)
            {
              // Print every time in devel/dbg modes
#ifndef NDEBUG
              libmesh_here();
              libMesh::err << "WARNING: Newton scheme has not converged in "
                           << cnt << " iterations:" << std::endl
                           << "   physical_point="
                           << physical_point
                           << "   physical_guess="
                           << physical_guess
                           << "   dp="
                           << dp
                           << "   p="
                           << p
                           << "   error=" << inverse_map_error
                           << "   in element " << elem->id()
                           << std::endl;

              elem->print_info(libMesh::err);
#else
              // In optimized mode, just print once that an inverse_map() call
              // had trouble converging its Newton iteration.
              libmesh_do_once(libMesh::err << "WARNING: At least one element took more than "
                              << max_cnt
                              << " iterations to converge in inverse_map()...\n"
                              << "Rerun in devel/dbg mode for more details."
                              << std::endl;);

#endif // NDEBUG

              if (cnt > 2*max_cnt)
                {
                  libMesh::err << "ERROR: Newton scheme FAILED to converge in "
                               << cnt
                               << " iterations in element "
                               << elem->id()
                               << " for physical point = "
                               << physical_point
                               << std::endl;

                  elem->print_info(libMesh::err);

                  libmesh_error_msg("Exiting...");
                }
            }
          //  Return a far off point when secure is false - this
          //  should only happen when we're trying to map a point
          //  that's outside the element
          else
            {
              for (unsigned int i=0; i != dim; ++i)
                p(i) = 1e6;

              return p;
            }
        }
    }
  while (inverse_map_error > tolerance);



  //  If we are in debug mode do two sanity checks.
#ifdef DEBUG

  if (secure)
    {
      // Make sure the point \p p on the reference element actually
      // does map to the point \p physical_point within a tolerance.

      const Point check = map (dim, elem, p);
      const Point diff  = physical_point - check;

      if (diff.norm() > tolerance)
        {
          libmesh_here();
          libMesh::err << "WARNING:  diff is "
                       << diff.norm()
                       << std::endl
                       << " point="
                       << physical_point;
          libMesh::err << " local=" << check;
          libMesh::err << " lref= " << p;

          elem->print_info(libMesh::err);
        }

      // Make sure the point \p p on the reference element actually
      // is

      if (!FEAbstract<RealType>::on_reference_element(p, elem->type(), 2*tolerance))
        {
          libmesh_here();
          libMesh::err << "WARNING:  inverse_map of physical point "
                       << physical_point
                       << " is not on element." << '\n';
          elem->print_info(libMesh::err);
        }
    }

#endif

  return p;
}



template <typename RealType>
void FEMapTempl<RealType>::inverse_map (const unsigned int dim,
                         const Elem * elem,
                         const std::vector<Point> & physical_points,
                         std::vector<Point> &       reference_points,
                         const Real tolerance,
                         const bool secure)
{
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
  if (elem->infinite())
    {
      InfFEMap::inverse_map(dim, elem, physical_points, reference_points, tolerance, secure);
      return;
      // libmesh_not_implemented();
    }
#endif

  // The number of points to find the
  // inverse map of
  const std::size_t n_points = physical_points.size();

  // Resize the vector to hold the points
  // on the reference element
  reference_points.resize(n_points);

  // Find the coordinates on the reference
  // element of each point in physical space
  for (std::size_t p=0; p<n_points; p++)
    reference_points[p] =
      inverse_map (dim, elem, physical_points[p], tolerance, secure);
}



template <typename RealType>
PointTempl<RealType> FEMapTempl<RealType>::map (const unsigned int dim,
                                 const Elem * elem,
                                 const Point & reference_point)
{
  libmesh_assert(elem);

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
  if (elem->infinite())
    return InfFEMap::map(dim, elem, reference_point);
#endif

  Point p;

  const FEFamily mapping_family = FEMap::map_fe_type(*elem);
  const ElemType type     = elem->type();
  const Order order       = elem->default_order();
  const FEType fe_type (order, mapping_family);

  const unsigned int n_sf = FEInterface::n_shape_functions(dim, fe_type, type);

  FEInterface::shape_ptr<RealType> shape_ptr =
    FEInterface::shape_function<RealType>(dim, fe_type);

  // Lagrange basis functions are used for mapping
  for (unsigned int i=0; i<n_sf; i++)
    p.add_scaled (elem->point(i),
                  shape_ptr(elem, order, i, reference_point, false));

  return p;
}



template <typename RealType>
PointTempl<RealType> FEMapTempl<RealType>::map_deriv (const unsigned int dim,
                                       const Elem * elem,
                                       const unsigned int j,
                                       const Point & reference_point)
{
  libmesh_assert(elem);

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
  if (elem->infinite())
    libmesh_not_implemented();
#endif

  Point p;

  const FEFamily mapping_family = FEMap::map_fe_type(*elem);
  const ElemType type     = elem->type();
  const Order order       = elem->default_order();
  const FEType fe_type (order, mapping_family);
  const unsigned int n_sf = FEInterface::n_shape_functions(dim, fe_type, type);

  FEInterface::shape_deriv_ptr<RealType> shape_deriv_ptr =
    FEInterface::shape_deriv_function<RealType>(dim, fe_type);

  // Lagrange basis functions are used for mapping
  for (unsigned int i=0; i<n_sf; i++)
    p.add_scaled (elem->point(i),
                  shape_deriv_ptr(elem, order, i, j, reference_point,
                                  false));

  return p;
}

} // namespace libMesh

#endif // LIBMESH_FE_MAP_IMPL_H
