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

#ifndef LIBMESH_HCURL_FE_TRANSFORMATION_IMPL_H
#define LIBMESH_HCURL_FE_TRANSFORMATION_IMPL_H

#include "libmesh/hcurl_fe_transformation.h"
#include "libmesh/fe_interface.h"
#include "libmesh/int_range.h"
#include "libmesh/fe_base.h"

namespace libMesh
{

namespace {
template <typename MathType, typename RealType>
struct HCurlShim;

template <typename RealType>
struct HCurlShim<RealVectorValue, RealType>
{
  typedef ElemTempl<RealType> Elem;
  typedef PointTempl<RealType> Point;
  typedef typename FEOutputTwoType<RealVectorValue,RealType>::OutputShape OutputShape;

  static void init_map_phi(const FEGenericBase<RealVectorValue,RealType> & fe)
    {
      fe.get_fe_map().get_dxidx();
    };

  static void init_map_dphi(const FEGenericBase<RealVectorValue,RealType> & fe)
    {
      fe.get_fe_map().get_dxidx();
    }

  static void init_map_d2phi(const FEGenericBase<RealVectorValue,RealType> & fe)
    {
      fe.get_fe_map().get_dxidx();
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
      fe.get_fe_map().get_d2xidxyz2();
#endif
    }

  static void map_phi(const unsigned int dim,
                      const Elem * const elem,
                      const std::vector<Point> & qp,
                      const FEGenericBase<RealVectorValue,RealType> & fe,
                      std::vector<std::vector<OutputShape>> & phi)
    {
      switch (dim)
      {
        case 0:
        case 1:
          libmesh_error_msg("These element transformations only make sense in 2D and 3D.");

        case 2:
        {
          const auto & dxidx_map = fe.get_fe_map().get_dxidx();
          const auto & dxidy_map = fe.get_fe_map().get_dxidy();

          const auto & detadx_map = fe.get_fe_map().get_detadx();
          const auto & detady_map = fe.get_fe_map().get_detady();

          // FIXME: Need to update for 2D elements in 3D space
          // phi = (dx/dxi)^-T * \hat{phi}
          // In 2D:
          // (dx/dxi)^{-1} = [  dxi/dx   dxi/dy
          // deta/dx  deta/dy ]
          //
          // so: dxi/dx^{-T} * \hat{phi} = [ dxi/dx  deta/dx   [ \hat{phi}_xi
          // dxi/dy  deta/dy ]   \hat{phi}_eta ]
          //
          // or in indicial notation:  phi_j = xi_{i,j}*\hat{phi}_i
          for (auto i : index_range(phi))
            for (auto p : index_range(phi[i]))
            {
              // Need to temporarily cache reference shape functions
              // TODO: PB: Might be worth trying to build phi_ref separately to see
              //           if we can get vectorization
              OutputShape phi_ref;
              FEInterface::shape(2, fe.get_fe_type(), elem, i, qp[p], phi_ref);

              phi[i][p](0) = dxidx_map[p]*phi_ref.slice(0) + detadx_map[p]*phi_ref.slice(1);

              phi[i][p](1) = dxidy_map[p]*phi_ref.slice(0) + detady_map[p]*phi_ref.slice(1);
            }

          break;
        }

        case 3:
        {
          const auto & dxidx_map = fe.get_fe_map().get_dxidx();
          const auto & dxidy_map = fe.get_fe_map().get_dxidy();
          const auto & dxidz_map = fe.get_fe_map().get_dxidz();

          const auto & detadx_map = fe.get_fe_map().get_detadx();
          const auto & detady_map = fe.get_fe_map().get_detady();
          const auto & detadz_map = fe.get_fe_map().get_detadz();

          const auto & dzetadx_map = fe.get_fe_map().get_dzetadx();
          const auto & dzetady_map = fe.get_fe_map().get_dzetady();
          const auto & dzetadz_map = fe.get_fe_map().get_dzetadz();

          // phi = (dx/dxi)^-T * \hat{phi}
          // In 3D:
          // dx/dxi^-1 = [  dxi/dx    dxi/dy    dxi/dz
          // deta/dx   deta/dy   deta/dz
          // dzeta/dx  dzeta/dy  dzeta/dz]
          //
          // so: dxi/dx^-T * \hat{phi} = [ dxi/dx  deta/dx  dzeta/dx   [ \hat{phi}_xi
          // dxi/dy  deta/dy  dzeta/dy     \hat{phi}_eta
          // dxi/dz  deta/dz  dzeta/dz ]   \hat{phi}_zeta ]
          //
          // or in indicial notation:  phi_j = xi_{i,j}*\hat{phi}_i
          for (auto i : index_range(phi))
            for (auto p : index_range(phi[i]))
            {
              // Need to temporarily cache reference shape functions
              // TODO: PB: Might be worth trying to build phi_ref separately to see
              //           if we can get vectorization
              OutputShape phi_ref;
              FEInterface::shape(3, fe.get_fe_type(), elem, i, qp[p], phi_ref);

              phi[i][p].slice(0) = dxidx_map[p]*phi_ref.slice(0) + detadx_map[p]*phi_ref.slice(1)
                + dzetadx_map[p]*phi_ref.slice(2);

              phi[i][p].slice(1) = dxidy_map[p]*phi_ref.slice(0) + detady_map[p]*phi_ref.slice(1)
                + dzetady_map[p]*phi_ref.slice(2);

              phi[i][p].slice(2) = dxidz_map[p]*phi_ref.slice(0) + detadz_map[p]*phi_ref.slice(1)
                + dzetadz_map[p]*phi_ref.slice(2);
            }
          break;
        }

        default:
          libmesh_error_msg("Invalid dim = " << dim);
      } // switch(dim)
    }

  static void map_curl(const unsigned int dim,
                       const Elem * const,
                       const std::vector<Point> &,
                       const FEGenericBase<RealVectorValue,RealType> & fe,
                       std::vector<std::vector<OutputShape>> & curl_phi)
    {
      switch (dim)
      {
        case 0:
        case 1:
          libmesh_error_msg("These element transformations only make sense in 2D and 3D.");

        case 2:
        {
          const auto & dphi_dxi = fe.get_dphidxi();
          const auto & dphi_deta = fe.get_dphideta();

          const auto & J = fe.get_fe_map().get_jacobian();

          // FIXME: I don't think this is valid for 2D elements in 3D space
          // In 2D: curl(phi) = J^{-1} * curl(\hat{phi})
          for (auto i : index_range(curl_phi))
            for (auto p : index_range(curl_phi[i]))
            {
              curl_phi[i][p].slice(0) = curl_phi[i][p].slice(1) = 0.0;
              curl_phi[i][p].slice(2) = ( dphi_dxi[i][p].slice(1) - dphi_deta[i][p].slice(0) )/J[p];
            }

          break;
        }

        case 3:
        {
          const auto & dphi_dxi = fe.get_dphidxi();
          const auto & dphi_deta = fe.get_dphideta();
          const auto & dphi_dzeta = fe.get_dphidzeta();

          const auto & dxyz_dxi   = fe.get_fe_map().get_dxyzdxi();
          const auto & dxyz_deta  = fe.get_fe_map().get_dxyzdeta();
          const auto & dxyz_dzeta = fe.get_fe_map().get_dxyzdzeta();

          const auto & J = fe.get_fe_map().get_jacobian();

          for (auto i : index_range(curl_phi))
            for (auto p : index_range(curl_phi[i]))
            {
              auto dx_dxi   = dxyz_dxi[p](0);
              auto dx_deta  = dxyz_deta[p](0);
              auto dx_dzeta = dxyz_dzeta[p](0);

              auto dy_dxi   = dxyz_dxi[p](1);
              auto dy_deta  = dxyz_deta[p](1);
              auto dy_dzeta = dxyz_dzeta[p](1);

              auto dz_dxi   = dxyz_dxi[p](2);
              auto dz_deta  = dxyz_deta[p](2);
              auto dz_dzeta = dxyz_dzeta[p](2);

              const auto inv_jac = 1.0/J[p];

              /* In 3D: curl(phi) = J^{-1} dx/dxi * curl(\hat{phi})

                 dx/dxi = [  dx/dxi  dx/deta  dx/dzeta
                 dy/dxi  dy/deta  dy/dzeta
                 dz/dxi  dz/deta  dz/dzeta ]

                 curl(u) = [ du_z/deta  - du_y/dzeta
                 du_x/dzeta - du_z/dxi
                 du_y/dxi   - du_x/deta ]
              */
              curl_phi[i][p].slice(0) = inv_jac*( dx_dxi*( dphi_deta[i][p].slice(2)  -
                                                           dphi_dzeta[i][p].slice(1)   ) +
                                                  dx_deta*( dphi_dzeta[i][p].slice(0) -
                                                            dphi_dxi[i][p].slice(2)     ) +
                                                  dx_dzeta*( dphi_dxi[i][p].slice(1) -
                                                             dphi_deta[i][p].slice(0)    ) );

              curl_phi[i][p].slice(1) = inv_jac*( dy_dxi*( dphi_deta[i][p].slice(2) -
                                                           dphi_dzeta[i][p].slice(1)  ) +
                                                  dy_deta*( dphi_dzeta[i][p].slice(0)-
                                                            dphi_dxi[i][p].slice(2)    ) +
                                                  dy_dzeta*( dphi_dxi[i][p].slice(1) -
                                                             dphi_deta[i][p].slice(0)   ) );

              curl_phi[i][p].slice(2) = inv_jac*( dz_dxi*( dphi_deta[i][p].slice(2) -
                                                           dphi_dzeta[i][p].slice(1)   ) +
                                                  dz_deta*( dphi_dzeta[i][p].slice(0) -
                                                            dphi_dxi[i][p].slice(2)     ) +
                                                  dz_dzeta*( dphi_dxi[i][p].slice(1) -
                                                             dphi_deta[i][p].slice(0)    ) );
            }

          break;
        }

        default:
          libmesh_error_msg("Invalid dim = " << dim);
      } // switch(dim)
    }
};

template <typename RealType>
struct HCurlShim<Real, RealType>
{
  typedef ElemTempl<RealType> Elem;
  typedef PointTempl<RealType> Point;
  typedef typename FEOutputTwoType<Real,RealType>::OutputShape OutputShape;

  static void init_map_phi(const FEGenericBase<Real,RealType> & )
    {
      libmesh_error_msg("HCurl transformations only make sense for vector-valued elements.");
    }

  static void init_map_dphi(const FEGenericBase<Real,RealType> & )
    {
      libmesh_error_msg("HCurl transformations only make sense for vector-valued elements.");
    }

  static void init_map_d2phi(const FEGenericBase<Real,RealType> & )
    {
      libmesh_error_msg("HCurl transformations only make sense for vector-valued elements.");
    }

  static void map_phi(const unsigned int,
                      const Elem * const,
                      const std::vector<Point> &,
                      const FEGenericBase<Real,RealType> &,
                      std::vector<std::vector<OutputShape>> &)
    {
      libmesh_error_msg("HCurl transformations only make sense for vector-valued elements.");
    }

  static void map_curl(const unsigned int,
                       const Elem * const,
                       const std::vector<Point> &,
                       const FEGenericBase<Real,RealType> &,
                       std::vector<std::vector<OutputShape>> &)
    {
      libmesh_error_msg("HCurl transformations only make sense for vector-valued elements.");
    }
};
}

template <typename MathType, typename RealType>
void HCurlFETransformation<MathType,RealType>::init_map_phi(const FEGenericBase<MathType,RealType> & fe) const
{
  HCurlShim<MathType,RealType>::init_map_phi(fe);
}



template <typename MathType, typename RealType>
void HCurlFETransformation<MathType,RealType>::init_map_dphi(const FEGenericBase<MathType,RealType> & fe) const
{
  HCurlShim<MathType,RealType>::init_map_dphi(fe);
}



template <typename MathType, typename RealType>
void HCurlFETransformation<MathType,RealType>::init_map_d2phi(const FEGenericBase<MathType,RealType> & fe) const
{
  HCurlShim<MathType,RealType>::init_map_d2phi(fe);
}



template <typename MathType, typename RealType>
void HCurlFETransformation<MathType,RealType>::map_phi(const unsigned int dim,
                                                       const Elem * const elem,
                                                       const std::vector<Point> & qp,
                                                       const FEGenericBase<MathType,RealType> & fe,
                                                       std::vector<std::vector<OutputShape>> & phi) const
{
  HCurlShim<MathType,RealType>::map_phi(dim, elem, qp, fe, phi);
}

template <typename MathType, typename RealType>
void HCurlFETransformation<MathType,RealType>::map_curl(const unsigned int dim,
                                                        const Elem * const elem,
                                                        const std::vector<Point> & qp,
                                                        const FEGenericBase<MathType,RealType> & fe,
                                                        std::vector<std::vector<OutputShape>> & curl_phi) const
{
   HCurlShim<MathType,RealType>::map_curl(dim, elem, qp, fe, curl_phi);
}

} // namespace libMesh

#endif // LIBMESH_HCURL_FE_TRANSFORMATION_IMPL_H
