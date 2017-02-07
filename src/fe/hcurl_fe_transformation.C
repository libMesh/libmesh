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

#include "libmesh/hcurl_fe_transformation.h"
#include "libmesh/fe_interface.h"

namespace libMesh
{

template< typename OutputShape >
void HCurlFETransformation<OutputShape>::init_map_phi(const FEGenericBase<OutputShape> & fe) const
{
  fe.get_fe_map().get_dxidx();
}



template< typename OutputShape >
void HCurlFETransformation<OutputShape>::init_map_dphi(const FEGenericBase<OutputShape> & fe) const
{
  fe.get_fe_map().get_dxidx();
}



template< typename OutputShape >
void HCurlFETransformation<OutputShape>::init_map_d2phi(const FEGenericBase<OutputShape> & fe) const
{
  fe.get_fe_map().get_dxidx();
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  fe.get_fe_map().get_d2xidxyz2();
#endif
}



template< typename OutputShape >
void HCurlFETransformation<OutputShape>::map_phi(const unsigned int dim,
                                                 const Elem * const elem,
                                                 const std::vector<Point> & qp,
                                                 const FEGenericBase<OutputShape> & fe,
                                                 std::vector<std::vector<OutputShape> > & phi) const
{
  switch (dim)
    {
    case 0:
    case 1:
      libmesh_error_msg("These element transformations only make sense in 2D and 3D.");

    case 2:
      {
        const std::vector<Real> & dxidx_map = fe.get_fe_map().get_dxidx();
        const std::vector<Real> & dxidy_map = fe.get_fe_map().get_dxidy();

        const std::vector<Real> & detadx_map = fe.get_fe_map().get_detadx();
        const std::vector<Real> & detady_map = fe.get_fe_map().get_detady();

        // FIXME: Need to update for 2D elements in 3D space
        /* phi = (dx/dxi)^-T * \hat{phi}
           In 2D:
           (dx/dxi)^{-1} = [  dxi/dx   dxi/dy
           deta/dx  deta/dy ]

           so: dxi/dx^{-T} * \hat{phi} = [ dxi/dx  deta/dx   [ \hat{phi}_xi
           dxi/dy  deta/dy ]   \hat{phi}_eta ]

           or in indicial notation:  phi_j = xi_{i,j}*\hat{phi}_i */

        for (std::size_t i=0; i<phi.size(); i++)
          for (std::size_t p=0; p<phi[i].size(); p++)
            {
              // Need to temporarily cache reference shape functions
              // TODO: PB: Might be worth trying to build phi_ref separately to see
              //           if we can get vectorization
              OutputShape phi_ref;
              FEInterface::shape<OutputShape>(2, fe.get_fe_type(), elem, i, qp[p], phi_ref);

              phi[i][p](0) = dxidx_map[p]*phi_ref.slice(0) + detadx_map[p]*phi_ref.slice(1);

              phi[i][p](1) = dxidy_map[p]*phi_ref.slice(0) + detady_map[p]*phi_ref.slice(1);
            }

        break;
      }

    case 3:
      {
        const std::vector<Real> & dxidx_map = fe.get_fe_map().get_dxidx();
        const std::vector<Real> & dxidy_map = fe.get_fe_map().get_dxidy();
        const std::vector<Real> & dxidz_map = fe.get_fe_map().get_dxidz();

        const std::vector<Real> & detadx_map = fe.get_fe_map().get_detadx();
        const std::vector<Real> & detady_map = fe.get_fe_map().get_detady();
        const std::vector<Real> & detadz_map = fe.get_fe_map().get_detadz();

        const std::vector<Real> & dzetadx_map = fe.get_fe_map().get_dzetadx();
        const std::vector<Real> & dzetady_map = fe.get_fe_map().get_dzetady();
        const std::vector<Real> & dzetadz_map = fe.get_fe_map().get_dzetadz();

        /* phi = (dx/dxi)^-T * \hat{phi}
           In 3D:
           dx/dxi^-1 = [  dxi/dx    dxi/dy    dxi/dz
           deta/dx   deta/dy   deta/dz
           dzeta/dx  dzeta/dy  dzeta/dz]

           so: dxi/dx^-T * \hat{phi} = [ dxi/dx  deta/dx  dzeta/dx   [ \hat{phi}_xi
           dxi/dy  deta/dy  dzeta/dy     \hat{phi}_eta
           dxi/dz  deta/dz  dzeta/dz ]   \hat{phi}_zeta ]

           or in indicial notation:  phi_j = xi_{i,j}*\hat{phi}_i */

        for (std::size_t i=0; i<phi.size(); i++)
          for (std::size_t p=0; p<phi[i].size(); p++)
            {
              // Need to temporarily cache reference shape functions
              // TODO: PB: Might be worth trying to build phi_ref separately to see
              //           if we can get vectorization
              OutputShape phi_ref;
              FEInterface::shape<OutputShape>(3, fe.get_fe_type(), elem, i, qp[p], phi_ref);

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

template< typename OutputShape >
void HCurlFETransformation<OutputShape>::map_curl(const unsigned int dim,
                                                  const Elem * const,
                                                  const std::vector<Point> &,
                                                  const FEGenericBase<OutputShape> & fe,
                                                  std::vector<std::vector<OutputShape> > & curl_phi) const
{
  switch (dim)
    {
    case 0:
    case 1:
      libmesh_error_msg("These element transformations only make sense in 2D and 3D.");

    case 2:
      {
        const std::vector<std::vector<OutputShape> > & dphi_dxi = fe.get_dphidxi();
        const std::vector<std::vector<OutputShape> > & dphi_deta = fe.get_dphideta();

        const std::vector<Real> & J = fe.get_fe_map().get_jacobian();

        // FIXME: I don't think this is valid for 2D elements in 3D space
        /* In 2D: curl(phi) = J^{-1} * curl(\hat{phi}) */
        for (std::size_t i=0; i<curl_phi.size(); i++)
          for (std::size_t p=0; p<curl_phi[i].size(); p++)
            {
              curl_phi[i][p].slice(0) = curl_phi[i][p].slice(1) = 0.0;

              curl_phi[i][p].slice(2) = ( dphi_dxi[i][p].slice(1) - dphi_deta[i][p].slice(0) )/J[p];
            }

        break;
      }

    case 3:
      {
        const std::vector<std::vector<OutputShape> > & dphi_dxi = fe.get_dphidxi();
        const std::vector<std::vector<OutputShape> > & dphi_deta = fe.get_dphideta();
        const std::vector<std::vector<OutputShape> > & dphi_dzeta = fe.get_dphidzeta();

        const std::vector<RealGradient> & dxyz_dxi   = fe.get_fe_map().get_dxyzdxi();
        const std::vector<RealGradient> & dxyz_deta  = fe.get_fe_map().get_dxyzdeta();
        const std::vector<RealGradient> & dxyz_dzeta = fe.get_fe_map().get_dxyzdzeta();

        const std::vector<Real> & J = fe.get_fe_map().get_jacobian();

        for (std::size_t i=0; i<curl_phi.size(); i++)
          for (std::size_t p=0; p<curl_phi[i].size(); p++)
            {
              Real dx_dxi   = dxyz_dxi[p](0);
              Real dx_deta  = dxyz_deta[p](0);
              Real dx_dzeta = dxyz_dzeta[p](0);

              Real dy_dxi   = dxyz_dxi[p](1);
              Real dy_deta  = dxyz_deta[p](1);
              Real dy_dzeta = dxyz_dzeta[p](1);

              Real dz_dxi   = dxyz_dxi[p](2);
              Real dz_deta  = dxyz_deta[p](2);
              Real dz_dzeta = dxyz_dzeta[p](2);

              const Real inv_jac = 1.0/J[p];

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

template class HCurlFETransformation<RealGradient>;

template<>
void HCurlFETransformation<Real>::init_map_phi(const FEGenericBase<Real> & ) const
{
  libmesh_error_msg("HCurl transformations only make sense for vector-valued elements.");
}

template<>
void HCurlFETransformation<Real>::init_map_dphi(const FEGenericBase<Real> & ) const
{
  libmesh_error_msg("HCurl transformations only make sense for vector-valued elements.");
}

template<>
void HCurlFETransformation<Real>::init_map_d2phi(const FEGenericBase<Real> & ) const
{
  libmesh_error_msg("HCurl transformations only make sense for vector-valued elements.");
}

template<>
void HCurlFETransformation<Real>::map_phi(const unsigned int,
                                          const Elem * const,
                                          const std::vector<Point> &,
                                          const FEGenericBase<Real> &,
                                          std::vector<std::vector<Real> > &) const
{
  libmesh_error_msg("HCurl transformations only make sense for vector-valued elements.");
}

template<>
void HCurlFETransformation<Real>::map_curl(const unsigned int,
                                           const Elem * const,
                                           const std::vector<Point> &,
                                           const FEGenericBase<Real> &,
                                           std::vector<std::vector<Real> > &) const
{
  libmesh_error_msg("HCurl transformations only make sense for vector-valued elements.");
}


} // namespace libMesh
