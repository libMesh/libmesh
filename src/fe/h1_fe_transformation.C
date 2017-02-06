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

#include "libmesh/fe_interface.h"
#include "libmesh/h1_fe_transformation.h"
#include "libmesh/tensor_value.h"
#include "libmesh/elem.h"

namespace libMesh
{
template< typename OutputShape >
void H1FETransformation<OutputShape>::init_map_phi(const FEGenericBase<OutputShape> & /* fe */ ) const
{
  // Nothing needed
}



template< typename OutputShape >
void H1FETransformation<OutputShape>::init_map_dphi(const FEGenericBase<OutputShape> & fe) const
{
  fe.get_fe_map().get_dxidx();
}



template< typename OutputShape >
void H1FETransformation<OutputShape>::init_map_d2phi(const FEGenericBase<OutputShape> & fe) const
{
  fe.get_fe_map().get_dxidx();
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  fe.get_fe_map().get_d2xidxyz2();
#endif
}

template< typename OutputShape >
void H1FETransformation<OutputShape>::map_phi( const unsigned int dim,
                                               const Elem * const elem,
                                               const std::vector<Point> & qp,
                                               const FEGenericBase<OutputShape> & fe,
                                               std::vector<std::vector<OutputShape> > & phi ) const
{
  switch(dim)
    {
    case 0:
      {
        for (std::size_t i=0; i<phi.size(); i++)
          {
            libmesh_assert_equal_to ( qp.size(), phi[i].size() );
            for (std::size_t p=0; p<phi[i].size(); p++)
              FEInterface::shape<OutputShape>(0, fe.get_fe_type(), elem, i, qp[p], phi[i][p]);
          }
        break;
      }
    case 1:
      {
        for (std::size_t i=0; i<phi.size(); i++)
          {
            libmesh_assert_equal_to ( qp.size(), phi[i].size() );
            for (std::size_t p=0; p<phi[i].size(); p++)
              FEInterface::shape<OutputShape>(1, fe.get_fe_type(), elem, i, qp[p], phi[i][p]);
          }
        break;
      }
    case 2:
      {
        for (std::size_t i=0; i<phi.size(); i++)
          {
            libmesh_assert_equal_to ( qp.size(), phi[i].size() );
            for (std::size_t p=0; p<phi[i].size(); p++)
              FEInterface::shape<OutputShape>(2, fe.get_fe_type(), elem, i, qp[p], phi[i][p]);
          }
        break;
      }
    case 3:
      {
        for (std::size_t i=0; i<phi.size(); i++)
          {
            libmesh_assert_equal_to ( qp.size(), phi[i].size() );
            for (std::size_t p=0; p<phi[i].size(); p++)
              FEInterface::shape<OutputShape>(3, fe.get_fe_type(), elem, i, qp[p], phi[i][p]);
          }
        break;
      }

    default:
      libmesh_error_msg("Invalid dim = " << dim);
    }
}


template< typename OutputShape >
void H1FETransformation<OutputShape>::map_dphi( const unsigned int dim,
                                                const Elem * const,
                                                const std::vector<Point> &,
                                                const FEGenericBase<OutputShape> & fe,
                                                std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputGradient> > & dphi,
                                                std::vector<std::vector<OutputShape> > & dphidx,
                                                std::vector<std::vector<OutputShape> > & dphidy,
                                                std::vector<std::vector<OutputShape> > & dphidz ) const
{
  switch(dim)
    {
    case 0: // No derivatives in 0D
      {
        for (std::size_t i=0; i<dphi.size(); i++)
          for (std::size_t p=0; p<dphi[i].size(); p++)
            dphi[i][p] = 0.;
        break;
      }

    case 1:
      {
        const std::vector<std::vector<OutputShape> > & dphidxi = fe.get_dphidxi();

        const std::vector<Real> & dxidx_map = fe.get_fe_map().get_dxidx();
#if LIBMESH_DIM>1
        const std::vector<Real> & dxidy_map = fe.get_fe_map().get_dxidy();
#endif
#if LIBMESH_DIM>2
        const std::vector<Real> & dxidz_map = fe.get_fe_map().get_dxidz();
#endif

        for (std::size_t i=0; i<dphi.size(); i++)
          for (std::size_t p=0; p<dphi[i].size(); p++)
            {
              // dphi/dx    = (dphi/dxi)*(dxi/dx)
              dphi[i][p].slice(0) = dphidx[i][p] = dphidxi[i][p]*dxidx_map[p];

#if LIBMESH_DIM>1
              dphi[i][p].slice(1)  = dphidy[i][p] = dphidxi[i][p]*dxidy_map[p];
#endif
#if LIBMESH_DIM>2
              dphi[i][p].slice(2) = dphidz[i][p] = dphidxi[i][p]*dxidz_map[p];
#endif
            }

        break;
      }

    case 2:
      {
        const std::vector<std::vector<OutputShape> > & dphidxi = fe.get_dphidxi();
        const std::vector<std::vector<OutputShape> > & dphideta = fe.get_dphideta();

        const std::vector<Real> & dxidx_map = fe.get_fe_map().get_dxidx();
        const std::vector<Real> & dxidy_map = fe.get_fe_map().get_dxidy();
#if LIBMESH_DIM > 2
        const std::vector<Real> & dxidz_map = fe.get_fe_map().get_dxidz();
#endif

        const std::vector<Real> & detadx_map = fe.get_fe_map().get_detadx();
        const std::vector<Real> & detady_map = fe.get_fe_map().get_detady();
#if LIBMESH_DIM > 2
        const std::vector<Real> & detadz_map = fe.get_fe_map().get_detadz();
#endif

        for (std::size_t i=0; i<dphi.size(); i++)
          for (std::size_t p=0; p<dphi[i].size(); p++)
            {
              // dphi/dx    = (dphi/dxi)*(dxi/dx) + (dphi/deta)*(deta/dx)
              dphi[i][p].slice(0) = dphidx[i][p] = (dphidxi[i][p]*dxidx_map[p] +
                                                    dphideta[i][p]*detadx_map[p]);

              // dphi/dy    = (dphi/dxi)*(dxi/dy) + (dphi/deta)*(deta/dy)
              dphi[i][p].slice(1) = dphidy[i][p] = (dphidxi[i][p]*dxidy_map[p] +
                                                    dphideta[i][p]*detady_map[p]);

#if LIBMESH_DIM > 2
              // dphi/dz    = (dphi/dxi)*(dxi/dz) + (dphi/deta)*(deta/dz)
              dphi[i][p].slice(2) = dphidz[i][p] = (dphidxi[i][p]*dxidz_map[p] +
                                                    dphideta[i][p]*detadz_map[p]);
#endif
            }

        break;
      }

    case 3:
      {
        const std::vector<std::vector<OutputShape> > & dphidxi = fe.get_dphidxi();
        const std::vector<std::vector<OutputShape> > & dphideta = fe.get_dphideta();
        const std::vector<std::vector<OutputShape> > & dphidzeta = fe.get_dphidzeta();

        const std::vector<Real> & dxidx_map = fe.get_fe_map().get_dxidx();
        const std::vector<Real> & dxidy_map = fe.get_fe_map().get_dxidy();
        const std::vector<Real> & dxidz_map = fe.get_fe_map().get_dxidz();

        const std::vector<Real> & detadx_map = fe.get_fe_map().get_detadx();
        const std::vector<Real> & detady_map = fe.get_fe_map().get_detady();
        const std::vector<Real> & detadz_map = fe.get_fe_map().get_detadz();

        const std::vector<Real> & dzetadx_map = fe.get_fe_map().get_dzetadx();
        const std::vector<Real> & dzetady_map = fe.get_fe_map().get_dzetady();
        const std::vector<Real> & dzetadz_map = fe.get_fe_map().get_dzetadz();

        for (std::size_t i=0; i<dphi.size(); i++)
          for (std::size_t p=0; p<dphi[i].size(); p++)
            {
              // dphi/dx    = (dphi/dxi)*(dxi/dx) + (dphi/deta)*(deta/dx) + (dphi/dzeta)*(dzeta/dx);
              dphi[i][p].slice(0) = dphidx[i][p] = (dphidxi[i][p]*dxidx_map[p] +
                                                    dphideta[i][p]*detadx_map[p] +
                                                    dphidzeta[i][p]*dzetadx_map[p]);

              // dphi/dy    = (dphi/dxi)*(dxi/dy) + (dphi/deta)*(deta/dy) + (dphi/dzeta)*(dzeta/dy);
              dphi[i][p].slice(1) = dphidy[i][p] = (dphidxi[i][p]*dxidy_map[p] +
                                                    dphideta[i][p]*detady_map[p] +
                                                    dphidzeta[i][p]*dzetady_map[p]);

              // dphi/dz    = (dphi/dxi)*(dxi/dz) + (dphi/deta)*(deta/dz) + (dphi/dzeta)*(dzeta/dz);
              dphi[i][p].slice(2) = dphidz[i][p] = (dphidxi[i][p]*dxidz_map[p] +
                                                    dphideta[i][p]*detadz_map[p] +
                                                    dphidzeta[i][p]*dzetadz_map[p]);
            }
        break;
      }

    default:
      libmesh_error_msg("Invalid dim = " << dim);
    } // switch(dim)
}

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
template< typename OutputShape >
void H1FETransformation<OutputShape>::map_d2phi(const unsigned int dim,
                                                const std::vector<Point> &,
                                                const FEGenericBase<OutputShape> & fe,
                                                std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputTensor> > & d2phi,
                                                std::vector<std::vector<OutputShape> > & d2phidx2,
                                                std::vector<std::vector<OutputShape> > & d2phidxdy,
                                                std::vector<std::vector<OutputShape> > & d2phidxdz,
                                                std::vector<std::vector<OutputShape> > & d2phidy2,
                                                std::vector<std::vector<OutputShape> > & d2phidydz,
                                                std::vector<std::vector<OutputShape> > & d2phidz2) const
{
  switch (dim)
    {
    case 0: // No derivatives in 0D
      {
        for (std::size_t i=0; i<d2phi.size(); i++)
          for (std::size_t p=0; p<d2phi[i].size(); p++)
            d2phi[i][p] = 0.;

        break;
      }

    case 1:
      {
        const std::vector<std::vector<OutputShape> > & d2phidxi2 = fe.get_d2phidxi2();

        const std::vector<Real> & dxidx_map = fe.get_fe_map().get_dxidx();
#if LIBMESH_DIM>1
        const std::vector<Real> & dxidy_map = fe.get_fe_map().get_dxidy();
#endif
#if LIBMESH_DIM>2
        const std::vector<Real> & dxidz_map = fe.get_fe_map().get_dxidz();
#endif

        // Shape function derivatives in reference space
        const std::vector<std::vector<OutputShape> > & dphidxi = fe.get_dphidxi();

        // Inverse map second derivatives
        const std::vector<std::vector<Real> > & d2xidxyz2 = fe.get_fe_map().get_d2xidxyz2();

        for (std::size_t i=0; i<d2phi.size(); i++)
          for (std::size_t p=0; p<d2phi[i].size(); p++)
            {
              // phi_{x x}
              d2phi[i][p].slice(0).slice(0) = d2phidx2[i][p] =
                d2phidxi2[i][p]*dxidx_map[p]*dxidx_map[p] + // (xi_x)^2 * phi_{xi xi}
                d2xidxyz2[p][0]*dphidxi[i][p];              // xi_{x x} * phi_{xi}

#if LIBMESH_DIM>1
              // phi_{x y}
              d2phi[i][p].slice(0).slice(1) = d2phi[i][p].slice(1).slice(0) = d2phidxdy[i][p] =
                d2phidxi2[i][p]*dxidx_map[p]*dxidy_map[p] + // xi_x * xi_y * phi_{xi xi}
                d2xidxyz2[p][1]*dphidxi[i][p];              // xi_{x y} * phi_{xi}
#endif

#if LIBMESH_DIM>2
              // phi_{x z}
              d2phi[i][p].slice(0).slice(2) = d2phi[i][p].slice(2).slice(0) = d2phidxdz[i][p] =
                d2phidxi2[i][p]*dxidx_map[p]*dxidz_map[p] + // xi_x * xi_z * phi_{xi xi}
                d2xidxyz2[p][2]*dphidxi[i][p];              // xi_{x z} * phi_{xi}
#endif


#if LIBMESH_DIM>1
              // phi_{y y}
              d2phi[i][p].slice(1).slice(1) = d2phidy2[i][p] =
                d2phidxi2[i][p]*dxidy_map[p]*dxidy_map[p] + // (xi_y)^2 * phi_{xi xi}
                d2xidxyz2[p][3]*dphidxi[i][p];              // xi_{y y} * phi_{xi}
#endif

#if LIBMESH_DIM>2
              // phi_{y z}
              d2phi[i][p].slice(1).slice(2) = d2phi[i][p].slice(2).slice(1) = d2phidydz[i][p] =
                d2phidxi2[i][p]*dxidy_map[p]*dxidz_map[p] + // xi_y * xi_z * phi_{xi xi}
                d2xidxyz2[p][4]*dphidxi[i][p];              // xi_{y z} * phi_{xi}

              // phi_{z z}
              d2phi[i][p].slice(2).slice(2) = d2phidz2[i][p] =
                d2phidxi2[i][p]*dxidz_map[p]*dxidz_map[p] + // (xi_z)^2 * phi_{xi xi}
                d2xidxyz2[p][5]*dphidxi[i][p];              // xi_{z z} * phi_{xi}
#endif
            }
        break;
      }

    case 2:
      {
        const std::vector<std::vector<OutputShape> > & d2phidxi2 = fe.get_d2phidxi2();
        const std::vector<std::vector<OutputShape> > & d2phidxideta = fe.get_d2phidxideta();
        const std::vector<std::vector<OutputShape> > & d2phideta2 = fe.get_d2phideta2();

        const std::vector<Real> & dxidx_map = fe.get_fe_map().get_dxidx();
        const std::vector<Real> & dxidy_map = fe.get_fe_map().get_dxidy();
#if LIBMESH_DIM > 2
        const std::vector<Real> & dxidz_map = fe.get_fe_map().get_dxidz();
#endif

        const std::vector<Real> & detadx_map = fe.get_fe_map().get_detadx();
        const std::vector<Real> & detady_map = fe.get_fe_map().get_detady();
#if LIBMESH_DIM > 2
        const std::vector<Real> & detadz_map = fe.get_fe_map().get_detadz();
#endif

        // Shape function derivatives in reference space
        const std::vector<std::vector<OutputShape> > & dphidxi = fe.get_dphidxi();
        const std::vector<std::vector<OutputShape> > & dphideta = fe.get_dphideta();

        // Inverse map second derivatives
        const std::vector<std::vector<Real> > & d2xidxyz2 = fe.get_fe_map().get_d2xidxyz2();
        const std::vector<std::vector<Real> > & d2etadxyz2 = fe.get_fe_map().get_d2etadxyz2();

        for (std::size_t i=0; i<d2phi.size(); i++)
          for (std::size_t p=0; p<d2phi[i].size(); p++)
            {
              // phi_{x x}
              d2phi[i][p].slice(0).slice(0) = d2phidx2[i][p] =
                d2phidxi2[i][p]*dxidx_map[p]*dxidx_map[p] +       // (xi_x)^2 * phi_{xi xi}
                d2phideta2[i][p]*detadx_map[p]*detadx_map[p] +    // (eta_x)^2 * phi_{eta eta}
                2*d2phidxideta[i][p]*dxidx_map[p]*detadx_map[p] + // 2 * xi_x * eta_x * phi_{xi eta}
                d2xidxyz2[p][0]*dphidxi[i][p] +                   // xi_{x x} * phi_{xi}
                d2etadxyz2[p][0]*dphideta[i][p];                  // eta_{x x} * phi_{eta}

              // phi_{x y}
              d2phi[i][p].slice(0).slice(1) = d2phi[i][p].slice(1).slice(0) = d2phidxdy[i][p] =
                d2phidxi2[i][p]*dxidx_map[p]*dxidy_map[p] +                                    // xi_x * xi_y * phi_{xi xi}
                d2phideta2[i][p]*detadx_map[p]*detady_map[p] +                                 // eta_x * eta_y * phi_{eta eta}
                d2phidxideta[i][p]*(dxidx_map[p]*detady_map[p] + detadx_map[p]*dxidy_map[p]) + // (xi_x*eta_y + eta_x*xi_y) * phi_{xi eta}
                d2xidxyz2[p][1]*dphidxi[i][p] +                                                // xi_{x y} * phi_{xi}
                d2etadxyz2[p][1]*dphideta[i][p];                                               // eta_{x y} * phi_{eta}

#if LIBMESH_DIM > 2
              // phi_{x z}
              d2phi[i][p].slice(0).slice(2) = d2phi[i][p].slice(2).slice(0) = d2phidxdz[i][p] =
                d2phidxi2[i][p]*dxidx_map[p]*dxidz_map[p] +                                    // xi_x * xi_z * phi_{xi xi}
                d2phideta2[i][p]*detadx_map[p]*detadz_map[p] +                                 // eta_x * eta_z * phi_{eta eta}
                d2phidxideta[i][p]*(dxidx_map[p]*detadz_map[p] + detadx_map[p]*dxidz_map[p]) + // (xi_x*eta_z + eta_x*xi_z) * phi_{xi eta}
                d2xidxyz2[p][2]*dphidxi[i][p] +                                                // xi_{x z} * phi_{xi}
                d2etadxyz2[p][2]*dphideta[i][p];                                               // eta_{x z} * phi_{eta}
#endif

              // phi_{y y}
              d2phi[i][p].slice(1).slice(1) = d2phidy2[i][p] =
                d2phidxi2[i][p]*dxidy_map[p]*dxidy_map[p] +       // (xi_y)^2 * phi_{xi xi}
                d2phideta2[i][p]*detady_map[p]*detady_map[p] +    // (eta_y)^2 * phi_{eta eta}
                2*d2phidxideta[i][p]*dxidy_map[p]*detady_map[p] + // 2 * xi_y * eta_y * phi_{xi eta}
                d2xidxyz2[p][3]*dphidxi[i][p] +                   // xi_{y y} * phi_{xi}
                d2etadxyz2[p][3]*dphideta[i][p];                  // eta_{y y} * phi_{eta}

#if LIBMESH_DIM > 2
              // phi_{y z}
              d2phi[i][p].slice(1).slice(2) = d2phi[i][p].slice(2).slice(1) = d2phidydz[i][p] =
                d2phidxi2[i][p]*dxidy_map[p]*dxidz_map[p] +                                    // xi_y * xi_z * phi_{xi xi}
                d2phideta2[i][p]*detady_map[p]*detadz_map[p] +                                 // eta_y * eta_z * phi_{eta eta}
                d2phidxideta[i][p]*(dxidy_map[p]*detadz_map[p] + detady_map[p]*dxidz_map[p]) + // (xi_y*eta_z + eta_y*xi_z) * phi_{xi eta}
                d2xidxyz2[p][4]*dphidxi[i][p] +                                                // xi_{y z} * phi_{xi}
                d2etadxyz2[p][4]*dphideta[i][p];                                               // eta_{y z} * phi_{eta}

              // phi_{z z}
              d2phi[i][p].slice(2).slice(2) = d2phidz2[i][p] =
                d2phidxi2[i][p]*dxidz_map[p]*dxidz_map[p] +       // (xi_z)^2 * phi_{xi xi}
                d2phideta2[i][p]*detadz_map[p]*detadz_map[p] +    // (eta_z)^2 * phi_{eta eta}
                2*d2phidxideta[i][p]*dxidz_map[p]*detadz_map[p] + // 2 * xi_z * eta_z * phi_{xi eta}
                d2xidxyz2[p][5]*dphidxi[i][p] +                   // xi_{z z} * phi_{xi}
                d2etadxyz2[p][5]*dphideta[i][p];                  // eta_{z z} * phi_{eta}
#endif
            }

        break;
      }

    case 3:
      {
        const std::vector<std::vector<OutputShape> > & d2phidxi2 = fe.get_d2phidxi2();
        const std::vector<std::vector<OutputShape> > & d2phidxideta = fe.get_d2phidxideta();
        const std::vector<std::vector<OutputShape> > & d2phideta2 = fe.get_d2phideta2();
        const std::vector<std::vector<OutputShape> > & d2phidxidzeta = fe.get_d2phidxidzeta();
        const std::vector<std::vector<OutputShape> > & d2phidetadzeta = fe.get_d2phidetadzeta();
        const std::vector<std::vector<OutputShape> > & d2phidzeta2 = fe.get_d2phidzeta2();

        const std::vector<Real> & dxidx_map = fe.get_fe_map().get_dxidx();
        const std::vector<Real> & dxidy_map = fe.get_fe_map().get_dxidy();
        const std::vector<Real> & dxidz_map = fe.get_fe_map().get_dxidz();

        const std::vector<Real> & detadx_map = fe.get_fe_map().get_detadx();
        const std::vector<Real> & detady_map = fe.get_fe_map().get_detady();
        const std::vector<Real> & detadz_map = fe.get_fe_map().get_detadz();

        const std::vector<Real> & dzetadx_map = fe.get_fe_map().get_dzetadx();
        const std::vector<Real> & dzetady_map = fe.get_fe_map().get_dzetady();
        const std::vector<Real> & dzetadz_map = fe.get_fe_map().get_dzetadz();

        // Shape function derivatives in reference space
        const std::vector<std::vector<OutputShape> > & dphidxi = fe.get_dphidxi();
        const std::vector<std::vector<OutputShape> > & dphideta = fe.get_dphideta();
        const std::vector<std::vector<OutputShape> > & dphidzeta = fe.get_dphidzeta();

        // Inverse map second derivatives
        const std::vector<std::vector<Real> > & d2xidxyz2 = fe.get_fe_map().get_d2xidxyz2();
        const std::vector<std::vector<Real> > & d2etadxyz2 = fe.get_fe_map().get_d2etadxyz2();
        const std::vector<std::vector<Real> > & d2zetadxyz2 = fe.get_fe_map().get_d2zetadxyz2();

        for (std::size_t i=0; i<d2phi.size(); i++)
          for (std::size_t p=0; p<d2phi[i].size(); p++)
            {
              // phi_{x x}
              d2phi[i][p].slice(0).slice(0) = d2phidx2[i][p] =
                d2phidxi2[i][p]*dxidx_map[p]*dxidx_map[p] +           // (xi_x)^2 * phi_{xi xi}
                d2phideta2[i][p]*detadx_map[p]*detadx_map[p] +        // (eta_x)^2 * phi_{eta eta}
                d2phidzeta2[i][p]*dzetadx_map[p]*dzetadx_map[p] +     // (zeta_x)^2 * phi_{zeta zeta}
                2*d2phidxideta[i][p]*dxidx_map[p]*detadx_map[p] +     // 2 * xi_x * eta_x * phi_{xi eta}
                2*d2phidxidzeta[i][p]*dxidx_map[p]*dzetadx_map[p] +   // 2 * xi_x * zeta_x * phi_{xi zeta}
                2*d2phidetadzeta[i][p]*detadx_map[p]*dzetadx_map[p] + // 2 * eta_x * zeta_x * phi_{eta zeta}
                d2xidxyz2[p][0]*dphidxi[i][p] +                       // xi_{x x} * phi_{xi}
                d2etadxyz2[p][0]*dphideta[i][p] +                     // eta_{x x} * phi_{eta}
                d2zetadxyz2[p][0]*dphidzeta[i][p];                    // zeta_{x x} * phi_{zeta}

              // phi_{x y}
              d2phi[i][p].slice(0).slice(1) = d2phi[i][p].slice(1).slice(0) = d2phidxdy[i][p] =
                d2phidxi2[i][p]*dxidx_map[p]*dxidy_map[p] +                                          // xi_x * xi_y * phi_{xi xi}
                d2phideta2[i][p]*detadx_map[p]*detady_map[p] +                                       // eta_x * eta_y * phi_{eta eta}
                d2phidzeta2[i][p]*dzetadx_map[p]*dzetady_map[p] +                                    // zeta_x * zeta_y * phi_{zeta zeta}
                d2phidxideta[i][p]*(dxidx_map[p]*detady_map[p] + detadx_map[p]*dxidy_map[p]) +       // (xi_x*eta_y + eta_x*xi_y) * phi_{xi eta}
                d2phidxidzeta[i][p]*(dxidx_map[p]*dzetady_map[p] + dzetadx_map[p]*dxidy_map[p]) +    // (zeta_x*xi_y + xi_x*zeta_y) * phi_{xi zeta}
                d2phidetadzeta[i][p]*(detadx_map[p]*dzetady_map[p] + dzetadx_map[p]*detady_map[p]) + // (zeta_x*eta_y + eta_x*zeta_y) * phi_{eta zeta}
                d2xidxyz2[p][1]*dphidxi[i][p] +                                                      // xi_{x y} * phi_{xi}
                d2etadxyz2[p][1]*dphideta[i][p] +                                                    // eta_{x y} * phi_{eta}
                d2zetadxyz2[p][1]*dphidzeta[i][p];                                                   // zeta_{x y} * phi_{zeta}

              // phi_{x z}
              d2phi[i][p].slice(0).slice(2) = d2phi[i][p].slice(2).slice(0) = d2phidy2[i][p] =
                d2phidxi2[i][p]*dxidx_map[p]*dxidz_map[p] +                                          // xi_x * xi_z * phi_{xi xi}
                d2phideta2[i][p]*detadx_map[p]*detadz_map[p] +                                       // eta_x * eta_z * phi_{eta eta}
                d2phidzeta2[i][p]*dzetadx_map[p]*dzetadz_map[p] +                                    // zeta_x * zeta_z * phi_{zeta zeta}
                d2phidxideta[i][p]*(dxidx_map[p]*detadz_map[p] + detadx_map[p]*dxidz_map[p]) +       // (xi_x*eta_z + eta_x*xi_z) * phi_{xi eta}
                d2phidxidzeta[i][p]*(dxidx_map[p]*dzetadz_map[p] + dzetadx_map[p]*dxidz_map[p]) +    // (zeta_x*xi_z + xi_x*zeta_z) * phi_{xi zeta}
                d2phidetadzeta[i][p]*(detadx_map[p]*dzetadz_map[p] + dzetadx_map[p]*detadz_map[p]) + // (zeta_x*eta_z + eta_x*zeta_z) * phi_{eta zeta}
                d2xidxyz2[p][2]*dphidxi[i][p] +                                                      // xi_{x z} * phi_{xi}
                d2etadxyz2[p][2]*dphideta[i][p] +                                                    // eta_{x z} * phi_{eta}
                d2zetadxyz2[p][2]*dphidzeta[i][p];                                                   // zeta_{x z} * phi_{zeta}

              // phi_{y y}
              d2phi[i][p].slice(1).slice(1) = d2phidxdz[i][p] =
                d2phidxi2[i][p]*dxidy_map[p]*dxidy_map[p] +           // (xi_y)^2 * phi_{xi xi}
                d2phideta2[i][p]*detady_map[p]*detady_map[p] +        // (eta_y)^2 * phi_{eta eta}
                d2phidzeta2[i][p]*dzetady_map[p]*dzetady_map[p] +     // (zeta_y)^2 * phi_{zeta zeta}
                2*d2phidxideta[i][p]*dxidy_map[p]*detady_map[p] +     // 2 * xi_y * eta_y * phi_{xi eta}
                2*d2phidxidzeta[i][p]*dxidy_map[p]*dzetady_map[p] +   // 2 * xi_y * zeta_y * phi_{xi zeta}
                2*d2phidetadzeta[i][p]*detady_map[p]*dzetady_map[p] + // 2 * eta_y * zeta_y * phi_{eta zeta}
                d2xidxyz2[p][3]*dphidxi[i][p] +                       // xi_{y y} * phi_{xi}
                d2etadxyz2[p][3]*dphideta[i][p] +                     // eta_{y y} * phi_{eta}
                d2zetadxyz2[p][3]*dphidzeta[i][p];                    // zeta_{y y} * phi_{zeta}

              // phi_{y z}
              d2phi[i][p].slice(1).slice(2) = d2phi[i][p].slice(2).slice(1) = d2phidydz[i][p] =
                d2phidxi2[i][p]*dxidy_map[p]*dxidz_map[p] +                                          // xi_y * xi_z * phi_{xi xi}
                d2phideta2[i][p]*detady_map[p]*detadz_map[p] +                                       // eta_y * eta_z * phi_{eta eta}
                d2phidzeta2[i][p]*dzetady_map[p]*dzetadz_map[p] +                                    // zeta_y * zeta_z * phi_{zeta zeta}
                d2phidxideta[i][p]*(dxidy_map[p]*detadz_map[p] + detady_map[p]*dxidz_map[p]) +       // (xi_y*eta_z + eta_y*xi_z) * phi_{xi eta}
                d2phidxidzeta[i][p]*(dxidy_map[p]*dzetadz_map[p] + dzetady_map[p]*dxidz_map[p]) +    // (zeta_y*xi_z + xi_y*zeta_z) * phi_{xi zeta}
                d2phidetadzeta[i][p]*(detady_map[p]*dzetadz_map[p] + dzetady_map[p]*detadz_map[p]) + // (zeta_y*eta_z + eta_y*zeta_z) * phi_{eta zeta}
                d2xidxyz2[p][4]*dphidxi[i][p] +                                                      // xi_{y z} * phi_{xi}
                d2etadxyz2[p][4]*dphideta[i][p] +                                                    // eta_{y z} * phi_{eta}
                d2zetadxyz2[p][4]*dphidzeta[i][p];                                                   // zeta_{y z} * phi_{zeta}

              // phi_{z z}
              d2phi[i][p].slice(2).slice(2) = d2phidz2[i][p] =
                d2phidxi2[i][p]*dxidz_map[p]*dxidz_map[p] +           // (xi_z)^2 * phi_{xi xi}
                d2phideta2[i][p]*detadz_map[p]*detadz_map[p] +        // (eta_z)^2 * phi_{eta eta}
                d2phidzeta2[i][p]*dzetadz_map[p]*dzetadz_map[p] +     // (zeta_z)^2 * phi_{zeta zeta}
                2*d2phidxideta[i][p]*dxidz_map[p]*detadz_map[p] +     // 2 * xi_z * eta_z * phi_{xi eta}
                2*d2phidxidzeta[i][p]*dxidz_map[p]*dzetadz_map[p] +   // 2 * xi_z * zeta_z * phi_{xi zeta}
                2*d2phidetadzeta[i][p]*detadz_map[p]*dzetadz_map[p] + // 2 * eta_z * zeta_z * phi_{eta zeta}
                d2xidxyz2[p][5]*dphidxi[i][p] +                       // xi_{z z} * phi_{xi}
                d2etadxyz2[p][5]*dphideta[i][p] +                     // eta_{z z} * phi_{eta}
                d2zetadxyz2[p][5]*dphidzeta[i][p];                    // zeta_{z z} * phi_{zeta}
            }

        break;
      }

    default:
      libmesh_error_msg("Invalid dim = " << dim);
    } // switch(dim)
}
#endif //LIBMESH_ENABLE_SECOND_DERIVATIVES

template< >
void H1FETransformation<Real>::map_curl(const unsigned int,
                                        const Elem * const,
                                        const std::vector<Point> &,
                                        const FEGenericBase<Real> &,
                                        std::vector<std::vector<Real> > &) const
{
  libmesh_error_msg("Computing the curl of a shape function only \nmakes sense for vector-valued elements.");
}

template< >
void H1FETransformation<RealGradient>::map_curl(const unsigned int dim,
                                                const Elem * const,
                                                const std::vector<Point> &,
                                                const FEGenericBase<RealGradient> & fe,
                                                std::vector<std::vector<RealGradient> > & curl_phi ) const
{
  switch (dim)
    {
    case 0:
    case 1:
      libmesh_error_msg("The curl only make sense in 2D and 3D");

    case 2:
      {
        const std::vector<std::vector<RealGradient> > & dphidxi = fe.get_dphidxi();
        const std::vector<std::vector<RealGradient> > & dphideta = fe.get_dphideta();

        const std::vector<Real> & dxidx_map = fe.get_fe_map().get_dxidx();
        const std::vector<Real> & dxidy_map = fe.get_fe_map().get_dxidy();
#if LIBMESH_DIM > 2
        const std::vector<Real> & dxidz_map = fe.get_fe_map().get_dxidz();
#endif

        const std::vector<Real> & detadx_map = fe.get_fe_map().get_detadx();
        const std::vector<Real> & detady_map = fe.get_fe_map().get_detady();
#if LIBMESH_DIM > 2
        const std::vector<Real> & detadz_map = fe.get_fe_map().get_detadz();
#endif

        /*
          For 2D elements in 3D space:

          curl = ( -dphi_y/dz, dphi_x/dz, dphi_y/dx - dphi_x/dy )
        */
        for (std::size_t i=0; i<curl_phi.size(); i++)
          for (std::size_t p=0; p<curl_phi[i].size(); p++)
            {

              Real dphiy_dx = (dphidxi[i][p].slice(1))*dxidx_map[p]
                + (dphideta[i][p].slice(1))*detadx_map[p];

              Real dphix_dy = (dphidxi[i][p].slice(0))*dxidy_map[p]
                + (dphideta[i][p].slice(0))*detady_map[p];

              curl_phi[i][p].slice(2) = dphiy_dx - dphix_dy;

#if LIBMESH_DIM > 2
              Real dphiy_dz = (dphidxi[i][p].slice(1))*dxidz_map[p]
                + (dphideta[i][p].slice(1))*detadz_map[p];

              Real dphix_dz = (dphidxi[i][p].slice(0))*dxidz_map[p]
                + (dphideta[i][p].slice(0))*detadz_map[p];

              curl_phi[i][p].slice(0) = -dphiy_dz;
              curl_phi[i][p].slice(1) = dphix_dz;
#endif
            }

        break;
      }
    case 3:
      {
        const std::vector<std::vector<RealGradient> > & dphidxi = fe.get_dphidxi();
        const std::vector<std::vector<RealGradient> > & dphideta = fe.get_dphideta();
        const std::vector<std::vector<RealGradient> > & dphidzeta = fe.get_dphidzeta();

        const std::vector<Real> & dxidx_map = fe.get_fe_map().get_dxidx();
        const std::vector<Real> & dxidy_map = fe.get_fe_map().get_dxidy();
        const std::vector<Real> & dxidz_map = fe.get_fe_map().get_dxidz();

        const std::vector<Real> & detadx_map = fe.get_fe_map().get_detadx();
        const std::vector<Real> & detady_map = fe.get_fe_map().get_detady();
        const std::vector<Real> & detadz_map = fe.get_fe_map().get_detadz();

        const std::vector<Real> & dzetadx_map = fe.get_fe_map().get_dzetadx();
        const std::vector<Real> & dzetady_map = fe.get_fe_map().get_dzetady();
        const std::vector<Real> & dzetadz_map = fe.get_fe_map().get_dzetadz();

        /*
          In 3D: curl = ( dphi_z/dy - dphi_y/dz, dphi_x/dz - dphi_z/dx, dphi_y/dx - dphi_x/dy )
        */
        for (std::size_t i=0; i<curl_phi.size(); i++)
          for (std::size_t p=0; p<curl_phi[i].size(); p++)
            {
              Real dphiz_dy = (dphidxi[i][p].slice(2))*dxidy_map[p]
                + (dphideta[i][p].slice(2))*detady_map[p]
                + (dphidzeta[i][p].slice(2))*dzetady_map[p];

              Real dphiy_dz = (dphidxi[i][p].slice(1))*dxidz_map[p]
                + (dphideta[i][p].slice(1))*detadz_map[p]
                + (dphidzeta[i][p].slice(1))*dzetadz_map[p];

              Real dphix_dz = (dphidxi[i][p].slice(0))*dxidz_map[p]
                + (dphideta[i][p].slice(0))*detadz_map[p]
                + (dphidzeta[i][p].slice(0))*dzetadz_map[p];

              Real dphiz_dx = (dphidxi[i][p].slice(2))*dxidx_map[p]
                + (dphideta[i][p].slice(2))*detadx_map[p]
                + (dphidzeta[i][p].slice(2))*dzetadx_map[p];

              Real dphiy_dx = (dphidxi[i][p].slice(1))*dxidx_map[p]
                + (dphideta[i][p].slice(1))*detadx_map[p]
                + (dphidzeta[i][p].slice(1))*dzetadx_map[p];

              Real dphix_dy = (dphidxi[i][p].slice(0))*dxidy_map[p]
                + (dphideta[i][p].slice(0))*detady_map[p]
                + (dphidzeta[i][p].slice(0))*dzetady_map[p];

              curl_phi[i][p].slice(0) = dphiz_dy - dphiy_dz;

              curl_phi[i][p].slice(1) = dphix_dz - dphiz_dx;

              curl_phi[i][p].slice(2) = dphiy_dx - dphix_dy;
            }

        break;
      }
    default:
      libmesh_error_msg("Invalid dim = " << dim);
    }
}


template< >
void H1FETransformation<Real>::map_div(const unsigned int,
                                       const Elem * const,
                                       const std::vector<Point> &,
                                       const FEGenericBase<Real> &,
                                       std::vector<std::vector<FEGenericBase<Real>::OutputDivergence> > & ) const
{
  libmesh_error_msg("Computing the divergence of a shape function only \nmakes sense for vector-valued elements.");
}


template<>
void H1FETransformation<RealGradient>::map_div(const unsigned int dim,
                                               const Elem * const,
                                               const std::vector<Point> &,
                                               const FEGenericBase<RealGradient> & fe,
                                               std::vector<std::vector<FEGenericBase<RealGradient>::OutputDivergence> > & div_phi) const
{
  switch (dim)
    {
    case 0:
    case 1:
      libmesh_error_msg("The divergence only make sense in 2D and 3D");

    case 2:
      {
        const std::vector<std::vector<RealGradient> > & dphidxi = fe.get_dphidxi();
        const std::vector<std::vector<RealGradient> > & dphideta = fe.get_dphideta();

        const std::vector<Real> & dxidx_map = fe.get_fe_map().get_dxidx();
        const std::vector<Real> & dxidy_map = fe.get_fe_map().get_dxidy();

        const std::vector<Real> & detadx_map = fe.get_fe_map().get_detadx();
        const std::vector<Real> & detady_map = fe.get_fe_map().get_detady();

        // In 2D: div = dphi_x/dx + dphi_y/dy
        for (std::size_t i=0; i<div_phi.size(); i++)
          for (std::size_t p=0; p<div_phi[i].size(); p++)
            {

              Real dphix_dx = (dphidxi[i][p].slice(0))*dxidx_map[p]
                + (dphideta[i][p].slice(0))*detadx_map[p];

              Real dphiy_dy = (dphidxi[i][p].slice(1))*dxidy_map[p]
                + (dphideta[i][p].slice(1))*detady_map[p];

              div_phi[i][p] = dphix_dx + dphiy_dy;
            }
        break;
      }
    case 3:
      {
        const std::vector<std::vector<RealGradient> > & dphidxi = fe.get_dphidxi();
        const std::vector<std::vector<RealGradient> > & dphideta = fe.get_dphideta();
        const std::vector<std::vector<RealGradient> > & dphidzeta = fe.get_dphidzeta();

        const std::vector<Real> & dxidx_map = fe.get_fe_map().get_dxidx();
        const std::vector<Real> & dxidy_map = fe.get_fe_map().get_dxidy();
        const std::vector<Real> & dxidz_map = fe.get_fe_map().get_dxidz();

        const std::vector<Real> & detadx_map = fe.get_fe_map().get_detadx();
        const std::vector<Real> & detady_map = fe.get_fe_map().get_detady();
        const std::vector<Real> & detadz_map = fe.get_fe_map().get_detadz();

        const std::vector<Real> & dzetadx_map = fe.get_fe_map().get_dzetadx();
        const std::vector<Real> & dzetady_map = fe.get_fe_map().get_dzetady();
        const std::vector<Real> & dzetadz_map = fe.get_fe_map().get_dzetadz();

        /*
          In 3D: div = dphi_x/dx + dphi_y/dy + dphi_z/dz
        */
        for (std::size_t i=0; i<div_phi.size(); i++)
          for (std::size_t p=0; p<div_phi[i].size(); p++)
            {
              Real dphix_dx = (dphidxi[i][p].slice(0))*dxidx_map[p]
                + (dphideta[i][p].slice(0))*detadx_map[p]
                + (dphidzeta[i][p].slice(0))*dzetadx_map[p];

              Real dphiy_dy = (dphidxi[i][p].slice(1))*dxidy_map[p]
                + (dphideta[i][p].slice(1))*detady_map[p]
                + (dphidzeta[i][p].slice(1))*dzetady_map[p];

              Real dphiz_dz = (dphidxi[i][p].slice(2))*dxidz_map[p]
                + (dphideta[i][p].slice(2))*detadz_map[p]
                + (dphidzeta[i][p].slice(2))*dzetadz_map[p];

              div_phi[i][p] = dphix_dx + dphiy_dy + dphiz_dz;
            }

        break;
      }

    default:                                                  \
      libmesh_error_msg("Invalid dim = " << dim);             \
    } // switch(dim)

  return;
}

// Explicit Instantations
template class H1FETransformation<Real>;
template class H1FETransformation<RealGradient>;

} //namespace libMesh
