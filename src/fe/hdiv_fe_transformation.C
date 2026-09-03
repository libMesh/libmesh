// The libMesh Finite Element Library.
// Copyright (C) 2002-2026 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#include "libmesh/hdiv_fe_transformation.h"
#include "libmesh/fe_interface.h"
#include "libmesh/int_range.h"

namespace {

using namespace libMesh;

// The contravariant Piola map phi = J^{-1} * (dx/dxi) * phihat, shared by
// map_phi and map_dphi (which also needs the physical shape function value
// for the term coming from the derivative of J^{-1}).
template<typename OutputShape>
OutputShape hdiv_piola_map(const RealGradient & dx_dxi,
                           const RealGradient & dx_deta,
                           Real J,
                           const OutputShape & phi_ref)
{
  return (dx_dxi*phi_ref(0) + dx_deta*phi_ref(1))/J;
}

template<typename OutputShape>
OutputShape hdiv_piola_map(const RealGradient & dx_dxi,
                           const RealGradient & dx_deta,
                           const RealGradient & dx_dzeta,
                           Real J,
                           const OutputShape & phi_ref)
{
  return (dx_dxi*phi_ref(0) + dx_deta*phi_ref(1) + dx_dzeta*phi_ref(2))/J;
}

} // anonymous namespace

namespace libMesh
{

template<typename OutputShape>
void HDivFETransformation<OutputShape>::init_map_phi(const FEGenericBase<OutputShape> & fe) const
{
  fe.get_fe_map().get_dxidx();
}



template<typename OutputShape>
void HDivFETransformation<OutputShape>::init_map_dphi(const FEGenericBase<OutputShape> & fe) const
{
  fe.get_fe_map().get_dxidx();
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  fe.get_fe_map().get_d2xyzdxi2();
#endif
}



template<typename OutputShape>
void HDivFETransformation<OutputShape>::init_map_d2phi(const FEGenericBase<OutputShape> & fe) const
{
  fe.get_fe_map().get_dxidx();
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  fe.get_fe_map().get_d2xidxyz2();
#endif
}



template<typename OutputShape>
void HDivFETransformation<OutputShape>::map_phi(const unsigned int dim,
                                                const Elem * const elem,
                                                const std::vector<Point> & qp,
                                                const FEGenericBase<OutputShape> & fe,
                                                std::vector<std::vector<OutputShape>> & phi,
                                                const bool /*add_p_level*/) const
{
  switch (dim)
    {
    case 0:
    case 1:
      libmesh_error_msg("These element transformations only make sense in 2D and 3D.");

    case 2:
      {
        const std::vector<RealGradient> & dxyz_dxi   = fe.get_fe_map().get_dxyzdxi();
        const std::vector<RealGradient> & dxyz_deta  = fe.get_fe_map().get_dxyzdeta();

        const std::vector<Real> & J = fe.get_fe_map().get_jacobian();

        // phi = J^{-1} * (dx/dxi) * \hat{phi}
        for (auto i : index_range(phi))
          for (auto p : index_range(phi[i]))
            {
              // Need to temporarily cache reference shape functions
              // We are computing mapping basis functions, so we explicitly ignore
              // any non-zero p_level() the Elem might have.
              OutputShape phi_ref;
              FEInterface::shape(fe.get_fe_type(), /*extra_order=*/0, elem, i, qp[p], phi_ref);

              phi[i][p] = hdiv_piola_map(dxyz_dxi[p], dxyz_deta[p], J[p], phi_ref);
            }

        break;
      }
    case 3:
      {
        const std::vector<RealGradient> & dxyz_dxi   = fe.get_fe_map().get_dxyzdxi();
        const std::vector<RealGradient> & dxyz_deta  = fe.get_fe_map().get_dxyzdeta();
        const std::vector<RealGradient> & dxyz_dzeta = fe.get_fe_map().get_dxyzdzeta();

        const std::vector<Real> & J = fe.get_fe_map().get_jacobian();

        // phi = J^{-1} * (dx/dxi) * \hat{phi}
        for (auto i : index_range(phi))
          for (auto p : index_range(phi[i]))
            {
              // Need to temporarily cache reference shape functions
              // We are computing mapping basis functions, so we explicitly ignore
              // any non-zero p_level() the Elem might have.
              OutputShape phi_ref;
              FEInterface::shape(fe.get_fe_type(), /*extra_order=*/0, elem, i, qp[p], phi_ref);

              phi[i][p] = hdiv_piola_map(dxyz_dxi[p], dxyz_deta[p], dxyz_dzeta[p], J[p], phi_ref);
            }

        break;
      }

    default:
      libmesh_error_msg("Invalid dim = " << dim);
    } // switch(dim)
}

template<typename OutputShape>
void HDivFETransformation<OutputShape>::map_dphi(const unsigned int dim,
                                                  const Elem * const elem,
                                                  const std::vector<Point> & qp,
                                                  const FEGenericBase<OutputShape> & fe,
                                                  std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputGradient>> & dphi,
                                                  std::vector<std::vector<OutputShape>> & dphidx,
                                                  std::vector<std::vector<OutputShape>> & dphidy,
                                                  std::vector<std::vector<OutputShape>> & dphidz) const
{
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  switch (dim)
    {
    case 0:
    case 1:
      libmesh_error_msg("These element transformations only make sense in 2D and 3D.");

    case 2:
      {
        const std::vector<RealGradient> & dxyz_dxi   = fe.get_fe_map().get_dxyzdxi();
        const std::vector<RealGradient> & dxyz_deta  = fe.get_fe_map().get_dxyzdeta();

        const std::vector<RealGradient> & d2xyz_dxi2    = fe.get_fe_map().get_d2xyzdxi2();
        const std::vector<RealGradient> & d2xyz_deta2   = fe.get_fe_map().get_d2xyzdeta2();
        const std::vector<RealGradient> & d2xyz_dxideta = fe.get_fe_map().get_d2xyzdxideta();

        const std::vector<Real> & J = fe.get_fe_map().get_jacobian();

        const std::vector<Real> & dxi_dx  = fe.get_fe_map().get_dxidx();
        const std::vector<Real> & dxi_dy  = fe.get_fe_map().get_dxidy();
        const std::vector<Real> & deta_dx = fe.get_fe_map().get_detadx();
        const std::vector<Real> & deta_dy = fe.get_fe_map().get_detady();
#if LIBMESH_DIM > 2
        const std::vector<Real> & dxi_dz  = fe.get_fe_map().get_dxidz();
        const std::vector<Real> & deta_dz = fe.get_fe_map().get_detadz();
#endif

        const std::vector<std::vector<OutputShape>> & dphi_dxi  = fe.get_dphidxi();
        const std::vector<std::vector<OutputShape>> & dphi_deta = fe.get_dphideta();

        for (auto i : index_range(dphi))
          for (auto p : index_range(dphi[i]))
            {
              // Need to temporarily cache reference shape functions
              // We are computing mapping basis functions, so we explicitly ignore
              // any non-zero p_level() the Elem might have.
              OutputShape phi_ref;
              FEInterface::shape(fe.get_fe_type(), /*extra_order=*/0, elem, i, qp[p], phi_ref);

              const OutputShape phi_val =
                hdiv_piola_map(dxyz_dxi[p], dxyz_deta[p], J[p], phi_ref);

              // Jacobi's formula: dJ/dxi_n = J * tr(F^{-1} * dF/dxi_n)
              Real dJ_dxi =
                J[p] * (dxi_dx[p]*d2xyz_dxi2[p](0)    + deta_dx[p]*d2xyz_dxideta[p](0) +
                        dxi_dy[p]*d2xyz_dxi2[p](1)    + deta_dy[p]*d2xyz_dxideta[p](1));
              Real dJ_deta =
                J[p] * (dxi_dx[p]*d2xyz_dxideta[p](0) + deta_dx[p]*d2xyz_deta2[p](0) +
                        dxi_dy[p]*d2xyz_dxideta[p](1) + deta_dy[p]*d2xyz_deta2[p](1));
#if LIBMESH_DIM > 2
              dJ_dxi  += J[p] * (dxi_dz[p]*d2xyz_dxi2[p](2)    + deta_dz[p]*d2xyz_dxideta[p](2));
              dJ_deta += J[p] * (dxi_dz[p]*d2xyz_dxideta[p](2) + deta_dz[p]*d2xyz_deta2[p](2));
#endif

              for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
                {
                  // dphi_k/dx_l = A + B + C, where (n, m summed over reference directions):
                  //   A = -phi_k(x) * SUM_n (dxi_n/dx_l) * (1/J) * (dJ/dxi_n)
                  //   B =  (1/J) * SUM_n (dxi_n/dx_l) * SUM_m (d^2 x_k/dxi_m dxi_n) * phihat_m
                  //   C =  (1/J) * SUM_n (dxi_n/dx_l) * SUM_m F_{km} * dphihat_m/dxi_n

                  // Term A: -phi_k(x) * SUM_n (dxi_n/dx_l)*(1/J)*(dJ/dxi_n)
                  const Real phik = (k == 0) ? phi_val(0) : (k == 1) ? phi_val(1) : phi_val(2);

                  // Term B: (1/J) * SUM_n (dxi_n/dx_l) * SUM_m (d^2 x_k/dxi_m dxi_n) * phihat_m
                  const Real d2x_k_dxi2    = (k == 0) ? d2xyz_dxi2[p](0)    : (k == 1) ? d2xyz_dxi2[p](1)    : d2xyz_dxi2[p](2);
                  const Real d2x_k_deta2   = (k == 0) ? d2xyz_deta2[p](0)   : (k == 1) ? d2xyz_deta2[p](1)   : d2xyz_deta2[p](2);
                  const Real d2x_k_dxideta = (k == 0) ? d2xyz_dxideta[p](0) : (k == 1) ? d2xyz_dxideta[p](1) : d2xyz_dxideta[p](2);

                  const Real dphik_dxi_B  = (d2x_k_dxi2*phi_ref(0)    + d2x_k_dxideta*phi_ref(1))/J[p];
                  const Real dphik_deta_B = (d2x_k_dxideta*phi_ref(0) + d2x_k_deta2*phi_ref(1))/J[p];

                  // Term C: (1/J) * SUM_n (dxi_n/dx_l) * SUM_m F_{km} * dphihat_m/dxi_n
                  const Real Fk_xi  = (k == 0) ? dxyz_dxi[p](0)  : (k == 1) ? dxyz_dxi[p](1)  : dxyz_dxi[p](2);
                  const Real Fk_eta = (k == 0) ? dxyz_deta[p](0) : (k == 1) ? dxyz_deta[p](1) : dxyz_deta[p](2);

                  const Real dphik_dxi_C  = (Fk_xi*dphi_dxi[i][p](0)  + Fk_eta*dphi_dxi[i][p](1))/J[p];
                  const Real dphik_deta_C = (Fk_xi*dphi_deta[i][p](0) + Fk_eta*dphi_deta[i][p](1))/J[p];

                  const Real dphik_dxi_total  = -phik*dJ_dxi/J[p]  + dphik_dxi_B  + dphik_dxi_C;
                  const Real dphik_deta_total = -phik*dJ_deta/J[p] + dphik_deta_B + dphik_deta_C;

                  dphidx[i][p](k) = dxi_dx[p]*dphik_dxi_total + deta_dx[p]*dphik_deta_total;
                  dphidy[i][p](k) = dxi_dy[p]*dphik_dxi_total + deta_dy[p]*dphik_deta_total;
#if LIBMESH_DIM > 2
                  dphidz[i][p](k) = dxi_dz[p]*dphik_dxi_total + deta_dz[p]*dphik_deta_total;
#endif
                }

              dphi[i][p].slice(0) = dphidx[i][p];
              dphi[i][p].slice(1) = dphidy[i][p];
#if LIBMESH_DIM > 2
              dphi[i][p].slice(2) = dphidz[i][p];
#endif
            }

        break;
      }
    case 3:
      {
        const std::vector<RealGradient> & dxyz_dxi   = fe.get_fe_map().get_dxyzdxi();
        const std::vector<RealGradient> & dxyz_deta  = fe.get_fe_map().get_dxyzdeta();
        const std::vector<RealGradient> & dxyz_dzeta = fe.get_fe_map().get_dxyzdzeta();

        const std::vector<RealGradient> & d2xyz_dxi2     = fe.get_fe_map().get_d2xyzdxi2();
        const std::vector<RealGradient> & d2xyz_deta2    = fe.get_fe_map().get_d2xyzdeta2();
        const std::vector<RealGradient> & d2xyz_dzeta2   = fe.get_fe_map().get_d2xyzdzeta2();
        const std::vector<RealGradient> & d2xyz_dxideta  = fe.get_fe_map().get_d2xyzdxideta();
        const std::vector<RealGradient> & d2xyz_dxidzeta = fe.get_fe_map().get_d2xyzdxidzeta();
        const std::vector<RealGradient> & d2xyz_detadzeta = fe.get_fe_map().get_d2xyzdetadzeta();

        const std::vector<Real> & J = fe.get_fe_map().get_jacobian();

        const std::vector<Real> & dxi_dx   = fe.get_fe_map().get_dxidx();
        const std::vector<Real> & dxi_dy   = fe.get_fe_map().get_dxidy();
        const std::vector<Real> & dxi_dz   = fe.get_fe_map().get_dxidz();
        const std::vector<Real> & deta_dx  = fe.get_fe_map().get_detadx();
        const std::vector<Real> & deta_dy  = fe.get_fe_map().get_detady();
        const std::vector<Real> & deta_dz  = fe.get_fe_map().get_detadz();
        const std::vector<Real> & dzeta_dx = fe.get_fe_map().get_dzetadx();
        const std::vector<Real> & dzeta_dy = fe.get_fe_map().get_dzetady();
        const std::vector<Real> & dzeta_dz = fe.get_fe_map().get_dzetadz();

        const std::vector<std::vector<OutputShape>> & dphi_dxi   = fe.get_dphidxi();
        const std::vector<std::vector<OutputShape>> & dphi_deta  = fe.get_dphideta();
        const std::vector<std::vector<OutputShape>> & dphi_dzeta = fe.get_dphidzeta();

        for (auto i : index_range(dphi))
          for (auto p : index_range(dphi[i]))
            {
              // Need to temporarily cache reference shape functions
              // We are computing mapping basis functions, so we explicitly ignore
              // any non-zero p_level() the Elem might have.
              OutputShape phi_ref;
              FEInterface::shape(fe.get_fe_type(), /*extra_order=*/0, elem, i, qp[p], phi_ref);

              const OutputShape phi_val =
                hdiv_piola_map(dxyz_dxi[p], dxyz_deta[p], dxyz_dzeta[p], J[p], phi_ref);

              // Jacobi's formula: dJ/dxi_n = J * tr(F^{-1} * dF/dxi_n)
              const Real dJ_dxi =
                J[p] * (dxi_dx[p]*d2xyz_dxi2[p](0)     + deta_dx[p]*d2xyz_dxideta[p](0)  + dzeta_dx[p]*d2xyz_dxidzeta[p](0) +
                        dxi_dy[p]*d2xyz_dxi2[p](1)     + deta_dy[p]*d2xyz_dxideta[p](1)  + dzeta_dy[p]*d2xyz_dxidzeta[p](1) +
                        dxi_dz[p]*d2xyz_dxi2[p](2)     + deta_dz[p]*d2xyz_dxideta[p](2)  + dzeta_dz[p]*d2xyz_dxidzeta[p](2));
              const Real dJ_deta =
                J[p] * (dxi_dx[p]*d2xyz_dxideta[p](0)  + deta_dx[p]*d2xyz_deta2[p](0)    + dzeta_dx[p]*d2xyz_detadzeta[p](0) +
                        dxi_dy[p]*d2xyz_dxideta[p](1)  + deta_dy[p]*d2xyz_deta2[p](1)    + dzeta_dy[p]*d2xyz_detadzeta[p](1) +
                        dxi_dz[p]*d2xyz_dxideta[p](2)  + deta_dz[p]*d2xyz_deta2[p](2)    + dzeta_dz[p]*d2xyz_detadzeta[p](2));
              const Real dJ_dzeta =
                J[p] * (dxi_dx[p]*d2xyz_dxidzeta[p](0) + deta_dx[p]*d2xyz_detadzeta[p](0) + dzeta_dx[p]*d2xyz_dzeta2[p](0) +
                        dxi_dy[p]*d2xyz_dxidzeta[p](1) + deta_dy[p]*d2xyz_detadzeta[p](1) + dzeta_dy[p]*d2xyz_dzeta2[p](1) +
                        dxi_dz[p]*d2xyz_dxidzeta[p](2) + deta_dz[p]*d2xyz_detadzeta[p](2) + dzeta_dz[p]*d2xyz_dzeta2[p](2));

              for (unsigned int k = 0; k < 3; ++k)
                {
                  // dphi_k/dx_l = A + B + C, where (n, m summed over reference directions):
                  //   A = -phi_k(x) * SUM_n (dxi_n/dx_l) * (1/J) * (dJ/dxi_n)
                  //   B =  (1/J) * SUM_n (dxi_n/dx_l) * SUM_m (d^2 x_k/dxi_m dxi_n) * phihat_m
                  //   C =  (1/J) * SUM_n (dxi_n/dx_l) * SUM_m F_{km} * dphihat_m/dxi_n

                  // Term A: -phi_k(x) * SUM_n (dxi_n/dx_l)*(1/J)*(dJ/dxi_n)
                  const Real phik = (k == 0) ? phi_val(0) : (k == 1) ? phi_val(1) : phi_val(2);

                  // Term B: (1/J) * SUM_n (dxi_n/dx_l) * SUM_m (d^2 x_k/dxi_m dxi_n) * phihat_m
                  const Real d2x_k_dxi2      = (k == 0) ? d2xyz_dxi2[p](0)      : (k == 1) ? d2xyz_dxi2[p](1)      : d2xyz_dxi2[p](2);
                  const Real d2x_k_deta2     = (k == 0) ? d2xyz_deta2[p](0)     : (k == 1) ? d2xyz_deta2[p](1)     : d2xyz_deta2[p](2);
                  const Real d2x_k_dzeta2    = (k == 0) ? d2xyz_dzeta2[p](0)    : (k == 1) ? d2xyz_dzeta2[p](1)    : d2xyz_dzeta2[p](2);
                  const Real d2x_k_dxideta   = (k == 0) ? d2xyz_dxideta[p](0)   : (k == 1) ? d2xyz_dxideta[p](1)   : d2xyz_dxideta[p](2);
                  const Real d2x_k_dxidzeta  = (k == 0) ? d2xyz_dxidzeta[p](0)  : (k == 1) ? d2xyz_dxidzeta[p](1)  : d2xyz_dxidzeta[p](2);
                  const Real d2x_k_detadzeta = (k == 0) ? d2xyz_detadzeta[p](0) : (k == 1) ? d2xyz_detadzeta[p](1) : d2xyz_detadzeta[p](2);

                  const Real dphik_dxi_B   = (d2x_k_dxi2*phi_ref(0)     + d2x_k_dxideta*phi_ref(1)   + d2x_k_dxidzeta*phi_ref(2))/J[p];
                  const Real dphik_deta_B  = (d2x_k_dxideta*phi_ref(0) + d2x_k_deta2*phi_ref(1)      + d2x_k_detadzeta*phi_ref(2))/J[p];
                  const Real dphik_dzeta_B = (d2x_k_dxidzeta*phi_ref(0) + d2x_k_detadzeta*phi_ref(1) + d2x_k_dzeta2*phi_ref(2))/J[p];

                  // Term C: (1/J) * SUM_n (dxi_n/dx_l) * SUM_m F_{km} * dphihat_m/dxi_n
                  const Real Fk_xi   = (k == 0) ? dxyz_dxi[p](0)   : (k == 1) ? dxyz_dxi[p](1)   : dxyz_dxi[p](2);
                  const Real Fk_eta  = (k == 0) ? dxyz_deta[p](0)  : (k == 1) ? dxyz_deta[p](1)  : dxyz_deta[p](2);
                  const Real Fk_zeta = (k == 0) ? dxyz_dzeta[p](0) : (k == 1) ? dxyz_dzeta[p](1) : dxyz_dzeta[p](2);

                  const Real dphik_dxi_C   = (Fk_xi*dphi_dxi[i][p](0)   + Fk_eta*dphi_dxi[i][p](1)   + Fk_zeta*dphi_dxi[i][p](2))/J[p];
                  const Real dphik_deta_C  = (Fk_xi*dphi_deta[i][p](0)  + Fk_eta*dphi_deta[i][p](1)  + Fk_zeta*dphi_deta[i][p](2))/J[p];
                  const Real dphik_dzeta_C = (Fk_xi*dphi_dzeta[i][p](0) + Fk_eta*dphi_dzeta[i][p](1) + Fk_zeta*dphi_dzeta[i][p](2))/J[p];

                  const Real dphik_dxi_total   = -phik*dJ_dxi/J[p]   + dphik_dxi_B   + dphik_dxi_C;
                  const Real dphik_deta_total  = -phik*dJ_deta/J[p]  + dphik_deta_B  + dphik_deta_C;
                  const Real dphik_dzeta_total = -phik*dJ_dzeta/J[p] + dphik_dzeta_B + dphik_dzeta_C;

                  dphidx[i][p](k) = dxi_dx[p]*dphik_dxi_total  + deta_dx[p]*dphik_deta_total  + dzeta_dx[p]*dphik_dzeta_total;
                  dphidy[i][p](k) = dxi_dy[p]*dphik_dxi_total  + deta_dy[p]*dphik_deta_total  + dzeta_dy[p]*dphik_dzeta_total;
                  dphidz[i][p](k) = dxi_dz[p]*dphik_dxi_total  + deta_dz[p]*dphik_deta_total  + dzeta_dz[p]*dphik_dzeta_total;
                }

              dphi[i][p].slice(0) = dphidx[i][p];
              dphi[i][p].slice(1) = dphidy[i][p];
              dphi[i][p].slice(2) = dphidz[i][p];
            }

        break;
      }

    default:
      libmesh_error_msg("Invalid dim = " << dim);
    } // switch(dim)
#else
  libmesh_ignore(dim, elem, qp, fe, dphi, dphidx, dphidy, dphidz);
  libmesh_error_msg("HDiv shape function gradients require the library to be configured "
                     "with --enable-second-derivatives (LIBMESH_ENABLE_SECOND_DERIVATIVES).");
#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES
}

template<typename OutputShape>
void HDivFETransformation<OutputShape>::map_div(const unsigned int dim,
                                                const Elem * const,
                                                const std::vector<Point> &,
                                                const FEGenericBase<OutputShape> & fe,
                                                std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputDivergence>> & div_phi) const
{
  switch (dim)
    {
    case 0:
    case 1:
      libmesh_error_msg("These element transformations only make sense in 2D and 3D.");

    case 2:
      {
        const std::vector<std::vector<OutputShape>> & dphi_dxi = fe.get_dphidxi();
        const std::vector<std::vector<OutputShape>> & dphi_deta = fe.get_dphideta();

        const std::vector<Real> & J = fe.get_fe_map().get_jacobian();

        // div(phi) = J^{-1} * div(\hat{phi})
        for (auto i : index_range(div_phi))
          for (auto p : index_range(div_phi[i]))
            {
              div_phi[i][p] = (dphi_dxi[i][p](0) + dphi_deta[i][p](1))/J[p];
            }

        break;
      }
    case 3:
      {
        const std::vector<std::vector<OutputShape>> & dphi_dxi = fe.get_dphidxi();
        const std::vector<std::vector<OutputShape>> & dphi_deta = fe.get_dphideta();
        const std::vector<std::vector<OutputShape>> & dphi_dzeta = fe.get_dphidzeta();

        const std::vector<Real> & J = fe.get_fe_map().get_jacobian();

        // div(phi) = J^{-1} * div(\hat{phi})
        for (auto i : index_range(div_phi))
          for (auto p : index_range(div_phi[i]))
            {
              div_phi[i][p] = (dphi_dxi[i][p](0) + dphi_deta[i][p](1) + dphi_dzeta[i][p](2))/J[p];
            }

        break;
      }

    default:
      libmesh_error_msg("Invalid dim = " << dim);
    } // switch(dim)
}

template class LIBMESH_EXPORT HDivFETransformation<RealGradient>;

template<>
void HDivFETransformation<Real>::init_map_phi(const FEGenericBase<Real> & ) const
{
  libmesh_error_msg("HDiv transformations only make sense for vector-valued elements.");
}

template<>
void HDivFETransformation<Real>::init_map_dphi(const FEGenericBase<Real> & ) const
{
  libmesh_error_msg("HDiv transformations only make sense for vector-valued elements.");
}

template<>
void HDivFETransformation<Real>::init_map_d2phi(const FEGenericBase<Real> & ) const
{
  libmesh_error_msg("HDiv transformations only make sense for vector-valued elements.");
}

template<>
void HDivFETransformation<Real>::map_phi(const unsigned int,
                                         const Elem * const,
                                         const std::vector<Point> &,
                                         const FEGenericBase<Real> &,
                                         std::vector<std::vector<Real>> &,
                                         bool) const
{
  libmesh_error_msg("HDiv transformations only make sense for vector-valued elements.");
}

template<>
void HDivFETransformation<Real>::map_dphi(const unsigned int,
                                          const Elem * const,
                                          const std::vector<Point> &,
                                          const FEGenericBase<Real> &,
                                          std::vector<std::vector<FEGenericBase<Real>::OutputGradient>> &,
                                          std::vector<std::vector<Real>> &,
                                          std::vector<std::vector<Real>> &,
                                          std::vector<std::vector<Real>> &) const
{
  libmesh_error_msg("HDiv transformations only make sense for vector-valued elements.");
}

template<>
void HDivFETransformation<Real>::map_div(const unsigned int,
                                         const Elem * const,
                                         const std::vector<Point> &,
                                         const FEGenericBase<Real> &,
                                         std::vector<std::vector<Real>> &) const
{
  libmesh_error_msg("HDiv transformations only make sense for vector-valued elements.");
}


} // namespace libMesh
