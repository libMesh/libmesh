// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
              Real dx_dxi   = dxyz_dxi[p](0);
              Real dx_deta  = dxyz_deta[p](0);

              Real dy_dxi   = dxyz_dxi[p](1);
              Real dy_deta  = dxyz_deta[p](1);

              // Need to temporarily cache reference shape functions
              // We are computing mapping basis functions, so we explicitly ignore
              // any non-zero p_level() the Elem might have.
              OutputShape phi_ref;
              FEInterface::shape(fe.get_fe_type(), /*extra_order=*/0, elem, i, qp[p], phi_ref);

              phi[i][p](0) = (dx_dxi*phi_ref(0) + dx_deta*phi_ref(1))/J[p];
              phi[i][p](1) = (dy_dxi*phi_ref(0) + dy_deta*phi_ref(1))/J[p];
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
              Real dx_dxi   = dxyz_dxi[p](0);
              Real dx_deta  = dxyz_deta[p](0);
              Real dx_dzeta = dxyz_dzeta[p](0);

              Real dy_dxi   = dxyz_dxi[p](1);
              Real dy_deta  = dxyz_deta[p](1);
              Real dy_dzeta = dxyz_dzeta[p](1);

              Real dz_dxi   = dxyz_dxi[p](2);
              Real dz_deta  = dxyz_deta[p](2);
              Real dz_dzeta = dxyz_dzeta[p](2);

              // Need to temporarily cache reference shape functions
              // We are computing mapping basis functions, so we explicitly ignore
              // any non-zero p_level() the Elem might have.
              OutputShape phi_ref;
              FEInterface::shape(fe.get_fe_type(), /*extra_order=*/0, elem, i, qp[p], phi_ref);

              phi[i][p](0) = (dx_dxi*phi_ref(0) + dx_deta*phi_ref(1) + dx_dzeta*phi_ref(2))/J[p];
              phi[i][p](1) = (dy_dxi*phi_ref(0) + dy_deta*phi_ref(1) + dy_dzeta*phi_ref(2))/J[p];
              phi[i][p](2) = (dz_dxi*phi_ref(0) + dz_deta*phi_ref(1) + dz_dzeta*phi_ref(2))/J[p];
            }

        break;
      }

    default:
      libmesh_error_msg("Invalid dim = " << dim);
    } // switch(dim)
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
void HDivFETransformation<Real>::map_div(const unsigned int,
                                         const Elem * const,
                                         const std::vector<Point> &,
                                         const FEGenericBase<Real> &,
                                         std::vector<std::vector<Real>> &) const
{
  libmesh_error_msg("HDiv transformations only make sense for vector-valued elements.");
}


} // namespace libMesh
