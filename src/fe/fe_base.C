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



// Local includes
#include "libmesh/fe.h"
#include "libmesh/inf_fe.h"
#include "libmesh/libmesh_logging.h"
// For projection code:
#include "libmesh/boundary_info.h"
#include "libmesh/mesh_base.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/fe_interface.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/periodic_boundary_base.h"
#include "libmesh/periodic_boundaries.h"
#include "libmesh/quadrature.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/tensor_value.h"
#include "libmesh/threads.h"
#include "libmesh/fe_type.h"

// Anonymous namespace, for a helper function for periodic boundary
// constraint calculations
namespace
{
using namespace libMesh;

// Find the "primary" element around a boundary point:
const Elem * primary_boundary_point_neighbor(const Elem * elem,
                                             const Point & p,
                                             const BoundaryInfo & boundary_info,
                                             const std::set<boundary_id_type> & boundary_ids)
{
  // If we don't find a better alternative, the user will have
  // provided the primary element
  const Elem * primary = elem;

  // Container to catch boundary IDs passed back by BoundaryInfo.
  std::vector<boundary_id_type> bc_ids;

  std::set<const Elem *> point_neighbors;
  elem->find_point_neighbors(p, point_neighbors);
  for (std::set<const Elem *>::const_iterator point_neighbors_iter =
         point_neighbors.begin();
       point_neighbors_iter != point_neighbors.end(); ++point_neighbors_iter)
    {
      const Elem * pt_neighbor = *point_neighbors_iter;

      // If this point neighbor isn't at least
      // as coarse as the current primary elem, or if it is at
      // the same level but has a lower id, then
      // we won't defer to it.
      if ((pt_neighbor->level() > primary->level()) ||
          (pt_neighbor->level() == primary->level() &&
           pt_neighbor->id() < primary->id()))
        continue;

      // Otherwise, we will defer to the point neighbor, but only if
      // one of its sides is on a relevant boundary and that side
      // contains this vertex
      bool vertex_on_periodic_side = false;
      for (unsigned short int ns = 0;
           ns != pt_neighbor->n_sides(); ++ns)
        {
          boundary_info.boundary_ids (pt_neighbor, ns, bc_ids);

          bool on_relevant_boundary = false;
          for (std::set<boundary_id_type>::const_iterator i =
                 boundary_ids.begin(); i != boundary_ids.end(); ++i)
            if (std::find(bc_ids.begin(), bc_ids.end(), *i) != bc_ids.end())
              on_relevant_boundary = true;

          if (!on_relevant_boundary)
            continue;

          if (!pt_neighbor->build_side_ptr(ns)->contains_point(p))
            continue;

          vertex_on_periodic_side = true;
          break;
        }

      if (vertex_on_periodic_side)
        primary = pt_neighbor;
    }

  return primary;
}

// Find the "primary" element around a boundary edge:
const Elem * primary_boundary_edge_neighbor(const Elem * elem,
                                            const Point & p1,
                                            const Point & p2,
                                            const BoundaryInfo & boundary_info,
                                            const std::set<boundary_id_type> & boundary_ids)
{
  // If we don't find a better alternative, the user will have
  // provided the primary element
  const Elem * primary = elem;

  std::set<const Elem *> edge_neighbors;
  elem->find_edge_neighbors(p1, p2, edge_neighbors);

  // Container to catch boundary IDs handed back by BoundaryInfo
  std::vector<boundary_id_type> bc_ids;

  for (std::set<const Elem *>::const_iterator edge_neighbors_iter =
         edge_neighbors.begin();
       edge_neighbors_iter != edge_neighbors.end(); ++edge_neighbors_iter)
    {
      const Elem * e_neighbor = *edge_neighbors_iter;

      // If this edge neighbor isn't at least
      // as coarse as the current primary elem, or if it is at
      // the same level but has a lower id, then
      // we won't defer to it.
      if ((e_neighbor->level() > primary->level()) ||
          (e_neighbor->level() == primary->level() &&
           e_neighbor->id() < primary->id()))
        continue;

      // Otherwise, we will defer to the edge neighbor, but only if
      // one of its sides is on this periodic boundary and that
      // side contains this edge
      bool vertex_on_periodic_side = false;
      for (unsigned short int ns = 0;
           ns != e_neighbor->n_sides(); ++ns)
        {
          boundary_info.boundary_ids (e_neighbor, ns, bc_ids);

          bool on_relevant_boundary = false;
          for (std::set<boundary_id_type>::const_iterator i =
                 boundary_ids.begin(); i != boundary_ids.end(); ++i)
            if (std::find(bc_ids.begin(), bc_ids.end(), *i) != bc_ids.end())
              on_relevant_boundary = true;

          if (!on_relevant_boundary)
            continue;

          UniquePtr<const Elem> periodic_side = e_neighbor->build_side_ptr(ns);
          if (!(periodic_side->contains_point(p1) &&
                periodic_side->contains_point(p2)))
            continue;

          vertex_on_periodic_side = true;
          break;
        }

      if (vertex_on_periodic_side)
        primary = e_neighbor;
    }

  return primary;
}

}

namespace libMesh
{



// ------------------------------------------------------------
// FEBase class members
template <>
UniquePtr<FEGenericBase<Real> >
FEGenericBase<Real>::build (const unsigned int dim,
                            const FEType & fet)
{
  switch (dim)
    {
      // 0D
    case 0:
      {
        switch (fet.family)
          {
          case CLOUGH:
            return UniquePtr<FEBase>(new FE<0,CLOUGH>(fet));

          case HERMITE:
            return UniquePtr<FEBase>(new FE<0,HERMITE>(fet));

          case LAGRANGE:
            return UniquePtr<FEBase>(new FE<0,LAGRANGE>(fet));

          case L2_LAGRANGE:
            return UniquePtr<FEBase>(new FE<0,L2_LAGRANGE>(fet));

          case HIERARCHIC:
            return UniquePtr<FEBase>(new FE<0,HIERARCHIC>(fet));

          case L2_HIERARCHIC:
            return UniquePtr<FEBase>(new FE<0,L2_HIERARCHIC>(fet));

          case MONOMIAL:
            return UniquePtr<FEBase>(new FE<0,MONOMIAL>(fet));

#ifdef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES
          case SZABAB:
            return UniquePtr<FEBase>(new FE<0,SZABAB>(fet));

          case BERNSTEIN:
            return UniquePtr<FEBase>(new FE<0,BERNSTEIN>(fet));
#endif

          case XYZ:
            return UniquePtr<FEBase>(new FEXYZ<0>(fet));

          case SCALAR:
            return UniquePtr<FEBase>(new FEScalar<0>(fet));

          default:
            libmesh_error_msg("ERROR: Bad FEType.family= " << fet.family);
          }
      }
      // 1D
    case 1:
      {
        switch (fet.family)
          {
          case CLOUGH:
            return UniquePtr<FEBase>(new FE<1,CLOUGH>(fet));

          case HERMITE:
            return UniquePtr<FEBase>(new FE<1,HERMITE>(fet));

          case LAGRANGE:
            return UniquePtr<FEBase>(new FE<1,LAGRANGE>(fet));

          case L2_LAGRANGE:
            return UniquePtr<FEBase>(new FE<1,L2_LAGRANGE>(fet));

          case HIERARCHIC:
            return UniquePtr<FEBase>(new FE<1,HIERARCHIC>(fet));

          case L2_HIERARCHIC:
            return UniquePtr<FEBase>(new FE<1,L2_HIERARCHIC>(fet));

          case MONOMIAL:
            return UniquePtr<FEBase>(new FE<1,MONOMIAL>(fet));

#ifdef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES
          case SZABAB:
            return UniquePtr<FEBase>(new FE<1,SZABAB>(fet));

          case BERNSTEIN:
            return UniquePtr<FEBase>(new FE<1,BERNSTEIN>(fet));
#endif

          case XYZ:
            return UniquePtr<FEBase>(new FEXYZ<1>(fet));

          case SCALAR:
            return UniquePtr<FEBase>(new FEScalar<1>(fet));

          default:
            libmesh_error_msg("ERROR: Bad FEType.family= " << fet.family);
          }
      }


      // 2D
    case 2:
      {
        switch (fet.family)
          {
          case CLOUGH:
            return UniquePtr<FEBase>(new FE<2,CLOUGH>(fet));

          case HERMITE:
            return UniquePtr<FEBase>(new FE<2,HERMITE>(fet));

          case LAGRANGE:
            return UniquePtr<FEBase>(new FE<2,LAGRANGE>(fet));

          case L2_LAGRANGE:
            return UniquePtr<FEBase>(new FE<2,L2_LAGRANGE>(fet));

          case HIERARCHIC:
            return UniquePtr<FEBase>(new FE<2,HIERARCHIC>(fet));

          case L2_HIERARCHIC:
            return UniquePtr<FEBase>(new FE<2,L2_HIERARCHIC>(fet));

          case MONOMIAL:
            return UniquePtr<FEBase>(new FE<2,MONOMIAL>(fet));

#ifdef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES
          case SZABAB:
            return UniquePtr<FEBase>(new FE<2,SZABAB>(fet));

          case BERNSTEIN:
            return UniquePtr<FEBase>(new FE<2,BERNSTEIN>(fet));
#endif

          case XYZ:
            return UniquePtr<FEBase>(new FEXYZ<2>(fet));

          case SCALAR:
            return UniquePtr<FEBase>(new FEScalar<2>(fet));

          case SUBDIVISION:
            return UniquePtr<FEBase>(new FESubdivision(fet));

          default:
            libmesh_error_msg("ERROR: Bad FEType.family= " << fet.family);
          }
      }


      // 3D
    case 3:
      {
        switch (fet.family)
          {
          case CLOUGH:
            libmesh_error_msg("ERROR: Clough-Tocher elements currently only support 1D and 2D");

          case HERMITE:
            return UniquePtr<FEBase>(new FE<3,HERMITE>(fet));

          case LAGRANGE:
            return UniquePtr<FEBase>(new FE<3,LAGRANGE>(fet));

          case L2_LAGRANGE:
            return UniquePtr<FEBase>(new FE<3,L2_LAGRANGE>(fet));

          case HIERARCHIC:
            return UniquePtr<FEBase>(new FE<3,HIERARCHIC>(fet));

          case L2_HIERARCHIC:
            return UniquePtr<FEBase>(new FE<3,L2_HIERARCHIC>(fet));

          case MONOMIAL:
            return UniquePtr<FEBase>(new FE<3,MONOMIAL>(fet));

#ifdef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES
          case SZABAB:
            return UniquePtr<FEBase>(new FE<3,SZABAB>(fet));

          case BERNSTEIN:
            return UniquePtr<FEBase>(new FE<3,BERNSTEIN>(fet));
#endif

          case XYZ:
            return UniquePtr<FEBase>(new FEXYZ<3>(fet));

          case SCALAR:
            return UniquePtr<FEBase>(new FEScalar<3>(fet));

          default:
            libmesh_error_msg("ERROR: Bad FEType.family= " << fet.family);
          }
      }

    default:
      libmesh_error_msg("Invalid dimension dim = " << dim);
    }

  libmesh_error_msg("We'll never get here!");
  return UniquePtr<FEBase>();
}



template <>
UniquePtr<FEGenericBase<RealGradient> >
FEGenericBase<RealGradient>::build (const unsigned int dim,
                                    const FEType & fet)
{
  switch (dim)
    {
      // 0D
    case 0:
      {
        switch (fet.family)
          {
          case LAGRANGE_VEC:
            return UniquePtr<FEVectorBase>(new FELagrangeVec<0>(fet));

          default:
            libmesh_error_msg("ERROR: Bad FEType.family= " << fet.family);
          }
      }
    case 1:
      {
        switch (fet.family)
          {
          case LAGRANGE_VEC:
            return UniquePtr<FEVectorBase>(new FELagrangeVec<1>(fet));

          default:
            libmesh_error_msg("ERROR: Bad FEType.family= " << fet.family);
          }
      }
    case 2:
      {
        switch (fet.family)
          {
          case LAGRANGE_VEC:
            return UniquePtr<FEVectorBase>(new FELagrangeVec<2>(fet));

          case NEDELEC_ONE:
            return UniquePtr<FEVectorBase>(new FENedelecOne<2>(fet));

          default:
            libmesh_error_msg("ERROR: Bad FEType.family= " << fet.family);
          }
      }
    case 3:
      {
        switch (fet.family)
          {
          case LAGRANGE_VEC:
            return UniquePtr<FEVectorBase>(new FELagrangeVec<3>(fet));

          case NEDELEC_ONE:
            return UniquePtr<FEVectorBase>(new FENedelecOne<3>(fet));

          default:
            libmesh_error_msg("ERROR: Bad FEType.family= " << fet.family);
          }
      }

    default:
      libmesh_error_msg("Invalid dimension dim = " << dim);
    } // switch(dim)

  libmesh_error_msg("We'll never get here!");
  return UniquePtr<FEVectorBase>();
}







#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS


template <>
UniquePtr<FEGenericBase<Real> >
FEGenericBase<Real>::build_InfFE (const unsigned int dim,
                                  const FEType & fet)
{
  switch (dim)
    {

      // 1D
    case 1:
      {
        switch (fet.radial_family)
          {
          case INFINITE_MAP:
            libmesh_error_msg("ERROR: Can't build an infinite element with FEFamily = " << fet.radial_family);

          case JACOBI_20_00:
            {
              switch (fet.inf_map)
                {
                case CARTESIAN:
                  return UniquePtr<FEBase>(new InfFE<1,JACOBI_20_00,CARTESIAN>(fet));

                default:
                  libmesh_error_msg("ERROR: Can't build an infinite element with InfMapType = " << fet.inf_map);
                }
            }

          case JACOBI_30_00:
            {
              switch (fet.inf_map)
                {
                case CARTESIAN:
                  return UniquePtr<FEBase>(new InfFE<1,JACOBI_30_00,CARTESIAN>(fet));

                default:
                  libmesh_error_msg("ERROR: Can't build an infinite element with InfMapType = " << fet.inf_map);
                }
            }

          case LEGENDRE:
            {
              switch (fet.inf_map)
                {
                case CARTESIAN:
                  return UniquePtr<FEBase>(new InfFE<1,LEGENDRE,CARTESIAN>(fet));

                default:
                  libmesh_error_msg("ERROR: Can't build an infinite element with InfMapType = " << fet.inf_map);
                }
            }

          case LAGRANGE:
            {
              switch (fet.inf_map)
                {
                case CARTESIAN:
                  return UniquePtr<FEBase>(new InfFE<1,LAGRANGE,CARTESIAN>(fet));

                default:
                  libmesh_error_msg("ERROR: Can't build an infinite element with InfMapType = " << fet.inf_map);
                }
            }

          default:
            libmesh_error_msg("ERROR: Bad FEType.radial_family= " << fet.radial_family);
          }
      }




      // 2D
    case 2:
      {
        switch (fet.radial_family)
          {
          case INFINITE_MAP:
            libmesh_error_msg("ERROR: Can't build an infinite element with FEFamily = " << fet.radial_family);

          case JACOBI_20_00:
            {
              switch (fet.inf_map)
                {
                case CARTESIAN:
                  return UniquePtr<FEBase>(new InfFE<2,JACOBI_20_00,CARTESIAN>(fet));

                default:
                  libmesh_error_msg("ERROR: Don't build an infinite element with InfMapType = " << fet.inf_map);
                }
            }

          case JACOBI_30_00:
            {
              switch (fet.inf_map)
                {
                case CARTESIAN:
                  return UniquePtr<FEBase>(new InfFE<2,JACOBI_30_00,CARTESIAN>(fet));

                default:
                  libmesh_error_msg("ERROR: Don't build an infinite element with InfMapType = " << fet.inf_map);
                }
            }

          case LEGENDRE:
            {
              switch (fet.inf_map)
                {
                case CARTESIAN:
                  return UniquePtr<FEBase>(new InfFE<2,LEGENDRE,CARTESIAN>(fet));

                default:
                  libmesh_error_msg("ERROR: Don't build an infinite element with InfMapType = " << fet.inf_map);
                }
            }

          case LAGRANGE:
            {
              switch (fet.inf_map)
                {
                case CARTESIAN:
                  return UniquePtr<FEBase>(new InfFE<2,LAGRANGE,CARTESIAN>(fet));

                default:
                  libmesh_error_msg("ERROR: Don't build an infinite element with InfMapType = " << fet.inf_map);
                }
            }

          default:
            libmesh_error_msg("ERROR: Bad FEType.radial_family= " << fet.radial_family);
          }
      }




      // 3D
    case 3:
      {
        switch (fet.radial_family)
          {
          case INFINITE_MAP:
            libmesh_error_msg("ERROR: Don't build an infinite element with FEFamily = " << fet.radial_family);

          case JACOBI_20_00:
            {
              switch (fet.inf_map)
                {
                case CARTESIAN:
                  return UniquePtr<FEBase>(new InfFE<3,JACOBI_20_00,CARTESIAN>(fet));

                default:
                  libmesh_error_msg("ERROR: Don't build an infinite element with InfMapType = " << fet.inf_map);
                }
            }

          case JACOBI_30_00:
            {
              switch (fet.inf_map)
                {
                case CARTESIAN:
                  return UniquePtr<FEBase>(new InfFE<3,JACOBI_30_00,CARTESIAN>(fet));

                default:
                  libmesh_error_msg("ERROR: Don't build an infinite element with InfMapType = " << fet.inf_map);
                }
            }

          case LEGENDRE:
            {
              switch (fet.inf_map)
                {
                case CARTESIAN:
                  return UniquePtr<FEBase>(new InfFE<3,LEGENDRE,CARTESIAN>(fet));

                default:
                  libmesh_error_msg("ERROR: Don't build an infinite element with InfMapType = " << fet.inf_map);
                }
            }

          case LAGRANGE:
            {
              switch (fet.inf_map)
                {
                case CARTESIAN:
                  return UniquePtr<FEBase>(new InfFE<3,LAGRANGE,CARTESIAN>(fet));

                default:
                  libmesh_error_msg("ERROR: Don't build an infinite element with InfMapType = " << fet.inf_map);
                }
            }

          default:
            libmesh_error_msg("ERROR: Bad FEType.radial_family= " << fet.radial_family);
          }
      }

    default:
      libmesh_error_msg("Invalid dimension dim = " << dim);
    }

  libmesh_error_msg("We'll never get here!");
  return UniquePtr<FEBase>();
}



template <>
UniquePtr<FEGenericBase<RealGradient> >
FEGenericBase<RealGradient>::build_InfFE (const unsigned int,
                                          const FEType & )
{
  // No vector types defined... YET.
  libmesh_not_implemented();
  return UniquePtr<FEVectorBase>();
}

#endif // ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS


template <typename OutputType>
void FEGenericBase<OutputType>::compute_shape_functions (const Elem * elem,
                                                         const std::vector<Point> & qp)
{
  //-------------------------------------------------------------------------
  // Compute the shape function values (and derivatives)
  // at the Quadrature points.  Note that the actual values
  // have already been computed via init_shape_functions

  // Start logging the shape function computation
  LOG_SCOPE("compute_shape_functions()", "FE");

  this->determine_calculations();

  if( calculate_phi )
    this->_fe_trans->map_phi( this->dim, elem, qp, (*this), this->phi );

  if( calculate_dphi )
    this->_fe_trans->map_dphi( this->dim, elem, qp, (*this), this->dphi,
                               this->dphidx, this->dphidy, this->dphidz );

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  if( calculate_d2phi )
    this->_fe_trans->map_d2phi( this->dim, qp, (*this), this->d2phi,
                                this->d2phidx2, this->d2phidxdy, this->d2phidxdz,
                                this->d2phidy2, this->d2phidydz, this->d2phidz2 );
#endif //LIBMESH_ENABLE_SECOND_DERIVATIVES

  // Only compute curl for vector-valued elements
  if( calculate_curl_phi && TypesEqual<OutputType,RealGradient>::value )
    this->_fe_trans->map_curl( this->dim, elem, qp, (*this), this->curl_phi );

  // Only compute div for vector-valued elements
  if( calculate_div_phi && TypesEqual<OutputType,RealGradient>::value )
    this->_fe_trans->map_div( this->dim, elem, qp, (*this), this->div_phi );
}


template <typename OutputType>
void FEGenericBase<OutputType>::print_phi(std::ostream & os) const
{
  for (std::size_t i=0; i<phi.size(); ++i)
    for (std::size_t j=0; j<phi[i].size(); ++j)
      os << " phi[" << i << "][" << j << "]=" << phi[i][j] << std::endl;
}




template <typename OutputType>
void FEGenericBase<OutputType>::print_dphi(std::ostream & os) const
{
  for (std::size_t i=0; i<dphi.size(); ++i)
    for (std::size_t j=0; j<dphi[i].size(); ++j)
      os << " dphi[" << i << "][" << j << "]=" << dphi[i][j];
}



template <typename OutputType>
void FEGenericBase<OutputType>::determine_calculations()
{
  this->calculations_started = true;

  // If the user forgot to request anything, we'll be safe and
  // calculate everything:
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  if (!this->calculate_phi && !this->calculate_dphi && !this->calculate_d2phi
      && !this->calculate_curl_phi && !this->calculate_div_phi)
    {
      this->calculate_phi = this->calculate_dphi = this->calculate_d2phi = this->calculate_dphiref = true;
      if( FEInterface::field_type(fe_type.family) == TYPE_VECTOR )
        {
          this->calculate_curl_phi = true;
          this->calculate_div_phi  = true;
        }
    }
#else
  if (!this->calculate_phi && !this->calculate_dphi && !this->calculate_curl_phi && !this->calculate_div_phi)
    {
      this->calculate_phi = this->calculate_dphi = this->calculate_dphiref = true;
      if( FEInterface::field_type(fe_type.family) == TYPE_VECTOR )
        {
          this->calculate_curl_phi = true;
          this->calculate_div_phi  = true;
        }
    }
#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES

  // Request whichever terms are necessary from the FEMap
  if( this->calculate_phi )
    this->_fe_trans->init_map_phi(*this);

  if( this->calculate_dphiref )
    this->_fe_trans->init_map_dphi(*this);

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  if( this->calculate_d2phi )
    this->_fe_trans->init_map_d2phi(*this);
#endif //LIBMESH_ENABLE_SECOND_DERIVATIVES
}



#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES


template <typename OutputType>
void FEGenericBase<OutputType>::print_d2phi(std::ostream & os) const
{
  for (std::size_t i=0; i<dphi.size(); ++i)
    for (std::size_t j=0; j<dphi[i].size(); ++j)
      os << " d2phi[" << i << "][" << j << "]=" << d2phi[i][j];
}

#endif



#ifdef LIBMESH_ENABLE_AMR

template <typename OutputType>
void
FEGenericBase<OutputType>::coarsened_dof_values(const NumericVector<Number> & old_vector,
                                                const DofMap & dof_map,
                                                const Elem * elem,
                                                DenseVector<Number> & Ue,
                                                const unsigned int var,
                                                const bool use_old_dof_indices)
{
  // Side/edge local DOF indices
  std::vector<unsigned int> new_side_dofs, old_side_dofs;

  // FIXME: what about 2D shells in 3D space?
  unsigned int dim = elem->dim();

  // We use local FE objects for now
  // FIXME: we should use more, external objects instead for efficiency
  const FEType & base_fe_type = dof_map.variable_type(var);
  UniquePtr<FEGenericBase<OutputShape> > fe
    (FEGenericBase<OutputShape>::build(dim, base_fe_type));
  UniquePtr<FEGenericBase<OutputShape> > fe_coarse
    (FEGenericBase<OutputShape>::build(dim, base_fe_type));

  UniquePtr<QBase> qrule     (base_fe_type.default_quadrature_rule(dim));
  UniquePtr<QBase> qedgerule (base_fe_type.default_quadrature_rule(1));
  UniquePtr<QBase> qsiderule (base_fe_type.default_quadrature_rule(dim-1));
  std::vector<Point> coarse_qpoints;

  // The values of the shape functions at the quadrature
  // points
  const std::vector<std::vector<OutputShape> > & phi_values =
    fe->get_phi();
  const std::vector<std::vector<OutputShape> > & phi_coarse =
    fe_coarse->get_phi();

  // The gradients of the shape functions at the quadrature
  // points on the child element.
  const std::vector<std::vector<OutputGradient> > * dphi_values =
    libmesh_nullptr;
  const std::vector<std::vector<OutputGradient> > * dphi_coarse =
    libmesh_nullptr;

  const FEContinuity cont = fe->get_continuity();

  if (cont == C_ONE)
    {
      const std::vector<std::vector<OutputGradient> > &
        ref_dphi_values = fe->get_dphi();
      dphi_values = &ref_dphi_values;
      const std::vector<std::vector<OutputGradient> > &
        ref_dphi_coarse = fe_coarse->get_dphi();
      dphi_coarse = &ref_dphi_coarse;
    }

  // The Jacobian * quadrature weight at the quadrature points
  const std::vector<Real> & JxW =
    fe->get_JxW();

  // The XYZ locations of the quadrature points on the
  // child element
  const std::vector<Point> & xyz_values =
    fe->get_xyz();



  FEType fe_type = base_fe_type, temp_fe_type;
  const ElemType elem_type = elem->type();
  fe_type.order = static_cast<Order>(fe_type.order +
                                     elem->max_descendant_p_level());

  // Number of nodes on parent element
  const unsigned int n_nodes = elem->n_nodes();

  // Number of dofs on parent element
  const unsigned int new_n_dofs =
    FEInterface::n_dofs(dim, fe_type, elem_type);

  // Fixed vs. free DoFs on edge/face projections
  std::vector<char> dof_is_fixed(new_n_dofs, false); // bools
  std::vector<int> free_dof(new_n_dofs, 0);

  DenseMatrix<Real> Ke;
  DenseVector<Number> Fe;
  Ue.resize(new_n_dofs); Ue.zero();


  // When coarsening, in general, we need a series of
  // projections to ensure a unique and continuous
  // solution.  We start by interpolating nodes, then
  // hold those fixed and project edges, then
  // hold those fixed and project faces, then
  // hold those fixed and project interiors

  // Copy node values first
  {
    std::vector<dof_id_type> node_dof_indices;
    if (use_old_dof_indices)
      dof_map.old_dof_indices (elem, node_dof_indices, var);
    else
      dof_map.dof_indices (elem, node_dof_indices, var);

    unsigned int current_dof = 0;
    for (unsigned int n=0; n!= n_nodes; ++n)
      {
        // FIXME: this should go through the DofMap,
        // not duplicate dof_indices code badly!
        const unsigned int my_nc =
          FEInterface::n_dofs_at_node (dim, fe_type,
                                       elem_type, n);
        if (!elem->is_vertex(n))
          {
            current_dof += my_nc;
            continue;
          }

        temp_fe_type = base_fe_type;
        // We're assuming here that child n shares vertex n,
        // which is wrong on non-simplices right now
        // ... but this code isn't necessary except on elements
        // where p refinement creates more vertex dofs; we have
        // no such elements yet.
        /*
          if (elem->child_ptr(n)->p_level() < elem->p_level())
          {
          temp_fe_type.order =
          static_cast<Order>(temp_fe_type.order +
          elem->child_ptr(n)->p_level());
          }
        */
        const unsigned int nc =
          FEInterface::n_dofs_at_node (dim, temp_fe_type,
                                       elem_type, n);
        for (unsigned int i=0; i!= nc; ++i)
          {
            Ue(current_dof) =
              old_vector(node_dof_indices[current_dof]);
            dof_is_fixed[current_dof] = true;
            current_dof++;
          }
      }
  }

  // In 3D, project any edge values next
  if (dim > 2 && cont != DISCONTINUOUS)
    for (unsigned int e=0; e != elem->n_edges(); ++e)
      {
        FEInterface::dofs_on_edge(elem, dim, fe_type,
                                  e, new_side_dofs);

        // Some edge dofs are on nodes and already
        // fixed, others are free to calculate
        unsigned int free_dofs = 0;
        for (std::size_t i=0; i != new_side_dofs.size(); ++i)
          if (!dof_is_fixed[new_side_dofs[i]])
            free_dof[free_dofs++] = i;
        Ke.resize (free_dofs, free_dofs); Ke.zero();
        Fe.resize (free_dofs); Fe.zero();
        // The new edge coefficients
        DenseVector<Number> Uedge(free_dofs);

        // Add projection terms from each child sharing
        // this edge
        for (unsigned int c=0; c != elem->n_children();
             ++c)
          {
            if (!elem->is_child_on_edge(c,e))
              continue;
            const Elem * child = elem->child_ptr(c);

            std::vector<dof_id_type> child_dof_indices;
            if (use_old_dof_indices)
              dof_map.old_dof_indices (child,
                                       child_dof_indices, var);
            else
              dof_map.dof_indices (child,
                                   child_dof_indices, var);
            const unsigned int child_n_dofs =
              cast_int<unsigned int>
              (child_dof_indices.size());

            temp_fe_type = base_fe_type;
            temp_fe_type.order =
              static_cast<Order>(temp_fe_type.order +
                                 child->p_level());

            FEInterface::dofs_on_edge(child, dim,
                                      temp_fe_type, e, old_side_dofs);

            // Initialize both child and parent FE data
            // on the child's edge
            fe->attach_quadrature_rule (qedgerule.get());
            fe->edge_reinit (child, e);
            const unsigned int n_qp = qedgerule->n_points();

            FEInterface::inverse_map (dim, fe_type, elem,
                                      xyz_values, coarse_qpoints);

            fe_coarse->reinit(elem, &coarse_qpoints);

            // Loop over the quadrature points
            for (unsigned int qp=0; qp<n_qp; qp++)
              {
                // solution value at the quadrature point
                OutputNumber fineval = libMesh::zero;
                // solution grad at the quadrature point
                OutputNumberGradient finegrad;

                // Sum the solution values * the DOF
                // values at the quadrature point to
                // get the solution value and gradient.
                for (unsigned int i=0; i<child_n_dofs;
                     i++)
                  {
                    fineval +=
                      (old_vector(child_dof_indices[i])*
                       phi_values[i][qp]);
                    if (cont == C_ONE)
                      finegrad += (*dphi_values)[i][qp] *
                        old_vector(child_dof_indices[i]);
                  }

                // Form edge projection matrix
                for (std::size_t sidei=0, freei=0; sidei != new_side_dofs.size(); ++sidei)
                  {
                    unsigned int i = new_side_dofs[sidei];
                    // fixed DoFs aren't test functions
                    if (dof_is_fixed[i])
                      continue;
                    for (std::size_t sidej=0, freej=0; sidej != new_side_dofs.size(); ++sidej)
                      {
                        unsigned int j =
                          new_side_dofs[sidej];
                        if (dof_is_fixed[j])
                          Fe(freei) -=
                            TensorTools::inner_product(phi_coarse[i][qp],
                                                       phi_coarse[j][qp]) *
                            JxW[qp] * Ue(j);
                        else
                          Ke(freei,freej) +=
                            TensorTools::inner_product(phi_coarse[i][qp],
                                                       phi_coarse[j][qp]) *
                            JxW[qp];
                        if (cont == C_ONE)
                          {
                            if (dof_is_fixed[j])
                              Fe(freei) -=
                                TensorTools::inner_product((*dphi_coarse)[i][qp],
                                                           (*dphi_coarse)[j][qp]) *
                                JxW[qp] * Ue(j);
                            else
                              Ke(freei,freej) +=
                                TensorTools::inner_product((*dphi_coarse)[i][qp],
                                                           (*dphi_coarse)[j][qp]) *
                                JxW[qp];
                          }
                        if (!dof_is_fixed[j])
                          freej++;
                      }
                    Fe(freei) += TensorTools::inner_product(phi_coarse[i][qp],
                                                            fineval) * JxW[qp];
                    if (cont == C_ONE)
                      Fe(freei) +=
                        TensorTools::inner_product(finegrad, (*dphi_coarse)[i][qp]) * JxW[qp];
                    freei++;
                  }
              }
          }
        Ke.cholesky_solve(Fe, Uedge);

        // Transfer new edge solutions to element
        for (unsigned int i=0; i != free_dofs; ++i)
          {
            Number & ui = Ue(new_side_dofs[free_dof[i]]);
            libmesh_assert(std::abs(ui) < TOLERANCE ||
                           std::abs(ui - Uedge(i)) < TOLERANCE);
            ui = Uedge(i);
            dof_is_fixed[new_side_dofs[free_dof[i]]] = true;
          }
      }

  // Project any side values (edges in 2D, faces in 3D)
  if (dim > 1 && cont != DISCONTINUOUS)
    for (unsigned int s=0; s != elem->n_sides(); ++s)
      {
        FEInterface::dofs_on_side(elem, dim, fe_type,
                                  s, new_side_dofs);

        // Some side dofs are on nodes/edges and already
        // fixed, others are free to calculate
        unsigned int free_dofs = 0;
        for (std::size_t i=0; i != new_side_dofs.size(); ++i)
          if (!dof_is_fixed[new_side_dofs[i]])
            free_dof[free_dofs++] = i;
        Ke.resize (free_dofs, free_dofs); Ke.zero();
        Fe.resize (free_dofs); Fe.zero();
        // The new side coefficients
        DenseVector<Number> Uside(free_dofs);

        // Add projection terms from each child sharing
        // this side
        for (unsigned int c=0; c != elem->n_children();
             ++c)
          {
            if (!elem->is_child_on_side(c,s))
              continue;
            const Elem * child = elem->child_ptr(c);

            std::vector<dof_id_type> child_dof_indices;
            if (use_old_dof_indices)
              dof_map.old_dof_indices (child,
                                       child_dof_indices, var);
            else
              dof_map.dof_indices (child,
                                   child_dof_indices, var);
            const unsigned int child_n_dofs =
              cast_int<unsigned int>
              (child_dof_indices.size());

            temp_fe_type = base_fe_type;
            temp_fe_type.order =
              static_cast<Order>(temp_fe_type.order +
                                 child->p_level());

            FEInterface::dofs_on_side(child, dim,
                                      temp_fe_type, s, old_side_dofs);

            // Initialize both child and parent FE data
            // on the child's side
            fe->attach_quadrature_rule (qsiderule.get());
            fe->reinit (child, s);
            const unsigned int n_qp = qsiderule->n_points();

            FEInterface::inverse_map (dim, fe_type, elem,
                                      xyz_values, coarse_qpoints);

            fe_coarse->reinit(elem, &coarse_qpoints);

            // Loop over the quadrature points
            for (unsigned int qp=0; qp<n_qp; qp++)
              {
                // solution value at the quadrature point
                OutputNumber fineval = libMesh::zero;
                // solution grad at the quadrature point
                OutputNumberGradient finegrad;

                // Sum the solution values * the DOF
                // values at the quadrature point to
                // get the solution value and gradient.
                for (unsigned int i=0; i<child_n_dofs;
                     i++)
                  {
                    fineval +=
                      old_vector(child_dof_indices[i]) *
                      phi_values[i][qp];
                    if (cont == C_ONE)
                      finegrad += (*dphi_values)[i][qp] *
                        old_vector(child_dof_indices[i]);
                  }

                // Form side projection matrix
                for (std::size_t sidei=0, freei=0; sidei != new_side_dofs.size(); ++sidei)
                  {
                    unsigned int i = new_side_dofs[sidei];
                    // fixed DoFs aren't test functions
                    if (dof_is_fixed[i])
                      continue;
                    for (std::size_t sidej=0, freej=0; sidej != new_side_dofs.size(); ++sidej)
                      {
                        unsigned int j =
                          new_side_dofs[sidej];
                        if (dof_is_fixed[j])
                          Fe(freei) -=
                            TensorTools::inner_product(phi_coarse[i][qp],
                                                       phi_coarse[j][qp]) *
                            JxW[qp] * Ue(j);
                        else
                          Ke(freei,freej) +=
                            TensorTools::inner_product(phi_coarse[i][qp],
                                                       phi_coarse[j][qp]) *
                            JxW[qp];
                        if (cont == C_ONE)
                          {
                            if (dof_is_fixed[j])
                              Fe(freei) -=
                                TensorTools::inner_product((*dphi_coarse)[i][qp],
                                                           (*dphi_coarse)[j][qp]) *
                                JxW[qp] * Ue(j);
                            else
                              Ke(freei,freej) +=
                                TensorTools::inner_product((*dphi_coarse)[i][qp],
                                                           (*dphi_coarse)[j][qp]) *
                                JxW[qp];
                          }
                        if (!dof_is_fixed[j])
                          freej++;
                      }
                    Fe(freei) += TensorTools::inner_product(fineval, phi_coarse[i][qp]) * JxW[qp];
                    if (cont == C_ONE)
                      Fe(freei) +=
                        TensorTools::inner_product(finegrad, (*dphi_coarse)[i][qp]) * JxW[qp];
                    freei++;
                  }
              }
          }
        Ke.cholesky_solve(Fe, Uside);

        // Transfer new side solutions to element
        for (unsigned int i=0; i != free_dofs; ++i)
          {
            Number & ui = Ue(new_side_dofs[free_dof[i]]);
            libmesh_assert(std::abs(ui) < TOLERANCE ||
                           std::abs(ui - Uside(i)) < TOLERANCE);
            ui = Uside(i);
            dof_is_fixed[new_side_dofs[free_dof[i]]] = true;
          }
      }

  // Project the interior values, finally

  // Some interior dofs are on nodes/edges/sides and
  // already fixed, others are free to calculate
  unsigned int free_dofs = 0;
  for (unsigned int i=0; i != new_n_dofs; ++i)
    if (!dof_is_fixed[i])
      free_dof[free_dofs++] = i;
  Ke.resize (free_dofs, free_dofs); Ke.zero();
  Fe.resize (free_dofs); Fe.zero();
  // The new interior coefficients
  DenseVector<Number> Uint(free_dofs);

  // Add projection terms from each child
  for (unsigned int c=0; c != elem->n_children(); ++c)
    {
      const Elem * child = elem->child_ptr(c);

      std::vector<dof_id_type> child_dof_indices;
      if (use_old_dof_indices)
        dof_map.old_dof_indices (child,
                                 child_dof_indices, var);
      else
        dof_map.dof_indices (child,
                             child_dof_indices, var);
      const unsigned int child_n_dofs =
        cast_int<unsigned int>
        (child_dof_indices.size());

      // Initialize both child and parent FE data
      // on the child's quadrature points
      fe->attach_quadrature_rule (qrule.get());
      fe->reinit (child);
      const unsigned int n_qp = qrule->n_points();

      FEInterface::inverse_map (dim, fe_type, elem,
                                xyz_values, coarse_qpoints);

      fe_coarse->reinit(elem, &coarse_qpoints);

      // Loop over the quadrature points
      for (unsigned int qp=0; qp<n_qp; qp++)
        {
          // solution value at the quadrature point
          OutputNumber fineval = libMesh::zero;
          // solution grad at the quadrature point
          OutputNumberGradient finegrad;

          // Sum the solution values * the DOF
          // values at the quadrature point to
          // get the solution value and gradient.
          for (unsigned int i=0; i<child_n_dofs; i++)
            {
              fineval +=
                (old_vector(child_dof_indices[i]) *
                 phi_values[i][qp]);
              if (cont == C_ONE)
                finegrad += (*dphi_values)[i][qp] *
                  old_vector(child_dof_indices[i]);
            }

          // Form interior projection matrix
          for (unsigned int i=0, freei=0;
               i != new_n_dofs; ++i)
            {
              // fixed DoFs aren't test functions
              if (dof_is_fixed[i])
                continue;
              for (unsigned int j=0, freej=0; j !=
                     new_n_dofs; ++j)
                {
                  if (dof_is_fixed[j])
                    Fe(freei) -=
                      TensorTools::inner_product(phi_coarse[i][qp],
                                                 phi_coarse[j][qp]) *
                      JxW[qp] * Ue(j);
                  else
                    Ke(freei,freej) +=
                      TensorTools::inner_product(phi_coarse[i][qp],
                                                 phi_coarse[j][qp]) *
                      JxW[qp];
                  if (cont == C_ONE)
                    {
                      if (dof_is_fixed[j])
                        Fe(freei) -=
                          TensorTools::inner_product((*dphi_coarse)[i][qp],
                                                     (*dphi_coarse)[j][qp]) *
                          JxW[qp] * Ue(j);
                      else
                        Ke(freei,freej) +=
                          TensorTools::inner_product((*dphi_coarse)[i][qp],
                                                     (*dphi_coarse)[j][qp]) *
                          JxW[qp];
                    }
                  if (!dof_is_fixed[j])
                    freej++;
                }
              Fe(freei) += TensorTools::inner_product(phi_coarse[i][qp], fineval) *
                JxW[qp];
              if (cont == C_ONE)
                Fe(freei) += TensorTools::inner_product(finegrad, (*dphi_coarse)[i][qp]) * JxW[qp];
              freei++;
            }
        }
    }
  Ke.cholesky_solve(Fe, Uint);

  // Transfer new interior solutions to element
  for (unsigned int i=0; i != free_dofs; ++i)
    {
      Number & ui = Ue(free_dof[i]);
      libmesh_assert(std::abs(ui) < TOLERANCE ||
                     std::abs(ui - Uint(i)) < TOLERANCE);
      ui = Uint(i);
      // We should be fixing all dofs by now; no need to keep track of
      // that unless we're debugging
#ifndef NDEBUG
      dof_is_fixed[free_dof[i]] = true;
#endif
    }

#ifndef NDEBUG
  // Make sure every DoF got reached!
  for (unsigned int i=0; i != new_n_dofs; ++i)
    libmesh_assert(dof_is_fixed[i]);
#endif
}



template <typename OutputType>
void
FEGenericBase<OutputType>::coarsened_dof_values(const NumericVector<Number> & old_vector,
                                                const DofMap & dof_map,
                                                const Elem * elem,
                                                DenseVector<Number> & Ue,
                                                const bool use_old_dof_indices)
{
  Ue.resize(0);

  for (unsigned int v=0; v != dof_map.n_variables(); ++v)
    {
      DenseVector<Number> Usub;

      coarsened_dof_values(old_vector, dof_map, elem, Usub,
                           use_old_dof_indices);

      Ue.append (Usub);
    }
}



template <typename OutputType>
void
FEGenericBase<OutputType>::compute_proj_constraints (DofConstraints & constraints,
                                                     DofMap & dof_map,
                                                     const unsigned int variable_number,
                                                     const Elem * elem)
{
  libmesh_assert(elem);

  const unsigned int Dim = elem->dim();

  // Only constrain elements in 2,3D.
  if (Dim == 1)
    return;

  // Only constrain active elements with this method
  if (!elem->active())
    return;

  const FEType & base_fe_type = dof_map.variable_type(variable_number);

  // Construct FE objects for this element and its neighbors.
  UniquePtr<FEGenericBase<OutputShape> > my_fe
    (FEGenericBase<OutputShape>::build(Dim, base_fe_type));
  const FEContinuity cont = my_fe->get_continuity();

  // We don't need to constrain discontinuous elements
  if (cont == DISCONTINUOUS)
    return;
  libmesh_assert (cont == C_ZERO || cont == C_ONE);

  UniquePtr<FEGenericBase<OutputShape> > neigh_fe
    (FEGenericBase<OutputShape>::build(Dim, base_fe_type));

  QGauss my_qface(Dim-1, base_fe_type.default_quadrature_order());
  my_fe->attach_quadrature_rule (&my_qface);
  std::vector<Point> neigh_qface;

  const std::vector<Real> & JxW = my_fe->get_JxW();
  const std::vector<Point> & q_point = my_fe->get_xyz();
  const std::vector<std::vector<OutputShape> > & phi = my_fe->get_phi();
  const std::vector<std::vector<OutputShape> > & neigh_phi =
    neigh_fe->get_phi();
  const std::vector<Point> * face_normals = libmesh_nullptr;
  const std::vector<std::vector<OutputGradient> > * dphi = libmesh_nullptr;
  const std::vector<std::vector<OutputGradient> > * neigh_dphi = libmesh_nullptr;

  std::vector<dof_id_type> my_dof_indices, neigh_dof_indices;
  std::vector<unsigned int> my_side_dofs, neigh_side_dofs;

  if (cont != C_ZERO)
    {
      const std::vector<Point> & ref_face_normals =
        my_fe->get_normals();
      face_normals = &ref_face_normals;
      const std::vector<std::vector<OutputGradient> > & ref_dphi =
        my_fe->get_dphi();
      dphi = &ref_dphi;
      const std::vector<std::vector<OutputGradient> > & ref_neigh_dphi =
        neigh_fe->get_dphi();
      neigh_dphi = &ref_neigh_dphi;
    }

  DenseMatrix<Real> Ke;
  DenseVector<Real> Fe;
  std::vector<DenseVector<Real> > Ue;

  // Look at the element faces.  Check to see if we need to
  // build constraints.
  for (unsigned int s=0; s<elem->n_sides(); s++)
    if (elem->neighbor_ptr(s) != libmesh_nullptr)
      {
        // Get pointers to the element's neighbor.
        const Elem * neigh = elem->neighbor_ptr(s);

        // h refinement constraints:
        // constrain dofs shared between
        // this element and ones coarser
        // than this element.
        if (neigh->level() < elem->level())
          {
            unsigned int s_neigh = neigh->which_neighbor_am_i(elem);
            libmesh_assert_less (s_neigh, neigh->n_neighbors());

            // Find the minimum p level; we build the h constraint
            // matrix with this and then constrain away all higher p
            // DoFs.
            libmesh_assert(neigh->active());
            const unsigned int min_p_level =
              std::min(elem->p_level(), neigh->p_level());

            // we may need to make the FE objects reinit with the
            // minimum shared p_level
            const unsigned int old_elem_level = elem->p_level();
            if (elem->p_level() != min_p_level)
              my_fe->set_fe_order(my_fe->get_fe_type().order.get_order() - old_elem_level + min_p_level);
            const unsigned int old_neigh_level = neigh->p_level();
            if (old_neigh_level != min_p_level)
              neigh_fe->set_fe_order(neigh_fe->get_fe_type().order.get_order() - old_neigh_level + min_p_level);

            my_fe->reinit(elem, s);

            // This function gets called element-by-element, so there
            // will be a lot of memory allocation going on.  We can
            // at least minimize this for the case of the dof indices
            // by efficiently preallocating the requisite storage.
            // n_nodes is not necessarily n_dofs, but it is better
            // than nothing!
            my_dof_indices.reserve    (elem->n_nodes());
            neigh_dof_indices.reserve (neigh->n_nodes());

            dof_map.dof_indices (elem, my_dof_indices,
                                 variable_number,
                                 min_p_level);
            dof_map.dof_indices (neigh, neigh_dof_indices,
                                 variable_number,
                                 min_p_level);

            const unsigned int n_qp = my_qface.n_points();

            FEInterface::inverse_map (Dim, base_fe_type, neigh,
                                      q_point, neigh_qface);

            neigh_fe->reinit(neigh, &neigh_qface);

            // We're only concerned with DOFs whose values (and/or first
            // derivatives for C1 elements) are supported on side nodes
            FEType elem_fe_type = base_fe_type;
            if (old_elem_level != min_p_level)
              elem_fe_type.order = base_fe_type.order.get_order() - old_elem_level + min_p_level;
            FEType neigh_fe_type = base_fe_type;
            if (old_neigh_level != min_p_level)
              neigh_fe_type.order = base_fe_type.order.get_order() - old_neigh_level + min_p_level;
            FEInterface::dofs_on_side(elem,  Dim, elem_fe_type,  s,       my_side_dofs);
            FEInterface::dofs_on_side(neigh, Dim, neigh_fe_type, s_neigh, neigh_side_dofs);

            const unsigned int n_side_dofs =
              cast_int<unsigned int>(my_side_dofs.size());
            libmesh_assert_equal_to (n_side_dofs, neigh_side_dofs.size());

            Ke.resize (n_side_dofs, n_side_dofs);
            Ue.resize(n_side_dofs);

            // Form the projection matrix, (inner product of fine basis
            // functions against fine test functions)
            for (unsigned int is = 0; is != n_side_dofs; ++is)
              {
                const unsigned int i = my_side_dofs[is];
                for (unsigned int js = 0; js != n_side_dofs; ++js)
                  {
                    const unsigned int j = my_side_dofs[js];
                    for (unsigned int qp = 0; qp != n_qp; ++qp)
                      {
                        Ke(is,js) += JxW[qp] * TensorTools::inner_product(phi[i][qp], phi[j][qp]);
                        if (cont != C_ZERO)
                          Ke(is,js) += JxW[qp] *
                            TensorTools::inner_product((*dphi)[i][qp] *
                                                       (*face_normals)[qp],
                                                       (*dphi)[j][qp] *
                                                       (*face_normals)[qp]);
                      }
                  }
              }

            // Form the right hand sides, (inner product of coarse basis
            // functions against fine test functions)
            for (unsigned int is = 0; is != n_side_dofs; ++is)
              {
                const unsigned int i = neigh_side_dofs[is];
                Fe.resize (n_side_dofs);
                for (unsigned int js = 0; js != n_side_dofs; ++js)
                  {
                    const unsigned int j = my_side_dofs[js];
                    for (unsigned int qp = 0; qp != n_qp; ++qp)
                      {
                        Fe(js) += JxW[qp] *
                          TensorTools::inner_product(neigh_phi[i][qp],
                                                     phi[j][qp]);
                        if (cont != C_ZERO)
                          Fe(js) += JxW[qp] *
                            TensorTools::inner_product((*neigh_dphi)[i][qp] *
                                                       (*face_normals)[qp],
                                                       (*dphi)[j][qp] *
                                                       (*face_normals)[qp]);
                      }
                  }
                Ke.cholesky_solve(Fe, Ue[is]);
              }

            for (unsigned int js = 0; js != n_side_dofs; ++js)
              {
                const unsigned int j = my_side_dofs[js];
                const dof_id_type my_dof_g = my_dof_indices[j];
                libmesh_assert_not_equal_to (my_dof_g, DofObject::invalid_id);

                // Hunt for "constraining against myself" cases before
                // we bother creating a constraint row
                bool self_constraint = false;
                for (unsigned int is = 0; is != n_side_dofs; ++is)
                  {
                    const unsigned int i = neigh_side_dofs[is];
                    const dof_id_type their_dof_g = neigh_dof_indices[i];
                    libmesh_assert_not_equal_to (their_dof_g, DofObject::invalid_id);

                    if (their_dof_g == my_dof_g)
                      {
#ifndef NDEBUG
                        const Real their_dof_value = Ue[is](js);
                        libmesh_assert_less (std::abs(their_dof_value-1.),
                                             10*TOLERANCE);

                        for (unsigned int k = 0; k != n_side_dofs; ++k)
                          libmesh_assert(k == is ||
                                         std::abs(Ue[k](js)) <
                                         10*TOLERANCE);
#endif

                        self_constraint = true;
                        break;
                      }
                  }

                if (self_constraint)
                  continue;

                DofConstraintRow * constraint_row;

                // we may be running constraint methods concurrently
                // on multiple threads, so we need a lock to
                // ensure that this constraint is "ours"
                {
                  Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);

                  if (dof_map.is_constrained_dof(my_dof_g))
                    continue;

                  constraint_row = &(constraints[my_dof_g]);
                  libmesh_assert(constraint_row->empty());
                }

                for (unsigned int is = 0; is != n_side_dofs; ++is)
                  {
                    const unsigned int i = neigh_side_dofs[is];
                    const dof_id_type their_dof_g = neigh_dof_indices[i];
                    libmesh_assert_not_equal_to (their_dof_g, DofObject::invalid_id);
                    libmesh_assert_not_equal_to (their_dof_g, my_dof_g);

                    const Real their_dof_value = Ue[is](js);

                    if (std::abs(their_dof_value) < 10*TOLERANCE)
                      continue;

                    constraint_row->insert(std::make_pair(their_dof_g,
                                                          their_dof_value));
                  }
              }

            my_fe->set_fe_order(my_fe->get_fe_type().order.get_order() + old_elem_level - min_p_level);
            neigh_fe->set_fe_order(neigh_fe->get_fe_type().order.get_order() + old_neigh_level - min_p_level);
          }

        // p refinement constraints:
        // constrain dofs shared between
        // active elements and neighbors with
        // lower polynomial degrees
        const unsigned int min_p_level =
          neigh->min_p_level_by_neighbor(elem, elem->p_level());
        if (min_p_level < elem->p_level())
          {
            // Adaptive p refinement of non-hierarchic bases will
            // require more coding
            libmesh_assert(my_fe->is_hierarchic());
            dof_map.constrain_p_dofs(variable_number, elem,
                                     s, min_p_level);
          }
      }
}

#endif // #ifdef LIBMESH_ENABLE_AMR



#ifdef LIBMESH_ENABLE_PERIODIC
template <typename OutputType>
void
FEGenericBase<OutputType>::
compute_periodic_constraints (DofConstraints & constraints,
                              DofMap & dof_map,
                              const PeriodicBoundaries & boundaries,
                              const MeshBase & mesh,
                              const PointLocatorBase * point_locator,
                              const unsigned int variable_number,
                              const Elem * elem)
{
  // Only bother if we truly have periodic boundaries
  if (boundaries.empty())
    return;

  libmesh_assert(elem);

  // Only constrain active elements with this method
  if (!elem->active())
    return;

  const unsigned int Dim = elem->dim();

  // We need sys_number and variable_number for DofObject methods
  // later
  const unsigned int sys_number = dof_map.sys_number();

  const FEType & base_fe_type = dof_map.variable_type(variable_number);

  // Construct FE objects for this element and its pseudo-neighbors.
  UniquePtr<FEGenericBase<OutputShape> > my_fe
    (FEGenericBase<OutputShape>::build(Dim, base_fe_type));
  const FEContinuity cont = my_fe->get_continuity();

  // We don't need to constrain discontinuous elements
  if (cont == DISCONTINUOUS)
    return;
  libmesh_assert (cont == C_ZERO || cont == C_ONE);

  // We'll use element size to generate relative tolerances later
  const Real primary_hmin = elem->hmin();

  UniquePtr<FEGenericBase<OutputShape> > neigh_fe
    (FEGenericBase<OutputShape>::build(Dim, base_fe_type));

  QGauss my_qface(Dim-1, base_fe_type.default_quadrature_order());
  my_fe->attach_quadrature_rule (&my_qface);
  std::vector<Point> neigh_qface;

  const std::vector<Real> & JxW = my_fe->get_JxW();
  const std::vector<Point> & q_point = my_fe->get_xyz();
  const std::vector<std::vector<OutputShape> > & phi = my_fe->get_phi();
  const std::vector<std::vector<OutputShape> > & neigh_phi =
    neigh_fe->get_phi();
  const std::vector<Point> * face_normals = libmesh_nullptr;
  const std::vector<std::vector<OutputGradient> > * dphi = libmesh_nullptr;
  const std::vector<std::vector<OutputGradient> > * neigh_dphi = libmesh_nullptr;
  std::vector<dof_id_type> my_dof_indices, neigh_dof_indices;
  std::vector<unsigned int> my_side_dofs, neigh_side_dofs;

  if (cont != C_ZERO)
    {
      const std::vector<Point> & ref_face_normals =
        my_fe->get_normals();
      face_normals = &ref_face_normals;
      const std::vector<std::vector<OutputGradient> > & ref_dphi =
        my_fe->get_dphi();
      dphi = &ref_dphi;
      const std::vector<std::vector<OutputGradient> > & ref_neigh_dphi =
        neigh_fe->get_dphi();
      neigh_dphi = &ref_neigh_dphi;
    }

  DenseMatrix<Real> Ke;
  DenseVector<Real> Fe;
  std::vector<DenseVector<Real> > Ue;

  // Container to catch the boundary ids that BoundaryInfo hands us.
  std::vector<boundary_id_type> bc_ids;

  // Look at the element faces.  Check to see if we need to
  // build constraints.
  for (unsigned short int s=0; s<elem->n_sides(); s++)
    {
      if (elem->neighbor_ptr(s))
        continue;

      mesh.get_boundary_info().boundary_ids (elem, s, bc_ids);

      for (std::vector<boundary_id_type>::const_iterator id_it=bc_ids.begin(); id_it!=bc_ids.end(); ++id_it)
        {
          const boundary_id_type boundary_id = *id_it;
          const PeriodicBoundaryBase * periodic = boundaries.boundary(boundary_id);
          if (periodic && periodic->is_my_variable(variable_number))
            {
              libmesh_assert(point_locator);

              // Get pointers to the element's neighbor.
              const Elem * neigh = boundaries.neighbor(boundary_id, *point_locator, elem, s);

              if (neigh == libmesh_nullptr)
                libmesh_error_msg("PeriodicBoundaries point locator object returned NULL!");

              // periodic (and possibly h refinement) constraints:
              // constrain dofs shared between
              // this element and ones as coarse
              // as or coarser than this element.
              if (neigh->level() <= elem->level())
                {
                  unsigned int s_neigh =
                    mesh.get_boundary_info().side_with_boundary_id(neigh, periodic->pairedboundary);
                  libmesh_assert_not_equal_to (s_neigh, libMesh::invalid_uint);

#ifdef LIBMESH_ENABLE_AMR
                  // Find the minimum p level; we build the h constraint
                  // matrix with this and then constrain away all higher p
                  // DoFs.
                  libmesh_assert(neigh->active());
                  const unsigned int min_p_level =
                    std::min(elem->p_level(), neigh->p_level());

                  // we may need to make the FE objects reinit with the
                  // minimum shared p_level
                  // FIXME - I hate using const_cast<> and avoiding
                  // accessor functions; there's got to be a
                  // better way to do this!
                  const unsigned int old_elem_level = elem->p_level();
                  if (old_elem_level != min_p_level)
                    (const_cast<Elem *>(elem))->hack_p_level(min_p_level);
                  const unsigned int old_neigh_level = neigh->p_level();
                  if (old_neigh_level != min_p_level)
                    (const_cast<Elem *>(neigh))->hack_p_level(min_p_level);
#endif // #ifdef LIBMESH_ENABLE_AMR

                  // We can do a projection with a single integration,
                  // due to the assumption of nested finite element
                  // subspaces.
                  // FIXME: it might be more efficient to do nodes,
                  // then edges, then side, to reduce the size of the
                  // Cholesky factorization(s)
                  my_fe->reinit(elem, s);

                  dof_map.dof_indices (elem, my_dof_indices,
                                       variable_number);
                  dof_map.dof_indices (neigh, neigh_dof_indices,
                                       variable_number);

                  const unsigned int n_qp = my_qface.n_points();

                  // Translate the quadrature points over to the
                  // neighbor's boundary
                  std::vector<Point> neigh_point(q_point.size());
                  for (std::size_t i=0; i != neigh_point.size(); ++i)
                    neigh_point[i] = periodic->get_corresponding_pos(q_point[i]);

                  FEInterface::inverse_map (Dim, base_fe_type, neigh,
                                            neigh_point, neigh_qface);

                  neigh_fe->reinit(neigh, &neigh_qface);

                  // We're only concerned with DOFs whose values (and/or first
                  // derivatives for C1 elements) are supported on side nodes
                  FEInterface::dofs_on_side(elem, Dim, base_fe_type, s, my_side_dofs);
                  FEInterface::dofs_on_side(neigh, Dim, base_fe_type, s_neigh, neigh_side_dofs);

                  // We're done with functions that examine Elem::p_level(),
                  // so let's unhack those levels
#ifdef LIBMESH_ENABLE_AMR
                  if (elem->p_level() != old_elem_level)
                    (const_cast<Elem *>(elem))->hack_p_level(old_elem_level);
                  if (neigh->p_level() != old_neigh_level)
                    (const_cast<Elem *>(neigh))->hack_p_level(old_neigh_level);
#endif // #ifdef LIBMESH_ENABLE_AMR

                  const unsigned int n_side_dofs =
                    cast_int<unsigned int>
                    (my_side_dofs.size());
                  libmesh_assert_equal_to (n_side_dofs, neigh_side_dofs.size());

                  Ke.resize (n_side_dofs, n_side_dofs);
                  Ue.resize(n_side_dofs);

                  // Form the projection matrix, (inner product of fine basis
                  // functions against fine test functions)
                  for (unsigned int is = 0; is != n_side_dofs; ++is)
                    {
                      const unsigned int i = my_side_dofs[is];
                      for (unsigned int js = 0; js != n_side_dofs; ++js)
                        {
                          const unsigned int j = my_side_dofs[js];
                          for (unsigned int qp = 0; qp != n_qp; ++qp)
                            {
                              Ke(is,js) += JxW[qp] *
                                TensorTools::inner_product(phi[i][qp],
                                                           phi[j][qp]);
                              if (cont != C_ZERO)
                                Ke(is,js) += JxW[qp] *
                                  TensorTools::inner_product((*dphi)[i][qp] *
                                                             (*face_normals)[qp],
                                                             (*dphi)[j][qp] *
                                                             (*face_normals)[qp]);
                            }
                        }
                    }

                  // Form the right hand sides, (inner product of coarse basis
                  // functions against fine test functions)
                  for (unsigned int is = 0; is != n_side_dofs; ++is)
                    {
                      const unsigned int i = neigh_side_dofs[is];
                      Fe.resize (n_side_dofs);
                      for (unsigned int js = 0; js != n_side_dofs; ++js)
                        {
                          const unsigned int j = my_side_dofs[js];
                          for (unsigned int qp = 0; qp != n_qp; ++qp)
                            {
                              Fe(js) += JxW[qp] *
                                TensorTools::inner_product(neigh_phi[i][qp],
                                                           phi[j][qp]);
                              if (cont != C_ZERO)
                                Fe(js) += JxW[qp] *
                                  TensorTools::inner_product((*neigh_dphi)[i][qp] *
                                                             (*face_normals)[qp],
                                                             (*dphi)[j][qp] *
                                                             (*face_normals)[qp]);
                            }
                        }
                      Ke.cholesky_solve(Fe, Ue[is]);
                    }

                  // Make sure we're not adding recursive constraints
                  // due to the redundancy in the way we add periodic
                  // boundary constraints
                  //
                  // In order for this to work while threaded or on
                  // distributed meshes, we need a rigorous way to
                  // avoid recursive constraints.  Here it is:
                  //
                  // For vertex DoFs, if there is a "prior" element
                  // (i.e. a coarser element or an equally refined
                  // element with a lower id) on this boundary which
                  // contains the vertex point, then we will avoid
                  // generating constraints; the prior element (or
                  // something prior to it) may do so.  If we are the
                  // most prior (or "primary") element on this
                  // boundary sharing this point, then we look at the
                  // boundary periodic to us, we find the primary
                  // element there, and if that primary is coarser or
                  // equal-but-lower-id, then our vertex dofs are
                  // constrained in terms of that element.
                  //
                  // For edge DoFs, if there is a coarser element
                  // on this boundary sharing this edge, then we will
                  // avoid generating constraints (we will be
                  // constrained indirectly via AMR constraints
                  // connecting us to the coarser element's DoFs).  If
                  // we are the coarsest element sharing this edge,
                  // then we generate constraints if and only if we
                  // are finer than the coarsest element on the
                  // boundary periodic to us sharing the corresponding
                  // periodic edge, or if we are at equal level but
                  // our edge nodes have higher ids than the periodic
                  // edge nodes (sorted from highest to lowest, then
                  // compared lexicographically)
                  //
                  // For face DoFs, we generate constraints if we are
                  // finer than our periodic neighbor, or if we are at
                  // equal level but our element id is higher than its
                  // element id.
                  //
                  // If the primary neighbor is also the current elem
                  // (a 1-element-thick mesh) then we choose which
                  // vertex dofs to constrain via lexicographic
                  // ordering on point locations

                  // FIXME: This code doesn't yet properly handle
                  // cases where multiple different periodic BCs
                  // intersect.
                  std::set<dof_id_type> my_constrained_dofs;

                  // Container to catch boundary IDs handed back by BoundaryInfo.
                  std::vector<boundary_id_type> new_bc_ids;

                  for (unsigned int n = 0; n != elem->n_nodes(); ++n)
                    {
                      if (!elem->is_node_on_side(n,s))
                        continue;

                      const Node & my_node = elem->node_ref(n);

                      if (elem->is_vertex(n))
                        {
                          // Find all boundary ids that include this
                          // point and have periodic boundary
                          // conditions for this variable
                          std::set<boundary_id_type> point_bcids;

                          for (unsigned int new_s = 0; new_s !=
                                 elem->n_sides(); ++new_s)
                            {
                              if (!elem->is_node_on_side(n,new_s))
                                continue;

                              mesh.get_boundary_info().boundary_ids (elem, s, new_bc_ids);

                              for (std::vector<boundary_id_type>::const_iterator
                                     new_id_it=new_bc_ids.begin(); new_id_it!=new_bc_ids.end(); ++new_id_it)
                                {
                                  const boundary_id_type new_boundary_id = *new_id_it;
                                  const PeriodicBoundaryBase * new_periodic = boundaries.boundary(new_boundary_id);
                                  if (new_periodic && new_periodic->is_my_variable(variable_number))
                                    {
                                      point_bcids.insert(new_boundary_id);
                                    }
                                }
                            }

                          // See if this vertex has point neighbors to
                          // defer to
                          if (primary_boundary_point_neighbor
                              (elem, my_node, mesh.get_boundary_info(), point_bcids)
                              != elem)
                            continue;

                          // Find the complementary boundary id set
                          std::set<boundary_id_type> point_pairedids;
                          for (std::set<boundary_id_type>::const_iterator i =
                                 point_bcids.begin(); i != point_bcids.end(); ++i)
                            {
                              const boundary_id_type new_boundary_id = *i;
                              const PeriodicBoundaryBase * new_periodic = boundaries.boundary(new_boundary_id);
                              point_pairedids.insert(new_periodic->pairedboundary);
                            }

                          // What do we want to constrain against?
                          const Elem * primary_elem = libmesh_nullptr;
                          const Elem * main_neigh = libmesh_nullptr;
                          Point main_pt = my_node,
                            primary_pt = my_node;

                          for (std::set<boundary_id_type>::const_iterator i =
                                 point_bcids.begin(); i != point_bcids.end(); ++i)
                            {
                              // Find the corresponding periodic point and
                              // its primary neighbor
                              const boundary_id_type new_boundary_id = *i;
                              const PeriodicBoundaryBase * new_periodic = boundaries.boundary(new_boundary_id);

                              const Point neigh_pt =
                                new_periodic->get_corresponding_pos(my_node);

                              // If the point is getting constrained
                              // to itself by this PBC then we don't
                              // generate any constraints
                              if (neigh_pt.absolute_fuzzy_equals
                                  (my_node, primary_hmin*TOLERANCE))
                                continue;

                              // Otherwise we'll have a constraint in
                              // one direction or another
                              if (!primary_elem)
                                primary_elem = elem;

                              const Elem * primary_neigh =
                                primary_boundary_point_neighbor(neigh, neigh_pt,
                                                                mesh.get_boundary_info(),
                                                                point_pairedids);

                              libmesh_assert(primary_neigh);

                              if (new_boundary_id == boundary_id)
                                {
                                  main_neigh = primary_neigh;
                                  main_pt = neigh_pt;
                                }

                              // Finer elements will get constrained in
                              // terms of coarser neighbors, not the
                              // other way around
                              if ((primary_neigh->level() > primary_elem->level()) ||

                                  // For equal-level elements, the one with
                                  // higher id gets constrained in terms of
                                  // the one with lower id
                                  (primary_neigh->level() == primary_elem->level() &&
                                   primary_neigh->id() > primary_elem->id()) ||

                                  // On a one-element-thick mesh, we compare
                                  // points to see what side gets constrained
                                  (primary_neigh == primary_elem &&
                                   (neigh_pt > primary_pt)))
                                continue;

                              primary_elem = primary_neigh;
                              primary_pt = neigh_pt;
                            }

                          if (!primary_elem ||
                              primary_elem != main_neigh ||
                              primary_pt != main_pt)
                            continue;
                        }
                      else if (elem->is_edge(n))
                        {
                          // Find which edge we're on
                          unsigned int e=0;
                          for (; e != elem->n_edges(); ++e)
                            {
                              if (elem->is_node_on_edge(n,e))
                                break;
                            }
                          libmesh_assert_less (e, elem->n_edges());

                          // Find the edge end nodes
                          const Node
                            * e1 = libmesh_nullptr,
                            * e2 = libmesh_nullptr;
                          for (unsigned int nn = 0; nn != elem->n_nodes(); ++nn)
                            {
                              if (nn == n)
                                continue;

                              if (elem->is_node_on_edge(nn, e))
                                {
                                  if (e1 == libmesh_nullptr)
                                    {
                                      e1 = elem->node_ptr(nn);
                                    }
                                  else
                                    {
                                      e2 = elem->node_ptr(nn);
                                      break;
                                    }
                                }
                            }
                          libmesh_assert (e1 && e2);

                          // Find all boundary ids that include this
                          // edge and have periodic boundary
                          // conditions for this variable
                          std::set<boundary_id_type> edge_bcids;

                          for (unsigned int new_s = 0; new_s !=
                                 elem->n_sides(); ++new_s)
                            {
                              if (!elem->is_node_on_side(n,new_s))
                                continue;

                              // We're reusing the new_bc_ids vector created outside the loop over nodes.
                              mesh.get_boundary_info().boundary_ids (elem, s, new_bc_ids);

                              for (std::vector<boundary_id_type>::const_iterator
                                     new_id_it=new_bc_ids.begin(); new_id_it!=new_bc_ids.end(); ++new_id_it)
                                {
                                  const boundary_id_type new_boundary_id = *new_id_it;
                                  const PeriodicBoundaryBase * new_periodic = boundaries.boundary(new_boundary_id);
                                  if (new_periodic && new_periodic->is_my_variable(variable_number))
                                    {
                                      edge_bcids.insert(new_boundary_id);
                                    }
                                }
                            }


                          // See if this edge has neighbors to defer to
                          if (primary_boundary_edge_neighbor
                              (elem, *e1, *e2, mesh.get_boundary_info(), edge_bcids)
                              != elem)
                            continue;

                          // Find the complementary boundary id set
                          std::set<boundary_id_type> edge_pairedids;
                          for (std::set<boundary_id_type>::const_iterator i =
                                 edge_bcids.begin(); i != edge_bcids.end(); ++i)
                            {
                              const boundary_id_type new_boundary_id = *i;
                              const PeriodicBoundaryBase * new_periodic = boundaries.boundary(new_boundary_id);
                              edge_pairedids.insert(new_periodic->pairedboundary);
                            }

                          // What do we want to constrain against?
                          const Elem * primary_elem = libmesh_nullptr;
                          const Elem * main_neigh = libmesh_nullptr;
                          Point main_pt1 = *e1,
                            main_pt2 = *e2,
                            primary_pt1 = *e1,
                            primary_pt2 = *e2;

                          for (std::set<boundary_id_type>::const_iterator i =
                                 edge_bcids.begin(); i != edge_bcids.end(); ++i)
                            {
                              // Find the corresponding periodic edge and
                              // its primary neighbor
                              const boundary_id_type new_boundary_id = *i;
                              const PeriodicBoundaryBase * new_periodic = boundaries.boundary(new_boundary_id);

                              Point neigh_pt1 = new_periodic->get_corresponding_pos(*e1),
                                neigh_pt2 = new_periodic->get_corresponding_pos(*e2);

                              // If the edge is getting constrained
                              // to itself by this PBC then we don't
                              // generate any constraints
                              if (neigh_pt1.absolute_fuzzy_equals
                                  (*e1, primary_hmin*TOLERANCE) &&
                                  neigh_pt2.absolute_fuzzy_equals
                                  (*e2, primary_hmin*TOLERANCE))
                                continue;

                              // Otherwise we'll have a constraint in
                              // one direction or another
                              if (!primary_elem)
                                primary_elem = elem;

                              const Elem * primary_neigh = primary_boundary_edge_neighbor
                                (neigh, neigh_pt1, neigh_pt2,
                                 mesh.get_boundary_info(), edge_pairedids);

                              libmesh_assert(primary_neigh);

                              if (new_boundary_id == boundary_id)
                                {
                                  main_neigh = primary_neigh;
                                  main_pt1 = neigh_pt1;
                                  main_pt2 = neigh_pt2;
                                }

                              // If we have a one-element thick mesh,
                              // we'll need to sort our points to get a
                              // consistent ordering rule
                              //
                              // Use >= in this test to make sure that,
                              // for angular constraints, no node gets
                              // constrained to itself.
                              if (primary_neigh == primary_elem)
                                {
                                  if (primary_pt1 > primary_pt2)
                                    std::swap(primary_pt1, primary_pt2);
                                  if (neigh_pt1 > neigh_pt2)
                                    std::swap(neigh_pt1, neigh_pt2);

                                  if (neigh_pt2 >= primary_pt2)
                                    continue;
                                }

                              // Otherwise:
                              // Finer elements will get constrained in
                              // terms of coarser ones, not the other way
                              // around
                              if ((primary_neigh->level() > primary_elem->level()) ||

                                  // For equal-level elements, the one with
                                  // higher id gets constrained in terms of
                                  // the one with lower id
                                  (primary_neigh->level() == primary_elem->level() &&
                                   primary_neigh->id() > primary_elem->id()))
                                continue;

                              primary_elem = primary_neigh;
                              primary_pt1 = neigh_pt1;
                              primary_pt2 = neigh_pt2;
                            }

                          if (!primary_elem ||
                              primary_elem != main_neigh ||
                              primary_pt1 != main_pt1 ||
                              primary_pt2 != main_pt2)
                            continue;
                        }
                      else if (elem->is_face(n))
                        {
                          // If we have a one-element thick mesh,
                          // use the ordering of the face node and its
                          // periodic counterpart to determine what
                          // gets constrained
                          if (neigh == elem)
                            {
                              const Point neigh_pt =
                                periodic->get_corresponding_pos(my_node);
                              if (neigh_pt > my_node)
                                continue;
                            }

                          // Otherwise:
                          // Finer elements will get constrained in
                          // terms of coarser ones, not the other way
                          // around
                          if ((neigh->level() > elem->level()) ||

                              // For equal-level elements, the one with
                              // higher id gets constrained in terms of
                              // the one with lower id
                              (neigh->level() == elem->level() &&
                               neigh->id() > elem->id()))
                            continue;
                        }

                      // If we made it here without hitting a continue
                      // statement, then we're at a node whose dofs
                      // should be constrained by this element's
                      // calculations.
                      const unsigned int n_comp =
                        my_node.n_comp(sys_number, variable_number);

                      for (unsigned int i=0; i != n_comp; ++i)
                        my_constrained_dofs.insert
                          (my_node.dof_number
                           (sys_number, variable_number, i));
                    }

                  // FIXME: old code for disambiguating periodic BCs:
                  // this is not threadsafe nor safe to run on a
                  // non-serialized mesh.
                  /*
                    std::vector<bool> recursive_constraint(n_side_dofs, false);

                    for (unsigned int is = 0; is != n_side_dofs; ++is)
                    {
                    const unsigned int i = neigh_side_dofs[is];
                    const dof_id_type their_dof_g = neigh_dof_indices[i];
                    libmesh_assert_not_equal_to (their_dof_g, DofObject::invalid_id);

                    {
                    Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);

                    if (!dof_map.is_constrained_dof(their_dof_g))
                    continue;
                    }

                    DofConstraintRow & their_constraint_row =
                    constraints[their_dof_g].first;

                    for (unsigned int js = 0; js != n_side_dofs; ++js)
                    {
                    const unsigned int j = my_side_dofs[js];
                    const dof_id_type my_dof_g = my_dof_indices[j];
                    libmesh_assert_not_equal_to (my_dof_g, DofObject::invalid_id);

                    if (their_constraint_row.count(my_dof_g))
                    recursive_constraint[js] = true;
                    }
                    }
                  */

                  for (unsigned int js = 0; js != n_side_dofs; ++js)
                    {
                      // FIXME: old code path
                      // if (recursive_constraint[js])
                      //  continue;

                      const unsigned int j = my_side_dofs[js];
                      const dof_id_type my_dof_g = my_dof_indices[j];
                      libmesh_assert_not_equal_to (my_dof_g, DofObject::invalid_id);

                      // FIXME: new code path
                      if (!my_constrained_dofs.count(my_dof_g))
                        continue;

                      DofConstraintRow * constraint_row;

                      // we may be running constraint methods concurretly
                      // on multiple threads, so we need a lock to
                      // ensure that this constraint is "ours"
                      {
                        Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);

                        if (dof_map.is_constrained_dof(my_dof_g))
                          continue;

                        constraint_row = &(constraints[my_dof_g]);
                        libmesh_assert(constraint_row->empty());
                      }

                      for (unsigned int is = 0; is != n_side_dofs; ++is)
                        {
                          const unsigned int i = neigh_side_dofs[is];
                          const dof_id_type their_dof_g = neigh_dof_indices[i];
                          libmesh_assert_not_equal_to (their_dof_g, DofObject::invalid_id);

                          // Periodic constraints should never be
                          // self-constraints
                          // libmesh_assert_not_equal_to (their_dof_g, my_dof_g);

                          const Real their_dof_value = Ue[is](js);

                          if (their_dof_g == my_dof_g)
                            {
                              libmesh_assert_less (std::abs(their_dof_value-1.), 1.e-5);
                              for (unsigned int k = 0; k != n_side_dofs; ++k)
                                libmesh_assert(k == is || std::abs(Ue[k](js)) < 1.e-5);
                              continue;
                            }

                          if (std::abs(their_dof_value) < 10*TOLERANCE)
                            continue;

                          constraint_row->insert(std::make_pair(their_dof_g,
                                                                their_dof_value));
                        }
                    }
                }
              // p refinement constraints:
              // constrain dofs shared between
              // active elements and neighbors with
              // lower polynomial degrees
#ifdef LIBMESH_ENABLE_AMR
              const unsigned int min_p_level =
                neigh->min_p_level_by_neighbor(elem, elem->p_level());
              if (min_p_level < elem->p_level())
                {
                  // Adaptive p refinement of non-hierarchic bases will
                  // require more coding
                  libmesh_assert(my_fe->is_hierarchic());
                  dof_map.constrain_p_dofs(variable_number, elem,
                                           s, min_p_level);
                }
#endif // #ifdef LIBMESH_ENABLE_AMR
            }
        }
    }
}

#endif // LIBMESH_ENABLE_PERIODIC

// ------------------------------------------------------------
// Explicit instantiations
template class FEGenericBase<Real>;
template class FEGenericBase<RealGradient>;

} // namespace libMesh
