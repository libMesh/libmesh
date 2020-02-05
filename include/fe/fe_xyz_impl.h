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

#ifndef LIBMESH_FE_XYZ_IMPL_H
#define LIBMESH_FE_XYZ_IMPL_H

// Local includes
#include "libmesh/libmesh_logging.h"
#include "libmesh/fe.h"
#include "libmesh/elem.h"
#include "libmesh/fe_interface.h"
#include "libmesh/enum_to_string.h"
#include "libmesh/int_range.h"

namespace libMesh
{

// ------------------------------------------------------------
// XYZ-specific implementations

// Anonymous namespace for local helper functions
namespace {

template <typename RealType>
void xyz_nodal_soln(const ElemTempl<RealType> * elem,
                    const Order order,
                    const std::vector<Number> & elem_soln,
                    std::vector<Number> &       nodal_soln,
                    unsigned Dim)
{
  const unsigned int n_nodes = elem->n_nodes();

  const ElemType elem_type = elem->type();

  nodal_soln.resize(n_nodes);

  const Order totalorder = static_cast<Order>(order + elem->p_level());

  switch (totalorder)
    {
      // Constant shape functions
    case CONSTANT:
      {
        libmesh_assert_equal_to (elem_soln.size(), 1);

        const Number val = elem_soln[0];

        for (unsigned int n=0; n<n_nodes; n++)
          nodal_soln[n] = val;

        return;
      }


      // For other orders do interpolation at the nodes
      // explicitly.
    default:
      {
        // FEType object to be passed to various FEInterface functions below.
        FEType fe_type(totalorder, XYZ);

        const unsigned int n_sf =
          // FE<Dim,T>::n_shape_functions(elem_type, totalorder);
          FEInterface::n_shape_functions(Dim, fe_type, elem_type);

        for (unsigned int n=0; n<n_nodes; n++)
          {
            libmesh_assert_equal_to (elem_soln.size(), n_sf);

            // Zero before summation
            nodal_soln[n] = 0;

            // u_i = Sum (alpha_i phi_i)
            for (unsigned int i=0; i<n_sf; i++)
              nodal_soln[n] += elem_soln[i] *
                // FE<Dim,T>::shape(elem, order, i, elem->point(n));
                FEInterface::shape(Dim, fe_type, elem, i, elem->point(n));
          }

        return;
      } // default
    } // switch
} // xyz_nodal_soln()





unsigned int xyz_n_dofs(const ElemType t, const Order o)
{
  switch (o)
    {

      // constant shape functions
      // no matter what shape there is only one DOF.
    case CONSTANT:
      return (t != INVALID_ELEM) ? 1 : 0;


      // Discontinuous linear shape functions
      // expressed in the XYZ monomials.
    case FIRST:
      {
        switch (t)
          {
          case NODEELEM:
            return 1;

          case EDGE2:
          case EDGE3:
          case EDGE4:
            return 2;

          case TRI3:
          case TRISHELL3:
          case TRI6:
          case QUAD4:
          case QUADSHELL4:
          case QUAD8:
          case QUADSHELL8:
          case QUAD9:
            return 3;

          case TET4:
          case TET10:
          case HEX8:
          case HEX20:
          case HEX27:
          case PRISM6:
          case PRISM15:
          case PRISM18:
          case PYRAMID5:
          case PYRAMID13:
          case PYRAMID14:
            return 4;

          case INVALID_ELEM:
            return 0;

          default:
            libmesh_error_msg("ERROR: Bad ElemType = " << Utility::enum_to_string(t) << " for " << Utility::enum_to_string(o) << " order approximation!");
          }
      }


      // Discontinuous quadratic shape functions
      // expressed in the XYZ monomials.
    case SECOND:
      {
        switch (t)
          {
          case NODEELEM:
            return 1;

          case EDGE2:
          case EDGE3:
          case EDGE4:
            return 3;

          case TRI3:
          case TRISHELL3:
          case TRI6:
          case QUAD4:
          case QUADSHELL4:
          case QUAD8:
          case QUADSHELL8:
          case QUAD9:
            return 6;

          case TET4:
          case TET10:
          case HEX8:
          case HEX20:
          case HEX27:
          case PRISM6:
          case PRISM15:
          case PRISM18:
          case PYRAMID5:
          case PYRAMID13:
          case PYRAMID14:
            return 10;

          case INVALID_ELEM:
            return 0;

          default:
            libmesh_error_msg("ERROR: Bad ElemType = " << Utility::enum_to_string(t) << " for " << Utility::enum_to_string(o) << " order approximation!");
          }
      }


      // Discontinuous cubic shape functions
      // expressed in the XYZ monomials.
    case THIRD:
      {
        switch (t)
          {
          case NODEELEM:
            return 1;

          case EDGE2:
          case EDGE3:
          case EDGE4:
            return 4;

          case TRI3:
          case TRISHELL3:
          case TRI6:
          case QUAD4:
          case QUADSHELL4:
          case QUAD8:
          case QUADSHELL8:
          case QUAD9:
            return 10;

          case TET4:
          case TET10:
          case HEX8:
          case HEX20:
          case HEX27:
          case PRISM6:
          case PRISM15:
          case PRISM18:
          case PYRAMID5:
          case PYRAMID13:
          case PYRAMID14:
            return 20;

          case INVALID_ELEM:
            return 0;

          default:
            libmesh_error_msg("ERROR: Bad ElemType = " << Utility::enum_to_string(t) << " for " << Utility::enum_to_string(o) << " order approximation!");
          }
      }


      // Discontinuous quartic shape functions
      // expressed in the XYZ monomials.
    case FOURTH:
      {
        switch (t)
          {
          case NODEELEM:
            return 1;

          case EDGE2:
          case EDGE3:
            return 5;

          case TRI3:
          case TRISHELL3:
          case TRI6:
          case QUAD4:
          case QUADSHELL4:
          case QUAD8:
          case QUADSHELL8:
          case QUAD9:
            return 15;

          case TET4:
          case TET10:
          case HEX8:
          case HEX20:
          case HEX27:
          case PRISM6:
          case PRISM15:
          case PRISM18:
          case PYRAMID5:
          case PYRAMID13:
          case PYRAMID14:
            return 35;

          case INVALID_ELEM:
            return 0;

          default:
            libmesh_error_msg("ERROR: Bad ElemType = " << Utility::enum_to_string(t) << " for " << Utility::enum_to_string(o) << " order approximation!");
          }
      }


    default:
      {
        const unsigned int order = static_cast<unsigned int>(o);
        switch (t)
          {
          case NODEELEM:
            return 1;

          case EDGE2:
          case EDGE3:
            return (order+1);

          case TRI3:
          case TRISHELL3:
          case TRI6:
          case QUAD4:
          case QUADSHELL4:
          case QUAD8:
          case QUADSHELL8:
          case QUAD9:
            return (order+1)*(order+2)/2;

          case TET4:
          case TET10:
          case HEX8:
          case HEX20:
          case HEX27:
          case PRISM6:
          case PRISM15:
          case PRISM18:
          case PYRAMID5:
          case PYRAMID13:
          case PYRAMID14:
            return (order+1)*(order+2)*(order+3)/6;

          case INVALID_ELEM:
            return 0;

          default:
            libmesh_error_msg("ERROR: Bad ElemType = " << Utility::enum_to_string(t) << " for " << Utility::enum_to_string(o) << " order approximation!");
          }
      }
    }
}




unsigned int xyz_n_dofs_per_elem(const ElemType t,
                                 const Order o)
{
  switch (o)
    {
      // constant shape functions always have 1 DOF per element
    case CONSTANT:
      return (t != INVALID_ELEM) ? 1 : 0;


      // Discontinuous linear shape functions
      // expressed in the XYZ monomials.
    case FIRST:
      {
        switch (t)
          {
          case NODEELEM:
            return 1;

            // 1D linears have 2 DOFs per element
          case EDGE2:
          case EDGE3:
          case EDGE4:
            return 2;

            // 2D linears have 3 DOFs per element
          case TRI3:
          case TRISHELL3:
          case TRI6:
          case QUAD4:
          case QUADSHELL4:
          case QUAD8:
          case QUADSHELL8:
          case QUAD9:
            return 3;

            // 3D linears have 4 DOFs per element
          case TET4:
          case TET10:
          case HEX8:
          case HEX20:
          case HEX27:
          case PRISM6:
          case PRISM15:
          case PRISM18:
          case PYRAMID5:
          case PYRAMID13:
          case PYRAMID14:
            return 4;

          case INVALID_ELEM:
            return 0;

          default:
            libmesh_error_msg("ERROR: Bad ElemType = " << Utility::enum_to_string(t) << " for " << Utility::enum_to_string(o) << " order approximation!");
          }
      }


      // Discontinuous quadratic shape functions
      // expressed in the XYZ monomials.
    case SECOND:
      {
        switch (t)
          {
          case NODEELEM:
            return 1;

            // 1D quadratics have 3 DOFs per element
          case EDGE2:
          case EDGE3:
          case EDGE4:
            return 3;

            // 2D quadratics have 6 DOFs per element
          case TRI3:
          case TRISHELL3:
          case TRI6:
          case QUAD4:
          case QUADSHELL4:
          case QUAD8:
          case QUADSHELL8:
          case QUAD9:
            return 6;

            // 3D quadratics have 10 DOFs per element
          case TET4:
          case TET10:
          case HEX8:
          case HEX20:
          case HEX27:
          case PRISM6:
          case PRISM15:
          case PRISM18:
          case PYRAMID5:
          case PYRAMID13:
          case PYRAMID14:
            return 10;

          case INVALID_ELEM:
            return 0;

          default:
            libmesh_error_msg("ERROR: Bad ElemType = " << Utility::enum_to_string(t) << " for " << Utility::enum_to_string(o) << " order approximation!");
          }
      }


      // Discontinuous cubic shape functions
      // expressed in the XYZ monomials.
    case THIRD:
      {
        switch (t)
          {
          case NODEELEM:
            return 1;

          case EDGE2:
          case EDGE3:
          case EDGE4:
            return 4;

          case TRI3:
          case TRISHELL3:
          case TRI6:
          case QUAD4:
          case QUADSHELL4:
          case QUAD8:
          case QUADSHELL8:
          case QUAD9:
            return 10;

          case TET4:
          case TET10:
          case HEX8:
          case HEX20:
          case HEX27:
          case PRISM6:
          case PRISM15:
          case PRISM18:
          case PYRAMID5:
          case PYRAMID13:
          case PYRAMID14:
            return 20;

          case INVALID_ELEM:
            return 0;

          default:
            libmesh_error_msg("ERROR: Bad ElemType = " << Utility::enum_to_string(t) << " for " << Utility::enum_to_string(o) << " order approximation!");
          }
      }


      // Discontinuous quartic shape functions
      // expressed in the XYZ monomials.
    case FOURTH:
      {
        switch (t)
          {
          case NODEELEM:
            return 1;

          case EDGE2:
          case EDGE3:
          case EDGE4:
            return 5;

          case TRI3:
          case TRISHELL3:
          case TRI6:
          case QUAD4:
          case QUADSHELL4:
          case QUAD8:
          case QUADSHELL8:
          case QUAD9:
            return 15;

          case TET4:
          case TET10:
          case HEX8:
          case HEX20:
          case HEX27:
          case PRISM6:
          case PRISM15:
          case PRISM18:
          case PYRAMID5:
          case PYRAMID13:
          case PYRAMID14:
            return 35;

          case INVALID_ELEM:
            return 0;

          default:
            libmesh_error_msg("ERROR: Bad ElemType = " << Utility::enum_to_string(t) << " for " << Utility::enum_to_string(o) << " order approximation!");
          }
      }

    default:
      {
        const unsigned int order = static_cast<unsigned int>(o);
        switch (t)
          {
          case NODEELEM:
            return 1;

          case EDGE2:
          case EDGE3:
            return (order+1);

          case TRI3:
          case TRISHELL3:
          case TRI6:
          case QUAD4:
          case QUADSHELL4:
          case QUAD8:
          case QUADSHELL8:
          case QUAD9:
            return (order+1)*(order+2)/2;

          case TET4:
          case TET10:
          case HEX8:
          case HEX20:
          case HEX27:
          case PRISM6:
          case PRISM15:
          case PRISM18:
          case PYRAMID5:
          case PYRAMID13:
          case PYRAMID14:
            return (order+1)*(order+2)*(order+3)/6;

          case INVALID_ELEM:
            return 0;

          default:
            libmesh_error_msg("ERROR: Bad ElemType = " << Utility::enum_to_string(t) << " for " << Utility::enum_to_string(o) << " order approximation!");
          }
      }
      return 0;
    }
}


} // anonymous namespace


template <unsigned int Dim, typename RealType>
void FEXYZ<Dim,RealType>::init_shape_functions(const std::vector<Point> & qp,
                                               const Elem * libmesh_dbg_var(elem))
{
  libmesh_assert(elem);
  this->calculations_started = true;

  // If the user forgot to request anything, we'll be safe and
  // calculate everything:
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  if (!this->calculate_phi && !this->calculate_dphi && !this->calculate_d2phi)
    this->calculate_phi = this->calculate_dphi = this->calculate_d2phi = true;
#else
  if (!this->calculate_phi && !this->calculate_dphi)
    this->calculate_phi = this->calculate_dphi = true;
#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES

  // Start logging the shape function initialization
  LOG_SCOPE("init_shape_functions()", "FE");

  // The number of quadrature points.
  const std::size_t n_qp = qp.size();

  // Number of shape functions in the finite element approximation
  // space.
  const unsigned int n_approx_shape_functions =
    this->n_shape_functions(this->get_type(),
                            this->get_order());

  // resize the vectors to hold current data
  // Phi are the shape functions used for the FE approximation
  {
    // (note: GCC 3.4.0 requires the use of this-> here)
    if (this->calculate_phi)
      this->phi.resize     (n_approx_shape_functions);
    if (this->calculate_dphi)
      {
        this->dphi.resize    (n_approx_shape_functions);
        this->dphidx.resize  (n_approx_shape_functions);
        this->dphidy.resize  (n_approx_shape_functions);
        this->dphidz.resize  (n_approx_shape_functions);
      }

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
    if (this->calculate_d2phi)
      {
        this->d2phi.resize     (n_approx_shape_functions);
        this->d2phidx2.resize  (n_approx_shape_functions);
        this->d2phidxdy.resize (n_approx_shape_functions);
        this->d2phidxdz.resize (n_approx_shape_functions);
        this->d2phidy2.resize  (n_approx_shape_functions);
        this->d2phidydz.resize (n_approx_shape_functions);
        this->d2phidz2.resize  (n_approx_shape_functions);
        this->d2phidxi2.resize (n_approx_shape_functions);
      }
#endif // ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

    for (unsigned int i=0; i<n_approx_shape_functions; i++)
      {
        if (this->calculate_phi)
          this->phi[i].resize           (n_qp);
        if (this->calculate_dphi)
          {
            this->dphi[i].resize        (n_qp);
            this->dphidx[i].resize      (n_qp);
            this->dphidy[i].resize      (n_qp);
            this->dphidz[i].resize      (n_qp);
          }
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
        if (this->calculate_d2phi)
          {
            this->d2phi[i].resize       (n_qp);
            this->d2phidx2[i].resize    (n_qp);
            this->d2phidxdy[i].resize   (n_qp);
            this->d2phidxdz[i].resize   (n_qp);
            this->d2phidy2[i].resize    (n_qp);
            this->d2phidydz[i].resize   (n_qp);
            this->d2phidz2[i].resize    (n_qp);
          }
#endif // ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
      }
  }



#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
  //------------------------------------------------------------
  // Initialize the data fields, which should only be used for infinite
  // elements, to some sensible values, so that using a FE with the
  // variational formulation of an InfFE, correct element matrices are
  // returned

  {
    this->weight.resize  (n_qp);
    this->dweight.resize (n_qp);
    this->dphase.resize  (n_qp);

    for (unsigned int p=0; p<n_qp; p++)
      {
        this->weight[p] = 1.;
        this->dweight[p].zero();
        this->dphase[p].zero();
      }

  }
#endif // ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
}




template <unsigned int Dim, typename RealType>
void FEXYZ<Dim,RealType>::compute_shape_functions (const Elem * elem,
                                                   const std::vector<Point> &)
{
  libmesh_assert(elem);

  //-------------------------------------------------------------------------
  // Compute the shape function values (and derivatives)
  // at the Quadrature points.  Note that the actual values
  // have already been computed via init_shape_functions

  // Start logging the shape function computation
  LOG_SCOPE("compute_shape_functions()", "FE");

  const std::vector<Point> & xyz_qp = this->get_xyz();

  // Compute the value of the derivative shape function i at quadrature point p
  switch (this->dim)
    {

    case 1:
      {
        if (this->calculate_phi)
          for (auto i : index_range(this->phi))
            for (auto p : index_range(this->phi[i]))
              this->phi[i][p] = FEShim<Dim,XYZ,RealType>::shape (elem, this->fe_type.order, i, xyz_qp[p]);

        if (this->calculate_dphi)
          for (auto i : index_range(this->dphi))
            for (auto p : index_range(this->dphi[i]))
              {
                this->dphi[i][p](0) =
                  this->dphidx[i][p] = FEShim<Dim,XYZ,RealType>::shape_deriv (elem, this->fe_type.order, i, 0, xyz_qp[p]);

                this->dphi[i][p](1) = this->dphidy[i][p] = 0.;
                this->dphi[i][p](2) = this->dphidz[i][p] = 0.;
              }
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
        if (this->calculate_d2phi)
          for (auto i : index_range(this->d2phi))
            for (auto p : index_range(this->d2phi[i]))
              {
                this->d2phi[i][p](0,0) =
                  this->d2phidx2[i][p] = FEShim<Dim,XYZ,RealType>::shape_second_deriv (elem, this->fe_type.order, i, 0, xyz_qp[p]);

#if LIBMESH_DIM>1
                this->d2phi[i][p](0,1) = this->d2phidxdy[i][p] =
                  this->d2phi[i][p](1,0) = 0.;
                this->d2phi[i][p](1,1) = this->d2phidy2[i][p] = 0.;
#if LIBMESH_DIM>2
                this->d2phi[i][p](0,2) = this->d2phidxdz[i][p] =
                  this->d2phi[i][p](2,0) = 0.;
                this->d2phi[i][p](1,2) = this->d2phidydz[i][p] =
                  this->d2phi[i][p](2,1) = 0.;
                this->d2phi[i][p](2,2) = this->d2phidz2[i][p] = 0.;
#endif
#endif
              }
#endif // ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

        // All done
        break;
      }

    case 2:
      {
        if (this->calculate_phi)
          for (auto i : index_range(this->phi))
            for (auto p : index_range(this->phi[i]))
              this->phi[i][p] = FEShim<Dim,XYZ,RealType>::shape (elem, this->fe_type.order, i, xyz_qp[p]);

        if (this->calculate_dphi)
          for (auto i : index_range(this->dphi))
            for (auto p : index_range(this->dphi[i]))
              {
                this->dphi[i][p](0) =
                  this->dphidx[i][p] = FEShim<Dim,XYZ,RealType>::shape_deriv (elem, this->fe_type.order, i, 0, xyz_qp[p]);

                this->dphi[i][p](1) =
                  this->dphidy[i][p] = FEShim<Dim,XYZ,RealType>::shape_deriv (elem, this->fe_type.order, i, 1, xyz_qp[p]);

#if LIBMESH_DIM == 3
                this->dphi[i][p](2) = // can only assign to the Z component if LIBMESH_DIM==3
#endif
                  this->dphidz[i][p] = 0.;
              }
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
        if (this->calculate_d2phi)
          for (auto i : index_range(this->d2phi))
            for (auto p : index_range(this->d2phi[i]))
              {
                this->d2phi[i][p](0,0) =
                  this->d2phidx2[i][p] = FEShim<Dim,XYZ,RealType>::shape_second_deriv (elem, this->fe_type.order, i, 0, xyz_qp[p]);

                this->d2phi[i][p](0,1) = this->d2phidxdy[i][p] =
                  this->d2phi[i][p](1,0) = FEShim<Dim,XYZ,RealType>::shape_second_deriv (elem, this->fe_type.order, i, 1, xyz_qp[p]);
                this->d2phi[i][p](1,1) =
                  this->d2phidy2[i][p] = FEShim<Dim,XYZ,RealType>::shape_second_deriv (elem, this->fe_type.order, i, 2, xyz_qp[p]);
#if LIBMESH_DIM>2
                this->d2phi[i][p](0,2) = this->d2phidxdz[i][p] =
                  this->d2phi[i][p](2,0) = 0.;
                this->d2phi[i][p](1,2) = this->d2phidydz[i][p] =
                  this->d2phi[i][p](2,1) = 0.;
                this->d2phi[i][p](2,2) = this->d2phidz2[i][p] = 0.;
#endif
              }
#endif // ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

        // All done
        break;
      }

    case 3:
      {
        if (this->calculate_phi)
          for (auto i : index_range(this->phi))
            for (auto p : index_range(this->phi[i]))
              this->phi[i][p] = FEShim<Dim,XYZ,RealType>::shape (elem, this->fe_type.order, i, xyz_qp[p]);

        if (this->calculate_dphi)
          for (auto i : index_range(this->dphi))
            for (auto p : index_range(this->dphi[i]))
              {
                this->dphi[i][p](0) =
                  this->dphidx[i][p] = FEShim<Dim,XYZ,RealType>::shape_deriv (elem, this->fe_type.order, i, 0, xyz_qp[p]);

                this->dphi[i][p](1) =
                  this->dphidy[i][p] = FEShim<Dim,XYZ,RealType>::shape_deriv (elem, this->fe_type.order, i, 1, xyz_qp[p]);

                this->dphi[i][p](2) =
                  this->dphidz[i][p] = FEShim<Dim,XYZ,RealType>::shape_deriv (elem, this->fe_type.order, i, 2, xyz_qp[p]);
              }
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
        if (this->calculate_d2phi)
          for (auto i : index_range(this->d2phi))
            for (auto p : index_range(this->d2phi[i]))
              {
                this->d2phi[i][p](0,0) =
                  this->d2phidx2[i][p] = FEShim<Dim,XYZ,RealType>::shape_second_deriv (elem, this->fe_type.order, i, 0, xyz_qp[p]);

                this->d2phi[i][p](0,1) = this->d2phidxdy[i][p] =
                  this->d2phi[i][p](1,0) = FEShim<Dim,XYZ,RealType>::shape_second_deriv (elem, this->fe_type.order, i, 1, xyz_qp[p]);
                this->d2phi[i][p](1,1) =
                  this->d2phidy2[i][p] = FEShim<Dim,XYZ,RealType>::shape_second_deriv (elem, this->fe_type.order, i, 2, xyz_qp[p]);
                this->d2phi[i][p](0,2) = this->d2phidxdz[i][p] =
                  this->d2phi[i][p](2,0) = FEShim<Dim,XYZ,RealType>::shape_second_deriv (elem, this->fe_type.order, i, 3, xyz_qp[p]);
                this->d2phi[i][p](1,2) = this->d2phidydz[i][p] =
                  this->d2phi[i][p](2,1) = FEShim<Dim,XYZ,RealType>::shape_second_deriv (elem, this->fe_type.order, i, 4, xyz_qp[p]);
                this->d2phi[i][p](2,2) = this->d2phidz2[i][p] = FEShim<Dim,XYZ,RealType>::shape_second_deriv (elem, this->fe_type.order, i, 5, xyz_qp[p]);
              }
#endif // ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

        // All done
        break;
      }

    default:
      libmesh_error_msg("ERROR: Invalid dimension " << this->dim);
    }
}




// Do full-specialization for every dimension, instead
// of explicit instantiation at the end of this file.
// This could be macro-ified so that it fits on one line...
template <typename RealType>
void FEShim<0,XYZ,RealType>::nodal_soln(const Elem * elem,
                           const Order order,
                           const std::vector<Number> & elem_soln,
                           std::vector<Number> & nodal_soln)
{ xyz_nodal_soln(elem, order, elem_soln, nodal_soln, /*Dim=*/0); }

template <typename RealType>
void FEShim<1,XYZ,RealType>::nodal_soln(const Elem * elem,
                           const Order order,
                           const std::vector<Number> & elem_soln,
                           std::vector<Number> & nodal_soln)
{ xyz_nodal_soln(elem, order, elem_soln, nodal_soln, /*Dim=*/1); }

template <typename RealType>
void FEShim<2,XYZ,RealType>::nodal_soln(const Elem * elem,
                           const Order order,
                           const std::vector<Number> & elem_soln,
                           std::vector<Number> & nodal_soln)
{ xyz_nodal_soln(elem, order, elem_soln, nodal_soln, /*Dim=*/2); }

template <typename RealType>
void FEShim<3,XYZ,RealType>::nodal_soln(const Elem * elem,
                           const Order order,
                           const std::vector<Number> & elem_soln,
                           std::vector<Number> & nodal_soln)
{ xyz_nodal_soln(elem, order, elem_soln, nodal_soln, /*Dim=*/3); }



// Full specialization of n_dofs() function for every dimension
template <typename RealType> unsigned int FEShim<0,XYZ,RealType>::n_dofs(const ElemType t, const Order o) { return xyz_n_dofs(t, o); }
template <typename RealType> unsigned int FEShim<1,XYZ,RealType>::n_dofs(const ElemType t, const Order o) { return xyz_n_dofs(t, o); }
template <typename RealType> unsigned int FEShim<2,XYZ,RealType>::n_dofs(const ElemType t, const Order o) { return xyz_n_dofs(t, o); }
template <typename RealType> unsigned int FEShim<3,XYZ,RealType>::n_dofs(const ElemType t, const Order o) { return xyz_n_dofs(t, o); }

// Full specialization of n_dofs_at_node() function for every dimension.
// XYZ FEMs have no dofs at nodes, only element dofs.
template <typename RealType> unsigned int FEShim<0,XYZ,RealType>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }
template <typename RealType> unsigned int FEShim<1,XYZ,RealType>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }
template <typename RealType> unsigned int FEShim<2,XYZ,RealType>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }
template <typename RealType> unsigned int FEShim<3,XYZ,RealType>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }

// Full specialization of n_dofs_per_elem() function for every dimension.
template <typename RealType> unsigned int FEShim<0,XYZ,RealType>::n_dofs_per_elem(const ElemType t, const Order o) { return xyz_n_dofs_per_elem(t, o); }
template <typename RealType> unsigned int FEShim<1,XYZ,RealType>::n_dofs_per_elem(const ElemType t, const Order o) { return xyz_n_dofs_per_elem(t, o); }
template <typename RealType> unsigned int FEShim<2,XYZ,RealType>::n_dofs_per_elem(const ElemType t, const Order o) { return xyz_n_dofs_per_elem(t, o); }
template <typename RealType> unsigned int FEShim<3,XYZ,RealType>::n_dofs_per_elem(const ElemType t, const Order o) { return xyz_n_dofs_per_elem(t, o); }

// Full specialization of get_continuity() function for every dimension.
template <typename RealType> FEContinuity FEShim<0,XYZ,RealType>::get_continuity() { return DISCONTINUOUS; }
template <typename RealType> FEContinuity FEShim<1,XYZ,RealType>::get_continuity() { return DISCONTINUOUS; }
template <typename RealType> FEContinuity FEShim<2,XYZ,RealType>::get_continuity() { return DISCONTINUOUS; }
template <typename RealType> FEContinuity FEShim<3,XYZ,RealType>::get_continuity() { return DISCONTINUOUS; }

// Full specialization of is_hierarchic() function for every dimension.
// The XYZ shape functions are hierarchic!
template <typename RealType> bool FEShim<0,XYZ,RealType>::is_hierarchic() { return true; }
template <typename RealType> bool FEShim<1,XYZ,RealType>::is_hierarchic() { return true; }
template <typename RealType> bool FEShim<2,XYZ,RealType>::is_hierarchic() { return true; }
template <typename RealType> bool FEShim<3,XYZ,RealType>::is_hierarchic() { return true; }

#ifdef LIBMESH_ENABLE_AMR

// Full specialization of compute_constraints() function for 2D and
// 3D only.  There are no constraints for discontinuous elements, so
// we do nothing.
template <typename RealType> void FEShim<2,XYZ,RealType>::compute_constraints (DofConstraints &, DofMap &, const unsigned int, const ElemTempl<Real> *) {}
template <typename RealType> void FEShim<3,XYZ,RealType>::compute_constraints (DofConstraints &, DofMap &, const unsigned int, const ElemTempl<Real> *) {}

#endif // #ifdef LIBMESH_ENABLE_AMR

// Full specialization of shapes_need_reinit() function for every dimension.
template <typename RealType> bool FEShim<0,XYZ,RealType>::shapes_need_reinit() { return true; }
template <typename RealType> bool FEShim<1,XYZ,RealType>::shapes_need_reinit() { return true; }
template <typename RealType> bool FEShim<2,XYZ,RealType>::shapes_need_reinit() { return true; }
template <typename RealType> bool FEShim<3,XYZ,RealType>::shapes_need_reinit() { return true; }

} // namespace libMesh

#endif // LIBMESH_FE_XYZ_IMPL_H
