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

#include "libmesh/fe_interface_impl.h"

namespace libMesh
{

//------------------------------------------------------------
//FEInterface class members
FEInterface::FEInterface()
{
  libmesh_error_msg("ERROR: Do not define an object of this type.");
}



#ifndef LIBMESH_ENABLE_INFINITE_ELEMENTS

bool
FEInterface::is_InfFE_elem(const ElemType)
{
  return false;
}

#else

bool
FEInterface::is_InfFE_elem(const ElemType et)
{

  switch (et)
    {
    case INFEDGE2:
    case INFQUAD4:
    case INFQUAD6:
    case INFHEX8:
    case INFHEX16:
    case INFHEX18:
    case INFPRISM6:
    case INFPRISM12:
      {
        return true;
      }

    default:
      {
        return false;
      }
    }
}

#endif //ifndef LIBMESH_ENABLE_INFINITE_ELEMENTS


unsigned int FEInterface::n_shape_functions(const unsigned int dim,
                                            const FEType & fe_t,
                                            const ElemType t)
{

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
  /*
   * Since the FEType, stored in DofMap/(some System child), has to
   * be the _same_ for InfFE and FE, we have to catch calls
   * to infinite elements through the element type.
   */

  if (is_InfFE_elem(t))
    return ifem_n_shape_functions(dim, fe_t, t);

#endif

  const Order o = fe_t.order;

  fe_with_vec_switch(n_shape_functions(t, o), Real);
}





unsigned int FEInterface::n_dofs(const unsigned int dim,
                                 const FEType & fe_t,
                                 const ElemType t)
{
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

  if (is_InfFE_elem(t))
    return ifem_n_dofs(dim, fe_t, t);

#endif

  const Order o = fe_t.order;

  fe_with_vec_switch(n_dofs(t, o), Real);
}


unsigned int FEInterface::n_dofs_at_node(const unsigned int dim,
                                         const FEType & fe_t,
                                         const ElemType t,
                                         const unsigned int n)
{
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

  if (is_InfFE_elem(t))
    return ifem_n_dofs_at_node(dim, fe_t, t, n);

#endif

  const Order o = fe_t.order;

  fe_with_vec_switch(n_dofs_at_node(t, o, n), Real);
}



FEInterface::n_dofs_at_node_ptr
FEInterface::n_dofs_at_node_function(const unsigned int dim,
                                     const FEType & fe_t)
{
  fe_with_vec_switch(n_dofs_at_node, Real);
}


unsigned int FEInterface::n_dofs_per_elem(const unsigned int dim,
                                          const FEType & fe_t,
                                          const ElemType t)
{
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

  if (is_InfFE_elem(t))
    return ifem_n_dofs_per_elem(dim, fe_t, t);

#endif

  const Order o = fe_t.order;

  fe_with_vec_switch(n_dofs_per_elem(t, o), Real);
}


unsigned int FEInterface::max_order(const FEType & fe_t,
                                    const ElemType & el_t)
{
  // Yeah, I know, infinity is much larger than 11, but our
  // solvers don't seem to like high degree polynomials, and our
  // quadrature rules and number_lookups tables
  // need to go up higher.
  const unsigned int unlimited = 11;

  // If we used 0 as a default, then elements missing from this
  // table (e.g. infinite elements) would be considered broken.
  const unsigned int unknown = unlimited;

  switch (fe_t.family)
    {
    case LAGRANGE:
    case L2_LAGRANGE: // TODO: L2_LAGRANGE can have higher "max_order" than LAGRANGE
    case LAGRANGE_VEC:
      switch (el_t)
        {
        case EDGE2:
        case EDGE3:
        case EDGE4:
          return 3;
        case TRI3:
        case TRISHELL3:
          return 1;
        case TRI6:
          return 2;
        case QUAD4:
        case QUADSHELL4:
          return 1;
        case QUAD8:
        case QUADSHELL8:
        case QUAD9:
          return 2;
        case TET4:
          return 1;
        case TET10:
          return 2;
        case HEX8:
          return 1;
        case HEX20:
        case HEX27:
          return 2;
        case PRISM6:
          return 1;
        case PRISM15:
        case PRISM18:
          return 2;
        case PYRAMID5:
          return 1;
        case PYRAMID13:
        case PYRAMID14:
          return 2;
        default:
          return unknown;
        }
      break;
    case MONOMIAL:
    case MONOMIAL_VEC:
      switch (el_t)
        {
        case EDGE2:
        case EDGE3:
        case EDGE4:
        case TRI3:
        case TRISHELL3:
        case TRI6:
        case QUAD4:
        case QUADSHELL4:
        case QUAD8:
        case QUADSHELL8:
        case QUAD9:
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
          return unlimited;
        default:
          return unknown;
        }
      break;
#ifdef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES
    case BERNSTEIN:
    case RATIONAL_BERNSTEIN:
      switch (el_t)
        {
        case EDGE2:
        case EDGE3:
        case EDGE4:
          return unlimited;
        case TRI3:
        case TRISHELL3:
          return 1;
        case TRI6:
          return 6;
        case QUAD4:
        case QUADSHELL4:
          return 1;
        case QUAD8:
        case QUADSHELL8:
        case QUAD9:
          return unlimited;
        case TET4:
          return 1;
        case TET10:
          return 2;
        case HEX8:
          return 1;
        case HEX20:
          return 2;
        case HEX27:
          return 4;
        case PRISM6:
        case PRISM15:
        case PRISM18:
        case PYRAMID5:
        case PYRAMID13:
        case PYRAMID14:
          return 0;
        default:
          return unknown;
        }
      break;
    case SZABAB:
      switch (el_t)
        {
        case EDGE2:
          return 1;
        case EDGE3:
        case EDGE4:
          return 7;
        case TRI3:
        case TRISHELL3:
          return 1;
        case TRI6:
          return 7;
        case QUAD4:
        case QUADSHELL4:
          return 1;
        case QUAD8:
        case QUADSHELL8:
        case QUAD9:
          return 7;
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
          return 0;
        default:
          return unknown;
        }
      break;
#endif
    case XYZ:
      switch (el_t)
        {
        case EDGE2:
        case EDGE3:
        case EDGE4:
        case TRI3:
        case TRISHELL3:
        case TRI6:
        case QUAD4:
        case QUADSHELL4:
        case QUAD8:
        case QUADSHELL8:
        case QUAD9:
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
          return unlimited;
        default:
          return unknown;
        }
      break;
    case CLOUGH:
      switch (el_t)
        {
        case EDGE2:
        case EDGE3:
          return 3;
        case EDGE4:
        case TRI3:
        case TRISHELL3:
          return 0;
        case TRI6:
          return 3;
        case QUAD4:
        case QUADSHELL4:
        case QUAD8:
        case QUADSHELL8:
        case QUAD9:
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
          return 0;
        default:
          return unknown;
        }
      break;
    case HERMITE:
      switch (el_t)
        {
        case EDGE2:
        case EDGE3:
          return unlimited;
        case EDGE4:
        case TRI3:
        case TRISHELL3:
        case TRI6:
          return 0;
        case QUAD4:
        case QUADSHELL4:
          return 3;
        case QUAD8:
        case QUADSHELL8:
        case QUAD9:
          return unlimited;
        case TET4:
        case TET10:
          return 0;
        case HEX8:
          return 3;
        case HEX20:
        case HEX27:
          return unlimited;
        case PRISM6:
        case PRISM15:
        case PRISM18:
        case PYRAMID5:
        case PYRAMID13:
        case PYRAMID14:
          return 0;
        default:
          return unknown;
        }
      break;
    case HIERARCHIC:
      switch (el_t)
        {
        case EDGE2:
        case EDGE3:
        case EDGE4:
          return unlimited;
        case TRI3:
        case TRISHELL3:
          return 1;
        case TRI6:
          return unlimited;
        case QUAD4:
        case QUADSHELL4:
          return 1;
        case QUAD8:
        case QUADSHELL8:
        case QUAD9:
          return unlimited;
        case TET4:
        case TET10:
          return 0;
        case HEX8:
        case HEX20:
          return 1;
        case HEX27:
          return unlimited;
        case PRISM6:
        case PRISM15:
        case PRISM18:
        case PYRAMID5:
        case PYRAMID13:
        case PYRAMID14:
          return 0;
        default:
          return unknown;
        }
      break;
    case L2_HIERARCHIC:
      switch (el_t)
        {
        case EDGE2:
        case EDGE3:
        case EDGE4:
          return unlimited;
        case TRI3:
        case TRISHELL3:
          return 1;
        case TRI6:
          return unlimited;
        case QUAD4:
        case QUADSHELL4:
          return 1;
        case QUAD8:
        case QUADSHELL8:
        case QUAD9:
          return unlimited;
        case TET4:
        case TET10:
          return 0;
        case HEX8:
        case HEX20:
          return 1;
        case HEX27:
          return unlimited;
        case PRISM6:
        case PRISM15:
        case PRISM18:
        case PYRAMID5:
        case PYRAMID13:
        case PYRAMID14:
          return 0;
        default:
          return unknown;
        }
      break;
    case SUBDIVISION:
      switch (el_t)
        {
        case TRI3SUBDIVISION:
          return unlimited;
        default:
          return unknown;
        }
      break;
    case NEDELEC_ONE:
      switch (el_t)
        {
        case TRI6:
        case QUAD8:
        case QUAD9:
        case HEX20:
        case HEX27:
          return 1;
        default:
          return 0;
        }
      break;
    default:
      return 0;
      break;
    }
}



bool FEInterface::extra_hanging_dofs(const FEType & fe_t)
{
  switch (fe_t.family)
    {
    case LAGRANGE:
    case L2_LAGRANGE:
    case MONOMIAL:
    case MONOMIAL_VEC:
    case L2_HIERARCHIC:
    case XYZ:
    case SUBDIVISION:
    case LAGRANGE_VEC:
    case NEDELEC_ONE:
      return false;
    case CLOUGH:
    case HERMITE:
    case HIERARCHIC:
#ifdef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES
    case BERNSTEIN:
    case SZABAB:
    case RATIONAL_BERNSTEIN:
#endif
    default:
      return true;
    }
}

FEFieldType FEInterface::field_type(const FEType & fe_type)
{
  return FEInterface::field_type(fe_type.family);
}

FEFieldType FEInterface::field_type (const FEFamily & fe_family)
{
  switch (fe_family)
    {
    case LAGRANGE_VEC:
    case NEDELEC_ONE:
    case MONOMIAL_VEC:
      return TYPE_VECTOR;
    default:
      return TYPE_SCALAR;
    }
}


FEContinuity FEInterface::get_continuity(const FEType & fe_type)
{
  switch (fe_type.family)
    {
      // Discontinuous elements
    case MONOMIAL:
    case MONOMIAL_VEC:
    case L2_HIERARCHIC:
    case L2_LAGRANGE:
    case XYZ:
    case SCALAR:
      return DISCONTINUOUS;

      // C0 elements
    case LAGRANGE:
    case HIERARCHIC:
    case BERNSTEIN:
    case SZABAB:
    case RATIONAL_BERNSTEIN:
    case INFINITE_MAP:
    case JACOBI_20_00:
    case JACOBI_30_00:
    case LEGENDRE:
    case LAGRANGE_VEC:
      return C_ZERO;

      // C1 elements
    case CLOUGH:
    case HERMITE:
    case SUBDIVISION:
      return C_ONE;

    case NEDELEC_ONE:
      return H_CURL;

    default:
      libmesh_error_msg("Unknown FE Family " << Utility::enum_to_string(fe_type.family));
    }
}

void FEInterface::compute_data(const unsigned int dim,
                               const FEType & fe_t,
                               const ElemTempl<Real> * elem,
                               FEComputeData & data)
{
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

  if (elem && is_InfFE_elem(elem->type()))
    {
      data.init();
      ifem_compute_data(dim, fe_t, elem, data);
      return;
    }

#endif

  const unsigned int n_dof = n_dofs (dim, fe_t, elem);
  const PointTempl<Real> & p = data.p;
  data.shape.resize(n_dof);

  if (data.need_derivative())
    {
      data.dshape.resize(n_dof);
      data.local_transform.resize(dim);

      for (unsigned int d=0; d<dim; d++)
        data.local_transform[d].resize(dim);

      UniquePtr<FEGenericBase<Real,Real>> fe (FEGenericBase<Real,Real>::build(dim, fe_t));
      std::vector<PointTempl<Real>> pt(1);
      pt[0]=p;
      fe->get_dphideta(); // to compute the map
      fe->reinit(elem, &pt);

      // compute the reference->physical map.
      data.local_transform[0][0] = fe->get_dxidx()[0];
      if (dim > 1)
        {
          data.local_transform[1][0] = fe->get_detadx()[0];
          data.local_transform[1][1] = fe->get_detady()[0];
          data.local_transform[0][1] = fe->get_dxidy()[0];
          if (dim > 2)
            {
              data.local_transform[2][0] = fe->get_dzetadx()[0];
              data.local_transform[2][1] = fe->get_dzetady()[0];
              data.local_transform[2][2] = fe->get_dzetadz()[0];
              data.local_transform[1][2] = fe->get_detadz()[0];
              data.local_transform[0][2] = fe->get_dxidz()[0];
            }
        }
    }

  // set default values for all the output fields
  data.init();

  for (unsigned int n=0; n<n_dof; n++)
    {
      // Here we pass the original fe_t object. Additional p-levels
      // (if any) are handled internally by the shape() and
      // shape_deriv() functions since they have access to the elem
      // pointer. Note that we are already using the n_dof value
      // appropriate to the elevated p-level.
      data.shape[n] = shape(dim, fe_t, elem, n, p);
      if (data.need_derivative())
        {
          for (unsigned int j=0; j<dim; j++)
            data.dshape[n](j) = shape_deriv(dim, fe_t, elem, n, j, p);
        }
    }
}



#ifdef LIBMESH_ENABLE_AMR

void FEInterface::compute_constraints (DofConstraints & constraints,
                                       DofMap & dof_map,
                                       const unsigned int variable_number,
                                       const ElemTempl<Real> * elem)
{
  libmesh_assert(elem);

  const FEType & fe_t = dof_map.variable_type(variable_number);

  switch (elem->dim())
    {
    case 0:
    case 1:
      {
        // No constraints in 0D/1D.
        return;
      }


    case 2:
      {
        switch (fe_t.family)
          {
          case CLOUGH:
            FE<2,CLOUGH>::compute_constraints (constraints,
                                               dof_map,
                                               variable_number,
                                               elem); return;

          case HERMITE:
            FE<2,HERMITE>::compute_constraints (constraints,
                                                dof_map,
                                                variable_number,
                                                elem); return;

          case LAGRANGE:
            FE<2,LAGRANGE>::compute_constraints (constraints,
                                                 dof_map,
                                                 variable_number,
                                                 elem); return;

          case HIERARCHIC:
            FE<2,HIERARCHIC>::compute_constraints (constraints,
                                                   dof_map,
                                                   variable_number,
                                                   elem); return;

          case L2_HIERARCHIC:
            FE<2,L2_HIERARCHIC>::compute_constraints (constraints,
                                                      dof_map,
                                                      variable_number,
                                                      elem); return;

          case LAGRANGE_VEC:
            FE<2,LAGRANGE_VEC>::compute_constraints (constraints,
                                                     dof_map,
                                                     variable_number,
                                                     elem); return;


#ifdef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES
          case SZABAB:
            FE<2,SZABAB>::compute_constraints (constraints,
                                               dof_map,
                                               variable_number,
                                               elem); return;

          case BERNSTEIN:
            FE<2,BERNSTEIN>::compute_constraints (constraints,
                                                  dof_map,
                                                  variable_number,
                                                  elem); return;

          case RATIONAL_BERNSTEIN:
            FE<2,RATIONAL_BERNSTEIN>::compute_constraints (constraints,
                                                           dof_map,
                                                           variable_number,
                                                           elem); return;

#endif
          default:
            return;
          }
      }


    case 3:
      {
        switch (fe_t.family)
          {
          case HERMITE:
            FE<3,HERMITE>::compute_constraints (constraints,
                                                dof_map,
                                                variable_number,
                                                elem); return;

          case LAGRANGE:
            FE<3,LAGRANGE>::compute_constraints (constraints,
                                                 dof_map,
                                                 variable_number,
                                                 elem); return;

          case HIERARCHIC:
            FE<3,HIERARCHIC>::compute_constraints (constraints,
                                                   dof_map,
                                                   variable_number,
                                                   elem); return;

          case L2_HIERARCHIC:
            FE<3,L2_HIERARCHIC>::compute_constraints (constraints,
                                                      dof_map,
                                                      variable_number,
                                                      elem); return;

          case LAGRANGE_VEC:
            FE<3,LAGRANGE_VEC>::compute_constraints (constraints,
                                                     dof_map,
                                                     variable_number,
                                                     elem); return;
#ifdef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES
          case SZABAB:
            FE<3,SZABAB>::compute_constraints (constraints,
                                               dof_map,
                                               variable_number,
                                               elem); return;

          case BERNSTEIN:
            FE<3,BERNSTEIN>::compute_constraints (constraints,
                                                  dof_map,
                                                  variable_number,
                                                  elem); return;

          case RATIONAL_BERNSTEIN:
            FE<3,RATIONAL_BERNSTEIN>::compute_constraints (constraints,
                                                           dof_map,
                                                           variable_number,
                                                           elem); return;

#endif
          default:
            return;
          }
      }


    default:
      libmesh_error_msg("Invalid dimension = " << elem->dim());
    }
}

#endif // #ifdef LIBMESH_ENABLE_AMR



#ifdef LIBMESH_ENABLE_PERIODIC

void FEInterface::compute_periodic_constraints (DofConstraints & constraints,
                                                DofMap & dof_map,
                                                const PeriodicBoundaries & boundaries,
                                                const MeshBaseTempl<Real> & mesh,
                                                const PointLocatorBase * point_locator,
                                                const unsigned int variable_number,
                                                const ElemTempl<Real> * elem)
{
  // No element-specific optimizations currently exist
  FEGenericBase<Real,Real>::compute_periodic_constraints (constraints,
                                                          dof_map,
                                                          boundaries,
                                                          mesh,
                                                          point_locator,
                                                          variable_number,
                                                          elem);
}

#endif // #ifdef LIBMESH_ENABLE_PERIODIC

INSTANTIATE_FE_INTERFACE_METHODS(Real);
} // namespace libMesh
