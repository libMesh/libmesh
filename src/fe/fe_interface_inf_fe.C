// The libMesh Finite Element Library.
// Copyright (C) 2002-2020 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

#include "libmesh/fe_interface.h"
#include "libmesh/fe_interface_macros.h"
#include "libmesh/inf_fe.h"
#include "libmesh/elem.h"

namespace libMesh
{

//------------------------------------------------------------
//FEInterface class members handling calls to InfFE



unsigned int FEInterface::ifem_n_shape_functions(const unsigned int dim,
                                                 const FEType & fe_t,
                                                 const ElemType t)
{
  switch (dim)
    {
      // 1D
    case 1:
      /*
       * Since InfFE<Dim,T_radial,T_map>::n_shape_functions(...)
       * is actually independent of T_radial and T_map, we can use
       * just any T_radial and T_map
       */
      return InfFE<1,JACOBI_20_00,CARTESIAN>::n_shape_functions(fe_t, t);

      // 2D
    case 2:
      return InfFE<2,JACOBI_20_00,CARTESIAN>::n_shape_functions(fe_t, t);

      // 3D
    case 3:
      return InfFE<3,JACOBI_20_00,CARTESIAN>::n_shape_functions(fe_t, t);

    default:
      libmesh_error_msg("Unsupported dim = " << dim);
    }
}





unsigned int FEInterface::ifem_n_dofs(const unsigned int dim,
                                      const FEType & fe_t,
                                      const ElemType t)
{
  libmesh_deprecated();

  switch (dim)
    {
      // 1D
    case 1:
      /*
       * Since InfFE<Dim,T_radial,T_map>::n_dofs(...)
       * is actually independent of T_radial and T_map, we can use
       * just any T_radial and T_map
       */
      return InfFE<1,JACOBI_20_00,CARTESIAN>::n_dofs(fe_t, t);

      // 2D
    case 2:
      return InfFE<2,JACOBI_20_00,CARTESIAN>::n_dofs(fe_t, t);

      // 3D
    case 3:
      return InfFE<3,JACOBI_20_00,CARTESIAN>::n_dofs(fe_t, t);

    default:
      libmesh_error_msg("Unsupported dim = " << dim);
    }
}



unsigned int
FEInterface::ifem_n_dofs(const FEType & fe_t,
                         const Elem * elem)
{
  switch (elem->dim())
    {
      // 1D
    case 1:
      /*
       * Since InfFE<Dim,T_radial,T_map>::n_dofs(...)
       * is actually independent of T_radial and T_map, we can use
       * just any T_radial and T_map
       */
      return InfFE<1,JACOBI_20_00,CARTESIAN>::n_dofs(fe_t, elem);

      // 2D
    case 2:
      return InfFE<2,JACOBI_20_00,CARTESIAN>::n_dofs(fe_t, elem);

      // 3D
    case 3:
      return InfFE<3,JACOBI_20_00,CARTESIAN>::n_dofs(fe_t, elem);

    default:
      libmesh_error_msg("Unsupported dim = " << elem->dim());
    }
}



unsigned int FEInterface::ifem_n_dofs_at_node(const unsigned int dim,
                                              const FEType & fe_t,
                                              const ElemType t,
                                              const unsigned int n)
{
  // TODO:
  // libmesh_deprecated();

  switch (dim)
    {
      // 1D
    case 1:
      /*
       * Since InfFE<Dim,T_radial,T_map>::n_dofs_at_node(...)
       * is actually independent of T_radial and T_map, we can use
       * just any T_radial and T_map
       */
      return InfFE<1,JACOBI_20_00,CARTESIAN>::n_dofs_at_node(fe_t, t, n);

      // 2D
    case 2:
      return InfFE<2,JACOBI_20_00,CARTESIAN>::n_dofs_at_node(fe_t, t, n);

      // 3D
    case 3:
      return InfFE<3,JACOBI_20_00,CARTESIAN>::n_dofs_at_node(fe_t, t, n);

    default:
      libmesh_error_msg("Unsupported dim = " << dim);
    }
}



unsigned int FEInterface::ifem_n_dofs_at_node(const FEType & fe_t,
                                              const Elem * elem,
                                              const unsigned int n)
{
  switch (elem->dim())
    {
      // 1D
    case 1:
      /*
       * Since InfFE<Dim,T_radial,T_map>::n_dofs_at_node(...)
       * is actually independent of T_radial and T_map, we can use
       * just any T_radial and T_map
       */
      return InfFE<1,JACOBI_20_00,CARTESIAN>::n_dofs_at_node(fe_t, elem, n);

      // 2D
    case 2:
      return InfFE<2,JACOBI_20_00,CARTESIAN>::n_dofs_at_node(fe_t, elem, n);

      // 3D
    case 3:
      return InfFE<3,JACOBI_20_00,CARTESIAN>::n_dofs_at_node(fe_t, elem, n);

    default:
      libmesh_error_msg("Unsupported dim = " << elem->dim());
    }
}





unsigned int FEInterface::ifem_n_dofs_per_elem(const unsigned int dim,
                                               const FEType & fe_t,
                                               const ElemType t)
{
  // TODO:
  // libmesh_deprecated();

  switch (dim)
    {
      // 1D
    case 1:
      /*
       * Since InfFE<Dim,T_radial,T_map>::n_dofs(...)
       * is actually independent of T_radial and T_map, we can use
       * just any T_radial and T_map
       */
      return InfFE<1,JACOBI_20_00,CARTESIAN>::n_dofs_per_elem(fe_t, t);

      // 2D
    case 2:
      return InfFE<2,JACOBI_20_00,CARTESIAN>::n_dofs_per_elem(fe_t, t);

      // 3D
    case 3:
      return InfFE<3,JACOBI_20_00,CARTESIAN>::n_dofs_per_elem(fe_t, t);

    default:
      libmesh_error_msg("Unsupported dim = " << dim);
    }
}



unsigned int FEInterface::ifem_n_dofs_per_elem(const FEType & fe_t,
                                               const Elem * elem)
{
  switch (elem->dim())
    {
      // 1D
    case 1:
      /*
       * Since InfFE<Dim,T_radial,T_map>::n_dofs(...)
       * is actually independent of T_radial and T_map, we can use
       * just any T_radial and T_map
       */
      return InfFE<1,JACOBI_20_00,CARTESIAN>::n_dofs_per_elem(fe_t, elem);

      // 2D
    case 2:
      return InfFE<2,JACOBI_20_00,CARTESIAN>::n_dofs_per_elem(fe_t, elem);

      // 3D
    case 3:
      return InfFE<3,JACOBI_20_00,CARTESIAN>::n_dofs_per_elem(fe_t, elem);

    default:
      libmesh_error_msg("Unsupported dim = " << elem->dim());
    }
}




void FEInterface::ifem_nodal_soln(const unsigned int dim,
                                  const FEType & fe_t,
                                  const Elem * elem,
                                  const std::vector<Number> & elem_soln,
                                  std::vector<Number> & nodal_soln)
{
  switch (dim)
    {

      // 1D
    case 1:
      {
        switch (fe_t.radial_family)
          {
          case INFINITE_MAP:
            libmesh_error_msg("ERROR: INFINITE_MAP is not a valid shape family for radial approximation.");

          case JACOBI_20_00:
            {
              switch (fe_t.inf_map)
                {
                case CARTESIAN:
                  {
                    InfFE<1,JACOBI_20_00,CARTESIAN>::nodal_soln(fe_t, elem, elem_soln, nodal_soln);
                    break;
                  }
                default:
                  libmesh_error_msg("ERROR: Spherical & Ellipsoidal IFEMs not implemented.");
                }
              break;
            }

          case JACOBI_30_00:
            {
              switch (fe_t.inf_map)
                {
                case CARTESIAN:
                  {
                    InfFE<1,JACOBI_30_00,CARTESIAN>::nodal_soln(fe_t, elem, elem_soln, nodal_soln);
                    break;
                  }
                default:
                  libmesh_error_msg("ERROR: Spherical & Ellipsoidal IFEMs not implemented.");
                }
              break;
            }

          case LEGENDRE:
            {
              switch (fe_t.inf_map)
                {
                case CARTESIAN:
                  {
                    InfFE<1,LEGENDRE,CARTESIAN>::nodal_soln(fe_t, elem, elem_soln, nodal_soln);
                    break;
                  }
                default:
                  libmesh_error_msg("ERROR: Spherical & Ellipsoidal IFEMs not implemented.");
                }
              break;
            }

          case LAGRANGE:
            {
              switch (fe_t.inf_map)
                {
                case CARTESIAN:
                  {
                    InfFE<1,LAGRANGE,CARTESIAN>::nodal_soln(fe_t, elem, elem_soln, nodal_soln);
                    break;
                  }
                default:
                  libmesh_error_msg("ERROR: Spherical & Ellipsoidal IFEMs not implemented.");
                }
              break;
            }

          default:
            libmesh_error_msg("ERROR: Bad FEType.radial_family= " << fe_t.radial_family);
          }

        break;
      }




      // 2D
    case 2:
      {
        switch (fe_t.radial_family)
          {
          case INFINITE_MAP:
            libmesh_error_msg("ERROR: INFINITE_MAP is not a valid shape family for radial approximation.");

          case JACOBI_20_00:
            {
              switch (fe_t.inf_map)
                {
                case CARTESIAN:
                  {
                    InfFE<2,JACOBI_20_00,CARTESIAN>::nodal_soln(fe_t, elem, elem_soln, nodal_soln);
                    break;
                  }
                default:
                  libmesh_error_msg("ERROR: Spherical & Ellipsoidal IFEMs not implemented.");
                }
              break;
            }

          case JACOBI_30_00:
            {
              switch (fe_t.inf_map)
                {
                case CARTESIAN:
                  {
                    InfFE<2,JACOBI_30_00,CARTESIAN>::nodal_soln(fe_t, elem, elem_soln, nodal_soln);
                    break;
                  }
                default:
                  libmesh_error_msg("ERROR: Spherical & Ellipsoidal IFEMs not implemented.");
                }
              break;
            }

          case LEGENDRE:
            {
              switch (fe_t.inf_map)
                {
                case CARTESIAN:
                  {
                    InfFE<2,LEGENDRE,CARTESIAN>::nodal_soln(fe_t, elem, elem_soln, nodal_soln);
                    break;
                  }
                default:
                  libmesh_error_msg("ERROR: Spherical & Ellipsoidal IFEMs not implemented.");
                }
              break;
            }

          case LAGRANGE:
            {
              switch (fe_t.inf_map)
                {
                case CARTESIAN:
                  {
                    InfFE<2,LAGRANGE,CARTESIAN>::nodal_soln(fe_t, elem, elem_soln, nodal_soln);
                    break;
                  }
                default:
                  libmesh_error_msg("ERROR: Spherical & Ellipsoidal IFEMs not implemented.");
                }
              break;
            }

          default:
            libmesh_error_msg("ERROR: Bad FEType.radial_family= " << fe_t.radial_family);
          }

        break;
      }




      // 3D
    case 3:
      {
        switch (fe_t.radial_family)
          {
          case INFINITE_MAP:
            libmesh_error_msg("ERROR: INFINITE_MAP is not a valid shape family for radial approximation.");

          case JACOBI_20_00:
            {
              switch (fe_t.inf_map)
                {
                case CARTESIAN:
                  {
                    InfFE<3,JACOBI_20_00,CARTESIAN>::nodal_soln(fe_t, elem, elem_soln, nodal_soln);
                    break;
                  }
                default:
                  libmesh_error_msg("ERROR: Spherical & Ellipsoidal IFEMs not implemented.");
                }
              break;
            }

          case JACOBI_30_00:
            {
              switch (fe_t.inf_map)
                {
                case CARTESIAN:
                  {
                    InfFE<3,JACOBI_30_00,CARTESIAN>::nodal_soln(fe_t, elem, elem_soln, nodal_soln);
                    break;
                  }
                default:
                  libmesh_error_msg("ERROR: Spherical & Ellipsoidal IFEMs not implemented.");
                }
              break;
            }

          case LEGENDRE:
            {
              switch (fe_t.inf_map)
                {
                case CARTESIAN:
                  {
                    InfFE<3,LEGENDRE,CARTESIAN>::nodal_soln(fe_t, elem, elem_soln, nodal_soln);
                    break;
                  }
                default:
                  libmesh_error_msg("ERROR: Spherical & Ellipsoidal IFEMs not implemented.");
                }
              break;
            }

          case LAGRANGE:
            {
              switch (fe_t.inf_map)
                {
                case CARTESIAN:
                  {
                    InfFE<3,LAGRANGE,CARTESIAN>::nodal_soln(fe_t, elem, elem_soln, nodal_soln);
                    break;
                  }
                default:
                  libmesh_error_msg("ERROR: Spherical & Ellipsoidal IFEMs not implemented.");
                }
              break;
            }



          default:
            libmesh_error_msg("ERROR: Bad FEType.radial_family= " << fe_t.radial_family);
          }

        break;
      }

    default:
      libmesh_error_msg("Invalid dim = " << dim);
    }
}




Point FEInterface::ifem_map (const unsigned int dim,
                             const FEType & fe_t,
                             const Elem * elem,
                             const Point & p)
{
  switch (fe_t.inf_map)
    {
    case CARTESIAN:
      {
        switch (dim)
          {
          case 1:
            return InfFE<1,JACOBI_20_00,CARTESIAN>::map(elem, p);
          case 2:
            return InfFE<2,JACOBI_20_00,CARTESIAN>::map(elem, p);
          case 3:
            return InfFE<3,JACOBI_20_00,CARTESIAN>::map(elem, p);
          default:
            libmesh_error_msg("Invalid dim = " << dim);
          }
      }
    case SPHERICAL:
    case ELLIPSOIDAL:
      libmesh_not_implemented_msg("ERROR: Spherical and Ellipsoidal IFEMs not (yet) implemented.");
    default:
      libmesh_error_msg("Invalid map = " << fe_t.inf_map);
    }
}



Point FEInterface::ifem_inverse_map (const unsigned int dim,
                                     const FEType & fe_t,
                                     const Elem * elem,
                                     const Point & p,
                                     const Real tolerance,
                                     const bool secure)
{
  switch (dim)
    {
      // 1D
    case 1:
      {
        switch (fe_t.inf_map)
          {
          case CARTESIAN:
            return InfFE<1,JACOBI_20_00,CARTESIAN>::inverse_map(elem, p, tolerance, secure);

          case SPHERICAL:
          case ELLIPSOIDAL:
            libmesh_not_implemented_msg("ERROR: Spherical and Ellipsoidal IFEMs not (yet) implemented.");

            /*
              case SPHERICAL:
              return InfFE<1,JACOBI_20_00,SPHERICAL>::inverse_map(elem, p, tolerance);

              case ELLIPSOIDAL:
              return InfFE<1,JACOBI_20_00,ELLIPSOIDAL>::inverse_map(elem, p, tolerance);
            */

          default:
            libmesh_error_msg("Invalid map = " << fe_t.inf_map);
          }
      }


      // 2D
    case 2:
      {
        switch (fe_t.inf_map)
          {
          case CARTESIAN:
            return InfFE<2,JACOBI_20_00,CARTESIAN>::inverse_map(elem, p, tolerance, secure);

          case SPHERICAL:
          case ELLIPSOIDAL:
            libmesh_not_implemented_msg("ERROR: Spherical and Ellipsoidal IFEMs not (yet) implemented.");

            /*
              case SPHERICAL:
              return InfFE<2,JACOBI_20_00,SPHERICAL>::inverse_map(elem, p, tolerance);

              case ELLIPSOIDAL:
              return InfFE<2,JACOBI_20_00,ELLIPSOIDAL>::inverse_map(elem, p, tolerance);
            */

          default:
            libmesh_error_msg("Invalid map = " << fe_t.inf_map);
          }
      }


      // 3D
    case 3:
      {
        switch (fe_t.inf_map)
          {
          case CARTESIAN:
            return InfFE<3,JACOBI_20_00,CARTESIAN>::inverse_map(elem, p, tolerance, secure);

          case SPHERICAL:
          case ELLIPSOIDAL:
            libmesh_not_implemented_msg("ERROR: Spherical and Ellipsoidal IFEMs not (yet) implemented.");

            /*
              case SPHERICAL:
              return InfFE<3,JACOBI_20_00,SPHERICAL>::inverse_map(elem, p, tolerance);

              case ELLIPSOIDAL:
              return InfFE<3,JACOBI_20_00,ELLIPSOIDAL>::inverse_map(elem, p, tolerance);
            */

          default:
            libmesh_error_msg("Invalid map = " << fe_t.inf_map);
          }
      }

    default:
      libmesh_error_msg("Invalid dim = " << dim);
    }
}



void FEInterface::ifem_inverse_map (const unsigned int dim,
                                    const FEType & fe_t,
                                    const Elem * elem,
                                    const std::vector<Point> & physical_points,
                                    std::vector<Point> &       reference_points,
                                    const Real tolerance,
                                    const bool secure)
{
  switch (dim)
    {
      // 1D
    case 1:
      {
        switch (fe_t.inf_map)
          {
          case CARTESIAN:
            InfFE<1,JACOBI_20_00,CARTESIAN>::inverse_map(elem, physical_points, reference_points, tolerance, secure);
            return;

          default:
            libmesh_error_msg("Invalid map = " << fe_t.inf_map);
          }
      }


      // 2D
    case 2:
      {
        switch (fe_t.inf_map)
          {
          case CARTESIAN:
            InfFE<2,JACOBI_20_00,CARTESIAN>::inverse_map(elem, physical_points, reference_points, tolerance, secure);
            return;

          default:
            libmesh_error_msg("Invalid map = " << fe_t.inf_map);
          }
      }


      // 3D
    case 3:
      {
        switch (fe_t.inf_map)
          {
          case CARTESIAN:
            InfFE<3,JACOBI_20_00,CARTESIAN>::inverse_map(elem, physical_points, reference_points, tolerance, secure);
            return;

          default:
            libmesh_error_msg("Invalid map = " << fe_t.inf_map);
          }
      }

    default:
      libmesh_error_msg("Invalid dim = " << dim);
    }
}




bool FEInterface::ifem_on_reference_element(const Point & p,
                                            const ElemType t,
                                            const Real eps)
{
  return FEBase::on_reference_element(p,t,eps);
}



Real FEInterface::ifem_shape(const unsigned int dim,
                             const FEType & fe_t,
                             const ElemType t,
                             const unsigned int i,
                             const Point & p)
{
  // TODO:
  // libmesh_deprecated();

  inf_fe_switch(shape(fe_t, t, i, p));
}




Real FEInterface::ifem_shape(const unsigned int dim,
                             const FEType & fe_t,
                             const Elem * elem,
                             const unsigned int i,
                             const Point & p)
{
  // TODO:
  // libmesh_deprecated();

  inf_fe_switch( shape(fe_t, elem, i, p));
}

Real FEInterface::ifem_shape(const FEType & fe_t,
                             const Elem * elem,
                             const unsigned int i,
                             const Point & p)
{
  // The inf_fe_switch macro requires a "dim" parameter.
  auto dim = elem->dim();

  inf_fe_switch( shape(fe_t, elem, i, p));
}

Real FEInterface::ifem_shape_deriv (const unsigned int dim,
                                    const FEType & fe_t,
                                    const Elem * elem,
                                    const unsigned int i,
                                    const unsigned int j,
                                    const Point & p)
{
  // TODO:
  // libmesh_deprecated();

  inf_fe_switch(shape_deriv(fe_t, elem, i, j, p));
}


Real FEInterface::ifem_shape_deriv(const unsigned int dim,
                                   const FEType & fe_t,
                                   const ElemType t,
                                   const unsigned int i,
                                   const unsigned int j,
                                   const Point & p)
{
  // TODO:
  // libmesh_deprecated();

  inf_fe_switch(shape_deriv(fe_t, t, i, j, p));
}



Real FEInterface::ifem_shape_deriv (const FEType & fe_t,
                                    const Elem * elem,
                                    const unsigned int i,
                                    const unsigned int j,
                                    const Point & p)
{
  // The inf_fe_switch macro requires a "dim" parameter.
  auto dim = elem->dim();

  inf_fe_switch(shape_deriv(fe_t, elem, i, j, p));
}


void FEInterface::ifem_compute_data(const unsigned int dim,
                                    const FEType & fe_t,
                                    const Elem * elem,
                                    FEComputeData & data)
{
  switch (dim)
    {
    case 1:
      {
        inf_fe_family_mapping_switch(1, compute_data(fe_t, elem,data), , ;break;);
        break;
      }
    case 2:
      {
        inf_fe_family_mapping_switch(2, compute_data(fe_t, elem,data), , ;break;);
        break;
      }
    case 3:
      {
        inf_fe_family_mapping_switch(3, compute_data(fe_t, elem,data), , ;break;);
        break;
      }


    default:
      libmesh_error_msg("Invalid dim = " << dim);
      break;
    }
}

} // namespace libMesh

#endif // ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
