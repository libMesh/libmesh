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
#include "libmesh/inf_fe.h"
#include "libmesh/inf_fe_macro.h"
#include "libmesh/fe.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fe_compute_data.h"
#include "libmesh/elem.h"

#include "libmesh/remote_elem.h"

// to fetch NodeConstraintRow:
#include "libmesh/dof_map.h"

namespace libMesh
{


// ------------------------------------------------------------
// InfFE class static member initialization
template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
ElemType InfFE<Dim,T_radial,T_map>::_compute_node_indices_fast_current_elem_type = INVALID_ELEM;

#ifdef DEBUG

template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
bool InfFE<Dim,T_radial,T_map>::_warned_for_nodal_soln = false;


template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
bool InfFE<Dim,T_radial,T_map>::_warned_for_shape      = false;

template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
bool InfFE<Dim,T_radial,T_map>::_warned_for_dshape     = false;

#endif




// ------------------------------------------------------------
// InfFE static class members
template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
unsigned int InfFE<Dim,T_radial,T_map>::n_dofs (const FEType & fet,
                                                const ElemType inf_elem_type)
{
  libmesh_deprecated();

  const ElemType base_et (InfFEBase::get_elem_type(inf_elem_type));

  if (Dim > 1)
    return FEInterface::n_dofs(Dim-1, fet, base_et) *
      InfFERadial::n_dofs(fet.radial_order);
  else
    return InfFERadial::n_dofs(fet.radial_order);
}



template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
unsigned int InfFE<Dim,T_radial,T_map>::n_dofs(const FEType & fet,
                                               const Elem * inf_elem)
{
  // The "base" Elem is a non-infinite Elem corresponding to side 0 of
  // the InfElem. This builds a "lightweight" proxy and so should be
  // relatively fast.
  auto base_elem = inf_elem->build_side_ptr(0);

  if (Dim > 1)
    return FEInterface::n_dofs(fet, base_elem.get()) *
      InfFERadial::n_dofs(fet.radial_order);
  else
    return InfFERadial::n_dofs(fet.radial_order);
}



template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
unsigned int InfFE<Dim,T_radial,T_map>::n_dofs_at_node (const FEType & fet,
                                                        const ElemType inf_elem_type,
                                                        const unsigned int n)
{
  // TODO:
  // libmesh_deprecated();

  const ElemType base_et (InfFEBase::get_elem_type(inf_elem_type));

  unsigned int n_base, n_radial;
  compute_node_indices(inf_elem_type, n, n_base, n_radial);

  //   libMesh::out << "elem_type=" << inf_elem_type
  //     << ",  fet.radial_order=" << fet.radial_order
  //     << ",  n=" << n
  //     << ",  n_radial=" << n_radial
  //     << ",  n_base=" << n_base
  //     << std::endl;

  if (Dim > 1)
    return FEInterface::n_dofs_at_node(Dim-1, fet, base_et, n_base)
      * InfFERadial::n_dofs_at_node(fet.radial_order, n_radial);
  else
    return InfFERadial::n_dofs_at_node(fet.radial_order, n_radial);
}



template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
unsigned int InfFE<Dim,T_radial,T_map>::n_dofs_at_node (const FEType & fet,
                                                        const Elem * inf_elem,
                                                        const unsigned int n)
{
  // The "base" Elem is a non-infinite Elem corresponding to side 0 of
  // the InfElem. This builds a "lightweight" proxy and so should be
  // relatively fast.
  auto base_elem = inf_elem->build_side_ptr(0);

  unsigned int n_base, n_radial;
  compute_node_indices(inf_elem->type(), n, n_base, n_radial);

  if (Dim > 1)
    return FEInterface::n_dofs_at_node(fet, base_elem.get(), n_base)
      * InfFERadial::n_dofs_at_node(fet.radial_order, n_radial);
  else
    return InfFERadial::n_dofs_at_node(fet.radial_order, n_radial);
}



template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
unsigned int InfFE<Dim,T_radial,T_map>::n_dofs_per_elem (const FEType & fet,
                                                         const ElemType inf_elem_type)
{
  // TODO:
  // libmesh_deprecated();

  const ElemType base_et (InfFEBase::get_elem_type(inf_elem_type));

  if (Dim > 1)
    return FEInterface::n_dofs_per_elem(Dim-1, fet, base_et)
      * InfFERadial::n_dofs_per_elem(fet.radial_order);
  else
    return InfFERadial::n_dofs_per_elem(fet.radial_order);
}



template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
unsigned int InfFE<Dim,T_radial,T_map>::n_dofs_per_elem (const FEType & fet,
                                                         const Elem * inf_elem)
{
  // The "base" Elem is a non-infinite Elem corresponding to side 0 of
  // the InfElem. This builds a "lightweight" proxy and so should be
  // relatively fast.
  auto base_elem = inf_elem->build_side_ptr(0);

  if (Dim > 1)
    return FEInterface::n_dofs_per_elem(fet, base_elem.get())
      * InfFERadial::n_dofs_per_elem(fet.radial_order);
  else
    return InfFERadial::n_dofs_per_elem(fet.radial_order);
}



template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
void InfFE<Dim,T_radial,T_map>::nodal_soln(const FEType & /* fet */,
                                           const Elem * /* elem */,
                                           const std::vector<Number> & /* elem_soln */,
                                           std::vector<Number> &       nodal_soln)
{
#ifdef DEBUG
  if (!_warned_for_nodal_soln)
    {
      libMesh::err << "WARNING: nodal_soln(...) does _not_ work for infinite elements." << std::endl
                   << " Will return an empty nodal solution.  Use " << std::endl
                   << " InfFE<Dim,T_radial,T_map>::compute_data(..) instead!" << std::endl;
      _warned_for_nodal_soln = true;
    }
#endif

  /*
   * In the base the infinite element couples to
   * conventional finite elements.  To not destroy
   * the values there, clear \p nodal_soln.  This
   * indicates to the user of \p nodal_soln to
   * not use this result.
   */
  nodal_soln.clear();
  libmesh_assert (nodal_soln.empty());
  return;
}








template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
Real InfFE<Dim,T_radial,T_map>::shape(const FEType & fet,
                                      const ElemType inf_elem_type,
                                      const unsigned int i,
                                      const Point & p)
{
  // TODO - if possible, not sure if we can easily fully remove this function.
  // libmesh_deprecated();

  libmesh_assert_not_equal_to (Dim, 0);

#ifdef DEBUG
  // this makes only sense when used for mapping
  if ((T_radial != INFINITE_MAP) && !_warned_for_shape)
    {
      libMesh::err << "WARNING: InfFE<Dim,T_radial,T_map>::shape(...) does _not_" << std::endl
                   << " return the correct trial function!  Use " << std::endl
                   << " InfFE<Dim,T_radial,T_map>::compute_data(..) instead!"
                   << std::endl;
      _warned_for_shape = true;
    }
#endif

  const ElemType     base_et  (InfFEBase::get_elem_type(inf_elem_type));
  const Order        o_radial (fet.radial_order);
  const Real         v        (p(Dim-1));

  unsigned int i_base, i_radial;
  compute_shape_indices(fet, inf_elem_type, i, i_base, i_radial);

  //TODO:[SP/DD]  exp(ikr) is still missing here!
  // but is it intended?  It would be probably somehow nice, but than it would be Number, not Real !
  // --> thus it would destroy the interface...
  if (Dim > 1)
    return FEInterface::shape(Dim-1, fet, base_et, i_base, p)
      * InfFE<Dim,T_radial,T_map>::eval(v, o_radial, i_radial)
      * InfFERadial::decay(Dim,v);
  else
    return InfFE<Dim,T_radial,T_map>::eval(v, o_radial, i_radial)
      * InfFERadial::decay(Dim,v);
}






template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
Real InfFE<Dim,T_radial,T_map>::shape(const FEType & fet,
                                      const Elem * inf_elem,
                                      const unsigned int i,
                                      const Point & p)
{
  libmesh_assert(inf_elem);
  libmesh_assert_not_equal_to (Dim, 0);

#ifdef DEBUG
  // this makes only sense when used for mapping
  if ((T_radial != INFINITE_MAP) && !_warned_for_shape)
    {
      libMesh::err << "WARNING: InfFE<Dim,T_radial,T_map>::shape(...) does _not_" << std::endl
                   << " return the correct trial function!  Use " << std::endl
                   << " InfFE<Dim,T_radial,T_map>::compute_data(..) instead!"
                   << std::endl;
      _warned_for_shape = true;
    }
#endif

  const Order o_radial (fet.radial_order);
  const Real v (p(Dim-1));
  std::unique_ptr<const Elem> base_el (inf_elem->build_side_ptr(0));

  unsigned int i_base, i_radial;
  compute_shape_indices(fet, inf_elem, i, i_base, i_radial);

  if (Dim > 1)
    return FEInterface::shape(fet, base_el.get(), i_base, p)
      * InfFE<Dim,T_radial,T_map>::eval(v, o_radial, i_radial)
      * InfFERadial::decay(Dim,v);
  else
    return InfFE<Dim,T_radial,T_map>::eval(v, o_radial, i_radial)
      * InfFERadial::decay(Dim,v);
}



template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
Real InfFE<Dim,T_radial,T_map>::shape(const FEType fet,
                                      const Elem * inf_elem,
                                      const unsigned int i,
                                      const Point & p,
                                      const bool add_p_level)
{
  if (add_p_level)
    {
      FEType tmp_fet=fet;
      tmp_fet = static_cast<Order>(fet.order + add_p_level * inf_elem->p_level());
      return InfFE<Dim,T_radial,T_map>::shape(tmp_fet, inf_elem, i, p);
    }
  return InfFE<Dim,T_radial,T_map>::shape(fet, inf_elem, i, p);
}


template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
Real InfFE<Dim,T_radial,T_map>::shape_deriv (const FEType & fe_t,
                                             const ElemType inf_elem_type,
                                             const unsigned int i,
                                             const unsigned int j,
                                             const Point & p)
{
  // TODO - if possible, not sure if we can easily fully remove this function.
  // libmesh_deprecated();

  libmesh_assert_not_equal_to (Dim, 0);
  libmesh_assert_greater (Dim,j);
#ifdef DEBUG
  // this makes only sense when used for mapping
  if ((T_radial != INFINITE_MAP) && !_warned_for_dshape)
    {
      libMesh::err << "WARNING: InfFE<Dim,T_radial,T_map>::shape_deriv(...) does _not_" << std::endl
                   << " return the correct trial function gradients!  Use " << std::endl
                   << " InfFE<Dim,T_radial,T_map>::compute_data(..) instead!"
                   << std::endl;
      _warned_for_dshape = true;
    }
#endif

  const ElemType     base_et  (InfFEBase::get_elem_type(inf_elem_type));
  const Order o_radial (fe_t.radial_order);
  const Real v (p(Dim-1));

  unsigned int i_base, i_radial;
  compute_shape_indices(fe_t, inf_elem_type, i, i_base, i_radial);

  if (j== Dim -1)
    {
      Real RadialDeriv = InfFE<Dim,T_radial,T_map>::eval_deriv(v, o_radial, i_radial)
        * InfFERadial::decay(Dim,v)
        + InfFE<Dim,T_radial,T_map>::eval(v, o_radial, i_radial)
        * InfFERadial::decay_deriv(Dim, v);

      return FEInterface::shape(Dim-1, fe_t, base_et, i_base, p)*RadialDeriv;
    }

  return FEInterface::shape_deriv(Dim-1, fe_t, base_et, i_base, j, p)
    * InfFE<Dim,T_radial,T_map>::eval(v, o_radial, i_radial)
    * InfFERadial::decay(Dim,v);
}


template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
Real InfFE<Dim,T_radial,T_map>::shape_deriv (const FEType & fe_t,
                                             const Elem * inf_elem,
                                             const unsigned int i,
                                             const unsigned int j,
                                             const Point & p)
{
  libmesh_assert_not_equal_to (Dim, 0);
  libmesh_assert_greater (Dim,j);
#ifdef DEBUG
  // this makes only sense when used for mapping
  if ((T_radial != INFINITE_MAP) && !_warned_for_dshape)
    {
      libMesh::err << "WARNING: InfFE<Dim,T_radial,T_map>::shape_deriv(...) does _not_" << std::endl
                   << " return the correct trial function gradients!  Use " << std::endl
                   << " InfFE<Dim,T_radial,T_map>::compute_data(..) instead!"
                   << std::endl;
      _warned_for_dshape = true;
    }
#endif
  const Order o_radial (fe_t.radial_order);
  const Real v (p(Dim-1));

  std::unique_ptr<const Elem> base_el (inf_elem->build_side_ptr(0));

  unsigned int i_base, i_radial;

  if ((-1. > v ) || (v  > 1.))
    {
      //TODO: This is for debugging. We should never come here.
      //   Therefore we can do very useless things then:
      i_base=0;
    }
  compute_shape_indices(fe_t, inf_elem, i, i_base, i_radial);

  if (j== Dim -1)
    {
      Real RadialDeriv = InfFE<Dim,T_radial,T_map>::eval_deriv(v, o_radial, i_radial)
        * InfFERadial::decay(Dim,v)
        + InfFE<Dim,T_radial,T_map>::eval(v, o_radial, i_radial)
        * InfFERadial::decay_deriv(Dim,v);

      return FEInterface::shape(fe_t, base_el.get(), i_base, p)*RadialDeriv;
    }
  return FEInterface::shape_deriv(fe_t, base_el.get(), i_base, j, p)
    * InfFE<Dim,T_radial,T_map>::eval(v, o_radial, i_radial)
    * InfFERadial::decay(Dim,v);
}





template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
void InfFE<Dim,T_radial,T_map>::compute_data(const FEType & fet,
                                             const Elem * inf_elem,
                                             FEComputeData & data)
{
  libmesh_assert(inf_elem);
  libmesh_assert_not_equal_to (Dim, 0);

  const Order        o_radial             (fet.radial_order);
  const Order        radial_mapping_order (InfFERadial::mapping_order());
  const Point &      p                    (data.p);
  const Real         v                    (p(Dim-1));
  std::unique_ptr<const Elem> base_el (inf_elem->build_side_ptr(0));

  /*
   * compute \p interpolated_dist containing the mapping-interpolated
   * distance of the base point to the origin.  This is the same
   * for all shape functions.  Set \p interpolated_dist to 0, it
   * is added to.
   */
  Real interpolated_dist = 0.;
  switch (Dim)
    {
    case 1:
      {
        libmesh_assert_equal_to (inf_elem->type(), INFEDGE2);
        interpolated_dist =  Point(inf_elem->point(0) - inf_elem->point(1)).norm();
        break;
      }

    case 2:
      {
        const unsigned int n_base_nodes = base_el->n_nodes();

        const Point    origin                 = inf_elem->origin();
        const Order    base_mapping_order     (base_el->default_order());
        const ElemType base_mapping_elem_type (base_el->type());

        // interpolate the base nodes' distances
        for (unsigned int n=0; n<n_base_nodes; n++)
          interpolated_dist += Point(base_el->point(n) - origin).norm()
            * FE<1,LAGRANGE>::shape (base_mapping_elem_type, base_mapping_order, n, p);
        break;
      }

    case 3:
      {
        const unsigned int n_base_nodes = base_el->n_nodes();

        const Point    origin                 = inf_elem->origin();
        const Order    base_mapping_order     (base_el->default_order());
        const ElemType base_mapping_elem_type (base_el->type());

        // interpolate the base nodes' distances
        for (unsigned int n=0; n<n_base_nodes; n++)
          interpolated_dist += Point(base_el->point(n) - origin).norm()
            * FE<2,LAGRANGE>::shape (base_mapping_elem_type, base_mapping_order, n, p);
        break;
      }

    default:
      libmesh_error_msg("Unknown Dim = " << Dim);
    }


  const Real speed = data.speed;

  //TODO: I find it inconvenient to have a quantity phase which is phase/speed.
  //    But it might be better than redefining a quantities meaning.
  data.phase = interpolated_dist                                                      /* together with next line:  */
    * InfFE<Dim,INFINITE_MAP,T_map>::eval(v, radial_mapping_order, 1)/speed;          /* phase(s,t,v)/c            */

  // We assume time-harmonic behavior in this function!

#ifdef LIBMESH_USE_COMPLEX_NUMBERS
  // the wave number
  const Number wavenumber = 2. * libMesh::pi * data.frequency / speed;

  // the exponent for time-harmonic behavior
  // \note: this form is much less general than the implementation of dphase, which can be easily extended to
  //     other forms than e^{i kr}.
  const Number exponent = imaginary                                                   /* imaginary unit            */
    * wavenumber                                                                      /* k  (can be complex)       */
    * data.phase*speed;

  const Number time_harmonic = exp(exponent);                                         /* e^(i*k*phase(s,t,v))      */
#else
  const Number time_harmonic = 1;
#endif //LIBMESH_USE_COMPLEX_NUMBERS

  /*
   * compute \p shape for all dof in the element
   */
  if (Dim > 1)
    {
      const unsigned int n_dof = n_dofs (fet, inf_elem->type());
      data.shape.resize(n_dof);
      if (data.need_derivative())
        {
          data.dshape.resize(n_dof);
          data.local_transform.resize(Dim);

          for (unsigned int d=0; d<Dim; d++)
            data.local_transform[d].resize(Dim);

          // compute the reference->physical map at the point \p p.
          // Use another fe_map to avoid interference with \p this->_fe_map
          // which is initialized at the quadrature points...
          UniquePtr<FEBase> fe (FEBase::build_InfFE(Dim, fet));
          std::vector<Point> pt(1);
          pt[0]=p;
          fe->get_dphideta(); // to compute the map
          fe->reinit(inf_elem, &pt);

          // compute the reference->physical map.
          data.local_transform[0][0] = fe->get_dxidx()[0];
          data.local_transform[1][0] = fe->get_detadx()[0];
          data.local_transform[1][1] = fe->get_detady()[0];
          data.local_transform[0][1] = fe->get_dxidy()[0];
          if (Dim > 2)
            {
              data.local_transform[2][0] = fe->get_dzetadx()[0];
              data.local_transform[2][1] = fe->get_dzetady()[0];
              data.local_transform[2][2] = fe->get_dzetadz()[0];
              data.local_transform[1][2] = fe->get_detadz()[0];
              data.local_transform[0][2] = fe->get_dxidz()[0];
            }
        } // endif data.need_derivative()

      for (unsigned int i=0; i<n_dof; i++)
        {
          // compute base and radial shape indices
          unsigned int i_base, i_radial;
          compute_shape_indices(fet, inf_elem, i, i_base, i_radial);

          data.shape[i] = (InfFERadial::decay(Dim,v)                                    /* (1.-v)/2. in 3D          */
                           * FEInterface::shape(fet, base_el.get(), i_base, p)          /* S_n(s,t)                 */
                           * InfFE<Dim,T_radial,T_map>::eval(v, o_radial, i_radial))    /* L_n(v)                   */
                           * time_harmonic;                                             /* e^(i*k*phase(s,t,v)      */

          // use differentiation of the above equation
          if (data.need_derivative())
            {
              data.dshape[i](0)     = (InfFERadial::decay(Dim,v)
                                       * FEInterface::shape_deriv(fet, base_el.get(), i_base, 0, p)
                                       * InfFE<Dim,T_radial,T_map>::eval(v, o_radial, i_radial))
                                       * time_harmonic;

              if (Dim > 2)
                {
                  data.dshape[i](1)   = (InfFERadial::decay(Dim,v)
                                         * FEInterface::shape_deriv(fet, base_el.get(), i_base, 1, p)
                                         * InfFE<Dim,T_radial,T_map>::eval(v, o_radial, i_radial))
                                         * time_harmonic;

                }
              data.dshape[i](Dim-1)  = (InfFERadial::decay_deriv(Dim, v) * InfFE<Dim,T_radial,T_map>::eval(v, o_radial, i_radial)
                                        +InfFERadial::decay(Dim,v) * InfFE<Dim,T_radial,T_map>::eval_deriv(v, o_radial, i_radial))
                                        * FEInterface::shape(fet, base_el.get(), i_base, p) * time_harmonic;

#ifdef LIBMESH_USE_COMPLEX_NUMBERS
              // derivative of time_harmonic (works for harmonic behavior only):
              data.dshape[i](Dim-1)+= data.shape[i]*imaginary*wavenumber
                                      *interpolated_dist*InfFE<Dim,INFINITE_MAP,T_map>::eval_deriv(v, radial_mapping_order, 1);

#else
            /*
             * The gradient in infinite elements is dominated by the contribution due to the oscillating phase.
             * Since this term is imaginary, I think there is no means to look at it without having complex numbers.
             */
            libmesh_not_implemented();
            // Maybe we can solve it with a warning as well, but I think one really should not do this...
#endif
            }
        }
    }

  else
    libmesh_error_msg("compute_data() for 1-dimensional InfFE not implemented.");
}


template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
Real InfFE<Dim,T_radial,T_map>::shape_deriv(const FEType fet,
                                            const Elem * inf_elem,
                                            const unsigned int i,
                                            const unsigned int j,
                                            const Point & p,
                                            const bool add_p_level)
{
  if (add_p_level)
    {
      FEType tmp_fet=fet;
      tmp_fet = static_cast<Order>(fet.order + add_p_level * inf_elem->p_level());
      return InfFE<Dim,T_radial,T_map>::shape_deriv(tmp_fet, inf_elem, i, j, p);
    }
  return InfFE<Dim,T_radial,T_map>::shape_deriv(fet, inf_elem, i, j, p);
}






template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
void InfFE<Dim,T_radial,T_map>::compute_node_indices (const ElemType inf_elem_type,
                                                      const unsigned int outer_node_index,
                                                      unsigned int & base_node,
                                                      unsigned int & radial_node)
{
  switch (inf_elem_type)
    {
    case INFEDGE2:
      {
        libmesh_assert_less (outer_node_index, 2);
        base_node   = 0;
        radial_node = outer_node_index;
        return;
      }


      // linear base approximation, easy to determine
    case INFQUAD4:
      {
        libmesh_assert_less (outer_node_index, 4);
        base_node   = outer_node_index % 2;
        radial_node = outer_node_index / 2;
        return;
      }

    case INFPRISM6:
      {
        libmesh_assert_less (outer_node_index, 6);
        base_node   = outer_node_index % 3;
        radial_node = outer_node_index / 3;
        return;
      }

    case INFHEX8:
      {
        libmesh_assert_less (outer_node_index, 8);
        base_node   = outer_node_index % 4;
        radial_node = outer_node_index / 4;
        return;
      }


      // higher order base approximation, more work necessary
    case INFQUAD6:
      {
        switch (outer_node_index)
          {
          case 0:
          case 1:
            {
              radial_node = 0;
              base_node   = outer_node_index;
              return;
            }

          case 2:
          case 3:
            {
              radial_node = 1;
              base_node   = outer_node_index-2;
              return;
            }

          case 4:
            {
              radial_node = 0;
              base_node   = 2;
              return;
            }

          case 5:
            {
              radial_node = 1;
              base_node   = 2;
              return;
            }

          default:
            libmesh_error_msg("Unrecognized outer_node_index = " << outer_node_index);
          }
      }


    case INFHEX16:
    case INFHEX18:
      {
        switch (outer_node_index)
          {
          case 0:
          case 1:
          case 2:
          case 3:
            {
              radial_node = 0;
              base_node   = outer_node_index;
              return;
            }

          case 4:
          case 5:
          case 6:
          case 7:
            {
              radial_node = 1;
              base_node   = outer_node_index-4;
              return;
            }

          case 8:
          case 9:
          case 10:
          case 11:
            {
              radial_node = 0;
              base_node   = outer_node_index-4;
              return;
            }

          case 12:
          case 13:
          case 14:
          case 15:
            {
              radial_node = 1;
              base_node   = outer_node_index-8;
              return;
            }

          case 16:
            {
              libmesh_assert_equal_to (inf_elem_type, INFHEX18);
              radial_node = 0;
              base_node   = 8;
              return;
            }

          case 17:
            {
              libmesh_assert_equal_to (inf_elem_type, INFHEX18);
              radial_node = 1;
              base_node   = 8;
              return;
            }

          default:
            libmesh_error_msg("Unrecognized outer_node_index = " << outer_node_index);
          }
      }


    case INFPRISM12:
      {
        switch (outer_node_index)
          {
          case 0:
          case 1:
          case 2:
            {
              radial_node = 0;
              base_node   = outer_node_index;
              return;
            }

          case 3:
          case 4:
          case 5:
            {
              radial_node = 1;
              base_node   = outer_node_index-3;
              return;
            }

          case 6:
          case 7:
          case 8:
            {
              radial_node = 0;
              base_node   = outer_node_index-3;
              return;
            }

          case 9:
          case 10:
          case 11:
            {
              radial_node = 1;
              base_node   = outer_node_index-6;
              return;
            }

          default:
            libmesh_error_msg("Unrecognized outer_node_index = " << outer_node_index);
          }
      }


    default:
      libmesh_error_msg("ERROR: Bad infinite element type=" << inf_elem_type << ", node=" << outer_node_index);
    }
}






template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
void InfFE<Dim,T_radial,T_map>::compute_node_indices_fast (const ElemType inf_elem_type,
                                                           const unsigned int outer_node_index,
                                                           unsigned int & base_node,
                                                           unsigned int & radial_node)
{
  libmesh_assert_not_equal_to (inf_elem_type, INVALID_ELEM);

  static std::vector<unsigned int> _static_base_node_index;
  static std::vector<unsigned int> _static_radial_node_index;

  /*
   * fast counterpart to compute_node_indices(), uses local static buffers
   * to store the index maps.  The class member
   * \p _compute_node_indices_fast_current_elem_type remembers
   * the current element type.
   *
   * Note that there exist non-static members storing the
   * same data.  However, you never know what element type
   * is currently used by the \p InfFE object, and what
   * request is currently directed to the static \p InfFE
   * members (which use \p compute_node_indices_fast()).
   * So separate these.
   *
   * check whether the work for this elemtype has already
   * been done.  If so, use this index.  Otherwise, refresh
   * the buffer to this element type.
   */
  if (inf_elem_type==_compute_node_indices_fast_current_elem_type)
    {
      base_node   = _static_base_node_index  [outer_node_index];
      radial_node = _static_radial_node_index[outer_node_index];
      return;
    }
  else
    {
      // store the map for _all_ nodes for this element type
      _compute_node_indices_fast_current_elem_type = inf_elem_type;

      unsigned int n_nodes = libMesh::invalid_uint;

      switch (inf_elem_type)
        {
        case INFEDGE2:
          {
            n_nodes = 2;
            break;
          }
        case INFQUAD4:
          {
            n_nodes = 4;
            break;
          }
        case INFQUAD6:
          {
            n_nodes = 6;
            break;
          }
        case INFHEX8:
          {
            n_nodes = 8;
            break;
          }
        case INFHEX16:
          {
            n_nodes = 16;
            break;
          }
        case INFHEX18:
          {
            n_nodes = 18;
            break;
          }
        case INFPRISM6:
          {
            n_nodes = 6;
            break;
          }
        case INFPRISM12:
          {
            n_nodes = 12;
            break;
          }
        default:
          libmesh_error_msg("ERROR: Bad infinite element type=" << inf_elem_type << ", node=" << outer_node_index);
        }


      _static_base_node_index.resize  (n_nodes);
      _static_radial_node_index.resize(n_nodes);

      for (unsigned int n=0; n<n_nodes; n++)
        compute_node_indices (inf_elem_type,
                              n,
                              _static_base_node_index  [outer_node_index],
                              _static_radial_node_index[outer_node_index]);

      // and return for the specified node
      base_node   = _static_base_node_index  [outer_node_index];
      radial_node = _static_radial_node_index[outer_node_index];
      return;
    }
}






template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
void InfFE<Dim,T_radial,T_map>::compute_shape_indices (const FEType & fet,
                                                       const Elem * inf_elem,
                                                       const unsigned int i,
                                                       unsigned int & base_shape,
                                                       unsigned int & radial_shape)
{
  // (Temporarily) call version of this function taking an
  // ElemType. Eventually there should only be one version of this
  // function that takes an Elem*.
  compute_shape_indices(fet, inf_elem->type(), i, base_shape, radial_shape);
}



template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
void InfFE<Dim,T_radial,T_map>::compute_shape_indices (const FEType & fet,
                                                       const ElemType inf_elem_type,
                                                       const unsigned int i,
                                                       unsigned int & base_shape,
                                                       unsigned int & radial_shape)
{
  // TODO: eventually figure out a way to deprecated this
  // function. Note that we can't go the other way around and have
  // this function call the Elem* version because there's not really a
  // clean way to create the required Elem object on the fly...
  // libmesh_deprecated();

  // An example is provided:  the numbers in comments refer to
  // a fictitious InfHex18.  The numbers are chosen as exemplary
  // values.  There is currently no base approximation that
  // requires this many dof's at nodes, sides, faces and in the element.
  //
  // the order of the shape functions is heavily related with the
  // order the dofs are assigned in \p DofMap::distributed_dofs().
  // Due to the infinite elements with higher-order base approximation,
  // some more effort is necessary.
  //
  // numbering scheme:
  // 1. all vertices in the base, assign node->n_comp() dofs to each vertex
  // 2. all vertices further out: innermost loop: radial shapes,
  //    then the base approximation shapes
  // 3. all side nodes in the base, assign node->n_comp() dofs to each side node
  // 4. all side nodes further out: innermost loop: radial shapes,
  //    then the base approximation shapes
  // 5. (all) face nodes in the base, assign node->n_comp() dofs to each face node
  // 6. (all) face nodes further out: innermost loop: radial shapes,
  //    then the base approximation shapes
  // 7. element-associated dof in the base
  // 8. element-associated dof further out

  const unsigned int radial_order       = static_cast<unsigned int>(fet.radial_order.get_order()); // 4
  const unsigned int radial_order_p_one = radial_order+1;                                          // 5

  const ElemType base_elem_type           (InfFEBase::get_elem_type(inf_elem_type));               // QUAD9

  // assume that the number of dof is the same for all vertices
  unsigned int n_base_vertices         = libMesh::invalid_uint;                                    // 4
  const unsigned int n_base_vertex_dof = FEInterface::n_dofs_at_node  (Dim-1, fet, base_elem_type, 0);// 2

  unsigned int n_base_side_nodes       = libMesh::invalid_uint;                                    // 4
  unsigned int n_base_side_dof         = libMesh::invalid_uint;                                    // 3

  unsigned int n_base_face_nodes       = libMesh::invalid_uint;                                    // 1
  unsigned int n_base_face_dof         = libMesh::invalid_uint;                                    // 5

  const unsigned int n_base_elem_dof   = FEInterface::n_dofs_per_elem (Dim-1, fet, base_elem_type);// 9


  switch (inf_elem_type)
    {
    case INFEDGE2:
      {
        n_base_vertices   = 1;
        n_base_side_nodes = 0;
        n_base_face_nodes = 0;
        n_base_side_dof   = 0;
        n_base_face_dof   = 0;
        break;
      }

    case INFQUAD4:
      {
        n_base_vertices   = 2;
        n_base_side_nodes = 0;
        n_base_face_nodes = 0;
        n_base_side_dof   = 0;
        n_base_face_dof   = 0;
        break;
      }

    case INFQUAD6:
      {
        n_base_vertices   = 2;
        n_base_side_nodes = 1;
        n_base_face_nodes = 0;
        n_base_side_dof   = FEInterface::n_dofs_at_node (Dim-1, fet,base_elem_type, n_base_vertices);
        n_base_face_dof   = 0;
        break;
      }

    case INFHEX8:
      {
        n_base_vertices   = 4;
        n_base_side_nodes = 0;
        n_base_face_nodes = 0;
        n_base_side_dof   = 0;
        n_base_face_dof   = 0;
        break;
      }

    case INFHEX16:
      {
        n_base_vertices   = 4;
        n_base_side_nodes = 4;
        n_base_face_nodes = 0;
        n_base_side_dof   = FEInterface::n_dofs_at_node (Dim-1, fet,base_elem_type, n_base_vertices);
        n_base_face_dof   = 0;
        break;
      }

    case INFHEX18:
      {
        n_base_vertices   = 4;
        n_base_side_nodes = 4;
        n_base_face_nodes = 1;
        n_base_side_dof   = FEInterface::n_dofs_at_node (Dim-1, fet,base_elem_type, n_base_vertices);
        n_base_face_dof   = FEInterface::n_dofs_at_node (Dim-1, fet,base_elem_type, 8);
        break;
      }


    case INFPRISM6:
      {
        n_base_vertices   = 3;
        n_base_side_nodes = 0;
        n_base_face_nodes = 0;
        n_base_side_dof   = 0;
        n_base_face_dof   = 0;
        break;
      }

    case INFPRISM12:
      {
        n_base_vertices   = 3;
        n_base_side_nodes = 3;
        n_base_face_nodes = 0;
        n_base_side_dof   = FEInterface::n_dofs_at_node (Dim-1, fet,base_elem_type, n_base_vertices);
        n_base_face_dof   = 0;
        break;
      }

    default:
      libmesh_error_msg("Unrecognized inf_elem_type = " << inf_elem_type);
    }


  {
    // these are the limits describing the intervals where the shape function lies
    const unsigned int n_dof_at_base_vertices = n_base_vertices*n_base_vertex_dof;                 // 8
    const unsigned int n_dof_at_all_vertices  = n_dof_at_base_vertices*radial_order_p_one;         // 40

    const unsigned int n_dof_at_base_sides    = n_base_side_nodes*n_base_side_dof;                 // 12
    const unsigned int n_dof_at_all_sides     = n_dof_at_base_sides*radial_order_p_one;            // 60

    const unsigned int n_dof_at_base_face     = n_base_face_nodes*n_base_face_dof;                 // 5
    const unsigned int n_dof_at_all_faces     = n_dof_at_base_face*radial_order_p_one;             // 25


    // start locating the shape function
    if (i < n_dof_at_base_vertices)                                              // range of i: 0..7
      {
        // belongs to vertex in the base
        radial_shape = 0;
        base_shape   = i;
      }

    else if (i < n_dof_at_all_vertices)                                          // range of i: 8..39
      {
        /* belongs to vertex in the outer shell
         *
         * subtract the number of dof already counted,
         * so that i_offset contains only the offset for the base
         */
        const unsigned int i_offset = i - n_dof_at_base_vertices;                // 0..31

        // first the radial dof are counted, then the base dof
        radial_shape = (i_offset % radial_order) + 1;
        base_shape   = i_offset / radial_order;
      }

    else if (i < n_dof_at_all_vertices+n_dof_at_base_sides)                      // range of i: 40..51
      {
        // belongs to base, is a side node
        radial_shape = 0;
        base_shape = i - radial_order * n_dof_at_base_vertices;                  //  8..19
      }

    else if (i < n_dof_at_all_vertices+n_dof_at_all_sides)                       // range of i: 52..99
      {
        // belongs to side node in the outer shell
        const unsigned int i_offset = i - (n_dof_at_all_vertices
                                           + n_dof_at_base_sides);               // 0..47
        radial_shape = (i_offset % radial_order) + 1;
        base_shape   = (i_offset / radial_order) + n_dof_at_base_vertices;
      }

    else if (i < n_dof_at_all_vertices+n_dof_at_all_sides+n_dof_at_base_face)    // range of i: 100..104
      {
        // belongs to the node in the base face
        radial_shape = 0;
        base_shape = i - radial_order*(n_dof_at_base_vertices
                                       + n_dof_at_base_sides);                   //  20..24
      }

    else if (i < n_dof_at_all_vertices+n_dof_at_all_sides+n_dof_at_all_faces)    // range of i: 105..124
      {
        // belongs to the node in the outer face
        const unsigned int i_offset = i - (n_dof_at_all_vertices
                                           + n_dof_at_all_sides
                                           + n_dof_at_base_face);                // 0..19
        radial_shape = (i_offset % radial_order) + 1;
        base_shape   = (i_offset / radial_order) + n_dof_at_base_vertices + n_dof_at_base_sides;
      }

    else if (i < n_dof_at_all_vertices+n_dof_at_all_sides+n_dof_at_all_faces+n_base_elem_dof)      // range of i: 125..133
      {
        // belongs to the base and is an element associated shape
        radial_shape = 0;
        base_shape = i - (n_dof_at_all_vertices
                          + n_dof_at_all_sides
                          + n_dof_at_all_faces);                                 // 0..8
      }

    else                                                                         // range of i: 134..169
      {
        libmesh_assert_less (i, n_dofs(fet, inf_elem_type));
        // belongs to the outer shell and is an element associated shape
        const unsigned int i_offset = i - (n_dof_at_all_vertices
                                           + n_dof_at_all_sides
                                           + n_dof_at_all_faces
                                           + n_base_elem_dof);                   // 0..19
        radial_shape = (i_offset % radial_order) + 1;
        base_shape   = (i_offset / radial_order) + n_dof_at_base_vertices + n_dof_at_base_sides + n_dof_at_base_face;
      }
  }

  return;
}


#ifdef LIBMESH_ENABLE_AMR
#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS

template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
void InfFE<Dim, T_radial, T_map>::inf_compute_node_constraints (NodeConstraints & constraints,const Elem * elem)
{
  // only constrain elements in 2,3d.
  if (Dim == 1)
    return;

  libmesh_assert(child_elem);

  // only constrain active and ancestor elements
  if (child_elem->subactive())
    return;

  // for infinite elements, the computation of constraints is somewhat different
  // than for Lagrange elements:
  // 1) Only the base element (i.e. side(0) ) may be refined.
  //    Thus, in radial direction no constraints must be considered.
  // 2) Due to the tensorial structure of shape functions (base_shape * radial_function),
  //    it must be ensured that all element DOFs inherit that constraint.
  // Consequently, the constraints are computed on the base (baseh_shape) but must
  // be applied to all DOFs with the respective base_shape index (i.e. for all radial_functions).
  //
  // FIXME: In the current form, this function does not work for infinite elements
  //        because constraining the non-base points requires knowledge of the T_map and T_radial
  //        parameters; but they are not accessible via the element and may differ between variables.
  //
  // For the moment being, we just check if this element can be skipped and fail otherwise.

  // if one of the sides needs a constraint, an error is thrown.
  // In other cases, we leave the function regularly.
  for (auto s : elem->side_index_range())
    {
      if (elem->neighbor_ptr(s) != nullptr &&
          elem->neighbor_ptr(s) != remote_elem)
        if (elem->neighbor_ptr(s)->level() < elem->level())
          {
            libmesh_not_implemented();
          }
    }
}

#endif //LIBMESH_ENABLE_NODE_CONSTRAINTS


template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
void InfFE<Dim, T_radial, T_map>::inf_compute_constraints (DofConstraints & constraints,
                                                           DofMap & dof_map,
                                                           const unsigned int variable_number,
                                                           const Elem * child_elem)
{

  // only constrain elements in 2,3d.
  if (Dim == 1)
    return;

  libmesh_assert(child_elem);

  // only constrain active and ancestor elements
  if (child_elem->subactive())
    return;

  // Before we start to compute anything, lets check if any confinement is needed:
  bool need_constraints=false;
  for (auto child_neighbor : child_elem->neighbor_ptr_range())
    if (child_neighbor->level() < child_elem->level())
      {
        need_constraints = true;
        break;
      }
  if (!need_constraints)
    return;

  // For infinite elements, the computation of constraints is somewhat different
  // than for Lagrange elements:
  // 1) When an infinite element is refined, only the base element (i.e. side(0) ) is refined.
  //
  // 2) Due to the tensorial structure of shape functions (base_shape * radial_function),
  //    it must be ensured that all element DOFs inherit that constraint.
  //    It is important here to distinguish the (total) DOF from base DOF and radial DOF contributions.
  //
  // 3) Due to the generality of radial polynomial (of type fe_type.radial_family and with order fe_type.radial_order)
  //    here basis functions cannot be mapped to nodes: Independent from the radial polynomial,
  //    infinite elements have one set of nodes at the base (side(0)) and a second set at twice the distance to their origin.
  //
  //    Independent from the polynomial and degree used, the first radial DOF is 1 at the base while all others are 0 there
  //
  //Constraining of DOFs is only needed when a DOF is nonzero at the elements face shared with a coarser element.
  // Thus, the following scheme is used here:
  //
  //  -If the coarser element is the neighbor(0) (i.e. we share only the base), we must constrain
  //   all DOFs that correspond to the first order radial contribution.
  //  -if an infinite neighbor is coarser (than 'child_elem'), all those DOFs must be constrained
  //   whose contribution from the base is non-zero at the interface.
  //   In this case, we lack a point-assignement between DOFs and nodes, but since there is no refinement in radial direction,
  //   the radial polynomials coincide on neighboring elements.
  //   Thus, if one constraines these DOFs at one (arbitrary) point correctly, they match for each point along the radial direction.
  //   Hence, we constrain them with the same values as those DOFs belonging to the first order polynomial, obtaining consistent
  //   constraints that mimick constraints that are computed at the support points for each radial polynomial contribution.

  FEType fe_type = dof_map.variable_type(variable_number);

  libmesh_assert(fe_type.family == LAGRANGE);

  std::vector<dof_id_type> child_base_dof_indices, parent_base_dof_indices;
  std::vector<dof_id_type> child_elem_dof_indices, parent_elem_dof_indices;

  const Elem * parent_elem = child_elem->parent();

  // This can't happen...  Only level-0 elements have nullptr
  // parents, and no level-0 elements can be at a higher
  // level than their neighbors!
  libmesh_assert(parent_elem);

  dof_map.dof_indices (child_elem, child_elem_dof_indices,
                       variable_number);
  dof_map.dof_indices (parent_elem, parent_elem_dof_indices,
                       variable_number);

  const unsigned int n_total_dofs = child_elem_dof_indices.size();
  // fill the elements shape index map: we will have to use it later
  // to find the elements dofs that correspond to certain base_elem_dofs.
  std::vector<unsigned int> radial_shape_index(n_total_dofs);
  std::vector<unsigned int> base_shape_index(n_total_dofs);
  // fill the shape index map
#ifdef DEBUG
  unsigned int max_base_id=0;
  unsigned int max_radial_id=0;
#endif
  for (unsigned int n=0; n<n_total_dofs; ++n)
    {
      compute_shape_indices (fe_type,
                             child_elem,
                             n,
                             base_shape_index[n],
                             radial_shape_index[n]);

#ifdef DEBUG
      if (base_shape_index[n] > max_base_id)
        max_base_id = base_shape_index[n];
      if (radial_shape_index[n] > max_radial_id)
        max_radial_id = radial_shape_index[n];
#endif
    }

#ifdef DEBUG
  libmesh_assert_equal_to( (max_base_id+1)*(max_radial_id+1), n_total_dofs );
#endif

  for (auto s : child_elem->side_index_range())
    if (child_elem->neighbor_ptr(s) != nullptr &&
        child_elem->neighbor_ptr(s) != remote_elem)
      if (child_elem->neighbor_ptr(s)->level() < child_elem->level())
        {
          // we ALWAYS take the base element for reference:
          // - For s=0, we refine all dofs with `radial_shape_index == 0
          // - for s>0, we refine all dofs whose corresponding base_shape has its support point shared with neighbor(s)
          std::unique_ptr<const Elem> child_base, parent_base;
          child_elem->build_side_ptr(child_base, 0);
          parent_elem->build_side_ptr(parent_base, 0);

          const unsigned int n_base_dofs =
            FEInterface::n_dofs(fe_type, child_base.get());

          // We need global DOF indices for both base and 'full' elements
          dof_map.dof_indices (child_base.get(), child_base_dof_indices,
                               variable_number);
          dof_map.dof_indices (parent_base.get(), parent_base_dof_indices,
                               variable_number);


          // First we loop over the childs base DOFs (nodes) and check which of them needs constraint
          // and which can be skipped.
          for (unsigned int child_base_dof=0; child_base_dof != n_base_dofs; ++child_base_dof)
            {
              libmesh_assert_less (child_base_dof, child_base->n_nodes());

              // Childs global dof index.
              const dof_id_type child_base_dof_g = child_base_dof_indices[child_base_dof];

              // Hunt for "constraining against myself" cases before
              // we bother creating a constraint row
              bool self_constraint = false;
              for (unsigned int parent_base_dof=0;
                   parent_base_dof != n_base_dofs; parent_base_dof++)
                {
                  libmesh_assert_less (parent_base_dof, parent_base->n_nodes());

                  // Their global dof index.
                  const dof_id_type parent_base_dof_g =
                    parent_base_dof_indices[parent_base_dof];

                  if (parent_base_dof_g == child_base_dof_g)
                    {
                      self_constraint = true;
                      break;
                    }
                }

              if (self_constraint)
                continue;

              // now we need to constrain all __child_elem__ DOFs whose base corresponds to
              // child_base_dof.
              //  --> loop over all child_elem dofs whose base_shape_index == child_base_dof
              unsigned int n_elem_dofs = FEInterface::n_dofs(fe_type, child_elem);
              libmesh_assert_equal_to(n_elem_dofs, n_total_dofs);
              for(unsigned int child_elem_dof=0; child_elem_dof != n_elem_dofs; ++child_elem_dof)
                {
                  if (base_shape_index[child_elem_dof] != child_base_dof)
                    continue;

                  // independent from the radial description, the first radial DOF is 1 at the base
                  // while all others start with 0.
                  // Thus, to confine for the bases neighbor, we only need to refine DOFs that correspond
                  // to the first radial DOF
                  if (s==0)
                    {
                      if (radial_shape_index[child_elem_dof] > 0)
                        continue;
                    }
                  else
                    {
                      // If the neighbor is not the base, we must check now if the support point of the dof
                      // is actually shared with that neighbor:
                      if ( !child_elem->neighbor_ptr(s)->contains_point(child_base->point(child_base_dof)) )
                        continue;
                    }


                  const dof_id_type child_elem_dof_g = child_elem_dof_indices[child_elem_dof];

                  DofConstraintRow * constraint_row;

                  // we may be running constraint methods concurrently
                  // on multiple threads, so we need a lock to
                  // ensure that this constraint is "ours"
                  {
                    Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);

                    if (dof_map.is_constrained_dof(child_elem_dof_g))
                      continue;

                    constraint_row = &(constraints[child_elem_dof_g]);
                    libmesh_assert(constraint_row->empty());
                  }

                  // The support point of the DOF
                  const Point & support_point = child_base->point(child_base_dof);

                  // Figure out where my (base) node lies on the parents reference element.
                  const Point mapped_point = FEMap::inverse_map(Dim-1,
                                                                parent_base.get(),
                                                                support_point);

                  // now we need the parents base DOFs, evaluated at the mapped_point for refinement:
                  for (unsigned int parent_base_dof=0;
                       parent_base_dof != n_base_dofs; parent_base_dof++)
                    {

                      const Real parent_base_dof_value = FEInterface::shape(Dim-1,
                                                                            fe_type,
                                                                            parent_base.get(),
                                                                            parent_base_dof,
                                                                            mapped_point);


                      // all parent elements DOFs whose base_index corresponds to parent_base_dof
                      //  must be constrained with the parent_base_dof_value.

                      // The value of the radial function does not play a role here:
                      // 1) only the function with radial_shape_index[] == 0 are 1 at the base,
                      //    the others are 0.
                      // 2) The radial basis is (usually) not a Lagrange polynomial.
                      //    Thus, constraining according to a support point doesn't work.
                      //    However, they reach '1' at a certain (radial) distance which is the same for parent and child.
                      for (unsigned int parent_elem_dof=0;
                           parent_elem_dof != n_elem_dofs; parent_elem_dof++)
                        {
                          if (base_shape_index[parent_elem_dof] != parent_base_dof)
                            continue;

                          // only constrain with coinciding radial DOFs.
                          // Otherwise, we start coupling all DOFs with each other and end up in a mess.
                          if (radial_shape_index[parent_elem_dof] != radial_shape_index[child_elem_dof])
                            continue;

                          // Their global dof index.
                          const dof_id_type parent_elem_dof_g =
                            parent_elem_dof_indices[parent_elem_dof];

                          // Only add non-zero and non-identity values
                          // for Lagrange basis functions. (parent_base is assumed to be of Lagrange-type).
                          if ((std::abs(parent_base_dof_value) > 1.e-5) &&
                              (std::abs(parent_base_dof_value) < .999))
                            {
                              constraint_row->emplace(parent_elem_dof_g, parent_base_dof_value);
                            }
#ifdef DEBUG
                          // Protect for the case u_i = 0.999 u_j,
                          // in which case i better equal j.
                          else if (parent_base_dof_value >= .999)
                            {
                              libmesh_assert_equal_to (child_base_dof_g, parent_base_dof_indices[parent_base_dof]);
                              libmesh_assert_equal_to (child_elem_dof_g, parent_elem_dof_g);
                            }
#endif
                        }

                    }
                }

            }
        }
}

#endif // LIBMESH_ENABLE_AMR


//--------------------------------------------------------------
// Explicit instantiations
//#include "libmesh/inf_fe_instantiate_1D.h"
//#include "libmesh/inf_fe_instantiate_2D.h"
//#include "libmesh/inf_fe_instantiate_3D.h"

INSTANTIATE_INF_FE_MBRF(1, CARTESIAN, unsigned int, n_dofs(const FEType &, const ElemType));
INSTANTIATE_INF_FE_MBRF(2, CARTESIAN, unsigned int, n_dofs(const FEType &, const ElemType));
INSTANTIATE_INF_FE_MBRF(3, CARTESIAN, unsigned int, n_dofs(const FEType &, const ElemType));
INSTANTIATE_INF_FE_MBRF(1, CARTESIAN, unsigned int, n_dofs(const FEType &, const Elem*));
INSTANTIATE_INF_FE_MBRF(2, CARTESIAN, unsigned int, n_dofs(const FEType &, const Elem*));
INSTANTIATE_INF_FE_MBRF(3, CARTESIAN, unsigned int, n_dofs(const FEType &, const Elem*));
INSTANTIATE_INF_FE_MBRF(1, CARTESIAN, unsigned int, n_dofs_per_elem(const FEType &, const ElemType));
INSTANTIATE_INF_FE_MBRF(2, CARTESIAN, unsigned int, n_dofs_per_elem(const FEType &, const ElemType));
INSTANTIATE_INF_FE_MBRF(3, CARTESIAN, unsigned int, n_dofs_per_elem(const FEType &, const ElemType));
INSTANTIATE_INF_FE_MBRF(1, CARTESIAN, unsigned int, n_dofs_per_elem(const FEType &, const Elem *));
INSTANTIATE_INF_FE_MBRF(2, CARTESIAN, unsigned int, n_dofs_per_elem(const FEType &, const Elem *));
INSTANTIATE_INF_FE_MBRF(3, CARTESIAN, unsigned int, n_dofs_per_elem(const FEType &, const Elem *));
INSTANTIATE_INF_FE_MBRF(1, CARTESIAN, unsigned int, n_dofs_at_node(const FEType &, const ElemType, const unsigned int));
INSTANTIATE_INF_FE_MBRF(2, CARTESIAN, unsigned int, n_dofs_at_node(const FEType &, const ElemType, const unsigned int));
INSTANTIATE_INF_FE_MBRF(3, CARTESIAN, unsigned int, n_dofs_at_node(const FEType &, const ElemType, const unsigned int));
INSTANTIATE_INF_FE_MBRF(1, CARTESIAN, unsigned int, n_dofs_at_node(const FEType &, const Elem *, const unsigned int));
INSTANTIATE_INF_FE_MBRF(2, CARTESIAN, unsigned int, n_dofs_at_node(const FEType &, const Elem *, const unsigned int));
INSTANTIATE_INF_FE_MBRF(3, CARTESIAN, unsigned int, n_dofs_at_node(const FEType &, const Elem *, const unsigned int));
INSTANTIATE_INF_FE_MBRF(1, CARTESIAN, void, compute_shape_indices(const FEType &, const ElemType, const unsigned int, unsigned int &, unsigned int &));
INSTANTIATE_INF_FE_MBRF(2, CARTESIAN, void, compute_shape_indices(const FEType &, const ElemType, const unsigned int, unsigned int &, unsigned int &));
INSTANTIATE_INF_FE_MBRF(3, CARTESIAN, void, compute_shape_indices(const FEType &, const ElemType, const unsigned int, unsigned int &, unsigned int &));
INSTANTIATE_INF_FE_MBRF(1, CARTESIAN, void, compute_shape_indices(const FEType &, const Elem *, const unsigned int, unsigned int &, unsigned int &));
INSTANTIATE_INF_FE_MBRF(2, CARTESIAN, void, compute_shape_indices(const FEType &, const Elem *, const unsigned int, unsigned int &, unsigned int &));
INSTANTIATE_INF_FE_MBRF(3, CARTESIAN, void, compute_shape_indices(const FEType &, const Elem *, const unsigned int, unsigned int &, unsigned int &));
INSTANTIATE_INF_FE_MBRF(1, CARTESIAN, void, compute_node_indices(const ElemType, const unsigned int, unsigned int &, unsigned int &));
INSTANTIATE_INF_FE_MBRF(2, CARTESIAN, void, compute_node_indices(const ElemType, const unsigned int, unsigned int &, unsigned int &));
INSTANTIATE_INF_FE_MBRF(3, CARTESIAN, void, compute_node_indices(const ElemType, const unsigned int, unsigned int &, unsigned int &));
INSTANTIATE_INF_FE_MBRF(1, CARTESIAN, Real, shape(const FEType &, const Elem *, const unsigned int, const Point & p));
INSTANTIATE_INF_FE_MBRF(2, CARTESIAN, Real, shape(const FEType &, const Elem *, const unsigned int, const Point & p));
INSTANTIATE_INF_FE_MBRF(3, CARTESIAN, Real, shape(const FEType &, const Elem *, const unsigned int, const Point & p));
INSTANTIATE_INF_FE_MBRF(1, CARTESIAN, Real, shape(const FEType, const Elem *, const unsigned int, const Point &, const bool));
INSTANTIATE_INF_FE_MBRF(2, CARTESIAN, Real, shape(const FEType, const Elem *, const unsigned int, const Point &, const bool));
INSTANTIATE_INF_FE_MBRF(3, CARTESIAN, Real, shape(const FEType, const Elem *, const unsigned int, const Point &, const bool));
INSTANTIATE_INF_FE_MBRF(1, CARTESIAN, Real, shape(const FEType &, const ElemType, const unsigned int, const Point &));
INSTANTIATE_INF_FE_MBRF(2, CARTESIAN, Real, shape(const FEType &, const ElemType, const unsigned int, const Point &));
INSTANTIATE_INF_FE_MBRF(3, CARTESIAN, Real, shape(const FEType &, const ElemType, const unsigned int, const Point &));
INSTANTIATE_INF_FE_MBRF(1, CARTESIAN, Real, shape_deriv(const FEType &, const Elem *, const unsigned int, const unsigned int, const Point &));
INSTANTIATE_INF_FE_MBRF(2, CARTESIAN, Real, shape_deriv(const FEType &, const Elem *, const unsigned int, const unsigned int, const Point &));
INSTANTIATE_INF_FE_MBRF(3, CARTESIAN, Real, shape_deriv(const FEType &, const Elem *, const unsigned int, const unsigned int, const Point &));
INSTANTIATE_INF_FE_MBRF(1, CARTESIAN, Real, shape_deriv(const FEType &, const ElemType, const unsigned int, const unsigned int, const Point &));
INSTANTIATE_INF_FE_MBRF(2, CARTESIAN, Real, shape_deriv(const FEType &, const ElemType, const unsigned int, const unsigned int, const Point &));
INSTANTIATE_INF_FE_MBRF(3, CARTESIAN, Real, shape_deriv(const FEType &, const ElemType, const unsigned int, const unsigned int, const Point &));
INSTANTIATE_INF_FE_MBRF(1, CARTESIAN, Real, shape_deriv(const FEType, const Elem *, const unsigned int, const unsigned int, const Point &, const bool));
INSTANTIATE_INF_FE_MBRF(2, CARTESIAN, Real, shape_deriv(const FEType, const Elem *, const unsigned int, const unsigned int, const Point &, const bool));
INSTANTIATE_INF_FE_MBRF(3, CARTESIAN, Real, shape_deriv(const FEType, const Elem *, const unsigned int, const unsigned int, const Point &, const bool));
INSTANTIATE_INF_FE_MBRF(1, CARTESIAN, void, compute_data(const FEType &, const Elem *, FEComputeData &));
INSTANTIATE_INF_FE_MBRF(2, CARTESIAN, void, compute_data(const FEType &, const Elem *, FEComputeData &));
INSTANTIATE_INF_FE_MBRF(3, CARTESIAN, void, compute_data(const FEType &, const Elem *, FEComputeData &));
INSTANTIATE_INF_FE_MBRF(1, CARTESIAN, void, nodal_soln(const FEType &, const Elem *, const std::vector<Number> &, std::vector<Number> &));
INSTANTIATE_INF_FE_MBRF(2, CARTESIAN, void, nodal_soln(const FEType &, const Elem *, const std::vector<Number> &, std::vector<Number> &));
INSTANTIATE_INF_FE_MBRF(3, CARTESIAN, void, nodal_soln(const FEType &, const Elem *, const std::vector<Number> &, std::vector<Number> &));
#ifdef LIBMESH_ENABLE_AMR
INSTANTIATE_INF_FE_MBRF(1, CARTESIAN, void, inf_compute_constraints(DofConstraints &, DofMap &, const unsigned int, const Elem *));
INSTANTIATE_INF_FE_MBRF(2, CARTESIAN, void, inf_compute_constraints(DofConstraints &, DofMap &, const unsigned int, const Elem *));
INSTANTIATE_INF_FE_MBRF(3, CARTESIAN, void, inf_compute_constraints(DofConstraints &, DofMap &, const unsigned int, const Elem *));
#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
INSTANTIATE_INF_FE_MBRF(1, CARTESIAN, void, inf_compute_node_constraints(NodeConstraints & constraints,const Elem * elem));
INSTANTIATE_INF_FE_MBRF(2, CARTESIAN, void, inf_compute_node_constraints(NodeConstraints & constraints,const Elem * elem));
INSTANTIATE_INF_FE_MBRF(3, CARTESIAN, void, inf_compute_node_constraints(NodeConstraints & constraints,const Elem * elem));
#endif
#endif

} // namespace libMesh

#endif //ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
