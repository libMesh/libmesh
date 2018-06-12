// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


#include "libmesh/libmesh_config.h"
#if defined(LIBMESH_HAVE_TRIANGLE) && defined(LIBMESH_HAVE_TETGEN)

#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/quadrature_trap.h"
#include "libmesh/quadrature_simpson.h"
#include "libmesh/quadrature_composite.h"
#include "libmesh/elem.h"
#include "libmesh/enum_quadrature_type.h"



namespace libMesh
{


template <class QSubCell>
QuadratureType QComposite<QSubCell>::type() const
{
  return QCOMPOSITE;
}



template <class QSubCell>
QComposite<QSubCell>::QComposite(unsigned int d,
                                 Order o) :
  QSubCell(d,o), // explicitly call base class constructor
  _q_subcell(d,o),
  _lagrange_fe(FEBase::build (d, FEType (FIRST, LAGRANGE)))
{
  // explicitly call the init function in 1D since the
  // other tensor-product rules require this one.
  // note that EDGE will not be used internally, however
  // if we called the function with INVALID_ELEM it would try to
  // be smart and return, thinking it had already done the work.
  if (_dim == 1)
    QSubCell::init(EDGE2);

  libmesh_assert (_lagrange_fe.get() != libmesh_nullptr);

  _lagrange_fe->attach_quadrature_rule (&_q_subcell);
}



template <class QSubCell>
void QComposite<QSubCell>::init (const Elem & elem,
                                 const std::vector<Real> & vertex_distance_func,
                                 unsigned int p_level)
{
  libmesh_assert_equal_to (vertex_distance_func.size(), elem.n_vertices());
  libmesh_assert_equal_to (_dim, elem.dim());

  // if we are not cut, revert to simple base class init() method.
  if (!_elem_cutter.is_cut (elem, vertex_distance_func))
    {
      _q_subcell.init (elem.type(), p_level);
      _points  = _q_subcell.get_points();
      _weights = _q_subcell.get_weights();

      //this->print_info();
      return;
    }

  // Get a pointer to the element's reference element.  We want to
  // perform cutting on the reference element such that the quadrature
  // point locations of the subelements live in the reference
  // coordinate system, thereby eliminating the need for inverse
  // mapping.
  const Elem * reference_elem = elem.reference_elem();

  libmesh_assert (reference_elem != libmesh_nullptr);

  _elem_cutter(*reference_elem, vertex_distance_func);
  //_elem_cutter(elem, vertex_distance_func);

  // clear our state & accumulate points from subelements
  _points.clear();
  _weights.clear();

  // inside subelem
  {
    const std::vector<Elem const *> & inside_elem (_elem_cutter.inside_elements());
    std::cout << inside_elem.size() << " elements inside\n";

    this->add_subelem_values(inside_elem);
  }

  // outside subelem
  {
    const std::vector<Elem const *> & outside_elem (_elem_cutter.outside_elements());
    std::cout << outside_elem.size() << " elements outside\n";

    this->add_subelem_values(outside_elem);
  }

  this->print_info();
}



template <class QSubCell>
void QComposite<QSubCell>::add_subelem_values (const std::vector<Elem const *> & subelem)

{
  const std::vector<Real>  & subelem_weights = _lagrange_fe->get_JxW();
  const std::vector<Point> & subelem_points  = _lagrange_fe->get_xyz();

  for (const auto & elem : subelem)
    {
      // tetgen seems to create 0-volume cells on occasion, but we *should*
      // be catching that appropriately now inside the ElemCutter class.
      // Just in case trap here, describe the error, and abort.
      libmesh_try
        {
          _lagrange_fe->reinit(elem);
          _weights.insert(_weights.end(),
                          subelem_weights.begin(), subelem_weights.end());

          _points.insert(_points.end(),
                         subelem_points.begin(), subelem_points.end());
        }
      libmesh_catch (...)
        {
          libMesh::err << "ERROR: found a bad cut cell!\n";

          for (unsigned int n=0; n<elem->n_nodes(); n++)
            libMesh::err << elem->point(n) << std::endl;

          libmesh_error_msg("Tetgen may have created a 0-volume cell during Cutcell integration.");
        }
    }
}



//--------------------------------------------------------------
// Explicit instantiations
template class QComposite<QGauss>;
template class QComposite<QTrap>;
template class QComposite<QSimpson>;

} // namespace libMesh

#endif // LIBMESH_HAVE_TRIANGLE && LIBMESH_HAVE_TETGEN
