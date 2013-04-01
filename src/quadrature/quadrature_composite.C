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


#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/quadrature_composite.h"
#include "libmesh/elem.h"



namespace libMesh
{


template <class QSubCell>
QComposite<QSubCell>::QComposite(const unsigned int d,
				 const Order o) :
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

  libmesh_assert (_lagrange_fe.get() != NULL);

  _lagrange_fe->attach_quadrature_rule (&_q_subcell);
}



template <class QSubCell>
QComposite<QSubCell>::~QComposite()
{}



template <class QSubCell>
void QComposite<QSubCell>::init (const Elem &elem,
				 const std::vector<Real> &vertex_distance_func,
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

  libmesh_here();

  _elem_cutter(elem, vertex_distance_func);


  // clear our state & accumulate points from subelements
  _points.clear();
  _weights.clear();


  const std::vector<Real>  &subelem_weights = _lagrange_fe->get_JxW();
  const std::vector<Point> &subelem_points  = _lagrange_fe->get_xyz();

  // inside subelem
  {
    const std::vector<Elem const*> &inside_elem (_elem_cutter.inside_elements());
    std::cout << inside_elem.size() << " elements inside\n";

    for (std::vector<Elem const*>::const_iterator it = inside_elem.begin();
	 it!=inside_elem.end(); ++it)
      {
	_lagrange_fe->reinit(*it);

	_weights.insert(_weights.end(),
			subelem_weights.begin(), subelem_weights.end());

	_points.insert(_points.end(),
		       subelem_points.begin(), subelem_points.end());
      }
  }

  // outside subelem
  {
    const std::vector<Elem const*> &outside_elem (_elem_cutter.outside_elements());
    std::cout << outside_elem.size() << " elements outside\n";

    for (std::vector<Elem const*>::const_iterator it = outside_elem.begin();
	 it!=outside_elem.end(); ++it)
      {
	_lagrange_fe->reinit(*it);

	_weights.insert(_weights.end(),
			subelem_weights.begin(), subelem_weights.end());

	_points.insert(_points.end(),
		       subelem_points.begin(), subelem_points.end());
      }
  }

  this->print_info();
}


//--------------------------------------------------------------
// Explicit instantiations
template class QComposite<QGauss>;

} // namespace libMesh
