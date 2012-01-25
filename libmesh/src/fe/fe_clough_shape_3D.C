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


// C++ inlcludes

// Local includes
#include "fe.h"
#include "elem.h"


// FIXME: 3D C1 finite elements are still a work in progress


// Anonymous namespace for persistant variables.
// This allows us to cache the global-to-local mapping transformation
// This should also screw up multithreading royally
namespace
{
  using namespace libMesh;

  static unsigned int old_elem_id = libMesh::invalid_uint;
  // Coefficient naming: d(1)d(2n) is the coefficient of the
  // global shape function corresponding to value 1 in terms of the
  // local shape function corresponding to normal derivative 2

Real clough_raw_shape_second_deriv(const unsigned int basis_num,
                                   const unsigned int deriv_type,
                                   const Point& p);
Real clough_raw_shape_deriv(const unsigned int basis_num,
                            const unsigned int deriv_type,
                            const Point& p);
Real clough_raw_shape(const unsigned int basis_num,
                      const Point& p);


// Compute the static coefficients for an element
void clough_compute_coefs(const Elem* elem)
{
  // Using static globals for old_elem_id, etc. will fail
  // horribly with more than one thread.
  libmesh_assert(libMesh::n_threads() == 1);

  // Coefficients are cached from old elements
  if (elem->id() == old_elem_id)
    return;

  old_elem_id = elem->id();

#if 0
  const Order mapping_order        (elem->default_order());
  const ElemType mapping_elem_type (elem->type());
  const int n_mapping_shape_functions =
    FE<3,LAGRANGE>::n_shape_functions(mapping_elem_type,
				      mapping_order);

  // Degrees of freedom are at vertices and edge midpoints
  std::vector<Point> dofpt;
#endif

}


unsigned char subtriangle_lookup(const Point&)
{
  return 0;
}

  // Return shape function second derivatives on the unit right
  // triangle
Real clough_raw_shape_second_deriv(const unsigned int basis_num,
                            const unsigned int deriv_type,
                            const Point& p)
{
  Real xi = p(0), eta = p(1), zeta = p(2);

  switch (deriv_type)
  {

  // second derivative in xi-xi direction
  case 0:
  switch (basis_num)
    {
      case 0:
        switch (subtriangle_lookup(p))
          {
            case 0:
              break;
          }
    }
  }

  libmesh_error();
  return xi + eta + zeta;
}



Real clough_raw_shape_deriv(const unsigned int basis_num,
                            const unsigned int deriv_type,
                            const Point& p)
{
  Real xi = p(0), eta = p(1), zeta = p(2);

  switch (deriv_type)
  {
  case 0:
  switch (basis_num)
    {
      case 0:
        switch (subtriangle_lookup(p))
          {
            case 0:
              break;
          }
    }
  }

  libmesh_error();
  return xi + eta + zeta;
}

Real clough_raw_shape(const unsigned int basis_num,
                      const Point& p)
{
  Real xi = p(0), eta = p(1), zeta = p(2);

  switch (basis_num)
    {
      case 0:
        switch (subtriangle_lookup(p))
          {
            case 0:
              break;
          }
    }

  libmesh_error();
  return xi + eta + zeta;
}


} // end anonymous namespace


namespace libMesh
{


template <>
Real FE<3,CLOUGH>::shape(const ElemType,
			     const Order,
			     const unsigned int,
			     const Point&)
{
  libMesh::err << "Clough-Tocher elements require the real element\n"
	        << "to construct gradient-based degrees of freedom."
	        << std::endl;

  libmesh_error();
  return 0.;
}



template <>
Real FE<3,CLOUGH>::shape(const Elem* elem,
			     const Order order,
			     const unsigned int,
			     const Point&)
{
  libmesh_assert (elem != NULL);

  libMesh::err << "3D Clough elements not yet implemented."
	        << std::endl;

  libmesh_error();

  clough_compute_coefs(elem);

  const ElemType type = elem->type();

  switch (order+elem->p_level())
    {
      // 3rd-order Clough-Tocher element
    case THIRD:
      {
	switch (type)
	  {
	  default:
            libMesh::err << "ERROR: Unsupported element type!" << std::endl;
	    libmesh_error();
	  }
      }
      // by default throw an error
    default:
      libMesh::err << "ERROR: Unsupported polynomial order!" << std::endl;
      libmesh_error();
    }

  libmesh_error();
  return 0.;
}



template <>
Real FE<3,CLOUGH>::shape_deriv(const ElemType,
				   const Order,
				   const unsigned int,
				   const unsigned int,
				   const Point&)
{
  libMesh::err << "Clough-Tocher elements require the real element\n"
	        << "to construct gradient-based degrees of freedom."
	        << std::endl;

  libmesh_error();
  return 0.;
}



template <>
Real FE<3,CLOUGH>::shape_deriv(const Elem* elem,
				   const Order order,
				   const unsigned int,
				   const unsigned int,
				   const Point&)
{
  libmesh_assert (elem != NULL);

  libMesh::err << "3D Clough elements not yet implemented."
	        << std::endl;

  libmesh_error();

  clough_compute_coefs(elem);

  const ElemType type = elem->type();

  switch (order+elem->p_level())
    {
      // 3rd-order Clough-Tocher element
    case THIRD:
      {
	switch (type)
	  {
	  default:
            libMesh::err << "ERROR: Unsupported element type!" << std::endl;
	    libmesh_error();
	  }
      }
      // by default throw an error
    default:
      libMesh::err << "ERROR: Unsupported polynomial order!" << std::endl;
      libmesh_error();
    }

  libmesh_error();
  return 0.;
}



template <>
Real FE<3,CLOUGH>::shape_second_deriv(const Elem* elem,
                                      const Order order,
                                      const unsigned int,
                                      const unsigned int,
                                      const Point&)
{
  libmesh_assert (elem != NULL);

  libMesh::err << "3D Clough elements not yet implemented."
	    << std::endl;

  libmesh_error();

  clough_compute_coefs(elem);

  const ElemType type = elem->type();

  switch (order+elem->p_level())
    {
      // 3rd-order Clough-Tocher element
    case THIRD:
      {
	switch (type)
	  {
	  default:
            libMesh::err << "ERROR: Unsupported element type!" << std::endl;
	    libmesh_error();
	  }
      }
      // by default throw an error
    default:
      libMesh::err << "ERROR: Unsupported polynomial order!" << std::endl;
      libmesh_error();
    }

  libmesh_error();
  return 0.;
}

} // namespace libMesh
