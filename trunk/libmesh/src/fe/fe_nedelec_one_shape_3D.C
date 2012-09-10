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

namespace libMesh
{

template <>
RealGradient FE<3,NEDELEC_ONE>::shape(const ElemType,
				      const Order,
				      const unsigned int,
				      const Point&)
{
#if LIBMESH_DIM == 3
  libMesh::err << "Nedelec elements require the element type\n"
	       << "because edge orientation is needed."
	       << std::endl;
  libmesh_error();
#endif

  libmesh_error();
  return RealGradient();
}



template <>
RealGradient FE<3,NEDELEC_ONE>::shape(const Elem* elem,
				      const Order order,
				      const unsigned int i,
				      const Point& /*p*/)
{
#if LIBMESH_DIM == 3
  libmesh_assert (elem != NULL);

  const Order totalorder = static_cast<Order>(order + elem->p_level());

  switch (totalorder)
    {
      // linear Lagrange shape functions
    case FIRST:
      {
	switch (elem->type())
	  {
	  case HEX20:
	  case HEX27:
	    {
	      libmesh_assert (i<12);

	      libmesh_not_implemented();
	      return RealGradient();
	    }

	  case TET10:
	    {
	      libmesh_assert(i<6);

	      libmesh_not_implemented();
	      return RealGradient();
	    }

	  default:
	    {
	      libMesh::err << "ERROR: Unsupported 3D element type!: " << elem->type()
			    << std::endl;
	      libmesh_error();
	    }
	  }
      }

      // unsupported order
    default:
      {
	libMesh::err << "ERROR: Unsupported 3D FE order!: " << totalorder
		      << std::endl;

	libmesh_error();
      }
    }
#endif

  libmesh_error();
  return RealGradient();
}




template <>
RealGradient FE<3,NEDELEC_ONE>::shape_deriv(const ElemType,
					    const Order,
					    const unsigned int,
					    const unsigned int,
					    const Point&)
{
#if LIBMESH_DIM == 3
  libMesh::err << "Nedelec elements require the element type\n"
	       << "because edge orientation is needed."
	       << std::endl;
  libmesh_error();
#endif

  libmesh_error();
  return RealGradient();
}

template <>
RealGradient FE<3,NEDELEC_ONE>::shape_deriv(const Elem* elem,
					    const Order order,
					    const unsigned int i,
					    const unsigned int j,
					    const Point& /*p*/)
{
#if LIBMESH_DIM == 3
  libmesh_assert (elem != NULL);
  libmesh_assert (j<3);

  const Order totalorder = static_cast<Order>(order + elem->p_level());

  switch (totalorder)
    {
    case FIRST:
      {
	switch (elem->type())
	  {
	  case HEX20:
	  case HEX27:
	    {
	      libmesh_assert (i<12);

	      libmesh_not_implemented();
	      return RealGradient();
	    }

	  case TET10:
	    {
	      libmesh_assert (i<6);

	      libmesh_not_implemented();
	      return RealGradient();
	    }

	  default:
	    {
	      libMesh::err << "ERROR: Unsupported 3D element type!: " << elem->type()
			    << std::endl;
	      libmesh_error();
	    }
	  }
      }

      // unsupported order
    default:
      {
	libMesh::err << "ERROR: Unsupported 3D FE order!: " << totalorder
		      << std::endl;
	libmesh_error();
      }
    }

#endif

  libmesh_error();
  return RealGradient();
}



template <>
RealGradient FE<3,NEDELEC_ONE>::shape_second_deriv(const ElemType,
						   const Order,
						   const unsigned int,
						   const unsigned int,
						   const Point&)
{
#if LIBMESH_DIM == 3
  libMesh::err << "Nedelec elements require the element type\n"
	       << "because edge orientation is needed."
	       << std::endl;
  libmesh_error();
#endif

  libmesh_error();
  return RealGradient();
}



template <>
RealGradient FE<3,NEDELEC_ONE>::shape_second_deriv(const Elem* elem,
						   const Order order,
						   const unsigned int i,
						   const unsigned int j,
						   const Point& /*p*/)
{
#if LIBMESH_DIM == 3

  libmesh_assert (elem != NULL);

  libmesh_assert (j<6);

  const Order totalorder = static_cast<Order>(order + elem->p_level());

  switch (totalorder)
    {
      // linear Lagrange shape functions
    case FIRST:
      {
	switch (elem->type())
	  {
	  case HEX20:
	  case HEX27:
	    {
	      libmesh_assert (i<12);
	      
	      libmesh_not_implemented();
	      return RealGradient();
	    }

	  case TET10:
	    {
	      libmesh_assert (i<6);

	      libmesh_not_implemented();
	      return RealGradient();
	    }

	  default:
	    {
	      libMesh::err << "ERROR: Unsupported 3D element type!: " << elem->type()
			    << std::endl;
	      libmesh_error();
	    }

	  } //switch(type)

      } // case FIRST:
      // unsupported order
    default:
      {
	libMesh::err << "ERROR: Unsupported 3D FE order!: " << totalorder
		      << std::endl;
	libmesh_error();
      }
    }

#endif

  libmesh_error();
  return RealGradient();
}

} // namespace libMesh
