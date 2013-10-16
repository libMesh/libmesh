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
#include "libmesh/fe.h"
#include "libmesh/elem.h"

namespace libMesh
{

template <>
RealGradient FE<2,NEDELEC_ONE>::shape(const ElemType,
				      const Order,
				      const unsigned int,
				      const Point&)
{
#if LIBMESH_DIM > 1
  libMesh::err << "Nedelec elements require the element type\n"
	       << "because edge orientation is needed."
	       << std::endl;
  libmesh_error();
#endif // LIBMESH_DIM > 1

  libmesh_error();
  return RealGradient();
}


// An excellent discussion of Nedelec shape functions is given in
// http://www.dealii.org/developer/reports/nedelec/nedelec.pdf
template <>
RealGradient FE<2,NEDELEC_ONE>::shape(const Elem* elem,
				      const Order order,
				      const unsigned int i,
				      const Point& p)
{
#if LIBMESH_DIM > 1
  libmesh_assert(elem);

  const Order total_order = static_cast<Order>(order + elem->p_level());

  switch (total_order)
    {
    case FIRST:
      {
	switch (elem->type())
	  {
	  case QUAD8:
	  case QUAD9:
	    {
	      libmesh_assert_less (i, 4);

	      const Real xi  = p(0);
	      const Real eta = p(1);

              // Even with a loose inverse_map tolerance we ought to
              // be nearly on the element interior in master
              // coordinates
	      libmesh_assert_less_equal ( std::fabs(xi), 1.0+10*TOLERANCE );
	      libmesh_assert_less_equal ( std::fabs(eta), 1.0+10*TOLERANCE );

	      switch(i)
		{
		case 0:
		  {
		    if( elem->point(0) > elem->point(1) )
		      return RealGradient( -0.25*(1.0-eta), 0.0 );
		    else
		      return RealGradient( 0.25*(1.0-eta), 0.0 );
		  }
		case 1:
		  {
		    if( elem->point(1) > elem->point(2) )
		      return RealGradient( 0.0, -0.25*(1.0+xi) );
		    else
		      return RealGradient( 0.0, 0.25*(1.0+xi) );
		  }

		case 2:
		  {
		    if( elem->point(2) > elem->point(3) )
		      return RealGradient( 0.25*(1.0+eta), 0.0 );
		    else
		      return RealGradient( -0.25*(1.0+eta), 0.0 );
		  }
		case 3:
		  {
		    if( elem->point(3) > elem->point(0) )
		      return RealGradient( 0.0, -0.25*(xi-1.0) );
		    else
		      return RealGradient( 0.0, 0.25*(xi-1.0) );
		  }

		default:
		  libmesh_error();

		}

	      return RealGradient();
	    }

	  case TRI6:
	    {
	      const Real xi  = p(0);
	      const Real eta = p(1);

	      libmesh_assert_less (i, 3);

	      switch(i)
		{
		case 0:
		  {
		    if( elem->point(0) > elem->point(1) )
		      return RealGradient( -1.0+eta, -xi );
		    else
		      return RealGradient( 1.0-eta, xi );
		  }
		case 1:
		  {
		    if( elem->point(1) > elem->point(2) )
		      return RealGradient( eta, -xi );
		    else
		      return RealGradient( -eta, xi );
		  }

		case 2:
		  {
		    if( elem->point(2) > elem->point(0) )
		      return RealGradient( eta, -xi+1.0 );
		    else
		      return RealGradient( -eta, xi-1.0 );
		  }

		default:
		  libmesh_error();

		}
	    }

	  default:
	    {
	      libMesh::err << "ERROR: Unsupported 2D element type!: " << elem->type()
			    << std::endl;
	      libmesh_error();
	    }
	  }
      }

      // unsupported order
    default:
      {
	libMesh::err << "ERROR: Unsupported 2D FE order!: " << total_order
		      << std::endl;
	libmesh_error();
      }
    }
#endif // LIBMESH_DIM > 1

  libmesh_error();
  return RealGradient();
}



template <>
RealGradient FE<2,NEDELEC_ONE>::shape_deriv(const ElemType,
					    const Order,
					    const unsigned int,
					    const unsigned int,
					    const Point&)
{
#if LIBMESH_DIM > 1
  libMesh::err << "Nedelec elements require the element type\n"
	       << "because edge orientation is needed."
	       << std::endl;
  libmesh_error();
#endif // LIBMESH_DIM > 1

  libmesh_error();
  return RealGradient();
}



template <>
RealGradient FE<2,NEDELEC_ONE>::shape_deriv(const Elem* elem,
					    const Order order,
					    const unsigned int i,
					    const unsigned int j,
					    const Point&)
{
#if LIBMESH_DIM > 1
  libmesh_assert(elem);
  libmesh_assert_less (j, 2);

  const Order total_order = static_cast<Order>(order + elem->p_level());

  switch (total_order)
    {
      // linear Lagrange shape functions
    case FIRST:
      {
	switch (elem->type())
	  {
	  case QUAD8:
	  case QUAD9:
	    {
	      libmesh_assert_less (i, 4);

	      switch (j)
		{
		  // d()/dxi
		case 0:
		  {
		    switch(i)
		      {
		      case 0:
		      case 2:
			return RealGradient();
		      case 1:
			{
			  if( elem->point(1) > elem->point(2) )
			    return RealGradient( 0.0, -0.25 );
			  else
			    return RealGradient( 0.0, 0.25 );
			}
		      case 3:
			{
			  if( elem->point(3) > elem->point(0) )
			    return RealGradient( 0.0, -0.25 );
			  else
			    return RealGradient( 0.0, 0.25 );
			}
		      default:
			libmesh_error();
		      }
		  } // j=0

		  // d()/deta
		case 1:
		  {
		    switch(i)
		      {
		      case 1:
		      case 3:
			return RealGradient();
		      case 0:
			{
			  if( elem->point(0) > elem->point(1) )
			    return RealGradient( 0.25 );
			  else
			    return RealGradient( -0.25 );
			}
		      case 2:
			{
			  if( elem->point(2) > elem->point(3) )
			    return RealGradient( 0.25 );
			  else
			    return RealGradient( -0.25 );
			}
		      default:
			libmesh_error();
		      }
		  } // j=1

		default:
		  libmesh_error();
		}

	      return RealGradient();
	    }

	  case TRI6:
	    {
	      libmesh_assert_less (i, 3);

	      // Account for edge flipping
	      Real f = 1.0;

	      switch(i)
		{
		case 0:
		  {
		    if( elem->point(0) > elem->point(1) )
		      f = -1.0;
		    break;
		  }
		case 1:
		  {
		    if( elem->point(1) > elem->point(2) )
		      f = -1.0;
		    break;
		  }
		case 2:
		  {
		    if( elem->point(2) > elem->point(0) )
		      f = -1.0;
		    break;
		  }
		default:
		  libmesh_error();
		}

	      switch (j)
		{
		  // d()/dxi
		case 0:
		  {
		    return RealGradient( 0.0, f*1.0);
		  }
		  // d()/deta
		case 1:
		  {
		    return RealGradient( f*(-1.0) );
		  }
		default:
		  libmesh_error();
		}
	    }

	  default:
	    {
	      libMesh::err << "ERROR: Unsupported 2D element type!: " << elem->type()
			    << std::endl;
	      libmesh_error();
	    }
	  }
      }
      // unsupported order
    default:
      {
	libMesh::err << "ERROR: Unsupported 2D FE order!: " << total_order
		      << std::endl;
	libmesh_error();
      }
    }
#endif // LIBMESH_DIM > 1

  libmesh_error();
  return RealGradient();
}




template <>
RealGradient FE<2,NEDELEC_ONE>::shape_second_deriv(const ElemType,
						   const Order,
						   const unsigned int,
						   const unsigned int,
						   const Point&)
{
#if LIBMESH_DIM > 1
  libMesh::err << "Nedelec elements require the element type\n"
	       << "because edge orientation is needed."
	       << std::endl;
  libmesh_error();
#endif // LIBMESH_DIM > 1

  libmesh_error();
  return RealGradient();
}



template <>
RealGradient FE<2,NEDELEC_ONE>::shape_second_deriv(const Elem* elem,
						   const Order order,
						   const unsigned int i,
						   const unsigned int j,
						   const Point&)
{
#if LIBMESH_DIM > 1
   libmesh_assert(elem);

  // j = 0 ==> d^2 phi / dxi^2
  // j = 1 ==> d^2 phi / dxi deta
  // j = 2 ==> d^2 phi / deta^2
  libmesh_assert_less (j, 3);

  const Order total_order = static_cast<Order>(order + elem->p_level());

  switch (total_order)
    {
      // linear Lagrange shape functions
    case FIRST:
      {
	switch (elem->type())
	  {
	  case QUAD8:
	  case QUAD9:
	    {
	      libmesh_assert_less (i, 4);
	      // All second derivatives for linear quads are zero.
	      return RealGradient();
	    }

	  case TRI6:
	    {
	      libmesh_assert_less (i, 3);
	      // All second derivatives for linear triangles are zero.
	      return RealGradient();
	    }

	  default:
	    {
	      libMesh::err << "ERROR: Unsupported 2D element type!: " << elem->type()
			    << std::endl;
	      libmesh_error();
	    }

	  } // end switch (type)
      } // end case FIRST

      // unsupported order
    default:
      {
	libMesh::err << "ERROR: Unsupported 2D FE order!: " << total_order
		      << std::endl;
	libmesh_error();
      }

    } // end switch (order)

#endif // LIBMESH_DIM > 1

  libmesh_error();
  return RealGradient();
}

} // namespace libMesh
