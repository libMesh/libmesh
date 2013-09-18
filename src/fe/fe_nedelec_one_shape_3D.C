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
				      const Point& p)
{
#if LIBMESH_DIM == 3
  libmesh_assert(elem);

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
	      libmesh_assert_less (i, 12);

              const Real xi   = p(0);
	      const Real eta  = p(1);
              const Real zeta = p(2);

	      libmesh_assert_less_equal ( std::fabs(xi),   1.0 );
	      libmesh_assert_less_equal ( std::fabs(eta),  1.0 );
              libmesh_assert_less_equal ( std::fabs(zeta), 1.0 );

              switch(i)
		{
                case 0:
                  {
                    if( elem->point(0) > elem->point(1) )
		      return RealGradient( -0.125*(1.0-eta-zeta+eta*zeta), 0.0, 0.0 );
		    else
		      return RealGradient(  0.125*(1.0-eta-zeta+eta*zeta), 0.0, 0.0 );
                  }
                case 1:
		  {
		    if( elem->point(1) > elem->point(2) )
		      return RealGradient( 0.0, -0.125*(1.0+xi-zeta-xi*zeta), 0.0 );
		    else
		      return RealGradient( 0.0,  0.125*(1.0+xi-zeta-xi*zeta), 0.0 );
		  }
		case 2:
		  {
		    if( elem->point(2) > elem->point(3) )
		      return RealGradient(  0.125*(1.0+eta-zeta-eta*zeta), 0.0, 0.0 );
		    else
		      return RealGradient( -0.125*(1.0+eta-zeta-eta*zeta), 0.0, 0.0 );
		  }
		case 3:
		  {
		    if( elem->point(3) > elem->point(0) )
		      return RealGradient( 0.0,  0.125*(1.0-xi-zeta+xi*zeta), 0.0 );
		    else
		      return RealGradient( 0.0, -0.125*(1.0-xi-zeta+xi*zeta), 0.0 );
		  }
                case 4:
		  {
		    if( elem->point(0) > elem->point(4) )
		      return RealGradient( 0.0, 0.0, -0.125*(1.0-xi-eta+xi*eta) );
		    else
		      return RealGradient( 0.0, 0.0,  0.125*(1.0-xi-eta+xi*eta) );
		  }
                case 5:
		  {
		    if( elem->point(1) > elem->point(5) )
		      return RealGradient( 0.0, 0.0, -0.125*(1.0+xi-eta-xi*eta) );
		    else
		      return RealGradient( 0.0, 0.0,  0.125*(1.0+xi-eta-xi*eta) );
		  }
                case 6:
		  {
		    if( elem->point(2) > elem->point(6) )
		      return RealGradient( 0.0, 0.0, -0.125*(1.0+xi+eta+xi*eta) );
		    else
		      return RealGradient( 0.0, 0.0,  0.125*(1.0+xi+eta+xi*eta) );
		  }
                case 7:
		  {
		    if( elem->point(3) > elem->point(7) )
		      return RealGradient( 0.0, 0.0, -0.125*(1.0-xi+eta-xi*eta) );
		    else
		      return RealGradient( 0.0, 0.0,  0.125*(1.0-xi+eta-xi*eta) );
		  }
                case 8:
                  {
                    if( elem->point(4) > elem->point(5) )
		      return RealGradient( -0.125*(1.0-eta+zeta-eta*zeta), 0.0, 0.0 );
		    else
		      return RealGradient(  0.125*(1.0-eta+zeta-eta*zeta), 0.0, 0.0 );
                  }
                case 9:
		  {
		    if( elem->point(5) > elem->point(6) )
		      return RealGradient( 0.0, -0.125*(1.0+xi+zeta+xi*zeta), 0.0 );
		    else
		      return RealGradient( 0.0,  0.125*(1.0+xi+zeta+xi*zeta), 0.0 );
		  }
		case 10:
		  {
		    if( elem->point(7) > elem->point(6) )
		      return RealGradient( -0.125*(1.0+eta+zeta+eta*zeta), 0.0, 0.0 );
		    else
		      return RealGradient(  0.125*(1.0+eta+zeta+eta*zeta), 0.0, 0.0 );
		  }
		case 11:
		  {
		    if( elem->point(4) > elem->point(7) )
		      return RealGradient( 0.0, -0.125*(1.0-xi+zeta-xi*zeta), 0.0 );
		    else
		      return RealGradient( 0.0,  0.125*(1.0-xi+zeta-xi*zeta), 0.0 );
		  }

                default:
		  libmesh_error();
                }

	      return RealGradient();
	    }

	  case TET10:
	    {
	      libmesh_assert_less (i, 6);

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
					    const Point& p)
{
#if LIBMESH_DIM == 3
  libmesh_assert(elem);
  libmesh_assert_less (j, 3);

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
	      libmesh_assert_less (i, 12);

	      const Real xi   = p(0);
	      const Real eta  = p(1);
              const Real zeta = p(2);

	      libmesh_assert_less_equal ( std::fabs(xi),   1.0 );
	      libmesh_assert_less_equal ( std::fabs(eta),  1.0 );
              libmesh_assert_less_equal ( std::fabs(zeta), 1.0 );

	      switch (j)
		{
		  // d()/dxi
		case 0:
		  {
		    switch(i)
		      {
		      case 0:
		      case 2:
                      case 8:
                      case 10:
			return RealGradient();
		      case 1:
			{
			  if( elem->point(1) > elem->point(2) )
			    return RealGradient( 0.0, -0.125*(1.0-zeta) );
			  else
			    return RealGradient( 0.0, 0.125*(1.0-zeta) );
			}
		      case 3:
			{
			  if( elem->point(3) > elem->point(0) )
			    return RealGradient( 0.0, 0.125*(-1.0+zeta) );
			  else
			    return RealGradient( 0.0, -0.125*(-1.0+zeta) );
			}
                      case 4:
                        {
                          if( elem->point(0) > elem->point(4) )
                            return RealGradient( 0.0, 0.0, -0.125*(-1.0+eta) );
                          else
                            return RealGradient( 0.0, 0.0,  0.125*(-1.0+eta) );
                        }
                      case 5:
                        {
                          if( elem->point(1) > elem->point(5) )
                            return RealGradient( 0.0, 0.0, -0.125*(1.0-eta) );
                          else
                            return RealGradient( 0.0, 0.0,  0.125*(1.0-eta) );
                        }
                      case 6:
                        {
                          if( elem->point(2) > elem->point(6) )
                            return RealGradient( 0.0, 0.0, -0.125*(1.0+eta) );
                          else
                            return RealGradient( 0.0, 0.0,  0.125*(1.0+eta) );
                        }
                      case 7:
                        {
                          if( elem->point(3) > elem->point(7) )
                            return RealGradient( 0.0, 0.0, -0.125*(-1.0-eta) );
                          else
                            return RealGradient( 0.0, 0.0,  0.125*(-1.0-eta) );
                        }
                      case 9:
                        {
                          if( elem->point(5) > elem->point(6) )
                            return RealGradient( 0.0, -0.125*(1.0+zeta), 0.0 );
                          else
                            return RealGradient( 0.0,  0.125*(1.0+zeta), 0.0 );
                        }
                      case 11:
                        {
                          if( elem->point(4) > elem->point(7) )
                            return RealGradient( 0.0, -0.125*(-1.0-zeta), 0.0 );
                          else
                            return RealGradient( 0.0,  0.125*(-1.0-zeta), 0.0 );
                        }
		      default:
			libmesh_error();
		      } // switch(i)

		  } // j=0

		  // d()/deta
		case 1:
		  {
		    switch(i)
		      {
		      case 1:
		      case 3:
                      case 9:
                      case 11:
			return RealGradient();
		      case 0:
                        {
                          if( elem->point(0) > elem->point(1) )
                            return RealGradient( -0.125*(-1.0+zeta), 0.0, 0.0 );
                          else
                            return RealGradient(  0.125*(-1.0+zeta), 0.0, 0.0 );
                        }
		      case 2:
                        {
                          if( elem->point(2) > elem->point(3) )
                            return RealGradient( 0.125*(1.0-zeta), 0.0, 0.0 );
                          else
                            return RealGradient( -0.125*(1.0-zeta), 0.0, 0.0 );
                        }
                      case 4:
                        {
                          if( elem->point(0) > elem->point(4) )
                            return RealGradient( 0.0, 0.0, -0.125*(-1.0+xi) );
                          else
                            return RealGradient( 0.0, 0.0,  0.125*(-1.0+xi) );
                        }
                      case 5:
                        {
                          if( elem->point(1) > elem->point(5) )
                            return RealGradient( 0.0, 0.0, -0.125*(-1.0-xi) );
                          else
                            return RealGradient( 0.0, 0.0,  0.125*(-1.0-xi) );
                        }
                      case 6:
                        {
                          if( elem->point(2) > elem->point(6) )
                            return RealGradient( 0.0, 0.0, -0.125*(1.0+xi) );
                          else
                            return RealGradient( 0.0, 0.0,  0.125*(1.0+xi) );
                        }
                      case 7:
                        {
                          if( elem->point(3) > elem->point(7) )
                            return RealGradient( 0.0, 0.0, -0.125*(1.0-xi) );
                          else
                            return RealGradient( 0.0, 0.0,  0.125*(1.0-xi) );
                        }
                      case 8:
                        {
                          if( elem->point(4) > elem->point(5) )
                            return RealGradient( -0.125*(-1.0-zeta), 0.0, 0.0 );
                          else
                            return RealGradient(  0.125*(-1.0-zeta), 0.0, 0.0 );
                        }
                      case 10:
                        {
                          if( elem->point(7) > elem->point(6) )
                            return RealGradient( -0.125*(1.0+zeta), 0.0, 0.0 );
                          else
                            return RealGradient(  0.125*(1.0+zeta), 0.0, 0.0 );
                        }
		      default:
			libmesh_error();
		      } // switch(i)

		  } // j=1

                  // d()/dzeta
		case 2:
		  {
		    switch(i)
		      {
                      case 4:
                      case 5:
                      case 6:
                      case 7:
			return RealGradient();

                      case 0:
                        {
                          if( elem->point(0) > elem->point(1) )
                            return RealGradient( -0.125*(-1.0+eta), 0.0, 0.0 );
                          else
                            return RealGradient(  0.125*(-1.0+eta), 0.0, 0.0 );
                        }
                      case 1:
                        {
                          if( elem->point(1) > elem->point(2) )
                            return RealGradient( 0.0, -0.125*(-1.0-xi), 0.0 );
                          else
                            return RealGradient( 0.0,  0.125*(-1.0-xi), 0.0 );
                        }
                      case 2:
                        {
                          if( elem->point(2) > elem->point(3) )
                            return RealGradient( 0.125*(-1.0-eta), 0.0, 0.0 );
                          else
                            return RealGradient( -0.125*(-1.0-eta), 0.0, 0.0 );
                        }
                      case 3:
                        {
                          if( elem->point(3) > elem->point(0) )
                            return RealGradient( 0.0, 0.125*(-1.0+xi), 0.0 );
                          else
                            return RealGradient( 0.0,  -0.125*(-1.0+xi), 0.0 );
                        }
                      case 8:
                        {
                          if( elem->point(4) > elem->point(5) )
                            return RealGradient( -0.125*(1.0-eta), 0.0, 0.0 );
                          else
                            return RealGradient(  0.125*(1.0-eta), 0.0, 0.0 );
                        }
                      case 9:
                        {
                          if( elem->point(5) > elem->point(6) )
                            return RealGradient( 0.0, -0.125*(1.0+xi), 0.0 );
                          else
                            return RealGradient( 0.0,  0.125*(1.0+xi), 0.0 );
                        }
                      case 10:
                        {
                          if( elem->point(7) > elem->point(6) )
                            return RealGradient( -0.125*(1.0+eta), 0.0, 0.0 );
                          else
                            return RealGradient(  0.125*(1.0+eta), 0.0, 0.0 );
                        }
                      case 11:
                        {
                          if( elem->point(4) > elem->point(7) )
                            return RealGradient( 0.0, -0.125*(1.0-xi), 0.0 );
                          else
                            return RealGradient( 0.0,  0.125*(1.0-xi), 0.0 );
                        }
                      default:
			libmesh_error();
		      } // switch(i)

                  } // j = 2

		default:
		  libmesh_error();
		}

	      return RealGradient();
	    }

	  case TET10:
	    {
	      libmesh_assert_less (i, 6);

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
						   const Point& p)
{
#if LIBMESH_DIM == 3

  libmesh_assert(elem);

  // j = 0 ==> d^2 phi / dxi^2
  // j = 1 ==> d^2 phi / dxi deta
  // j = 2 ==> d^2 phi / deta^2
  // j = 3 ==> d^2 phi / dxi dzeta
  // j = 4 ==> d^2 phi / deta dzeta
  // j = 5 ==> d^2 phi / dzeta^2
  libmesh_assert_less (j, 6);

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
	      libmesh_assert_less (i, 12);

              const Real xi   = p(0);
	      const Real eta  = p(1);
              const Real zeta = p(2);

	      libmesh_assert_less_equal ( std::fabs(xi),   1.0 );
	      libmesh_assert_less_equal ( std::fabs(eta),  1.0 );
              libmesh_assert_less_equal ( std::fabs(zeta), 1.0 );

	      switch (j)
		{
		  // d^2()/dxi^2
		case 0:
		  {
                    // All d^2()/dxi^2 derivatives for linear hexes are zero.
                    return RealGradient();
                  } // j=0

		  // d^2()/dxideta
		case 1:
		  {
		    switch(i)
		      {
                      case 0:
		      case 1:
                      case 2:
		      case 3:
                      case 8:
                      case 9:
                      case 10:
                      case 11:
			return RealGradient();
                      case 4:
                        {
                          if( elem->point(0) > elem->point(4) )
                            return RealGradient( 0.0, 0.0, -0.125 );
                          else
                            return RealGradient( 0.0, 0.0,  0.125 );
                        }
                      case 5:
                        {
                          if( elem->point(1) > elem->point(5) )
                            return RealGradient( 0.0, 0.0,  0.125 );
                          else
                            return RealGradient( 0.0, 0.0, -0.125 );
                        }
                      case 6:
                        {
                          if( elem->point(2) > elem->point(6) )
                            return RealGradient( 0.0, 0.0, -0.125 );
                          else
                            return RealGradient( 0.0, 0.0,  0.125 );
                        }
                      case 7:
                        {
                          if( elem->point(3) > elem->point(7) )
                            return RealGradient( 0.0, 0.0,  0.125 );
                          else
                            return RealGradient( 0.0, 0.0, -0.125 );
                        }
		      default:
			libmesh_error();
		      } // switch(i)

		  } // j=1

                  // d^2()/deta^2
		case 2:
		  {
		    // All d^2()/deta^2 derivatives for linear hexes are zero.
                    return RealGradient();
                  } // j = 2

                  // d^2()/dxidzeta
		case 3:
                  {
                    switch(i)
		      {
                      case 0:
		      case 2:
                      case 4:
                      case 5:
                      case 6:
                      case 7:
                      case 8:
                      case 10:
			return RealGradient();

                      case 1:
			{
			  if( elem->point(1) > elem->point(2) )
			    return RealGradient( 0.0,  0.125 );
			  else
			    return RealGradient( 0.0, -0.125 );
			}
		      case 3:
			{
			  if( elem->point(3) > elem->point(0) )
			    return RealGradient( 0.0, -0.125 );
			  else
			    return RealGradient( 0.0,  0.125 );
			}
                      case 9:
                        {
                          if( elem->point(5) > elem->point(6) )
                            return RealGradient( 0.0, -0.125, 0.0 );
                          else
                            return RealGradient( 0.0,  0.125, 0.0 );
                        }
                      case 11:
                        {
                          if( elem->point(4) > elem->point(7) )
                            return RealGradient( 0.0,  0.125, 0.0 );
                          else
                            return RealGradient( 0.0, -0.125, 0.0 );
                        }
                      default:
			libmesh_error();
		      } // switch(i)

                  } // j = 3

                  // d^2()/detadzeta
		case 4:
                  {
                    switch(i)
		      {
		      case 1:
		      case 3:
                      case 4:
                      case 5:
                      case 6:
                      case 7:
                      case 9:
                      case 11:
			return RealGradient();

                      case 0:
                        {
                          if( elem->point(0) > elem->point(1) )
                            return RealGradient( -0.125, 0.0, 0.0 );
                          else
                            return RealGradient(  0.125, 0.0, 0.0 );
                        }
		      case 2:
                        {
                          if( elem->point(2) > elem->point(3) )
                            return RealGradient(  0.125, 0.0, 0.0 );
                          else
                            return RealGradient( -0.125, 0.0, 0.0 );
                        }
                      case 8:
                        {
                          if( elem->point(4) > elem->point(5) )
                            return RealGradient(  0.125, 0.0, 0.0 );
                          else
                            return RealGradient( -0.125, 0.0, 0.0 );
                        }
                      case 10:
                        {
                          if( elem->point(7) > elem->point(6) )
                            return RealGradient( -0.125, 0.0, 0.0 );
                          else
                            return RealGradient(  0.125, 0.0, 0.0 );
                        }
                      default:
			libmesh_error();
		      } // switch(i)

                  } // j = 4

                  // d^2()/dzeta^2
		case 5:
                  {
                    // All d^2()/dzeta^2 derivatives for linear hexes are zero.
                    return RealGradient();
                  } // j = 5

		default:
		  libmesh_error();
		}

	      return RealGradient();
	    }

	  case TET10:
	    {
	      libmesh_assert_less (i, 6);

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
