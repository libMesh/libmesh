// $Id: fe_interface_inf_fe.C,v 1.6 2003-04-03 14:17:24 ddreyer Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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
#include "fe_interface.h"
#ifdef ENABLE_INFINITE_ELEMENTS

#include "elem.h"
#include "fe.h"
#include "inf_fe.h"
#include "fe_compute_data.h"




//------------------------------------------------------------
//FEInterface class members handling calls to InfFE



unsigned int FEInterface::ifem_n_shape_functions(const unsigned int dim,
						 const FEType& fe_t,
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
      error();
    }

  
  error();
  return 0;
}





unsigned int FEInterface::ifem_n_dofs(const unsigned int dim,
				      const FEType& fe_t,
				      const ElemType t)
{
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
      error();
    }

  
  error();
  return 0;
}

		


unsigned int FEInterface::ifem_n_dofs_at_node(const unsigned int dim,
					      const FEType& fe_t,
					      const ElemType t,
					      const unsigned int n)
{
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
      error();
    }

  
  error();
  return 0;
}





unsigned int FEInterface::ifem_n_dofs_per_elem(const unsigned int dim,
					       const FEType& fe_t,
					       const ElemType t)
{
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
      error();
    }

  
  error();
  return 0;
}




void FEInterface::ifem_nodal_soln(const unsigned int dim,
				  const FEType& fe_t,
				  const Elem* elem,
				  const std::vector<Number>& elem_soln,
				  std::vector<Number>& nodal_soln)
{
  switch (dim)
    {

      // 1D
    case 1:
      {
	switch (fe_t.radial_family)
	  {
	  case INFINITE_MAP:
	    {
	      std::cerr << "ERROR: INFINTE_MAP is not a valid shape family for radial approximation." << std::endl;
	      error();
	      break;
	    }

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
		    std::cerr << "ERROR: Spherical & Ellipsoidal IFEMs not implemented." << std::endl;
		    error();
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
		    std::cerr << "ERROR: Spherical & Ellipsoidal IFEMs not implemented." << std::endl;
		    error();
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
		    std::cerr << "ERROR: Spherical & Ellipsoidal IFEMs not implemented." << std::endl;
		    error();
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
		    std::cerr << "ERROR: Spherical & Ellipsoidal IFEMs not implemented." << std::endl;
		    error();
		}
	      break;
	    }


	    
	  default:
	    std::cerr << "ERROR: Bad FEType.radial_family= " << fe_t.radial_family << std::endl;
	    error();
	    break;
	  }

	break;
      }

      


      // 2D
    case 2:
      {
	switch (fe_t.radial_family)
	  {
	  case INFINITE_MAP:
	    {
	      std::cerr << "ERROR: INFINTE_MAP is not a valid shape family for radial approximation." << std::endl;
	      error();
	      break;
	    }

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
		    std::cerr << "ERROR: Spherical & Ellipsoidal IFEMs not implemented." << std::endl;			      
		    error();
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
		    std::cerr << "ERROR: Spherical & Ellipsoidal IFEMs not implemented." << std::endl;			      
		    error();
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
		    std::cerr << "ERROR: Spherical & Ellipsoidal IFEMs not implemented." << std::endl;			      
		    error();
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
		    std::cerr << "ERROR: Spherical & Ellipsoidal IFEMs not implemented." << std::endl;			      
		    error();
		}
	      break;
	    }


	    
	  default:
	    std::cerr << "ERROR: Bad FEType.radial_family= " << fe_t.radial_family << std::endl;
	    error();
	    break;
	  }

	break;
      }

      


      // 3D
    case 3:
      {
	switch (fe_t.radial_family)
	  {
	  case INFINITE_MAP:
	    {
	      std::cerr << "ERROR: INFINTE_MAP is not a valid shape family for radial approximation." << std::endl;
	      error();
	      break;
	    }

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
		    std::cerr << "ERROR: Spherical & Ellipsoidal IFEMs not implemented." << std::endl;			      
		    error();
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
		    std::cerr << "ERROR: Spherical & Ellipsoidal IFEMs not implemented." << std::endl;			      
		    error();
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
		    std::cerr << "ERROR: Spherical & Ellipsoidal IFEMs not implemented." << std::endl;		      
		    error();
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
		    std::cerr << "ERROR: Spherical & Ellipsoidal IFEMs not implemented." << std::endl;			      
		    error();
		}
	      break;
	    }


	    
	  default:
	    std::cerr << "ERROR: Bad FEType.radial_family= " << fe_t.radial_family << std::endl;
	    error();
	    break;
	  }

	break;
      }

    default:
      error();
    }
  return;
}







Point FEInterface::ifem_inverse_map (const unsigned int dim,
				     const FEType& fe_t,
				     const Elem* elem,
				     const Point& p,
				     const Real tolerance)
{
  switch (dim)
    {
      // 1D
    case 1:
      {
	switch (fe_t.inf_map)
	  {
	  case CARTESIAN:
	    return InfFE<1,JACOBI_20_00,CARTESIAN>::inverse_map(elem, p, tolerance);

	  case SPHERICAL:
	  case ELLIPSOIDAL:
	    {
	      std::cerr << "ERROR: Spherical and Ellipsoidal IFEMs not (yet) " << std::endl
			<< "implemented." << std::endl;
	      error();
	    }

/*
	  case SPHERICAL:
	    return InfFE<1,JACOBI_20_00,SPHERICAL>::inverse_map(elem, p, tolerance);

	  case ELLIPSOIDAL:
	    return InfFE<1,JACOBI_20_00,ELLIPSOIDAL>::inverse_map(elem, p, tolerance);
*/

	  default:
	    error();
	  }
      }

      
      // 2D
    case 2:
      {
	switch (fe_t.inf_map)
	  {
	  case CARTESIAN:
	    return InfFE<2,JACOBI_20_00,CARTESIAN>::inverse_map(elem, p, tolerance);

	  case SPHERICAL:
	  case ELLIPSOIDAL:
	    {
	      std::cerr << "ERROR: Spherical and Ellipsoidal IFEMs not (yet) " << std::endl
			<< "implemented." << std::endl;
	      error();
	    }

/*
	  case SPHERICAL:
	    return InfFE<2,JACOBI_20_00,SPHERICAL>::inverse_map(elem, p, tolerance);

	  case ELLIPSOIDAL:
	    return InfFE<2,JACOBI_20_00,ELLIPSOIDAL>::inverse_map(elem, p, tolerance);
*/

	  default:
	    error();
	  }

      }

      
      // 3D
    case 3:
      {
	switch (fe_t.inf_map)
	  {
	  case CARTESIAN:
	    return InfFE<3,JACOBI_20_00,CARTESIAN>::inverse_map(elem, p, tolerance);

	  case SPHERICAL:
	  case ELLIPSOIDAL:
	    {
	      std::cerr << "ERROR: Spherical and Ellipsoidal IFEMs not (yet) " << std::endl
			<< "implemented." << std::endl;
	      error();
	    }

/*
	  case SPHERICAL:
	    return InfFE<3,JACOBI_20_00,SPHERICAL>::inverse_map(elem, p, tolerance);

	  case ELLIPSOIDAL:
	    return InfFE<3,JACOBI_20_00,ELLIPSOIDAL>::inverse_map(elem, p, tolerance);
*/

	  default:
	    error();
	  }

      }


    default:
      error();
    }

  
  error();
  Point pt;
  return pt;
}



bool FEInterface::ifem_on_reference_element(const Point& p,
					    const ElemType t,
					    const Real eps)
{
  return FEBase::on_reference_element(p,t,eps);
}




Real FEInterface::ifem_shape(const unsigned int dim,
			     const FEType& fe_t,
			     const ElemType t,
			     const unsigned int i,
			     const Point& p)
{
  switch (dim)
    {
      // 1D
    case 1:
      {
	switch (fe_t.radial_family)
	  {
	    /*
	     * For no derivatives (and local coordinates, as
	     * given in \p p) the infinite element shapes
	     * are independent of mapping type
	     */
	  case INFINITE_MAP:
	    return InfFE<1,INFINITE_MAP,CARTESIAN>::shape(fe_t, t, i, p);

	  case JACOBI_20_00:
	    return InfFE<1,JACOBI_20_00,CARTESIAN>::shape(fe_t, t, i, p);

	  case JACOBI_30_00:
	    return InfFE<1,JACOBI_30_00,CARTESIAN>::shape(fe_t, t, i, p);

	  case LEGENDRE:   
	    return InfFE<1,LEGENDRE,CARTESIAN>::shape(fe_t, t, i, p);

	  case LAGRANGE:
	    return InfFE<1,LAGRANGE,CARTESIAN>::shape(fe_t, t, i, p);

	  default:
	    error();
	  }
      }

      
      // 2D
    case 2:
      {
	switch (fe_t.radial_family)
	  {
	  case INFINITE_MAP:
	    return InfFE<2,INFINITE_MAP,CARTESIAN>::shape(fe_t, t, i, p);

	  case JACOBI_20_00:
	    return InfFE<2,JACOBI_20_00,CARTESIAN>::shape(fe_t, t, i, p);

	  case JACOBI_30_00:
	    return InfFE<2,JACOBI_30_00,CARTESIAN>::shape(fe_t, t, i, p);

	  case LEGENDRE:   
	    return InfFE<2,LEGENDRE,CARTESIAN>::shape(fe_t, t, i, p);

	  case LAGRANGE:
	    return InfFE<2,LAGRANGE,CARTESIAN>::shape(fe_t, t, i, p);

	  default:
	    error();
	  }

      }

      
      // 3D
    case 3:
      {
	switch (fe_t.radial_family)
	  {
	  case INFINITE_MAP:
	    return InfFE<3,INFINITE_MAP,CARTESIAN>::shape(fe_t, t, i, p);

	  case JACOBI_20_00:
	    return InfFE<3,JACOBI_20_00,CARTESIAN>::shape(fe_t, t, i, p);

	  case JACOBI_30_00:
	    return InfFE<3,JACOBI_30_00,CARTESIAN>::shape(fe_t, t, i, p);

	  case LEGENDRE:   
	    return InfFE<3,LEGENDRE,CARTESIAN>::shape(fe_t, t, i, p);

	  case LAGRANGE:
	    return InfFE<3,LAGRANGE,CARTESIAN>::shape(fe_t, t, i, p);

	  default:
	    error();
	  }

      }


    default:
      error();
    }

  
  error();
  return 0.;
}




Real FEInterface::ifem_shape(const unsigned int dim,
			     const FEType& fe_t,
			     const Elem* elem,
			     const unsigned int i,
			     const Point& p)
{
  switch (dim)
    {
      // 1D
    case 1:
      {
	switch (fe_t.radial_family)
	  {
	    /*
	     * For no derivatives (and local coordinates, as
	     * given in \p p) the infinite element shapes
	     * are independent of mapping type
	     */
	  case INFINITE_MAP:
	    return InfFE<1,INFINITE_MAP,CARTESIAN>::shape(fe_t, elem, i, p);

	  case JACOBI_20_00:
	    return InfFE<1,JACOBI_20_00,CARTESIAN>::shape(fe_t, elem, i, p);

	  case JACOBI_30_00:
	    return InfFE<1,JACOBI_30_00,CARTESIAN>::shape(fe_t, elem, i, p);

	  case LEGENDRE:   
	    return InfFE<1,LEGENDRE,CARTESIAN>::shape(fe_t, elem, i, p);

	  case LAGRANGE:
	    return InfFE<1,LAGRANGE,CARTESIAN>::shape(fe_t, elem, i, p);

	  default:
	    error();
	  }
      }


      // 2D
    case 2:
      {
	switch (fe_t.radial_family)
	  {
	  case INFINITE_MAP:
	    return InfFE<2,INFINITE_MAP,CARTESIAN>::shape(fe_t, elem, i, p);

	  case JACOBI_20_00:
	    return InfFE<2,JACOBI_20_00,CARTESIAN>::shape(fe_t, elem, i, p);

	  case JACOBI_30_00:
	    return InfFE<2,JACOBI_30_00,CARTESIAN>::shape(fe_t, elem, i, p);

	  case LEGENDRE:   
	    return InfFE<2,LEGENDRE,CARTESIAN>::shape(fe_t, elem, i, p);

	  case LAGRANGE:
	    return InfFE<2,LAGRANGE,CARTESIAN>::shape(fe_t, elem, i, p);

	  default:
	    error();
	  }

      }

            
      // 3D
    case 3:
      {
	switch (fe_t.radial_family)
	  {
	  case INFINITE_MAP:
	    return InfFE<3,INFINITE_MAP,CARTESIAN>::shape(fe_t, elem, i, p);

	  case JACOBI_20_00:
	    return InfFE<3,JACOBI_20_00,CARTESIAN>::shape(fe_t, elem, i, p);

	  case JACOBI_30_00:
	    return InfFE<3,JACOBI_30_00,CARTESIAN>::shape(fe_t, elem, i, p);

	  case LEGENDRE:   
	    return InfFE<3,LEGENDRE,CARTESIAN>::shape(fe_t, elem, i, p);

	  case LAGRANGE:
	    return InfFE<3,LAGRANGE,CARTESIAN>::shape(fe_t, elem, i, p);

	  default:
	    error();
	  }

      }


    default:
      error();
    }

  
  error();
  return 0.;
}




void FEInterface::ifem_compute_data(const unsigned int dim,
				    const FEType& fe_t,
				    const Elem* elem,
				    FEComputeData& data)
{
  switch (dim)
    {
      // 1D
    case 1:
      {
	switch (fe_t.radial_family)
	  {
	    /*
	     * For no derivatives (and local coordinates, as
	     * given in \p p) the infinite element shapes
	     * are independent of mapping type
	     */
	  case INFINITE_MAP:
	    InfFE<1,INFINITE_MAP,CARTESIAN>::compute_data(fe_t, elem, data);
	    break;

	  case JACOBI_20_00:
	    InfFE<1,JACOBI_20_00,CARTESIAN>::compute_data(fe_t, elem, data);
	    break;

	  case JACOBI_30_00:
	    InfFE<1,JACOBI_30_00,CARTESIAN>::compute_data(fe_t, elem, data);
	    break;

	  case LEGENDRE:   
	    InfFE<1,LEGENDRE,CARTESIAN>::compute_data(fe_t, elem, data);
	    break;

	  case LAGRANGE:
	    InfFE<1,LAGRANGE,CARTESIAN>::compute_data(fe_t, elem, data);
	    break;

	  default:
	    error();
	  }

	break;
      }


      // 2D
    case 2:
      {
	switch (fe_t.radial_family)
	  {
	  case INFINITE_MAP:
	    InfFE<2,INFINITE_MAP,CARTESIAN>::compute_data(fe_t, elem, data);
	    break;

	  case JACOBI_20_00:
	    InfFE<2,JACOBI_20_00,CARTESIAN>::compute_data(fe_t, elem, data);
	    break;

	  case JACOBI_30_00:
	    InfFE<2,JACOBI_30_00,CARTESIAN>::compute_data(fe_t, elem, data);
	    break;

	  case LEGENDRE:   
	    InfFE<2,LEGENDRE,CARTESIAN>::compute_data(fe_t, elem, data);
	    break;

	  case LAGRANGE:
	    InfFE<2,LAGRANGE,CARTESIAN>::compute_data(fe_t, elem, data);
	    break;

	  default:
	    error();
	  }

	break;
      }

            
      // 3D
    case 3:
      {
	switch (fe_t.radial_family)
	  {
	  case INFINITE_MAP:
	    InfFE<3,INFINITE_MAP,CARTESIAN>::compute_data(fe_t, elem, data);
	    break;

	  case JACOBI_20_00:
	    InfFE<3,JACOBI_20_00,CARTESIAN>::compute_data(fe_t, elem, data);
	    break;

	  case JACOBI_30_00:
	    InfFE<3,JACOBI_30_00,CARTESIAN>::compute_data(fe_t, elem, data);
	    break;

	  case LEGENDRE:   
	    InfFE<3,LEGENDRE,CARTESIAN>::compute_data(fe_t, elem, data);
	    break;

	  case LAGRANGE:
	    InfFE<3,LAGRANGE,CARTESIAN>::compute_data(fe_t, elem, data);
	    break;

	  default:
	    error();
	  }

	break;
      }


    default:
      error();
      break;
    }

  return;
}

#endif // ifdef ENABLE_INFINITE_ELEMENTS

