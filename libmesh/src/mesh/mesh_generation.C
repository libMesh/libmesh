// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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



// C++ includes
#include <cmath> // for std::sqrt


// Local includes
#include "mesh_generation.h"
#include "unstructured_mesh.h"
// #include "elem.h"
#include "mesh_refinement.h"
#include "edge_edge2.h"
#include "edge_edge3.h"
#include "edge_edge4.h"
#include "face_tri3.h"
#include "face_tri6.h"
#include "face_quad4.h"
#include "face_quad8.h"
#include "face_quad9.h"
#include "cell_hex8.h"
#include "cell_hex20.h"
#include "cell_hex27.h"
#include "cell_prism6.h"
#include "cell_prism15.h"
#include "cell_prism18.h"
#include "cell_tet4.h"
#include "libmesh_logging.h"
#include "boundary_info.h"
#include "sphere.h"
#include "mesh_modification.h"
#include "mesh_smoother_laplace.h"

// ------------------------------------------------------------
// MeshTools::Generation function for mesh generation
void MeshTools::Generation::build_cube(UnstructuredMesh& mesh,
				       const unsigned int nx,
				       const unsigned int ny,
				       const unsigned int nz,
				       const Real xmin, const Real xmax,
				       const Real ymin, const Real ymax,
				       const Real zmin, const Real zmax,
				       const ElemType type,
				       const bool gauss_lobatto_grid)
{
  START_LOG("build_cube()", "MeshTools::Generation");

  // Declare that we are using the indexing utility routine
  // in the "Private" part of our current namespace.  If this doesn't
  // work in GCC 2.95.3 we can either remove it or stop supporting
  // 2.95.3 altogether.
  using MeshTools::Generation::Private::idx;
  
  // Clear the mesh and start from scratch
  mesh.clear();
  
  switch (mesh.mesh_dimension())
    {

      
      //---------------------------------------------------------------------
      // Build a 1D line
    case 1:      
      {
	libmesh_assert (nx != 0);

        // Reserve elements
        switch (type)
          {
          case INVALID_ELEM:
          case EDGE2:
          case EDGE3:
          case EDGE4:
            {
	      mesh.reserve_elem (nx);
              break;
            }

          default:
            {
	      std::cerr << "ERROR: Unrecognized 1D element type." << std::endl;
	      libmesh_error();
	    }
	  }

        // Reserve nodes
        switch (type)
          {
          case INVALID_ELEM:
          case EDGE2:
            {
              mesh.reserve_nodes(nx+1);
              break;
            }

          case EDGE3:
            {
              mesh.reserve_nodes(2*nx+1);
              break;
            }

          case EDGE4:
            {
              mesh.reserve_nodes(3*nx+1);
              break;
            }

          default:
            {
              std::cerr << "ERROR: Unrecognized 1D element type." << std::endl;
              libmesh_error();
            }
          }


        // Build the nodes, depends on whether we're using linears, 
        // quadratics or cubics and whether using uniform grid or Gauss-Lobatto
        unsigned int node_id = 0;
        switch(type)
        {
          case INVALID_ELEM:
          case EDGE2:
            {
              for (unsigned int i=0; i<=nx; i++)
              {
                if (gauss_lobatto_grid)
                  mesh.add_point (Point(0.5*(std::cos(libMesh::pi*static_cast<Real>(nx-i)/static_cast<Real>(nx))+1.0),
                        0, 
                        0), node_id++);
                else
                  mesh.add_point (Point(static_cast<Real>(i)/static_cast<Real>(nx), 
                        0, 
                        0), node_id++);
              }
              break;
            }

          case EDGE3:
            {
              for (unsigned int i=0; i<=2*nx; i++)
              {
                if (gauss_lobatto_grid)
                {
                  // The x location of the point.
                  Real x=0.;

                  // Shortcut variable
                  const Real pi = libMesh::pi;

                  // Shortcut quantities (do not depend on i)
                  const Real c = std::cos( pi*i / static_cast<Real>(2*nx) );

                  // If i is even, compute a normal Gauss-Lobatto point
                  if (i%2 == 0)
                    x = 0.5*(1.0 - c);

                  // Otherwise, it is the average of the previous and next points
                  else
                  {
                    Real cmin = std::cos( pi*(i-1) / static_cast<Real>(2*nx) );
                    Real cmax = std::cos( pi*(i+1) / static_cast<Real>(2*nx) );
                    
                    Real xmin = 0.5*(1.0 - cmin);
                    Real xmax = 0.5*(1.0 - cmax);
                    x = 0.5*(xmin + xmax);
                  }

                  mesh.add_point (Point(x,0.,0.), node_id++);
                }
                else
                  mesh.add_point (Point(static_cast<Real>(i)/static_cast<Real>(2*nx),
                        0,
                        0), node_id++);
              }
              break;
            }

          case EDGE4:
            {
              for (unsigned int i=0; i<=3*nx; i++)
              {
                if (gauss_lobatto_grid)
                {
                  // The x location of the point
                  Real x=0.;
                  
                  const Real pi = libMesh::pi;

                  // Shortcut quantities
                  const Real c = std::cos( pi*i / static_cast<Real>(3*nx) );

                  // If i is multiple of 3, compute a normal Gauss-Lobatto point
                  if (i%3 == 0)
                    x = 0.5*(1.0 - c);

                  // Otherwise, distribute points evenly within the element
                  else
                  {
                    if(i%3 == 1)
                    {
                      Real cmin = std::cos( pi*(i-1) / static_cast<Real>(3*nx) );
                      Real cmax = std::cos( pi*(i+2) / static_cast<Real>(3*nx) );

                      Real xmin = 0.5*(1.0 - cmin);
                      Real xmax = 0.5*(1.0 - cmax);

                      x = (2.*xmin + xmax)/3.;
                    }
                    else
                    if(i%3 == 2)
                    {
                      Real cmin = std::cos( pi*(i-2) / static_cast<Real>(3*nx) );
                      Real cmax = std::cos( pi*(i+1) / static_cast<Real>(3*nx) );

                      Real xmin = 0.5*(1.0 - cmin);
                      Real xmax = 0.5*(1.0 - cmax);

                      x = (xmin + 2.*xmax)/3.;
                    }

                  }

                  mesh.add_point (Point(x,0.,0.), node_id++);
                }
                else
                mesh.add_point (Point(static_cast<Real>(i)/static_cast<Real>(3*nx),
                        0,
                        0), node_id++);
              }


            
              break;
            }
            
          default:
            {
              std::cerr << "ERROR: Unrecognized 1D element type." << std::endl;
              libmesh_error();
            }
              
        }

        // Build the elements of the mesh
        switch(type)
          {
            case INVALID_ELEM:
            case EDGE2:
              {
                for (unsigned int i=0; i<nx; i++)
                {
                  Elem* elem = mesh.add_elem (new Edge2);
                  elem->set_node(0) = mesh.node_ptr(i);
                  elem->set_node(1) = mesh.node_ptr(i+1);
                }
              break;
              }

            case EDGE3:
              {
                for (unsigned int i=0; i<nx; i++)
                {
                  Elem* elem = mesh.add_elem (new Edge3);
                  elem->set_node(0) = mesh.node_ptr(2*i);
                  elem->set_node(2) = mesh.node_ptr(2*i+1);
                  elem->set_node(1) = mesh.node_ptr(2*i+2);
                }
              break;
              }

            case EDGE4:
              {
                for (unsigned int i=0; i<nx; i++)
                {
                  Elem* elem = mesh.add_elem (new Edge4);
                  elem->set_node(0) = mesh.node_ptr(3*i);
                  elem->set_node(2) = mesh.node_ptr(3*i+1);
                  elem->set_node(3) = mesh.node_ptr(3*i+2);
                  elem->set_node(1) = mesh.node_ptr(3*i+3);
                }
              break;
              }

            default:
              {
                std::cerr << "ERROR: Unrecognized 1D element type." << std::endl;
                libmesh_error();                
              }
          }


	// Scale the nodal positions
	for (unsigned int p=0; p<mesh.n_nodes(); p++)
	  mesh.node(p)(0) = (mesh.node(p)(0))*(xmax-xmin) + xmin;
	
	break;
      }







      
      
      
      //---------------------------------------------------------------------
      // Build a 2D quadrilateral
    case 2:
      {
	libmesh_assert (nx != 0);
	libmesh_assert (ny != 0);

	// Reserve elements.  The TRI3 and TRI6 meshes
	// have twice as many elements...
	switch (type)
	  {
	  case INVALID_ELEM:
	  case QUAD4:
	  case QUAD8:
	  case QUAD9:
	    {
	      mesh.reserve_elem (nx*ny);
	      break;
	    }

	  case TRI3:
	  case TRI6:
	    {
	      mesh.reserve_elem (2*nx*ny);
	      break;
	    }

	  default:
	    {
	      std::cerr << "ERROR: Unrecognized 2D element type." << std::endl;
	      libmesh_error();
	    }
	  }


	
	// Reserve nodes.  The quadratic element types
	// need to reserve more nodes than the linear types.
	switch (type)
	  {
	  case INVALID_ELEM:
	  case QUAD4:
	  case TRI3:
	    {
	      mesh.reserve_nodes( (nx+1)*(ny+1) );
	      break;
	    }
	    
	  case QUAD8:
	  case QUAD9:
	  case TRI6:
	    {
	      mesh.reserve_nodes( (2*nx+1)*(2*ny+1) );
	      break;
	    }

	    
	  default:
	    {
	      std::cerr << "ERROR: Unrecognized 2D element type." << std::endl;
	      libmesh_error();
	    }
	  }



	// Build the nodes. Depends on whether you are using a linear
	// or quadratic element, and whether you are using a uniform
	// grid or the Gauss-Lobatto grid points.
        unsigned int node_id = 0;
	switch (type)
	  {
	  case INVALID_ELEM:
	  case QUAD4:
	  case TRI3:
	    {
	      for (unsigned int j=0; j<=ny; j++)
		for (unsigned int i=0; i<=nx; i++)
		  {
		    if (gauss_lobatto_grid)
		      {
			// Shortcut variable
			const Real pi = libMesh::pi;

			mesh.add_point (Point(0.5*(1.0 - std::cos(pi*static_cast<Real>(i)/static_cast<Real>(nx))),
					      0.5*(1.0 - std::cos(pi*static_cast<Real>(j)/static_cast<Real>(ny))),
					      0.), node_id++);
		      }
		  
		    else
		      mesh.add_point (Point(static_cast<Real>(i)/static_cast<Real>(nx),
					    static_cast<Real>(j)/static_cast<Real>(ny),
					    0.), node_id++);
		  }

	      break;
	    }
	    
	  case QUAD8:
	  case QUAD9:
	  case TRI6:
	    {
	      for (unsigned int j=0; j<=(2*ny); j++)
		for (unsigned int i=0; i<=(2*nx); i++)
		  {
		    if (gauss_lobatto_grid)
		      {
			// The x,y locations of the point.
			Real x=0., y=0.;

			// Shortcut variable
			const Real pi = libMesh::pi;
			
			// Shortcut quantities (do not depend on i,j)
			const Real a = std::cos( pi / static_cast<Real>(2*nx) );
			const Real b = std::cos( pi / static_cast<Real>(2*ny) );

			// Shortcut quantities (depend on i,j)
			const Real c = std::cos( pi*i / static_cast<Real>(2*nx) );
			const Real d = std::cos( pi*j / static_cast<Real>(2*ny) );
			
			// If i is even, compute a normal Gauss-Lobatto point
			if (i%2 == 0)
			  x = 0.5*(1.0 - c);

			// Otherwise, it is the average of the previous and next points
			else
			  x = 0.5*(1.0 - a*c);
		      
			// If j is even, compute a normal Gauss-Lobatto point
			if (j%2 == 0)
			  y = 0.5*(1.0 - d);
		      
			// Otherwise, it is the average of the previous and next points
			else
			  y = 0.5*(1.0 - b*d);
		      

			mesh.add_point (Point(x,y,0.), node_id++);
		      } 

		  
		    else
		      mesh.add_point (Point(static_cast<Real>(i)/static_cast<Real>(2*nx),
					    static_cast<Real>(j)/static_cast<Real>(2*ny),
					    0), node_id++);
		}

	      break;
	    }

	    
	  default:
	    {
	      std::cerr << "ERROR: Unrecognized 2D element type." << std::endl;
	      libmesh_error();
	    }
	  }






	// Build the elements.  Each one is a bit different.
	switch (type)
	  {

	  case INVALID_ELEM:
	  case QUAD4:
	    {
	      for (unsigned int j=0; j<ny; j++)
		for (unsigned int i=0; i<nx; i++)
		  {
		    Elem* elem = mesh.add_elem(new Quad4);

		    elem->set_node(0) = mesh.node_ptr(idx(type,nx,i,j)    );
		    elem->set_node(1) = mesh.node_ptr(idx(type,nx,i+1,j)  );
		    elem->set_node(2) = mesh.node_ptr(idx(type,nx,i+1,j+1));
		    elem->set_node(3) = mesh.node_ptr(idx(type,nx,i,j+1)  );

		    if (j == 0)
		      mesh.boundary_info->add_side(elem, 0, 0);
		  
		    if (j == (ny-1))
		      mesh.boundary_info->add_side(elem, 2, 2);
		  
		    if (i == 0)
		      mesh.boundary_info->add_side(elem, 3, 3);
		  
		    if (i == (nx-1))
		      mesh.boundary_info->add_side(elem, 1, 1);
		  }
	      break;
	    }


	  case TRI3:
	    {
	      for (unsigned int j=0; j<ny; j++)
		for (unsigned int i=0; i<nx; i++)
		  {
		    Elem* elem = NULL;

		    // Add first Tri3
		    elem = mesh.add_elem(new Tri3);

		    elem->set_node(0) = mesh.node_ptr(idx(type,nx,i,j)    );
		    elem->set_node(1) = mesh.node_ptr(idx(type,nx,i+1,j)  );
		    elem->set_node(2) = mesh.node_ptr(idx(type,nx,i+1,j+1));

		    if (j == 0)
		      mesh.boundary_info->add_side(elem, 0, 0);

		    if (i == (nx-1))
		      mesh.boundary_info->add_side(elem, 1, 1);

		    // Add second Tri3
		    elem = mesh.add_elem(new Tri3);

		    elem->set_node(0) = mesh.node_ptr(idx(type,nx,i,j)    );
		    elem->set_node(1) = mesh.node_ptr(idx(type,nx,i+1,j+1));
		    elem->set_node(2) = mesh.node_ptr(idx(type,nx,i,j+1)  );

		    if (j == (ny-1))
		      mesh.boundary_info->add_side(elem, 1, 2);

		    if (i == 0)
		      mesh.boundary_info->add_side(elem, 2, 3);
		  }
	      break;
	    }



	  case QUAD8:
	  case QUAD9:
	    {
	      for (unsigned int j=0; j<(2*ny); j += 2)
		for (unsigned int i=0; i<(2*nx); i += 2)
		  {
		    Elem* elem = (type == QUAD8) ?
		      mesh.add_elem(new Quad8) :
		      mesh.add_elem(new Quad9);
		  

		    elem->set_node(0) = mesh.node_ptr(idx(type,nx,i,j)    );
		    elem->set_node(1) = mesh.node_ptr(idx(type,nx,i+2,j)  );
		    elem->set_node(2) = mesh.node_ptr(idx(type,nx,i+2,j+2));
		    elem->set_node(3) = mesh.node_ptr(idx(type,nx,i,j+2)  );
		    elem->set_node(4) = mesh.node_ptr(idx(type,nx,i+1,j)  );
		    elem->set_node(5) = mesh.node_ptr(idx(type,nx,i+2,j+1));
		    elem->set_node(6) = mesh.node_ptr(idx(type,nx,i+1,j+2));
		    elem->set_node(7) = mesh.node_ptr(idx(type,nx,i,j+1)  );
		    if (type == QUAD9)
		      elem->set_node(8) = mesh.node_ptr(idx(type,nx,i+1,j+1));
		  

		    if (j == 0)
		      mesh.boundary_info->add_side(elem, 0, 0);
		  
		    if (j == 2*(ny-1))
		      mesh.boundary_info->add_side(elem, 2, 2);
		  
		    if (i == 0)
		      mesh.boundary_info->add_side(elem, 3, 3);
		  
		    if (i == 2*(nx-1))
		      mesh.boundary_info->add_side(elem, 1, 1);
		  }
	      break;
	    }


	  case TRI6:
	    {
	      for (unsigned int j=0; j<(2*ny); j += 2)
		for (unsigned int i=0; i<(2*nx); i += 2)
		  {
		    Elem* elem = NULL;

		    // Add first Tri6
		    elem = mesh.add_elem(new Tri6);

		    elem->set_node(0) = mesh.node_ptr(idx(type,nx,i,j)    );
		    elem->set_node(1) = mesh.node_ptr(idx(type,nx,i+2,j)  );
		    elem->set_node(2) = mesh.node_ptr(idx(type,nx,i+2,j+2));
		    elem->set_node(3) = mesh.node_ptr(idx(type,nx,i+1,j)  );
		    elem->set_node(4) = mesh.node_ptr(idx(type,nx,i+2,j+1));
		    elem->set_node(5) = mesh.node_ptr(idx(type,nx,i+1,j+1));

		    if (j == 0)
		      mesh.boundary_info->add_side(elem, 0, 0);

		    if (i == 2*(nx-1))
		      mesh.boundary_info->add_side(elem, 1, 1);

		    // Add second Tri6
		    elem = mesh.add_elem(new Tri6);

		    elem->set_node(0) = mesh.node_ptr(idx(type,nx,i,j)    );
		    elem->set_node(1) = mesh.node_ptr(idx(type,nx,i+2,j+2));
		    elem->set_node(2) = mesh.node_ptr(idx(type,nx,i,j+2)  );
		    elem->set_node(3) = mesh.node_ptr(idx(type,nx,i+1,j+1));
		    elem->set_node(4) = mesh.node_ptr(idx(type,nx,i+1,j+2));
		    elem->set_node(5) = mesh.node_ptr(idx(type,nx,i,j+1)  );

		    if (j == 2*(ny-1))
		      mesh.boundary_info->add_side(elem, 1, 2);

		    if (i == 0)
		      mesh.boundary_info->add_side(elem, 2, 3);

		  }
	      break;
	    };

	    
	  default:
	    {
	      std::cerr << "ERROR: Unrecognized 2D element type." << std::endl;
	      libmesh_error();
	    }
	  }



	
	// Scale the nodal positions
	for (unsigned int p=0; p<mesh.n_nodes(); p++)
	  {
	    mesh.node(p)(0) = (mesh.node(p)(0))*(xmax-xmin) + xmin;
	    mesh.node(p)(1) = (mesh.node(p)(1))*(ymax-ymin) + ymin;
	  }

	
	break;
      }








      


      //---------------------------------------------------------------------
      // Build a 3D mesh using hexahedral or prismatic elements.
    case 3:
      {
	libmesh_assert (nx != 0);
	libmesh_assert (ny != 0);
	libmesh_assert (nz != 0);


	// Reserve elements.  Meshes with prismatic elements require
	// twice as many elements.
	switch (type)
	  {
	  case INVALID_ELEM:
	  case HEX8:
	  case HEX20:
	  case HEX27:
	  case TET4:  // TET4's are created from an initial HEX27 discretization
	  case TET10: // TET10's are created from an initial HEX27 discretization
	    {
	      mesh.reserve_elem(nx*ny*nz);
	      break;
	    }

	  case PRISM6:
	  case PRISM15:
	  case PRISM18:
	    {
	      mesh.reserve_elem(2*nx*ny*nz);
	      break;
	    }

	  default:
	    {
	      std::cerr << "ERROR: Unrecognized 3D element type." << std::endl;
	      libmesh_error();
	    }
	  }




	
	// Reserve nodes.  Quadratic elements need twice as many nodes as linear elements.
	switch (type)
	  {
	  case INVALID_ELEM:
	  case HEX8:
	  case PRISM6:
	    {
	      mesh.reserve_nodes( (nx+1)*(ny+1)*(nz+1) );
	      break;
	    }

	  case HEX20:
	  case HEX27:
	  case TET4: // TET4's are created from an initial HEX27 discretization
	  case TET10: // TET10's are created from an initial HEX27 discretization
	  case PRISM15:
	  case PRISM18:
	    {
	      mesh.reserve_nodes( (2*nx+1)*(2*ny+1)*(2*nz+1) );
	      break;
	    }

	  default:
	    {
	      std::cerr << "ERROR: Unrecognized 3D element type." << std::endl;
	      libmesh_error();
	    }
	  }




	// Build the nodes.
        unsigned int node_id = 0;
	switch (type)
	  {
	  case INVALID_ELEM:
	  case HEX8:
	  case PRISM6:
	    {
	      for (unsigned int k=0; k<=nz; k++)
		for (unsigned int j=0; j<=ny; j++)
		  for (unsigned int i=0; i<=nx; i++)
		    {
		      if (gauss_lobatto_grid)
			{
			  // Shortcut variable
			  const Real pi = libMesh::pi;

			  mesh.add_point (Point(0.5*(1.0 - std::cos(pi*static_cast<Real>(i)/static_cast<Real>(nx))),
						0.5*(1.0 - std::cos(pi*static_cast<Real>(j)/static_cast<Real>(ny))),
						0.5*(1.0 - std::cos(pi*static_cast<Real>(k)/static_cast<Real>(nz)))), node_id++);
			}
		      
		      else
			mesh.add_point(Point(static_cast<Real>(i)/static_cast<Real>(nx),
					     static_cast<Real>(j)/static_cast<Real>(ny),
					     static_cast<Real>(k)/static_cast<Real>(nz)), node_id++);
		    }
	      break;
	    }
	    
	  case HEX20:
	  case HEX27:
	  case TET4: // TET4's are created from an initial HEX27 discretization
	  case TET10: // TET10's are created from an initial HEX27 discretization
	  case PRISM15:
	  case PRISM18:
	    {
	      for (unsigned int k=0; k<=(2*nz); k++)
		for (unsigned int j=0; j<=(2*ny); j++)
		  for (unsigned int i=0; i<=(2*nx); i++)
		    {
		      if (gauss_lobatto_grid)
			{
			  // Shortcut variable
			  const Real pi = libMesh::pi;
			  
			  // The x,y locations of the point.
			  Real x=0., y=0., z=0.;

			  // Shortcut quantities (do not depend on i,j)
			  const Real a = std::cos( pi / static_cast<Real>(2*nx) );
			  const Real b = std::cos( pi / static_cast<Real>(2*ny) );

			  // Shortcut quantities (depend on i,j)
			  const Real c = std::cos( pi*i / static_cast<Real>(2*nx) );
			  const Real d = std::cos( pi*j / static_cast<Real>(2*ny) );
			  
			  // Additional shortcut quantities (for 3D)
			  const Real e = std::cos( pi / static_cast<Real>(2*nz) );
			  const Real f = std::cos( pi*k / static_cast<Real>(2*nz) );
			  
			  // If i is even, compute a normal Gauss-Lobatto point
			  if (i%2 == 0)
			    x = 0.5*(1.0 - c);

			  // Otherwise, it is the average of the previous and next points
			  else
			    x = 0.5*(1.0 - a*c);
		      
			  // If j is even, compute a normal Gauss-Lobatto point
			  if (j%2 == 0)
			    y = 0.5*(1.0 - d);
		      
			  // Otherwise, it is the average of the previous and next points
			  else
			    y = 0.5*(1.0 - b*d);

			  // If k is even, compute a normal Gauss-Lobatto point
			  if (k%2 == 0)
			    z = 0.5*(1.0 - f);
		      
			  // Otherwise, it is the average of the previous and next points
			  else
			    z = 0.5*(1.0 - e*f);
		      
			  
			  mesh.add_point (Point(x,y,z), node_id++);
			}

		      else
			mesh.add_point(Point(static_cast<Real>(i)/static_cast<Real>(2*nx),
					     static_cast<Real>(j)/static_cast<Real>(2*ny),
					     static_cast<Real>(k)/static_cast<Real>(2*nz)), node_id++);
		    }
	      break;
	    }


	  default:
	    {
	      std::cerr << "ERROR: Unrecognized 3D element type." << std::endl;
	      libmesh_error();
	    }
	  }

	    

	
	// Build the elements.
	switch (type)
	  {
	  case INVALID_ELEM:
	  case HEX8:
	    {
	      for (unsigned int k=0; k<nz; k++)
		for (unsigned int j=0; j<ny; j++)
		  for (unsigned int i=0; i<nx; i++)
		    {
		      Elem* elem = mesh.add_elem(new Hex8);
			
		      elem->set_node(0) = mesh.node_ptr(idx(type,nx,ny,i,j,k)      );
		      elem->set_node(1) = mesh.node_ptr(idx(type,nx,ny,i+1,j,k)    );
		      elem->set_node(2) = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k)  );
		      elem->set_node(3) = mesh.node_ptr(idx(type,nx,ny,i,j+1,k)    );
		      elem->set_node(4) = mesh.node_ptr(idx(type,nx,ny,i,j,k+1)    );
		      elem->set_node(5) = mesh.node_ptr(idx(type,nx,ny,i+1,j,k+1)  );
		      elem->set_node(6) = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k+1));
		      elem->set_node(7) = mesh.node_ptr(idx(type,nx,ny,i,j+1,k+1)  );
			
		      if (k == 0)
			mesh.boundary_info->add_side(elem, 0, 0);
			
		      if (k == (nz-1))
			mesh.boundary_info->add_side(elem, 5, 5);
			
		      if (j == 0)
			mesh.boundary_info->add_side(elem, 1, 1);

		      if (j == (ny-1))
			mesh.boundary_info->add_side(elem, 3, 3);
 			
		      if (i == 0)
			mesh.boundary_info->add_side(elem, 4, 4);
 			
		      if (i == (nx-1))
			mesh.boundary_info->add_side(elem, 2, 2);
		    }
	      break;
	    }




	  case PRISM6:
	    {
	      for (unsigned int k=0; k<nz; k++)
		for (unsigned int j=0; j<ny; j++)
		  for (unsigned int i=0; i<nx; i++)
		    {
		      // First Prism
		      Elem* elem = NULL;
		      elem = mesh.add_elem(new Prism6);

		      elem->set_node(0) = mesh.node_ptr(idx(type,nx,ny,i,j,k)      );
		      elem->set_node(1) = mesh.node_ptr(idx(type,nx,ny,i+1,j,k)    );
		      elem->set_node(2) = mesh.node_ptr(idx(type,nx,ny,i,j+1,k)    );
		      elem->set_node(3) = mesh.node_ptr(idx(type,nx,ny,i,j,k+1)    );
		      elem->set_node(4) = mesh.node_ptr(idx(type,nx,ny,i+1,j,k+1)  );
		      elem->set_node(5) = mesh.node_ptr(idx(type,nx,ny,i,j+1,k+1)  );

		      // Add sides for first prism to boundary info object
		      if (i==0)
			mesh.boundary_info->add_side(elem, 3, 4);

		      if (j==0)
			mesh.boundary_info->add_side(elem, 1, 1);

		      if (k==0)
			mesh.boundary_info->add_side(elem, 0, 0);

		      if (k == (nz-1))
			mesh.boundary_info->add_side(elem, 4, 5);

		      // Second Prism
		      elem = mesh.add_elem(new Prism6);

		      elem->set_node(0) = mesh.node_ptr(idx(type,nx,ny,i+1,j,k)    );
		      elem->set_node(1) = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k)  );
		      elem->set_node(2) = mesh.node_ptr(idx(type,nx,ny,i,j+1,k)    );
		      elem->set_node(3) = mesh.node_ptr(idx(type,nx,ny,i+1,j,k+1)  );
		      elem->set_node(4) = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k+1));
		      elem->set_node(5) = mesh.node_ptr(idx(type,nx,ny,i,j+1,k+1)  );

		      // Add sides for second prism to boundary info object
		      if (i == (nx-1))
			mesh.boundary_info->add_side(elem, 1, 2);

		      if (j == (ny-1))
			mesh.boundary_info->add_side(elem, 2, 3);

		      if (k==0)
			mesh.boundary_info->add_side(elem, 0, 0);

		      if (k == (nz-1))
			mesh.boundary_info->add_side(elem, 4, 5);
		    }
	      break;
	    }




	    
	    
	  case HEX20:
	  case HEX27:
	  case TET4: // TET4's are created from an initial HEX27 discretization
	  case TET10: // TET10's are created from an initial HEX27 discretization
	    {
	      for (unsigned int k=0; k<(2*nz); k += 2)
		for (unsigned int j=0; j<(2*ny); j += 2)
		  for (unsigned int i=0; i<(2*nx); i += 2)
		    {
		      Elem* elem = (type == HEX20) ?
			mesh.add_elem(new Hex20) :
			mesh.add_elem(new Hex27);
		    
		      elem->set_node(0)  = mesh.node_ptr(idx(type,nx,ny,i,  j,  k)  );
		      elem->set_node(1)  = mesh.node_ptr(idx(type,nx,ny,i+2,j,  k)  );
		      elem->set_node(2)  = mesh.node_ptr(idx(type,nx,ny,i+2,j+2,k)  );
		      elem->set_node(3)  = mesh.node_ptr(idx(type,nx,ny,i,  j+2,k)  );
		      elem->set_node(4)  = mesh.node_ptr(idx(type,nx,ny,i,  j,  k+2));
		      elem->set_node(5)  = mesh.node_ptr(idx(type,nx,ny,i+2,j,  k+2));
		      elem->set_node(6)  = mesh.node_ptr(idx(type,nx,ny,i+2,j+2,k+2));
		      elem->set_node(7)  = mesh.node_ptr(idx(type,nx,ny,i,  j+2,k+2));
		      elem->set_node(8)  = mesh.node_ptr(idx(type,nx,ny,i+1,j,  k)  );
		      elem->set_node(9)  = mesh.node_ptr(idx(type,nx,ny,i+2,j+1,k)  );
		      elem->set_node(10) = mesh.node_ptr(idx(type,nx,ny,i+1,j+2,k)  );
		      elem->set_node(11) = mesh.node_ptr(idx(type,nx,ny,i,  j+1,k)  );
		      elem->set_node(12) = mesh.node_ptr(idx(type,nx,ny,i,  j,  k+1));
		      elem->set_node(13) = mesh.node_ptr(idx(type,nx,ny,i+2,j,  k+1));
		      elem->set_node(14) = mesh.node_ptr(idx(type,nx,ny,i+2,j+2,k+1));
		      elem->set_node(15) = mesh.node_ptr(idx(type,nx,ny,i,  j+2,k+1));
		      elem->set_node(16) = mesh.node_ptr(idx(type,nx,ny,i+1,j,  k+2));
		      elem->set_node(17) = mesh.node_ptr(idx(type,nx,ny,i+2,j+1,k+2));
		      elem->set_node(18) = mesh.node_ptr(idx(type,nx,ny,i+1,j+2,k+2));
		      elem->set_node(19) = mesh.node_ptr(idx(type,nx,ny,i,  j+1,k+2));
		      if ((type == HEX27) || (type == TET4) || (type == TET10))
			{
			  elem->set_node(20) = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k)  );
			  elem->set_node(21) = mesh.node_ptr(idx(type,nx,ny,i+1,j,  k+1));
			  elem->set_node(22) = mesh.node_ptr(idx(type,nx,ny,i+2,j+1,k+1));
			  elem->set_node(23) = mesh.node_ptr(idx(type,nx,ny,i+1,j+2,k+1));
			  elem->set_node(24) = mesh.node_ptr(idx(type,nx,ny,i,  j+1,k+1));
			  elem->set_node(25) = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k+2));
			  elem->set_node(26) = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k+1));
			}
			
			
		      if (k == 0)
			mesh.boundary_info->add_side(elem, 0, 0);
			
		      if (k == 2*(nz-1))
			mesh.boundary_info->add_side(elem, 5, 5);
			
		      if (j == 0)
			mesh.boundary_info->add_side(elem, 1, 1);

		      if (j == 2*(ny-1))
			mesh.boundary_info->add_side(elem, 3, 3);
 			
		      if (i == 0)
			mesh.boundary_info->add_side(elem, 4, 4);
 			
		      if (i == 2*(nx-1))
			mesh.boundary_info->add_side(elem, 2, 2);
		    }
	      break;
	    }



	    
	  case PRISM15:
	  case PRISM18:
	    {
	      for (unsigned int k=0; k<(2*nz); k += 2)
		for (unsigned int j=0; j<(2*ny); j += 2)
		  for (unsigned int i=0; i<(2*nx); i += 2)
		    {
		      // First Prism
		      Elem* elem = NULL;
		      elem = ((type == PRISM15) ?
			      mesh.add_elem(new Prism15) :
			      mesh.add_elem(new Prism18));

		      elem->set_node(0)  = mesh.node_ptr(idx(type,nx,ny,i,  j,  k)  );
		      elem->set_node(1)  = mesh.node_ptr(idx(type,nx,ny,i+2,j,  k)  );
		      elem->set_node(2)  = mesh.node_ptr(idx(type,nx,ny,i,  j+2,k)  );
		      elem->set_node(3)  = mesh.node_ptr(idx(type,nx,ny,i,  j,  k+2));
		      elem->set_node(4)  = mesh.node_ptr(idx(type,nx,ny,i+2,j,  k+2));
		      elem->set_node(5)  = mesh.node_ptr(idx(type,nx,ny,i,  j+2,k+2));
		      elem->set_node(6)  = mesh.node_ptr(idx(type,nx,ny,i+1,j,  k)  );
		      elem->set_node(7)  = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k)  );
		      elem->set_node(8)  = mesh.node_ptr(idx(type,nx,ny,i,  j+1,k)  );
		      elem->set_node(9)  = mesh.node_ptr(idx(type,nx,ny,i,  j,  k+1));
		      elem->set_node(10) = mesh.node_ptr(idx(type,nx,ny,i+2,j,  k+1));
		      elem->set_node(11) = mesh.node_ptr(idx(type,nx,ny,i,  j+2,k+1));
		      elem->set_node(12) = mesh.node_ptr(idx(type,nx,ny,i+1,j,  k+2));
		      elem->set_node(13) = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k+2));
		      elem->set_node(14) = mesh.node_ptr(idx(type,nx,ny,i,  j+1,k+2));
		      if (type == PRISM18)
			{
			  elem->set_node(15) = mesh.node_ptr(idx(type,nx,ny,i+1,j,  k+1));
			  elem->set_node(16) = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k+1));
			  elem->set_node(17) = mesh.node_ptr(idx(type,nx,ny,i,  j+1,k+1));
			}

		      // Add sides for first prism to boundary info object
		      if (i==0)
			mesh.boundary_info->add_side(elem, 3, 4);

		      if (j==0)
			mesh.boundary_info->add_side(elem, 1, 1);

		      if (k==0)
			mesh.boundary_info->add_side(elem, 0, 0);

		      if (k == 2*(nz-1))
			mesh.boundary_info->add_side(elem, 4, 5);

		      
		      // Second Prism
		      elem = ((type == PRISM15) ?
			      mesh.add_elem(new Prism15) :
			      mesh.add_elem(new Prism18));
			
		      elem->set_node(0)  = mesh.node_ptr(idx(type,nx,ny,i+2,j,k)     );
		      elem->set_node(1)  = mesh.node_ptr(idx(type,nx,ny,i+2,j+2,k)   );
		      elem->set_node(2)  = mesh.node_ptr(idx(type,nx,ny,i,j+2,k)     );
		      elem->set_node(3)  = mesh.node_ptr(idx(type,nx,ny,i+2,j,k+2)   );
		      elem->set_node(4)  = mesh.node_ptr(idx(type,nx,ny,i+2,j+2,k+2) );
		      elem->set_node(5)  = mesh.node_ptr(idx(type,nx,ny,i,j+2,k+2)   );
		      elem->set_node(6)  = mesh.node_ptr(idx(type,nx,ny,i+2,j+1,k)  );
		      elem->set_node(7)  = mesh.node_ptr(idx(type,nx,ny,i+1,j+2,k)  );
		      elem->set_node(8)  = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k)  );
		      elem->set_node(9)  = mesh.node_ptr(idx(type,nx,ny,i+2,j,k+1)  );
		      elem->set_node(10) = mesh.node_ptr(idx(type,nx,ny,i+2,j+2,k+1));
		      elem->set_node(11) = mesh.node_ptr(idx(type,nx,ny,i,j+2,k+1)  );
		      elem->set_node(12) = mesh.node_ptr(idx(type,nx,ny,i+2,j+1,k+2));
		      elem->set_node(13) = mesh.node_ptr(idx(type,nx,ny,i+1,j+2,k+2));
		      elem->set_node(14) = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k+2));
		      if (type == PRISM18)
			{
			  elem->set_node(15)  = mesh.node_ptr(idx(type,nx,ny,i+2,j+1,k+1));
			  elem->set_node(16)  = mesh.node_ptr(idx(type,nx,ny,i+1,j+2,k+1));
			  elem->set_node(17)  = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k+1));
			}
		      
		      // Add sides for second prism to boundary info object
		      if (i == 2*(nx-1))
			mesh.boundary_info->add_side(elem, 1, 2);

		      if (j == 2*(ny-1))
			mesh.boundary_info->add_side(elem, 2, 3);

		      if (k==0)
			mesh.boundary_info->add_side(elem, 0, 0);

		      if (k == 2*(nz-1))
			mesh.boundary_info->add_side(elem, 4, 5);

		    }
	      break;
	    }




	    
	  default:
	    {
	      std::cerr << "ERROR: Unrecognized 3D element type." << std::endl;
	      libmesh_error();
	    }
	  }

	    


	//.......................................
	// Scale the nodal positions
	for (unsigned int p=0; p<mesh.n_nodes(); p++)
	  {
	    mesh.node(p)(0) = (mesh.node(p)(0))*(xmax-xmin) + xmin;
	    mesh.node(p)(1) = (mesh.node(p)(1))*(ymax-ymin) + ymin;
	    mesh.node(p)(2) = (mesh.node(p)(2))*(zmax-zmin) + zmin;
	  }




	// Additional work for TET4 and TET10: we take the existing HEX27 discretization
	// and split each element into 24 sub tets.  This isn't the minimum-possible
	// number of tets, but it obviates any concerns about the edge orientations
	// between the various elements.
	if ((type == TET4) || (type == TET10)) 
	  {
	    // Temporary storage for new elements. (24 per hex)
	    std::vector<Elem*> new_elements;
	    new_elements.reserve(24*mesh.n_elem());

     
	    // Create tetrahedra
	    {
	      MeshBase::element_iterator       el     = mesh.elements_begin();
	      const MeshBase::element_iterator end_el = mesh.elements_end();

	      for ( ; el != end_el;  ++el)
		{
		  // Get a pointer to the HEX27 element.
		  Elem* base_hex = *el;

		  // Get a pointer to the node located at the HEX27 centroid
		  Node* apex_node = base_hex->get_node(26);
	 
		  for (unsigned int s=0; s<base_hex->n_sides(); ++s)
		    {
		      // Get the boundary ID for this side
		      short int b_id = mesh.boundary_info->boundary_id(*el, s);
		      
		      // Need to build the full-ordered side!
		      AutoPtr<Elem> side = base_hex->build_side(s);
	     
		      for (unsigned int sub_tet=0; sub_tet<4; ++sub_tet)
			{
			  new_elements.push_back( new Tet4 );
			  Elem* sub_elem = new_elements.back();
			  sub_elem->set_node(0) = side->get_node(sub_tet);
			  sub_elem->set_node(1) = side->get_node(8);                           // centroid of the face
			  sub_elem->set_node(2) = side->get_node(sub_tet==3 ? 0 : sub_tet+1 ); // wrap-around
			  sub_elem->set_node(3) = apex_node;                                   // apex node always used!

			  // If the original hex was a boundary hex, add the new sub_tet's side
			  // 0 with the same b_id.  Note: the tets are all aligned so that their
			  // side 0 is on the boundary.
			  if (b_id != BoundaryInfo::invalid_id)
			    mesh.boundary_info->add_side(sub_elem, 0, b_id);
			}
		    }
		}
	    }
	    

	    // Delete the original HEX27 elements from the mesh, and the boundary info structure.
	    {
	      MeshBase::element_iterator       el     = mesh.elements_begin();
	      const MeshBase::element_iterator end_el = mesh.elements_end();
       
	      for ( ; el != end_el;  ++el)
		{
		  mesh.boundary_info->remove(*el); // Safe even if *el has no boundary info.
		  mesh.delete_elem(*el);
		}
	    }

	    // Add the new elements
	    for (unsigned int i=0; i<new_elements.size(); ++i)
	      mesh.add_elem(new_elements[i]);
	    
	  } // end if (type == TET4) || (type == TET10)


	// Use all_second_order to convert the TET4's to TET10's
	if (type == TET10)
	  {
	    mesh.all_second_order();
	  }

	
	break;
      } // end case dim==3

    default:
      {
	libmesh_error();
      }
    }  

  STOP_LOG("build_cube()", "MeshTools::Generation");


  
  // Done building the mesh.  Now prepare it for use.
  mesh.prepare_for_use ();  
}



void MeshTools::Generation::build_line (UnstructuredMesh& mesh,
                                        const unsigned int nx,
                                        const Real xmin, const Real xmax,
                                        const ElemType type,
                                        const bool gauss_lobatto_grid)
{
    // This method only makes sense in 1D!
    libmesh_assert(mesh.mesh_dimension() == 1);

    build_cube(mesh,
               nx, 0, 0,
               xmin, xmax,
               0., 0.,
               0., 0.,
               type,
               gauss_lobatto_grid);
}
                                        


void MeshTools::Generation::build_square (UnstructuredMesh& mesh,
					  const unsigned int nx,
					  const unsigned int ny,
					  const Real xmin, const Real xmax,
					  const Real ymin, const Real ymax,
					  const ElemType type,
					  const bool gauss_lobatto_grid)
{
  // This method only makes sense in 2D!
  libmesh_assert (mesh.mesh_dimension() == 2);

  // Call the build_cube() member to actually do the work for us.
  build_cube (mesh,
	      nx, ny, 0,
	      xmin, xmax,
	      ymin, ymax,
	      0., 0.,
	      type,
	      gauss_lobatto_grid);
}









#ifndef LIBMESH_ENABLE_AMR
void MeshTools::Generation::build_sphere (UnstructuredMesh&,
					  const Real,
					  const unsigned int,
					  const ElemType)
{
 	std::cout << "Building a circle/sphere only works with AMR." << std::endl;
 	libmesh_error();
}

#else
	
void MeshTools::Generation::build_sphere (UnstructuredMesh& mesh,
					  const Real rad,
					  const unsigned int nr,
					  const ElemType type)
{
  libmesh_assert (mesh.mesh_dimension() != 1);
  libmesh_assert (rad > 0.);
  //libmesh_assert (nr > 0); // must refine at least once otherwise will end up with a square/cube
  
  START_LOG("build_sphere()", "MeshTools::Generation");

  // Clear the mesh and start from scratch
  mesh.clear();
  
  // Sphere is centered at origin by default
  const Point cent;

  const Sphere sphere (cent, rad);
  
  switch (mesh.mesh_dimension())
    {
      //-----------------------------------------------------------------
      // Build a circle in two dimensions
    case 2:
      {
	const Real sqrt_2     = std::sqrt(2.);
        const Real rad_2      = .25*rad;
        const Real rad_sqrt_2 = rad/sqrt_2;

	// (Temporary) convenient storage for node pointers
	std::vector<Node*> nodes(8);
	    
	// Point 0
	nodes[0] = mesh.add_point (Point(-rad_2,-rad_2, 0.));
	    
	// Point 1
	nodes[1] = mesh.add_point (Point( rad_2,-rad_2, 0.));
	    
	// Point 2
	nodes[2] = mesh.add_point (Point( rad_2, rad_2, 0.));

	// Point 3
	nodes[3] = mesh.add_point (Point(-rad_2, rad_2, 0.));
	    
	// Point 4
	nodes[4] = mesh.add_point (Point(-rad_sqrt_2,-rad_sqrt_2, 0.));
	    
	// Point 5
	nodes[5] = mesh.add_point (Point( rad_sqrt_2,-rad_sqrt_2, 0.));
	    
	// Point 6
	nodes[6] = mesh.add_point (Point( rad_sqrt_2, rad_sqrt_2, 0.));
	    
	// Point 7
	nodes[7] = mesh.add_point (Point(-rad_sqrt_2, rad_sqrt_2, 0.));

	// Build the elements & set node pointers
	    
	// Element 0
	Elem* elem0 = mesh.add_elem (new Quad4);
	elem0->set_node(0) = nodes[0];
	elem0->set_node(1) = nodes[1];
	elem0->set_node(2) = nodes[2];
	elem0->set_node(3) = nodes[3];
	    
	// Element 1
	Elem* elem1 = mesh.add_elem (new Quad4);
	elem1->set_node(0) = nodes[4];
	elem1->set_node(1) = nodes[0];
	elem1->set_node(2) = nodes[3];
	elem1->set_node(3) = nodes[7];
	    
	// Element 2
	Elem* elem2 = mesh.add_elem (new Quad4);
	elem2->set_node(0) = nodes[4];
	elem2->set_node(1) = nodes[5];
	elem2->set_node(2) = nodes[1];
	elem2->set_node(3) = nodes[0];
	    
	// Element 3
	Elem* elem3 = mesh.add_elem (new Quad4);
	elem3->set_node(0) = nodes[1];
	elem3->set_node(1) = nodes[5];
	elem3->set_node(2) = nodes[6];
	elem3->set_node(3) = nodes[2];
	    
	// Element 4
	Elem* elem4 = mesh.add_elem (new Quad4);
	elem4->set_node(0) = nodes[3];
	elem4->set_node(1) = nodes[2];
	elem4->set_node(2) = nodes[6];
	elem4->set_node(3) = nodes[7];

	break;
      } // end case 2


      

      
      //-----------------------------------------------------------------
      // Build a sphere in three dimensions
    case 3:
      {
	// (Currently) supported types
	if (!((type == HEX8) || (type == HEX27)))
	  {
	    // FIXME: We'd need an all_tet() routine (which could also be used by
	    // build_square()) to do Tets, or Prisms for that matter.
	    std::cerr << "Error: Only HEX8/27 currently supported."
		      << std::endl;
	    libmesh_error();
	  }

	
	// 3D analog of 2D initial grid:
	const Real
	  r_small = 0.25*rad,                      //  0.25 *radius
	  r_med   = (0.125*std::sqrt(2.)+0.5)*rad; // .67677*radius 
	
	// (Temporary) convenient storage for node pointers
	std::vector<Node*> nodes(16);
	
	// Points 0-7 are the initial HEX8
	nodes[0] = mesh.add_point (Point(-r_small,-r_small, -r_small));
	nodes[1] = mesh.add_point (Point( r_small,-r_small, -r_small));
	nodes[2] = mesh.add_point (Point( r_small, r_small, -r_small));
	nodes[3] = mesh.add_point (Point(-r_small, r_small, -r_small));
	nodes[4] = mesh.add_point (Point(-r_small,-r_small,  r_small));
	nodes[5] = mesh.add_point (Point( r_small,-r_small,  r_small));
	nodes[6] = mesh.add_point (Point( r_small, r_small,  r_small));
	nodes[7] = mesh.add_point (Point(-r_small, r_small,  r_small));

	//  Points 8-15 are for the outer hexes, we number them in the same way
	nodes[8]  = mesh.add_point (Point(-r_med,-r_med, -r_med));
	nodes[9]  = mesh.add_point (Point( r_med,-r_med, -r_med));
	nodes[10] = mesh.add_point (Point( r_med, r_med, -r_med));
	nodes[11] = mesh.add_point (Point(-r_med, r_med, -r_med));
	nodes[12] = mesh.add_point (Point(-r_med,-r_med,  r_med));
	nodes[13] = mesh.add_point (Point( r_med,-r_med,  r_med));
	nodes[14] = mesh.add_point (Point( r_med, r_med,  r_med));
	nodes[15] = mesh.add_point (Point(-r_med, r_med,  r_med));

	// Now create the elements and add them to the mesh
	// Element 0 - center element
	{
	  Elem* elem0 = mesh.add_elem (new Hex8);
	  elem0->set_node(0) = nodes[0];
	  elem0->set_node(1) = nodes[1];
	  elem0->set_node(2) = nodes[2];
	  elem0->set_node(3) = nodes[3];
	  elem0->set_node(4) = nodes[4];
	  elem0->set_node(5) = nodes[5];
	  elem0->set_node(6) = nodes[6];
	  elem0->set_node(7) = nodes[7];
	}
	
	// Element 1 - "bottom"
	{
	  Elem* elem1 = mesh.add_elem (new Hex8);
	  elem1->set_node(0) = nodes[8];
	  elem1->set_node(1) = nodes[9];
	  elem1->set_node(2) = nodes[10];
	  elem1->set_node(3) = nodes[11];
	  elem1->set_node(4) = nodes[0];
	  elem1->set_node(5) = nodes[1];
	  elem1->set_node(6) = nodes[2];
	  elem1->set_node(7) = nodes[3];
	}
	
	// Element 2 - "front"
	{
	  Elem* elem2 = mesh.add_elem (new Hex8);
	  elem2->set_node(0) = nodes[8];
	  elem2->set_node(1) = nodes[9];
	  elem2->set_node(2) = nodes[1];
	  elem2->set_node(3) = nodes[0];
	  elem2->set_node(4) = nodes[12];
	  elem2->set_node(5) = nodes[13];
	  elem2->set_node(6) = nodes[5];
	  elem2->set_node(7) = nodes[4];
	}

	// Element 3 - "right"
	{
	  Elem* elem3 = mesh.add_elem (new Hex8);
	  elem3->set_node(0) = nodes[1];
	  elem3->set_node(1) = nodes[9];
	  elem3->set_node(2) = nodes[10];
	  elem3->set_node(3) = nodes[2];
	  elem3->set_node(4) = nodes[5];
	  elem3->set_node(5) = nodes[13];
	  elem3->set_node(6) = nodes[14];
	  elem3->set_node(7) = nodes[6];
	}

	// Element 4 - "back"
	{
	  Elem* elem4 = mesh.add_elem (new Hex8);
	  elem4->set_node(0) = nodes[3];
	  elem4->set_node(1) = nodes[2];
	  elem4->set_node(2) = nodes[10];
	  elem4->set_node(3) = nodes[11];
	  elem4->set_node(4) = nodes[7];
	  elem4->set_node(5) = nodes[6];
	  elem4->set_node(6) = nodes[14];
	  elem4->set_node(7) = nodes[15];
	}

	// Element 5 - "left"
	{
	  Elem* elem5 = mesh.add_elem (new Hex8);
	  elem5->set_node(0) = nodes[8];
	  elem5->set_node(1) = nodes[0];
	  elem5->set_node(2) = nodes[3];
	  elem5->set_node(3) = nodes[11];
	  elem5->set_node(4) = nodes[12];
	  elem5->set_node(5) = nodes[4];
	  elem5->set_node(6) = nodes[7];
	  elem5->set_node(7) = nodes[15];
	}

	// Element 6 - "top"
	{
	  Elem* elem6 = mesh.add_elem (new Hex8);
	  elem6->set_node(0) = nodes[4];
	  elem6->set_node(1) = nodes[5];
	  elem6->set_node(2) = nodes[6];
	  elem6->set_node(3) = nodes[7];
	  elem6->set_node(4) = nodes[12];
	  elem6->set_node(5) = nodes[13];
	  elem6->set_node(6) = nodes[14];
	  elem6->set_node(7) = nodes[15];
	}
	
	break;
      } // end case 3
      
    default:
      libmesh_error();
      
    } // end switch (dim)

  
	
  // Now we have the beginnings of a sphere.
  // Add some more elements by doing uniform refinements and
  // popping nodes to the boundary.
  MeshRefinement mesh_refinement (mesh);

  // Loop over the elements, refine, pop nodes to boundary.
  for (unsigned int r=0; r<nr; r++)
    {
      mesh_refinement.uniformly_refine(1);

      MeshBase::element_iterator       it  = mesh.active_elements_begin();
      const MeshBase::element_iterator end = mesh.active_elements_end(); 

      for (; it != end; ++it)
	{
	  Elem* elem = *it;
		
	  for (unsigned int s=0; s<elem->n_sides(); s++)
	    if (elem->neighbor(s) == NULL)
	      {
		AutoPtr<Elem> side(elem->build_side(s));

		// Pop each point to the sphere boundary
		for (unsigned int n=0; n<side->n_nodes(); n++)
		  side->point(n) =
		    sphere.closest_point(side->point(n));
	      }
	}
      }

  // The mesh now contains a refinement hierarchy due to the refinements
  // used to generate the grid.  In order to call other support functions
  // like all_tri() and all_second_order, you need a "flat" mesh file (with no
  // refinement trees) so
  MeshTools::Modification::flatten(mesh);

  // In 2D, convert all the quads to triangles if requested
  if (mesh.mesh_dimension()==2)
    {
      if ((type == TRI6) || (type == TRI3))
	{
	  MeshTools::Modification::all_tri(mesh);
	}
    }

  
  // Convert to second-order elements if the user requested it.
  if (Elem::second_order_equivalent_type(type) == INVALID_ELEM)
    {
      // type is already a second-order, determine if it is the
      // "full-ordered" second-order element, or the "serendipity"
      // second order element.  Note also that all_second_order
      // can't be called once the mesh has been refined.
      bool full_ordered = !((type==QUAD8) || (type==HEX20));
      mesh.all_second_order(full_ordered);

      // And pop to the boundary again...
      MeshBase::element_iterator       it  = mesh.active_elements_begin();
      const MeshBase::element_iterator end = mesh.active_elements_end(); 

      for (; it != end; ++it)
	{
	  Elem* elem = *it;
		
	  for (unsigned int s=0; s<elem->n_sides(); s++)
	    if (elem->neighbor(s) == NULL)
	      {
		AutoPtr<Elem> side(elem->build_side(s));

		// Pop each point to the sphere boundary
		for (unsigned int n=0; n<side->n_nodes(); n++)
		  side->point(n) =
		    sphere.closest_point(side->point(n));
	      }
	}
    }
  

  // The meshes could probably use some smoothing...
  LaplaceMeshSmoother smoother(mesh);
  smoother.smooth(2);
  
  STOP_LOG("build_sphere()", "MeshTools::Generation");

  
  // Done building the mesh.  Now prepare it for use.
  mesh.prepare_for_use();
}

#endif // #ifndef LIBMESH_ENABLE_AMR

