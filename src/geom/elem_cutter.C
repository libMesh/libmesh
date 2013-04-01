// The libMesh Finite Element Library.
// Copyright (C) 2002-2013 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/elem_cutter.h"
#include "libmesh/elem.h"
#include "libmesh/serial_mesh.h"
#include "libmesh/mesh_triangle_interface.h"
#include "libmesh/mesh_tetgen_interface.h"

// C++ includes

namespace
{
  unsigned int cut_cntr;
}
namespace libMesh
{


  // ------------------------------------------------------------
  // ElemCutter implementation
  ElemCutter::ElemCutter()
  {
    _inside_mesh_2D.reset  (new SerialMesh); /**/ _triangle_inside.reset  (new TriangleInterface (*_inside_mesh_2D));   
    _outside_mesh_2D.reset (new SerialMesh); /**/ _triangle_outside.reset (new TriangleInterface (*_outside_mesh_2D));
    
    _inside_mesh_3D.reset  (new SerialMesh); /**/ _tetgen_inside.reset  (new TetGenMeshInterface (*_inside_mesh_3D));   
    _outside_mesh_3D.reset (new SerialMesh); /**/ _tetgen_outside.reset (new TetGenMeshInterface (*_outside_mesh_3D));

    cut_cntr = 0;
  }


  
  ElemCutter::~ElemCutter()
  {}


  
  void ElemCutter::operator()(const Elem &elem,
			      const std::vector<Real> &vertex_distance_func)

  {
    libmesh_assert_equal_to (vertex_distance_func.size(), elem.n_vertices());

    _inside_elem.clear();
    _outside_elem.clear();
  
    // check for quick return?
    {
      Real
	nmin = vertex_distance_func.front(),
	nmax = nmin;

      for (std::vector<Real>::const_iterator it=vertex_distance_func.begin();
	   it!=vertex_distance_func.end(); ++it)
	{
	  nmin = std::min (nmin, *it);
	  nmax = std::max (nmax, *it);
	}

      // completely outside?
      if (nmin >= 0.)
	{
	  //std::cout << "element completely outside\n";
	  _outside_elem.push_back(&elem);
	  return;
	}

      // completely inside?
      else if (nmax <= 0.)
	{
	  //std::cout << "element completely inside\n";
	  _inside_elem.push_back(&elem);
	  return;
	}

      libmesh_assert_greater (nmax, 0.);
      libmesh_assert_less    (nmin, 0.);
    }


    
    // we now know we are in a cut element, find the intersecting points.
    this->find_intersection_points (elem, vertex_distance_func);
    
    // and then dispatch the proper method
    switch (elem.dim())
      {
      case 1: this->cut_1D(); break;
      case 2: this->cut_2D(elem, vertex_distance_func); break;
      case 3: this->cut_3D(); break;
      default: libmesh_error();
      }
  }

  

  void ElemCutter::find_intersection_points(const Elem &elem,
					    const std::vector<Real> &vertex_distance_func)
  {
    _intersection_pts.clear();

    for (unsigned int e=0; e<elem.n_edges(); e++)
      {
	AutoPtr<Elem> edge (elem.build_edge(e));

	// find the element nodes el0, el1 that map
	unsigned int
	  el0 = elem.get_node_index(edge->get_node(0)),
	  el1 = elem.get_node_index(edge->get_node(1));

	libmesh_assert (elem.is_vertex(el0));
	libmesh_assert (elem.is_vertex(el1));
	libmesh_assert_less (el0, vertex_distance_func.size());
	libmesh_assert_less (el1, vertex_distance_func.size());

	const Real
	  d0 = vertex_distance_func[el0],
	  d1 = vertex_distance_func[el1];

	// if this egde has a 0 crossing
	if (d0*d1 < 0.)
	  {
	    libmesh_assert_not_equal_to (d0, d1);
	    
	    // then find d_star in [0,1], the
	    // distance from el0 to el1 where the 0 lives.
	    const Real d_star = d0 / (d0 - d1);


	    // Prevent adding nodes trivially close to existing
	    // nodes.
	    const Real endpoint_tol = 0.01;
	    
	    if ( (d_star > endpoint_tol) &&
		 (d_star < (1.-endpoint_tol)) )
	      {
		const Point x_star = (edge->point(0)*(1.-d_star) +
				      edge->point(1)*d_star);
		
		std::cout << "adding cut point (d_star, x_star) = "
			  << d_star << " , " << x_star << std::endl;
		
		_intersection_pts.push_back (x_star);
	      }
	  }
      }
  }



  
  void ElemCutter::cut_1D ()
  {
    libmesh_not_implemented();
  }

  

  void ElemCutter::cut_2D (const Elem &elem,
			   const std::vector<Real> &vertex_distance_func)
  {
    //libmesh_not_implemented();
    std::cout << "Inside cut element!\n";

#ifdef LIBMESH_HAVE_TRIANGLE
    
    libmesh_assert (_inside_mesh_2D.get() != NULL);
    
    _inside_mesh_2D->clear();
    _outside_mesh_2D->clear();

    for (unsigned int v=0; v<elem.n_vertices(); v++)
      {
	if (vertex_distance_func[v] >= 0.)
	  _outside_mesh_2D->add_point (elem.point(v));
	
	if (vertex_distance_func[v] <= 0.)
	  _inside_mesh_2D->add_point (elem.point(v));
      }

    for (std::vector<Point>::const_iterator it=_intersection_pts.begin();
	 it != _intersection_pts.end(); ++it)
      {	
	_inside_mesh_2D->add_point(*it);
	_outside_mesh_2D->add_point(*it);
      }

    
    // Customize the variables for the triangulation
    // we will be cutting reference cell, and want as few triangles
    // as possible, so jack this up larger than the area we will be
    // triangulating so we are governed only by accurately defining
    // the boundaries.
    _triangle_inside->desired_area()  = 100.;
    _triangle_outside->desired_area() = 100.;

    // allow for small angles
    _triangle_inside->minimum_angle()  = 5.;
    _triangle_outside->minimum_angle() = 5.;
    
    // Turn off Laplacian mesh smoothing after generation.
    _triangle_inside->smooth_after_generating()  = false;
    _triangle_outside->smooth_after_generating() = false;
    
    // Triangulate!
    _triangle_inside->triangulate();
    _triangle_outside->triangulate();

    std::ostringstream name;
    
    name << "cut_face_"
	 << cut_cntr++
	 << ".dat";
    _inside_mesh_2D->write  ("in_"  + name.str());
    _outside_mesh_2D->write ("out_" + name.str());
#endif
  }

  

  void ElemCutter::cut_3D ()
  {
    libmesh_not_implemented();
  }

  

} // namespace libMesh
