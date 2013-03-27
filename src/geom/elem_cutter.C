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

// C++ includes


namespace libMesh
{


  // ------------------------------------------------------------
  // ElemCutter implementation
  ElemCutter::ElemCutter()
  {
    _working_mesh_2D.reset (new SerialMesh);
    _working_mesh_3D.reset (new SerialMesh);
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
	  _outside_elem.push_back(&elem);
	  return;
	}

      // completely inside?
      else if (nmax <= 0.)
	{
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
      case 2: this->cut_2D(); break;
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
	   
	    const Point x_star = (edge->point(0)*(1.-d_star) +
				  edge->point(1)*d_star);
	    
	    _intersection_pts.push_back (x_star);
	  }
      }
  }



  
  void ElemCutter::cut_1D ()
  {
    libmesh_not_implemented();
  }

  

  void ElemCutter::cut_2D ()
  {
    libmesh_not_implemented();
  }

  

  void ElemCutter::cut_3D ()
  {
    libmesh_not_implemented();
  }

  

} // namespace libMesh
