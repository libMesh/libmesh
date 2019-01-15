// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// Local includes
#include "libmesh/elem_cutter.h"
#include "libmesh/elem.h"
#include "libmesh/replicated_mesh.h"
#include "libmesh/mesh_triangle_interface.h"
#include "libmesh/mesh_tetgen_interface.h"
#include "libmesh/auto_ptr.h" // libmesh_make_unique

namespace
{
unsigned int cut_cntr;
}

namespace libMesh
{

ElemCutter::ElemCutter() :
  _inside_mesh_2D(libmesh_make_unique<ReplicatedMesh>(_comm_self,2)),
  _triangle_inside(libmesh_make_unique<TriangleInterface>(*_inside_mesh_2D)),
  _outside_mesh_2D(libmesh_make_unique<ReplicatedMesh>(_comm_self,2)),
  _triangle_outside(libmesh_make_unique<TriangleInterface>(*_outside_mesh_2D)),
  _inside_mesh_3D(libmesh_make_unique<ReplicatedMesh>(_comm_self,3)),
  _tetgen_inside(libmesh_make_unique<TetGenMeshInterface>(*_inside_mesh_3D)),
  _outside_mesh_3D(libmesh_make_unique<ReplicatedMesh>(_comm_self,3)),
  _tetgen_outside(libmesh_make_unique<TetGenMeshInterface>(*_outside_mesh_3D))
{
  cut_cntr = 0;
}



ElemCutter::~ElemCutter()
{}



bool ElemCutter::is_inside (const Elem & libmesh_dbg_var(elem),
                            const std::vector<Real> & vertex_distance_func) const
{
  libmesh_assert_equal_to (elem.n_vertices(), vertex_distance_func.size());

  for (const auto & val : vertex_distance_func)
    if (val > 0.)
      return false;

  // if the distance function is nonpositive, we are outside
  return true;
}



bool ElemCutter::is_outside (const Elem & libmesh_dbg_var(elem),
                             const std::vector<Real> & vertex_distance_func) const
{
  libmesh_assert_equal_to (elem.n_vertices(), vertex_distance_func.size());

  for (const auto & val : vertex_distance_func)
    if (val < 0.)
      return false;

  // if the distance function is nonnegative, we are outside
  return true;
}



bool ElemCutter::is_cut (const Elem & libmesh_dbg_var(elem),
                         const std::vector<Real> & vertex_distance_func) const
{
  libmesh_assert_equal_to (elem.n_vertices(), vertex_distance_func.size());

  Real
    vmin = vertex_distance_func.front(),
    vmax = vmin;

  for (const auto & val : vertex_distance_func)
    {
      vmin = std::min (vmin, val);
      vmax = std::max (vmax, val);
    }

  // if the distance function changes sign, we're cut.
  return (vmin*vmax < 0.);
}



void ElemCutter::operator()(const Elem & elem,
                            const std::vector<Real> & vertex_distance_func)

{
  libmesh_assert_equal_to (vertex_distance_func.size(), elem.n_vertices());

  _inside_elem.clear();
  _outside_elem.clear();

  // check for quick return?
  {
    // completely outside?
    if (this->is_outside(elem, vertex_distance_func))
      {
        //std::cout << "element completely outside\n";
        _outside_elem.push_back(& elem);
        return;
      }

    // completely inside?
    else if (this->is_inside(elem, vertex_distance_func))
      {
        //std::cout << "element completely inside\n";
        _inside_elem.push_back(&elem);
        return;
      }

    libmesh_assert (this->is_cut (elem, vertex_distance_func));
  }

  // we now know we are in a cut element, find the intersecting points.
  this->find_intersection_points (elem, vertex_distance_func);

  // and then dispatch the proper method
  switch (elem.dim())
    {
    case 1: this->cut_1D(elem, vertex_distance_func); break;
    case 2: this->cut_2D(elem, vertex_distance_func); break;
    case 3: this->cut_3D(elem, vertex_distance_func); break;
    default: libmesh_error_msg("Invalid element dimension: " << elem.dim());
    }
}



void ElemCutter::find_intersection_points(const Elem & elem,
                                          const std::vector<Real> & vertex_distance_func)
{
  _intersection_pts.clear();

  for (unsigned int e=0; e<elem.n_edges(); e++)
    {
      std::unique_ptr<const Elem> edge (elem.build_edge_ptr(e));

      // find the element nodes el0, el1 that map
      unsigned int
        el0 = elem.get_node_index(edge->node_ptr(0)),
        el1 = elem.get_node_index(edge->node_ptr(1));

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
              const Point x_star = (edge->point(0)*(1-d_star) +
                                    edge->point(1)*d_star);

              std::cout << "adding cut point (d_star, x_star) = "
                        << d_star << " , " << x_star << std::endl;

              _intersection_pts.push_back (x_star);
            }
        }
    }
}




void ElemCutter::cut_1D (const Elem & /*elem*/,
                         const std::vector<Real> &/*vertex_distance_func*/)
{
  libmesh_not_implemented();
}



void ElemCutter::cut_2D (const Elem & elem,
                         const std::vector<Real> & vertex_distance_func)
{
#ifndef LIBMESH_HAVE_TRIANGLE

  // current implementation requires triangle!
  libMesh::err << "ERROR: current libMesh ElemCutter 2D implementation requires\n"
               << "       the \"triangle\" library!\n"
               << std::endl;
  libmesh_not_implemented();

#else // OK, LIBMESH_HAVE_TRIANGLE

  std::cout << "Inside cut face element!\n";

  libmesh_assert (_inside_mesh_2D.get()  != nullptr);
  libmesh_assert (_outside_mesh_2D.get() != nullptr);

  _inside_mesh_2D->clear();
  _outside_mesh_2D->clear();

  for (unsigned int v=0; v<elem.n_vertices(); v++)
    {
      if (vertex_distance_func[v] >= 0.)
        _outside_mesh_2D->add_point (elem.point(v));

      if (vertex_distance_func[v] <= 0.)
        _inside_mesh_2D->add_point (elem.point(v));
    }

  for (const auto & pt : _intersection_pts)
    {
      _inside_mesh_2D->add_point(pt);
      _outside_mesh_2D->add_point(pt);
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

  // std::ostringstream name;

  // name << "cut_face_"
  //  << cut_cntr++
  //  << ".dat";
  // _inside_mesh_2D->write  ("in_"  + name.str());
  // _outside_mesh_2D->write ("out_" + name.str());

  // finally, add the elements to our lists.
  _inside_elem.clear();
  _outside_elem.clear();

  for (const auto & el : _inside_mesh_2D->element_ptr_range())
    _inside_elem.push_back (el);

  for (const auto & el : _outside_mesh_2D->element_ptr_range())
    _outside_elem.push_back (el);

#endif
}



void ElemCutter::cut_3D (const Elem & elem,
                         const std::vector<Real> & vertex_distance_func)
{
#ifndef LIBMESH_HAVE_TETGEN

  // current implementation requires tetgen!
  libMesh::err << "ERROR: current libMesh ElemCutter 3D implementation requires\n"
               << "       the \"tetgen\" library!\n"
               << std::endl;
  libmesh_not_implemented();

#else // OK, LIBMESH_HAVE_TETGEN

  std::cout << "Inside cut cell element!\n";

  libmesh_assert (_inside_mesh_3D.get()  != nullptr);
  libmesh_assert (_outside_mesh_3D.get() != nullptr);

  _inside_mesh_3D->clear();
  _outside_mesh_3D->clear();

  for (unsigned int v=0; v<elem.n_vertices(); v++)
    {
      if (vertex_distance_func[v] >= 0.)
        _outside_mesh_3D->add_point (elem.point(v));

      if (vertex_distance_func[v] <= 0.)
        _inside_mesh_3D->add_point (elem.point(v));
    }

  for (const auto & pt : _intersection_pts)
    {
      _inside_mesh_3D->add_point(pt);
      _outside_mesh_3D->add_point(pt);
    }


  // Triangulate!
  _tetgen_inside->triangulate_pointset();
  //_inside_mesh_3D->print_info();
  _tetgen_outside->triangulate_pointset();
  //_outside_mesh_3D->print_info();


  // (below generates some horribly expensive meshes,
  //  but seems immune to the 0 volume problem).
  // _tetgen_inside->pointset_convexhull();
  // _inside_mesh_3D->find_neighbors();
  // _inside_mesh_3D->print_info();
  // _tetgen_inside->triangulate_conformingDelaunayMesh (1.e3, 100.);
  // _inside_mesh_3D->print_info();

  // _tetgen_outside->pointset_convexhull();
  // _outside_mesh_3D->find_neighbors();
  // _outside_mesh_3D->print_info();
  // _tetgen_outside->triangulate_conformingDelaunayMesh (1.e3, 100.);
  // _outside_mesh_3D->print_info();

  std::ostringstream name;

  name << "cut_cell_"
       << cut_cntr++
       << ".dat";
  _inside_mesh_3D->write  ("in_"  + name.str());
  _outside_mesh_3D->write ("out_" + name.str());

  // finally, add the elements to our lists.
  _inside_elem.clear();
  _outside_elem.clear();

  for (const auto & el : _inside_mesh_3D->element_ptr_range())
    if (el->volume() > std::numeric_limits<Real>::epsilon())
      _inside_elem.push_back (el);

  for (const auto & el : _outside_mesh_3D->element_ptr_range())
    if (el->volume() > std::numeric_limits<Real>::epsilon())
      _outside_elem.push_back (el);

#endif
}



} // namespace libMesh

#endif // LIBMESH_HAVE_TRIANGLE && LIBMESH_HAVE_TETGEN
