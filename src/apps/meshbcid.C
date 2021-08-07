// The libMesh Finite Element Library.
// Copyright (C) 2002-2021 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// Open the mesh named in command line arguments,
// apply the named tests on every boundary side and normal vector,
// and give the specified boundary condition id to each side
// that passes the tests.

#include <limits>
#include <string>

#include "libmesh/libmesh.h"

#include "libmesh/boundary_info.h"
#include "libmesh/bounding_box.h"
#include "libmesh/elem.h"
#include "libmesh/fe.h"
#include "libmesh/getpot.h"
#include "libmesh/mesh.h"
#include "libmesh/point.h"
#include "libmesh/quadrature_gauss.h"

using namespace libMesh;

void usage_error(const char * progname)
{
  libMesh::out << "Usage: " << progname
               << " --dim d --input inputmesh --output outputmesh --newbcid idnum --tests --moretests"
               << std::endl;

  exit(1);
}

int main(int argc, char ** argv)
{
  LibMeshInit init(argc, argv);

  GetPot cl(argc, argv);

  unsigned char dim = -1;
  if (!cl.search("--dim"))
    {
      libMesh::err << "No --dim argument found!" << std::endl;
      usage_error(argv[0]);
    }
  dim = cl.next(dim);

  Mesh mesh(init.comm(), dim);

  if (!cl.search("--input"))
    {
      libMesh::err << "No --input argument found!" << std::endl;
      usage_error(argv[0]);
    }
  const char * meshname = cl.next("mesh.xda");

  mesh.read(meshname);
  libMesh::out << "Loaded mesh " << meshname << std::endl;

  if (!cl.search("--newbcid"))
    {
      libMesh::err << "No --bcid argument found!" << std::endl;
      usage_error(argv[0]);
    }
  boundary_id_type bcid = 0;
  bcid = cl.next(bcid);

  Point minpt(-std::numeric_limits<Real>::max());
#if LIBMESH_DIM > 1
  minpt(1) = -std::numeric_limits<Real>::max();
#endif
#if LIBMESH_DIM > 2
  minpt(2) = -std::numeric_limits<Real>::max();
#endif
  Point maxpt = -minpt;

  BoundingBox normals(minpt, maxpt),
              points(minpt, maxpt);

  if (cl.search("--minnormalx"))
    normals.min()(0) = cl.next(normals.min()(0));
  if (cl.search("--maxnormalx"))
    normals.max()(0) = cl.next(normals.max()(0));

  if (cl.search("--minpointx"))
    points.min()(0) = cl.next(points.min()(0));
  if (cl.search("--maxpointx"))
    points.max()(0) = cl.next(points.max()(0));

#if LIBMESH_DIM > 1
  if (cl.search("--minnormaly"))
    normals.min()(1) = cl.next(normals.min()(1));
  if (cl.search("--maxnormaly"))
    normals.max()(1) = cl.next(normals.max()(1));

  if (cl.search("--minpointy"))
    points.min()(1) = cl.next(points.min()(1));
  if (cl.search("--maxpointy"))
    points.max()(1) = cl.next(points.max()(1));
#endif

#if LIBMESH_DIM > 2
  if (cl.search("--minnormalz"))
    normals.min()(2) = cl.next(normals.min()(2));
  if (cl.search("--maxnormalz"))
    normals.max()(2) = cl.next(normals.max()(2));

  if (cl.search("--minpointz"))
    points.min()(2) = cl.next(points.min()(2));
  if (cl.search("--maxpointz"))
    points.max()(2) = cl.next(points.max()(2));
#endif

  libMesh::out << "min point = " << points.min() << std::endl;
  libMesh::out << "max point = " << points.max() << std::endl;
  libMesh::out << "min normal = " << normals.min() << std::endl;
  libMesh::out << "max normal = " << normals.max() << std::endl;

  bool matcholdbcid = false;
  boundary_id_type oldbcid = 0;
  if (cl.search("--oldbcid"))
    {
      matcholdbcid = true;
      oldbcid = cl.next(oldbcid);
      if (oldbcid < 0)
        oldbcid = BoundaryInfo::invalid_id;
    }

  std::unique_ptr<FEBase> fe = FEBase::build(dim, FEType(FIRST,LAGRANGE));
  QGauss qface(dim-1, CONSTANT);
  fe->attach_quadrature_rule(&qface);
  const std::vector<Point> & face_points = fe->get_xyz();
  const std::vector<Point> & face_normals = fe->get_normals();

  for (auto & elem : mesh.element_ptr_range())
    {
      unsigned int n_sides = elem->n_sides();

      // Container to catch ids handed back from BoundaryInfo
      std::vector<boundary_id_type> ids;

      for (unsigned short s=0; s != n_sides; ++s)
        {
          if (elem->neighbor_ptr(s))
            continue;

          fe->reinit(elem,s);
          const Point & p = face_points[0];
          const Point & n = face_normals[0];

          //libMesh::out << "elem = " << elem->id() << std::endl;
          //libMesh::out << "vertex average = " << elem->vertex_average() << std::endl;
          //libMesh::out << "p = " << p << std::endl;
          //libMesh::out << "n = " << n << std::endl;

          if (points.contains_point(p) &&
              normals.contains_point(n))
            {
              // Get the list of boundary ids for this side
              mesh.get_boundary_info().boundary_ids(elem, s, ids);

              // There should be at most one value present, otherwise the
              // logic here won't work.
              libmesh_assert(ids.size() <= 1);

              // A convenient name for the side's ID.
              boundary_id_type b_id = ids.empty() ? BoundaryInfo::invalid_id : ids[0];

              if (matcholdbcid && b_id != oldbcid)
                continue;

              mesh.get_boundary_info().remove_side(elem, s);
              mesh.get_boundary_info().add_side(elem, s, bcid);
              //libMesh::out << "Set element " << elem->id() << " side " << s <<
              //                " to boundary " << bcid << std::endl;
            }
        }
    }

  // We might have removed *every* instance of a given id, and if that
  // happened then we should make sure that file formats which write
  // out id sets do not write out the removed id.
  mesh.get_boundary_info().regenerate_id_sets();

  std::string outputname;
  if (cl.search("--output"))
    {
      outputname = cl.next("mesh.xda");
    }
  else
    {
      outputname = "new.";
      outputname += meshname;
    }


  mesh.write(outputname.c_str());
  libMesh::out << "Wrote mesh " << outputname << std::endl;

  return 0;
}
