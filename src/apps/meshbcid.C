// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

  if(!cl.search("--input"))
    {
      libMesh::err << "No --input argument found!" << std::endl;
      usage_error(argv[0]);
    }
  const char * meshname = cl.next("mesh.xda");

  mesh.read(meshname);
  libMesh::out << "Loaded mesh " << meshname << std::endl;

  if(!cl.search("--newbcid"))
    {
      libMesh::err << "No --bcid argument found!" << std::endl;
      usage_error(argv[0]);
    }
  boundary_id_type bcid = 0;
  bcid = cl.next(bcid);

  Point minnormal(-std::numeric_limits<Real>::max(),
                  -std::numeric_limits<Real>::max(),
                  -std::numeric_limits<Real>::max());
  Point maxnormal(std::numeric_limits<Real>::max(),
                  std::numeric_limits<Real>::max(),
                  std::numeric_limits<Real>::max());
  Point minpoint(-std::numeric_limits<Real>::max(),
                 -std::numeric_limits<Real>::max(),
                 -std::numeric_limits<Real>::max());
  Point maxpoint(std::numeric_limits<Real>::max(),
                 std::numeric_limits<Real>::max(),
                 std::numeric_limits<Real>::max());

  if (cl.search("--minnormalx"))
    minnormal(0) = cl.next(minnormal(0));
  if (cl.search("--minnormalx"))
    minnormal(0) = cl.next(minnormal(0));
  if (cl.search("--maxnormalx"))
    maxnormal(0) = cl.next(maxnormal(0));
  if (cl.search("--minnormaly"))
    minnormal(1) = cl.next(minnormal(1));
  if (cl.search("--maxnormaly"))
    maxnormal(1) = cl.next(maxnormal(1));
  if (cl.search("--minnormalz"))
    minnormal(2) = cl.next(minnormal(2));
  if (cl.search("--maxnormalz"))
    maxnormal(2) = cl.next(maxnormal(2));

  if (cl.search("--minpointx"))
    minpoint(0) = cl.next(minpoint(0));
  if (cl.search("--maxpointx"))
    maxpoint(0) = cl.next(maxpoint(0));
  if (cl.search("--minpointy"))
    minpoint(1) = cl.next(minpoint(1));
  if (cl.search("--maxpointy"))
    maxpoint(1) = cl.next(maxpoint(1));
  if (cl.search("--minpointz"))
    minpoint(2) = cl.next(minpoint(2));
  if (cl.search("--maxpointz"))
    maxpoint(2) = cl.next(maxpoint(2));

  libMesh::out << "min point = " << minpoint << std::endl;
  libMesh::out << "max point = " << maxpoint << std::endl;
  libMesh::out << "min normal = " << minnormal << std::endl;
  libMesh::out << "max normal = " << maxnormal << std::endl;

  bool matcholdbcid = false;
  boundary_id_type oldbcid = 0;
  if (cl.search("--oldbcid"))
    {
      matcholdbcid = true;
      oldbcid = cl.next(oldbcid);
      if (oldbcid < 0)
        oldbcid = BoundaryInfo::invalid_id;
    }

  UniquePtr<FEBase> fe = FEBase::build(dim, FEType(FIRST,LAGRANGE));
  QGauss qface(dim-1, CONSTANT);
  fe->attach_quadrature_rule(&qface);
  const std::vector<Point> & face_points = fe->get_xyz();
  const std::vector<Point> & face_normals = fe->get_normals();

  MeshBase::element_iterator           el = mesh.elements_begin();
  const MeshBase::element_iterator end_el = mesh.elements_end();
  for (; el != end_el; ++el)
    {
      Elem * elem = *el;
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
          //libMesh::out << "centroid = " << elem->centroid() << std::endl;
          //libMesh::out << "p = " << p << std::endl;
          //libMesh::out << "n = " << n << std::endl;

          if (p(0) > minpoint(0) && p(0) < maxpoint(0) &&
              p(1) > minpoint(1) && p(1) < maxpoint(1) &&
              p(2) > minpoint(2) && p(2) < maxpoint(2) &&
              n(0) > minnormal(0) && n(0) < maxnormal(0) &&
              n(1) > minnormal(1) && n(1) < maxnormal(1) &&
              n(2) > minnormal(2) && n(2) < maxnormal(2))
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
  if(cl.search("--output"))
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
