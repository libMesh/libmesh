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



// <h1>Subdomains Example 3 - Integrating discontinuous data that cuts the mesh</h1>
// \author Benjamin S. Kirk
// \date 2013


// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>

// Basic include file needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/quadrature_composite.h"
#include "libmesh/fe.h"
#include "libmesh/elem.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// declare the functions we will use
void integrate_function (const MeshBase & mesh);

// signed distance function
const Real radius = 0.5;

Real distance (const Point & p)
{
  Point cent(0.8, 0.9);
  return ((p-cent).norm() - radius);
}

Real integrand (const Point & p)
{
  return (distance(p) < 0) ? 10. : 1.;
}



// Begin the main program.
int main (int argc, char ** argv)
{
  // Initialize libMesh and any dependent libaries, like in example 2.
  LibMeshInit init (argc, argv);

  // This example requires Adaptive Mesh Refinement support - although
  // it only refines uniformly, the refinement code used is the same
  // underneath.  It also requires libmesh support for Triangle and
  // Tetgen, which means that libmesh must be configured with
  // --disable-strict-lgpl for this example to run.
#if !defined(LIBMESH_HAVE_TRIANGLE) || !defined(LIBMESH_HAVE_TETGEN) || !defined(LIBMESH_ENABLE_AMR)
  libmesh_example_requires(false, "--disable-strict-lgpl --enable-amr");
#else

  // Skip this 2D example if libMesh was compiled as 1D-only.
  libmesh_example_requires(2 <= LIBMESH_DIM, "2D support");

  // Read the mesh from file.  This is the coarse mesh that will be used
  // in example 10 to demonstrate adaptive mesh refinement.  Here we will
  // simply read it in and uniformly refine it 5 times before we compute
  // with it.
  Mesh mesh(init.comm());

  {
    unsigned int dim=2;

    if (argc == 3 && std::atoi(argv[2]) == 3)
      {
        libmesh_here();
        dim=3;
      }

    mesh.read ((dim==2) ? "mesh.xda" : "hybrid_3d.xda");
  }

  // Create a MeshRefinement object to handle refinement of our mesh.
  // This class handles all the details of mesh refinement and coarsening.
  MeshRefinement mesh_refinement (mesh);

  // Uniformly refine the mesh 4 times.  This is the
  // first time we use the mesh refinement capabilities
  // of the library.
  mesh_refinement.uniformly_refine (4);

  // Print information about the mesh to the screen.
  mesh.print_info();

  // integrate the desired function
  integrate_function (mesh);

  // All done.
  return 0;
#endif
}



void integrate_function (const MeshBase & mesh)
{
#if defined(LIBMESH_HAVE_TRIANGLE) && defined(LIBMESH_HAVE_TETGEN)
  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  std::vector<Real> vertex_distance;

  QComposite<QGauss> qrule (mesh.mesh_dimension(), FIRST);
  //QGauss qrule (mesh.mesh_dimension(), FIRST);

  UniquePtr<FEBase> fe (FEBase::build (mesh.mesh_dimension(), FEType (FIRST, LAGRANGE)));

  Real int_val=0.;

  const std::vector<Point> & q_points = fe->get_xyz();
  const std::vector<Real>  & JxW      = fe->get_JxW();

  for (; el!=end_el; ++el)
    {
      const Elem * elem = *el;

      vertex_distance.clear();

      for (unsigned int v=0; v<elem->n_vertices(); v++)
        vertex_distance.push_back (distance(elem->point(v)));

      qrule.init (*elem, vertex_distance);

      fe->reinit (elem,
                  &(qrule.get_points()),
                  &(qrule.get_weights()));


      // TODO:  would it be valuable to have the composite quadrature rule sort
      // from smallest to largest JxW value to help prevent
      // ... large + small + large + large + small ...
      // type truncation errors?
      for (std::size_t qp=0; qp<q_points.size(); qp++)
        int_val += JxW[qp] * integrand(q_points[qp]);
    }

  mesh.comm().sum (int_val);

  libMesh::out << "\n***********************************\n"
               << " int_val   = " << int_val << std::endl
               << " exact_val = " <<  1*(2*2 - radius*radius*pi) + 10.*(radius*radius*pi)
               << "\n***********************************\n"
               << std::endl;
#else
  libmesh_ignore(mesh);
#endif
}
