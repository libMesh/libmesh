/* The libMesh Finite Element Library. */
/* Copyright (C) 2013  Benjamin S. Kirk */

/* This library is free software; you can redistribute it and/or */
/* modify it under the terms of the GNU Lesser General Public */
/* License as published by the Free Software Foundation; either */
/* version 2.1 of the License, or (at your option) any later version. */

/* This library is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU */
/* Lesser General Public License for more details. */

/* You should have received a copy of the GNU Lesser General Public */
/* License along with this library; if not, write to the Free Software */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */



 // <h1>Subdomains Example 3 - Integrating discontinuous data that cuts the mesh</h1>
 //


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

// Bring in everything from the libMesh namespace
using namespace libMesh;

// declare the functions we will use
void integrate_function (const MeshBase &mesh);

// signed distance function
const Real radius = 0.5;
Real distance (const Point &p)
{
  Point cent(0.8, 0.9);
  return ((p-cent).size() - radius);
}

Real integrand (const Point &p)
{
  return (distance(p) < 0) ? 1. : 2.;
}



// Begin the main program.
int main (int argc, char** argv)
{
  // Initialize libMesh and any dependent libaries, like in example 2.
  LibMeshInit init (argc, argv);

  // This example requires Adaptive Mesh Refinement support - although
  // it only refines uniformly, the refinement code used is the same
  // underneath
#ifndef LIBMESH_ENABLE_AMR
  libmesh_example_assert(false, "--enable-amr");
#else

  // Skip this 2D example if libMesh was compiled as 1D-only.
  libmesh_example_assert(2 <= LIBMESH_DIM, "2D support");

  // Read the mesh from file.  This is the coarse mesh that will be used
  // in example 10 to demonstrate adaptive mesh refinement.  Here we will
  // simply read it in and uniformly refine it 5 times before we compute
  // with it.
  Mesh mesh;

  mesh.read ("mesh.xda");

  // Create a MeshRefinement object to handle refinement of our mesh.
  // This class handles all the details of mesh refinement and coarsening.
  MeshRefinement mesh_refinement (mesh);

  // Uniformly refine the mesh 5 times.  This is the
  // first time we use the mesh refinement capabilities
  // of the library.
  mesh_refinement.uniformly_refine (3);

  // Print information about the mesh to the screen.
  mesh.print_info();

  // integrate the desired function
  integrate_function (mesh);

  // All done.
  return 0;
#endif
}



void integrate_function (const MeshBase &mesh)
{
  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  std::vector<Real> vertex_distance;

  QComposite<QGauss> qrule (mesh.mesh_dimension(), THIRD);

  for (; el!=end_el; ++el)
    {
      const Elem *elem = *el;

      vertex_distance.clear();

      for (unsigned int v=0; v<elem->n_vertices(); v++)
	vertex_distance.push_back (distance(elem->point(v)));

      qrule.init (*elem, vertex_distance);
    }
}
