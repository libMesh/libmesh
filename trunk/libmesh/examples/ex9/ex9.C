/* $Id$ */

/* The Next Great Finite Element Library. */
/* Copyright (C) 2003  Benjamin S. Kirk */

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



 // <h1>Example 9 - Solving a Transient Linear System in Parallel</h1>
 //
 // This example shows how a simple, linear transient
 // system can be solved in parallel.  The system is simple
 // scalar convection-diffusion with a specified external
 // velocity.  The initial condition is given, and the
 // solution is advanced in time with a standard Crank-Nicolson
 // time-stepping strategy.

// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>

// Basic include file needed for the mesh functionality.
#include "libmesh.h"
#include "mesh.h"
#include "serial_mesh.h"
#include "mesh_refinement.h"
#include "gmv_io.h"
#include "tecplot_io.h"
#include "equation_systems.h"
#include "fe.h"
#include "quadrature_gauss.h"
#include "dof_map.h"
#include "sparse_matrix.h"
#include "numeric_vector.h"
#include "dense_matrix.h"
#include "dense_vector.h"
#include "xdr_io.h"
#include "parallel.h"

// Some (older) compilers do not offer full stream 
// functionality, OStringStream works around this.
#include "o_string_stream.h"

// This example will solve a linear transient system,
// so we need to include the \p TransientLinearImplicitSystem definition.
#include "linear_implicit_system.h"
#include "transient_system.h"
#include "vector_value.h"

// The definition of a geometric element
#include "elem.h"
#include "gmv_io.h"


// We can now begin the main program.  Note that this
// example will fail if you are using complex numbers
// since it was designed to be run only with real numbers.
int main (int argc, char** argv)
{
  // Initialize libMesh.
  libMesh::init (argc, argv);

#ifndef ENABLE_AMR
  if (libMesh::processor_id() == 0)
    std::cerr << "ERROR: This example requires libMesh to be\n"
              << "compiled with AMR support!"
              << std::endl;
  return 0;
#else

  {    
    // Create a two-dimensional mesh.
    Mesh mesh (3);
        
    // Read the mesh from file.  This is the coarse mesh that will be used
    // in example 10 to demonstrate adaptive mesh refinement.  Here we will
    // simply read it in and uniformly refine it 5 times before we compute
    // with it.
    //mesh.read ("../ex10/mesh.xda");
    //mesh.read ("../../reference_elements/3D/one_hex.xda");
    mesh.read ("/pr2/benkirk/fins/apps/sphere/grid/3D/sphere.111.xdr");
    
    // Create a MeshRefinement object to handle refinement of our mesh.
    // This class handles all the details of mesh refinement and coarsening.
    //MeshRefinement mesh_refinement (mesh);
    
    // Uniformly refine the mesh 5 times.  This is the
    // first time we use the mesh refinement capabilities
    // of the library.
    //mesh_refinement.uniformly_refine (2);
    
    // Print information about the mesh to the screen.
    //mesh.print_info();

    XdrIO io(mesh);
    io.legacy() = false;
    io.binary() = false;
    //io.partition_map_file_name() = ".";
    //io.subdomain_map_file_name() = ".";
    //io.polynomial_level_file_name() = ".";
    io.write("foo.xda");
    mesh.clear();
    Parallel::barrier();
    
    io.read("foo.xda");
    mesh.print_info();
    mesh.prepare_for_use();
    //io.write("foo.xda");

    GMVIO(mesh).write("foo.gmv"); 
    TecplotIO(mesh).write("foo.plt"); 
    
  }
#endif // #ifdef ENABLE_AMR

  // All done.  
  return libMesh::close ();
}
