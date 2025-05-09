// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



// <h1>Miscellaneous Example 8 - Meshfree Interpolation Utilities.</h1>
// \author Benjamin S. Kirk
// \date 2012
//
// LibMesh provides some utilities for pointcloud-type interpolation, as
// demonstrated in this example.

// Example include files
#include "libmesh/libmesh.h"
#include "libmesh/meshfree_interpolation.h"
#include "libmesh/radial_basis_interpolation.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/equation_systems.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/tecplot_io.h"
#include "libmesh/threads.h"
#include "libmesh/node.h"
#include "libmesh/getpot.h"
#include "libmesh/meshfree_interpolation_function.h"

// C++ includes
#include <cstdlib>


// Bring in everything from the libMesh namespace
using namespace libMesh;


void create_random_point_cloud (const unsigned int Npts,
                                std::vector<Point> & pts,
                                const Real max_range = 10)
{
  libMesh::out << "Generating "<< Npts << " local point cloud...";
  pts.resize(Npts);

  for (size_t i=0;i<Npts;i++)
    {
      pts[i](0) = max_range * (std::rand() % 1000) / Real(1000);
      pts[i](1) = max_range * (std::rand() % 1000) / Real(1000);
      pts[i](2) = max_range * (std::rand() % 1000) / Real(1000);
    }
  libMesh::out << "done\n";
}



Real exact_solution_u (const Point & p)
{
  const Real
    x = p(0),
    y = p(1),
    z = p(2);

  return (x*x*x   +
          y*y*y*y +
          z*z*z*z*z);
}



Real exact_solution_v (const Point & p)
{
  const Real
    x = p(0),
    y = p(1),
    z = p(2);

  return (x*x +
          y*y +
          z*z*z);
}

Number exact_value (const Point & p,
                    const Parameters &,
                    const std::string &,
                    const std::string &)
{
  return exact_solution_v(p);
}

// We now define the function which provides the
// initialization routines for the "Convection-Diffusion"
// system.  This handles things like setting initial
// conditions and boundary conditions.
void init_sys(EquationSystems & es,
              const std::string & system_name)
{
  // Get a reference to the Convection-Diffusion system object.
  System & system =
    es.get_system<System>(system_name);

  system.project_solution(exact_value, nullptr, es.parameters);
}




int main(int argc, char ** argv)
{
  // Skip this example if we do not meet certain requirements
  libmesh_example_requires(3 <= LIBMESH_DIM, "3D support");
#ifndef LIBMESH_HAVE_EIGEN
  libmesh_example_requires(false, "--enable-eigen");
#endif
#ifndef LIBMESH_HAVE_ZLIB_H
  libmesh_example_requires(false, "--enable-zlib");
#endif
  // Nanoflann is giving me trouble with extended Reals
#ifdef LIBMESH_DEFAULT_TRIPLE_PRECISION
  libmesh_example_requires(false, "--disable-triple-precision");
#endif
#ifdef LIBMESH_DEFAULT_QUADRUPLE_PRECISION
  libmesh_example_requires(false, "--disable-quadruple-precision");
#endif

  // Initialize libMesh.
  LibMeshInit init (argc, argv);
  {
    GetPot input(argc, argv);

    // Demonstration case 1
    {
      std::vector<Point>       tgt_pts;
      std::vector<Number>      tgt_data_idi, tgt_data_rbi;
      std::vector<std::string> field_vars {"u", "v"};

      InverseDistanceInterpolation<3> idi (init.comm(),
                                           /* n_interp_pts = */ 8,
                                           /* power =        */ 2);

      RadialBasisInterpolation<3> rbi (init.comm());

      idi.set_field_variables (field_vars);
      rbi.set_field_variables (field_vars);

      const int n_source_points = input("n_source_points", 100);

      // Source points are independent on each processor
      const int my_n_source_points =
        n_source_points / init.comm().size() +
        (init.comm().rank() < n_source_points % init.comm().size());

      create_random_point_cloud (my_n_source_points,
                                 idi.get_source_points());


      // Explicitly set the data values we will interpolate from
      {
        const std::vector<Point> & src_pts  (idi.get_source_points());
        std::vector<Number>      & src_vals (idi.get_source_vals());

        src_vals.clear(); src_vals.reserve(2*src_pts.size());

        for (std::vector<Point>::const_iterator pt_it=src_pts.begin();
             pt_it != src_pts.end(); ++pt_it)
          {
            src_vals.push_back (exact_solution_u (*pt_it));
            src_vals.push_back (exact_solution_v (*pt_it));
          }
      }

      // give rbi the same info as idi
      rbi.get_source_points() = idi.get_source_points();
      rbi.get_source_vals()   = idi.get_source_vals();

      idi.prepare_for_use();
      rbi.prepare_for_use();

      libMesh::out << idi;

      // Interpolate to some other random points, and evaluate the result
      {
        const int n_target_points = input("n_target_points", 100);

        create_random_point_cloud (n_target_points,
                                   tgt_pts);

        //tgt_pts = rbi.get_source_points();

        idi.interpolate_field_data (field_vars,
                                    tgt_pts,
                                    tgt_data_idi);

        rbi.interpolate_field_data (field_vars,
                                    tgt_pts,
                                    tgt_data_rbi);

        std::vector<Number>::const_iterator
          v_idi = tgt_data_idi.begin(),
          v_rbi = tgt_data_rbi.begin();

        for (std::vector<Point>::const_iterator p_it=tgt_pts.begin();
             p_it!=tgt_pts.end(); ++p_it)
          {
            libMesh::out << "\nAt target point " << *p_it
                         << "\n u_interp_idi="   << *v_idi
                         << ", u_interp_rbi="    << *v_rbi
                         << ", u_exact="         << exact_solution_u(*p_it);
            ++v_idi;
            ++v_rbi;
            libMesh::out << "\n v_interp_idi=" << *v_idi
                         << ", v_interp_rbi="  << *v_rbi
                         << ", v_exact="       << exact_solution_v(*p_it)
                         << std::endl;
            ++v_idi;
            ++v_rbi;
          }
      }
    }


    // Demonstration case 2
    {
      Mesh mesh_a(init.comm()), mesh_b(init.comm());

      mesh_a.read("struct.ucd.gz");
      mesh_b.read("unstruct.ucd.gz");

      // Refine the meshes if requested
      const int n_refinements = input("n_refinements", 0);

      // Skip adaptive runs on a non-adaptive libMesh build
#ifndef LIBMESH_ENABLE_AMR
      libmesh_example_requires(n_refinements==0, "--enable-amr");
#else
      MeshRefinement mesh_refinement_a(mesh_a);
      mesh_refinement_a.uniformly_refine(n_refinements);
      MeshRefinement mesh_refinement_b(mesh_b);
      mesh_refinement_b.uniformly_refine(n_refinements);
#endif

      // Create equation systems objects.
      EquationSystems
        es_a(mesh_a), es_b(mesh_b);

      System & sys_a = es_a.add_system<System>("src_system");
      System & sys_b = es_b.add_system<System>("dest_system");

      sys_a.add_variable ("Cp", FIRST);
      sys_b.add_variable ("Cp", FIRST);

      sys_a.attach_init_function (init_sys);
      es_a.init();

      // Write out the initial conditions.
      TecplotIO(mesh_a).write_equation_systems ("src.dat",
                                                es_a);

      InverseDistanceInterpolation<3> idi (init.comm(),
                                           /* n_interp_pts = */ 4,
                                           /* power =        */ 2);
      RadialBasisInterpolation<3> rbi (init.comm());

      std::vector<Point> & src_pts  (idi.get_source_points());
      std::vector<Number> & src_vals (idi.get_source_vals());
      idi.set_field_variables({"Cp"});

      // We now will loop over every node in the source mesh
      // and add it to a source point list, along with the solution
      for (const auto & node : mesh_a.local_node_ptr_range())
        {
          src_pts.push_back(*node);
          src_vals.push_back(sys_a.current_solution(node->dof_number(0, 0, 0)));
        }

      rbi.set_field_variables({"Cp"});
      rbi.get_source_points() = idi.get_source_points();
      rbi.get_source_vals()   = idi.get_source_vals();

      // We have only set local values - prepare for use by gathering remote data
      idi.prepare_for_use();
      rbi.prepare_for_use();

      // Create a MeshfreeInterpolationFunction that uses our InverseDistanceInterpolation
      // object.  Since each MeshfreeInterpolationFunction shares the same InverseDistanceInterpolation
      // object in a threaded environment we must also provide a locking mechanism.
      {
        Threads::spin_mutex mutex;
        MeshfreeInterpolationFunction mif(idi, mutex);

        // project the solution onto system b
        es_b.init();
        sys_b.project_solution (&mif);

        // Write the result
        TecplotIO(mesh_b).write_equation_systems ("dest_idi.dat",
                                                  es_b);
      }

      // Create a MeshfreeInterpolationFunction that uses our RadialBasisInterpolation
      // object.  Since each MeshfreeInterpolationFunction shares the same RadialBasisInterpolation
      // object in a threaded environment we must also provide a locking mechanism.
      {
        Threads::spin_mutex mutex;
        MeshfreeInterpolationFunction mif(rbi, mutex);

        // project the solution onto system b
        sys_b.project_solution (&mif);

        // Write the result
        TecplotIO(mesh_b).write_equation_systems ("dest_rbi.dat",
                                                  es_b);
      }
    }
  }
  return 0;
}
