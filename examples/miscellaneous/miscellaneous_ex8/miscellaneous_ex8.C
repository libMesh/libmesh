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



// <h1>Miscellaneous Example 7 - Can use the PetscDMNonlinearSolver (available in PETSc-3.3.0 or above) to solve a VI version of the problem.</h1>
//
// LibMesh interfaces directly with PETSc's variational inequality solver,
// this example shows how to do it.

// Example include files
#include "libmesh/libmesh.h"
#include "libmesh/meshfree_interpolation.h"
#include "libmesh/mesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/tecplot_io.h"
#include "meshless_interpolation_function.h"

// C++ includes
#include <cstdlib>


// Bring in everything from the libMesh namespace
using namespace libMesh;


void create_random_point_cloud (const unsigned int Npts,
				std::vector<Point> &pts,
				const Real max_range = 10)
{
  std::cout << "Generating "<< Npts << " point cloud...";
  pts.resize(Npts);

  for (size_t i=0;i<Npts;i++)
    {
      pts[i](0) = max_range * (std::rand() % 1000) / Real(1000);
      pts[i](1) = max_range * (std::rand() % 1000) / Real(1000);
      pts[i](2) = max_range * (std::rand() % 1000) / Real(1000);
    }
  std::cout << "done\n";
}



Real exact_solution_u (const Point &p)
{
  const Real
    x = p(0),
    y = p(1),
    z = p(2);
  
  return (x*x*x   +
	  y*y*y*y +
	  z*z*z*z*z);
}



Real exact_solution_v (const Point &p)
{
  const Real
    x = p(0),
    y = p(1),
    z = p(2);
  
  return (x*x   +
	  y*y +
	  z*z*z);
}

Number exact_value (const Point& p,
                    const Parameters&,
                    const std::string&,
                    const std::string&)
{
  return exact_solution_v(p);
}

// We now define the function which provides the
// initialization routines for the "Convection-Diffusion"
// system.  This handles things like setting initial
// conditions and boundary conditions.
void init_sys(EquationSystems& es,
              const std::string& system_name)
{
  // Get a reference to the Convection-Diffusion system object.
  System & system =
    es.get_system<System>(system_name);
  
  system.project_solution(exact_value, NULL, es.parameters);
}




int main(int argc, char** argv)
{
  // Initialize libMesh.
  LibMeshInit init (argc, argv);
  {
    // Demonstration case 1
    {
      std::vector<Point>       tgt_pts;
      std::vector<Number>      tgt_data;
      std::vector<std::string> field_vars;
      
      field_vars.push_back("u");
      field_vars.push_back("v");
      
      InverseDistanceInterpolation<3> idi (/* n_interp_pts = */ 8,
					   /* power =        */ 2);
      
      idi.set_field_variables (field_vars);
      
      create_random_point_cloud (1e5,
				 idi.get_source_points());
      
      // Explicitly set the data values we will interpolate from
      {
	const std::vector<Point> &src_pts  (idi.get_source_points());
	std::vector<Real>        &src_vals (idi.get_source_vals());
	
	src_vals.clear(); src_vals.reserve(src_pts.size());
	
	for (std::vector<Point>::const_iterator pt_it=src_pts.begin();
	     pt_it != src_pts.end(); ++pt_it)
	  {
	    src_vals.push_back (exact_solution_u (*pt_it));
	    src_vals.push_back (exact_solution_v (*pt_it));
	  }	  
      }

      std::cout << idi;
      
      // Interpolate to some other random points, and evaluate the result
      {
	create_random_point_cloud (10,
				   tgt_pts);
	
	idi.interpolate_field_data (field_vars,
				    tgt_pts,
				    tgt_data);
	
      
	std::vector<Number>::const_iterator v_it=tgt_data.begin();
	
	for (std::vector<Point>::const_iterator  p_it=tgt_pts.begin();
	     p_it!=tgt_pts.end(); ++p_it)
	  {
	    std::cout << "\nAt target point " << *p_it
		      << "\n u_interp=" << *v_it
		      << ", u_exact="  << exact_solution_u(*p_it);
	    ++v_it;
	    std::cout << "\n v_interp=" << *v_it
		      << ", v_exact="  << exact_solution_v(*p_it)
		      << std::endl;
	    ++v_it;
	  }
      }
    }


    // Demonstration case 2
    {
      Mesh mesh_a, mesh_b;

      mesh_a.read("struct.ucd.gz"); mesh_b.read("unstruct.ucd.gz");

      // Create equation systems objects.
      EquationSystems
	es_a(mesh_a), es_b(mesh_b);

      System
	&sys_a = es_a.add_system<System>("src_system"),
	&sys_b = es_b.add_system<System>("dest_system");

      sys_a.add_variable ("Cp", FIRST);
      sys_b.add_variable ("Cp", FIRST);

      sys_a.attach_init_function (init_sys);
      es_a.init();
      
      // Write out the initial conditions.
      TecplotIO(mesh_a).write_equation_systems ("src.dat",
						es_a);

      InverseDistanceInterpolation<3> idi (/* n_interp_pts = */ 4,
					   /* power =        */ 2);

      std::vector<Point>  &src_pts  (idi.get_source_points());
      std::vector<Number> &src_vals (idi.get_source_vals());
      std::vector<std::string> field_vars;      
      field_vars.push_back("Cp");
      idi.set_field_variables(field_vars);

      // We now will loop over every node in the source mesh
      // and add it to a source point list, along with the solution
      {
	MeshBase::const_node_iterator nd  = mesh_a.local_nodes_begin();
	MeshBase::const_node_iterator end = mesh_a.local_nodes_end();

	for (; nd!=end; ++nd)
	  {
	    const Node *node(*nd);
	    src_pts.push_back(*node);
	    src_vals.push_back(sys_a.current_solution(node->dof_number(0,0,0)));
	  }	  	
      }

      // We have only set local values - prepare for use by gathering remote gata
      idi.prepare_for_use();

      // Create a MeshlessInterpolationFunction that uses our InverseDistanceInterpolation
      // object
      MeshlessInterpolationFunction mif(idi);

      // project the solution onto system b
      es_b.init();
      sys_b.project_solution (&mif);
      
      // Write the result
      TecplotIO(mesh_b).write_equation_systems ("dest.dat",
						es_b);
    }


    
  }
  return 0;
}
