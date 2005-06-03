/* $Id: ex15.C,v 1.4 2005-06-03 15:49:57 jwpeterson Exp $ */

/* The Next Great Finite Element Library. */
/* Copyright (C) 2004  Benjamin S. Kirk, John W. Peterson */

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

 // <h1>Example 15 - Biharmonic Equation</h1>
 //
 // This example solves the Biharmonic equation on a square domain
 // using a Galerkin formulation with C1 elements approximating the
 // H^2_0 function space.
 // The initial mesh contains two TRI6 elements.
 // The mesh is provided in the standard libMesh ASCII format file
 // named "domain.xda".  In addition, an input file named "ex15.in"
 // is provided which allows the user to set several parameters for
 // the solution so that the problem can be re-run without a
 // re-compile.  The solution technique employed is to have a
 // refinement loop with a linear solve inside followed by a
 // refinement of the grid and projection of the solution to the new grid
 // In the final loop iteration, there is no additional
 // refinement after the solve.  In the input file "ex15.in", the variable
 // "max_r_steps" controls the number of refinement steps, and
 // "max_r_level" controls the maximum element refinement level.

// LibMesh include files.
#include "mesh.h"
#include "equation_systems.h"
#include "linear_implicit_system.h"
#include "gmv_io.h"
#include "tecplot_io.h"
#include "fe.h"
#include "quadrature_clough.h"
#include "dense_matrix.h"
#include "dense_vector.h"
#include "sparse_matrix.h"
#include "mesh_refinement.h"
#include "error_vector.h"
#include "kelly_error_estimator.h"
#include "getpot.h"
#include "exact_solution.h"
#include "dof_map.h"
#include "numeric_vector.h"

// Function prototype.  This is the function that will assemble
// the linear system for our Biharmonic problem.  Note that the
// function will take the \p EquationSystems object and the
// name of the system we are assembling as input.  From the
// \p EquationSystems object we have acess to the \p Mesh and
// other objects we might need.
void assemble_biharmonic(EquationSystems& es,
                      const std::string& system_name);


// Prototype for calculation of the exact solution.  Useful
// for setting boundary conditions.
Number exact_solution(const Point& p,
		      const Parameters&,   // parameters, not needed
		      const std::string&,  // sys_name, not needed
		      const std::string&); // unk_name, not needed);

// Prototype for calculation of the gradient of the exact solution.  
// Necessary for setting boundary conditions in H^2_0 and testing
// H^1 convergence of the solution
Gradient exact_derivative(const Point& p,
			  const Parameters&,   // parameters, not needed
			  const std::string&,  // sys_name, not needed
			  const std::string&); // unk_name, not needed);

Number forcing_function(const Point& p);



int main(int argc, char** argv)
{
  // Initialize libMesh.
  libMesh::init (argc, argv);

#ifndef ENABLE_SECOND_DERIVATIVES

  std::cerr << "ERROR: This example requires the library to be "
	    << "compiled with second derivatives support!"
	    << std::endl;
  here();

  return 0;

#else

  {
    // Set the dimensionality of the mesh = 2
    const unsigned int dim = 2;

    // Create a two-dimensional mesh.
    Mesh mesh (dim);
    
    // Parse the input file
    GetPot input_file("ex15.in");

    // Read in parameters from the input file
    const unsigned int max_r_steps = input_file("max_r_steps", 3);

    std::cerr.setf(std::ios::scientific);
    std::cerr.precision(3);
    
    // Output file for plotting the error 
    std::string output_file = "clough_uniform.m";
    
    std::ofstream out (output_file.c_str());
    out << "% dofs     L2-error     H1-error" << std::endl;
    out << "e = [" << std::endl;
    
    // Read in the mesh
    mesh.read("domain.xda");

    // Mesh Refinement object
    MeshRefinement mesh_refinement(mesh);
    
    // Create an equation systems object.
    EquationSystems equation_systems (mesh);

    // Declare the system and its variables.
    {
      // Creates a system named "Biharmonic"
      equation_systems.add_system<LinearImplicitSystem> ("Biharmonic");

      // Adds the variable "u" to "Biharmonic".  "u"
      // will be approximated using Clough-Tocher cubic C1 triangles
      equation_systems("Biharmonic").add_variable("u", THIRD, CLOUGH);
      
      // Give the system a pointer to the matrix assembly
      // function.
      equation_systems("Biharmonic").attach_assemble_function
		      (assemble_biharmonic);
      
      // Initialize the data structures for the equation system.
      equation_systems.init();

      // Set linear solver max iterations
      equation_systems.parameters.set<unsigned int>
		      ("linear solver maximum iterations") = 1000;

      // Linear solver tolerance.
      equation_systems.parameters.set<Real>
		      ("linear solver tolerance") = 1.e-13;
      
      // Prints information about the system to the screen.
      equation_systems.print_info();
    }

    // Construct ExactSolution object and attach function to compute exact solution
    ExactSolution exact_sol(equation_systems);
    exact_sol.attach_exact_value(exact_solution);
    exact_sol.attach_exact_deriv(exact_derivative);

    // Convenient reference to the system
    LinearImplicitSystem& system = 
      equation_systems.get_system<LinearImplicitSystem>("Biharmonic");

    // A refinement loop.
    for (unsigned int r_step=0; r_step<max_r_steps; r_step++)
      {
	std::cout << "Beginning Solve " << r_step << std::endl;
	
	// Solve the system "Biharmonic", just like example 2.
	equation_systems("Biharmonic").solve();

	std::cout << "System has: " << equation_systems.n_active_dofs()
		  << " degrees of freedom."
		  << std::endl;
	
	std::cout << "Linear solver converged at step: "
		  << system.n_linear_iterations()
		  << ", final residual: "
		  << system.final_linear_residual()
		  << std::endl;
	
	// Compute the error.
	exact_sol.compute_error("Biharmonic", "u");

	// Print out the error values
	std::cout << "L2-Error is: "
		  << exact_sol.l2_error("Biharmonic", "u")
		  << std::endl;
	std::cout << "H1-Error is: "
		  << exact_sol.h1_error("Biharmonic", "u")
		  << std::endl;

	// Print to output file
	out << equation_systems.n_active_dofs() << " "
	    << exact_sol.l2_error("Biharmonic", "u") << " "
	    << exact_sol.h1_error("Biharmonic", "u") << std::endl;

	// Possibly refine the mesh - Clough Tocher elements currently
	// do not support hanging nodes, so we use uniform refinement
	if (r_step+1 != max_r_steps)
	  {
	    std::cout << "  Refining the mesh..." << std::endl;

            mesh_refinement.uniformly_refine(1);
	    
	    // This call reinitializes the \p EquationSystems object for
	    // the newly refined mesh.  One of the steps in the
	    // reinitialization is projecting the \p solution,
	    // \p old_solution, etc... vectors from the old mesh to
	    // the current one.
	    equation_systems.reinit ();
	  }
      }	    
    
    

    
    // Write out the solution
    // After solving the system write the solution
    // to a GMV-formatted plot file.
    GMVIO (mesh).write_equation_systems ("domain.gmv",
    					 equation_systems);

    // Close up the output file.
    out << "];" << std::endl;
    out << "hold on" << std::endl;
    out << "plot(e(:,1), e(:,2), 'bo-');" << std::endl;
    out << "plot(e(:,1), e(:,3), 'ro-');" << std::endl;
    //    out << "set(gca,'XScale', 'Log');" << std::endl;
    //    out << "set(gca,'YScale', 'Log');" << std::endl;
    out << "xlabel('dofs');" << std::endl;
    out  << "title('Clough-Tocher elements');" << std::endl;
    out << "legend('L2-error', 'H1-error');" << std::endl;
    //     out << "disp('L2-error linear fit');" << std::endl;
    //     out << "polyfit(log10(e(:,1)), log10(e(:,2)), 1)" << std::endl;
    //     out << "disp('H1-error linear fit');" << std::endl;
    //     out << "polyfit(log10(e(:,1)), log10(e(:,3)), 1)" << std::endl;
  }

  
  // All done.  
  return libMesh::close ();
#endif
}




// We now define the exact solution
Number exact_solution(const Point& p,
		      const Parameters&,  // parameters, not needed
		      const std::string&, // sys_name, not needed
		      const std::string&) // unk_name, not needed
{
  const Real x = p(0);
  const Real y = p(1);
  
  return sin(x*y);
}


Number forcing_function(const Point& p)
{
  const Real x = p(0);
  const Real y = p(1);

  return (x*x + y*y) * (x*x + y*y) * sin(x*y) -
		  8*x*y*cos(x*y) - 4*sin(x*y);
}




// We now define the gradient of the exact solution
Gradient exact_derivative(const Point& p,
			  const Parameters&,  // parameters, not needed
			  const std::string&, // sys_name, not needed
			  const std::string&) // unk_name, not needed
{
  // Gradient value to be returned.
  Gradient gradu;
  
  // x and y coordinates in space
  const Real x = p(0);
  const Real y = p(1);

  gradu(0) = y * cos(x*y);
  gradu(1) = x * cos(x*y);

  return gradu;
}






// We now define the matrix assembly function for the
// Biharmonic system.  We need to first compute element
// matrices and right-hand sides, and then take into
// account the boundary conditions, which will be handled
// via a penalty method.
void assemble_biharmonic(EquationSystems& es,
                      const std::string& system_name)
{
  // It is a good idea to make sure we are assembling
  // the proper system.
  assert (system_name == "Biharmonic");

#ifdef ENABLE_SECOND_DERIVATIVES

  // Declare a performance log.  Give it a descriptive
  // string to identify what part of the code we are
  // logging, since there may be many PerfLogs in an
  // application.
  PerfLog perf_log ("Matrix Assembly",false);
  
    // Get a constant reference to the mesh object.
  const Mesh& mesh = es.get_mesh();

  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to the LinearImplicitSystem we are solving
  LinearImplicitSystem& system = es.get_system<LinearImplicitSystem>("Biharmonic");
  
  // A reference to the \p DofMap object for this system.  The \p DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.  We will talk more about the \p DofMap
  // in future examples.
  const DofMap& dof_map = system.get_dof_map();

  // Get a constant reference to the Finite Element type
  // for the first (and only) variable in the system.
  FEType fe_type = dof_map.variable_type(0);

  // Build a Finite Element object of the specified type.  Since the
  // \p FEBase::build() member dynamically creates memory we will
  // store the object as an \p AutoPtr<FEBase>.  This can be thought
  // of as a pointer that will clean up after itself.
  AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));
  
  // A 7th order Clough quadrature rule for numerical integration.
  // With 2D triangles, the Clough quadrature rule puts a Gaussian
  // quadrature rule on each of the 3 subelements
  QClough qrule (dim, SEVENTH);

  // Tell the finite element object to use our quadrature rule.
  fe->attach_quadrature_rule (&qrule);

  // Declare a special finite element object for
  // boundary integration.
  AutoPtr<FEBase> fe_face (FEBase::build(dim, fe_type));
	      
  // Boundary integration requires another quadraure rule,
  // with dimensionality one less than the dimensionality
  // of the element.
  // In 1D, the Clough and Gauss quadrature rules are identical.
  QClough qface(dim-1, SIXTH);
  
  // Tell the finte element object to use our
  // quadrature rule.
  fe_face->attach_quadrature_rule (&qface);

  // Here we define some references to cell-specific data that
  // will be used to assemble the linear system.
  // We begin with the element Jacobian * quadrature weight at each
  // integration point.   
  const std::vector<Real>& JxW = fe->get_JxW();

  // The physical XY locations of the quadrature points on the element.
  // These might be useful for evaluating spatially varying material
  // properties at the quadrature points.
  const std::vector<Point>& q_point = fe->get_xyz();

  // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector<Real> >& phi = fe->get_phi();

  // The element shape function second derivatives evaluated at the
  // quadrature points.  Note that for the simple biharmonic, shape
  // function first derivatives are unnecessary.
  const std::vector<std::vector<RealTensor> >& d2phi = fe->get_d2phi();

  // Define data structures to contain the element matrix
  // and right-hand-side vector contribution.  Following
  // basic finite element terminology we will denote these
  // "Ke" and "Fe". More detail is in example 3.
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<unsigned int> dof_indices;

  // Now we will loop over all the elements in the mesh.  We will
  // compute the element matrix and right-hand-side contribution.  See
  // example 3 for a discussion of the element iterators.

  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end(); 
  
  for ( ; el != end_el; ++el)
    {
      // Start logging the shape function initialization.
      // This is done through a simple function call with
      // the name of the event to log.
      perf_log.start_event("elem init");      

      // Store a pointer to the element we are currently
      // working on.  This allows for nicer syntax later.
      const Elem* elem = *el;

      // Get the degree of freedom indices for the
      // current element.  These define where in the global
      // matrix and right-hand-side this element will
      // contribute to.
      dof_map.dof_indices (elem, dof_indices);

      // Compute the element-specific data for the current
      // element.  This involves computing the location of the
      // quadrature points (q_point) and the shape functions
      // (phi, dphi) for the current element.
      fe->reinit (elem);

      // Zero the element matrix and right-hand side before
      // summing them.
      Ke.resize (dof_indices.size(),
		 dof_indices.size());

      Fe.resize (dof_indices.size());

      // Stop logging the shape function initialization.
      // If you forget to stop logging an event the PerfLog
      // object will probably catch the error and abort.
      perf_log.stop_event("elem init");      

      // Now we will build the element matrix.  This involves
      // a double loop to integrate laplacians of the test funcions
      // (i) against laplacians of the trial functions (j).
      //
      // This step is why we need the Clough-Tocher elements -
      // these C1 differentiable elements have square-integrable
      // second derivatives.
      //
      // Now start logging the element matrix computation
      perf_log.start_event ("Ke");

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
	for (unsigned int i=0; i<phi.size(); i++)
	  for (unsigned int j=0; j<phi.size(); j++)
	    Ke(i,j) += JxW[qp]*(d2phi[i][qp](0,0)+d2phi[i][qp](1,1))
			    *(d2phi[j][qp](0,0)+d2phi[j][qp](1,1));
	    

      // Stop logging the matrix computation
      perf_log.stop_event ("Ke");


      // At this point the interior element integration has
      // been completed.  However, we have not yet addressed
      // boundary conditions.  For this example we will only
      // consider simple Dirichlet boundary conditions imposed
      // via the penalty method.  Note that this is a fourth-order
      // problem: Dirichlet boundary conditions include *both*
      // boundary values and boundary normal fluxes.
      {
	// Start logging the boundary condition computation
	perf_log.start_event ("BCs");

	// The penalty value.  
	const Real penalty = 1e10;
	const Real penalty2 = 1e10;

	// The following loops over the sides of the element.
	// If the element has no neighbor on a side then that
	// side MUST live on a boundary of the domain.
	for (unsigned int s=0; s<elem->n_sides(); s++)
	  if (elem->neighbor(s) == NULL)
	    {
              // The value of the shape functions at the quadrature
              // points.
	      const std::vector<std::vector<Real> >&  phi_face =
			      fe_face->get_phi();

	      // The value of the shape function derivatives at the
	      // quadrature points.
              const std::vector<std::vector<RealGradient> >& dphi_face =
			      fe_face->get_dphi();

              // The Jacobian * Quadrature Weight at the quadrature
              // points on the face.
              const std::vector<Real>& JxW_face = fe_face->get_JxW();
                                                                               
              // The XYZ locations (in physical space) of the
              // quadrature points on the face.  This is where
              // we will interpolate the boundary value function.
              const std::vector<Point >& qface_point = fe_face->get_xyz();

	      const std::vector<Point>& face_normals =
			      fe_face->get_normals();

              // Compute the shape function values on the element
              // face.
              fe_face->reinit(elem, s);
                                                                                
              // Loop over the face quagrature points for integration.
              for (unsigned int qp=0; qp<qface.n_points(); qp++)
                {
                  // The boundary value.
		  const Number value = exact_solution(qface_point[qp],
						      es.parameters, "null",
						      "void");
		  const Gradient flux =
				  exact_derivative(qface_point[qp], es.parameters,
						   "null", "void");

                  // Matrix contribution of the L2 projection.
		  // Note that the basis function values are
		  // integrated against test function values while
		  // basis fluxes are integrated against test function
		  // fluxes.
                  for (unsigned int i=0; i<phi_face.size(); i++)
                    for (unsigned int j=0; j<phi_face.size(); j++)
		      Ke(i,j) += JxW_face[qp] *
				 (penalty * phi_face[i][qp] *
				  phi_face[j][qp] + penalty2
				  * (dphi_face[i][qp] *
				  face_normals[qp]) *
				  (dphi_face[j][qp] *
				   face_normals[qp]));

                  // Right-hand-side contribution of the L2
                  // projection.
                  for (unsigned int i=0; i<phi_face.size(); i++)
                    Fe(i) += JxW_face[qp] *
				    (penalty * value * phi_face[i][qp]
				     + penalty2 * 
				     (flux * face_normals[qp])
				    * (dphi_face[i][qp]
				       * face_normals[qp]));
                }
	    } 
	
	// Stop logging the boundary condition computation
	perf_log.stop_event ("BCs");
      } 

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
	for (unsigned int i=0; i<phi.size(); i++)
	  Fe(i) += JxW[qp]*phi[i][qp]*forcing_function(q_point[qp]);
	    
      

      // The element matrix and right-hand-side are now built
      // for this element.  Add them to the global matrix and
      // right-hand-side vector.  The \p PetscMatrix::add_matrix()
      // and \p PetscVector::add_vector() members do this for us.
      // Start logging the insertion of the local (element)
      // matrix and vector into the global matrix and vector
      perf_log.start_event ("matrix insertion");

      dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);
      system.matrix->add_matrix (Ke, dof_indices);
      system.rhs->add_vector    (Fe, dof_indices);

      // Stop logging the insertion of the local (element)
      // matrix and vector into the global matrix and vector
      perf_log.stop_event ("matrix insertion");
    }

  // That's it.  We don't need to do anything else to the
  // PerfLog.  When it goes out of scope (at this function return)
  // it will print its log to the screen. Pretty easy, huh?

  #else

#endif

}
