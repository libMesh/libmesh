/* $Id: ex14.C,v 1.3 2004-05-24 19:58:36 jwpeterson Exp $ */

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

 // <h1>Example 14 - Laplace Equation in the L-Shaped Domain</h1>
 //
 // This example solves the Laplace equation on the classic "L-shaped"
 // domain with adaptive mesh refinement.  In this case, the exact
 // solution is u(r,\theta) = r^{2/3} * \sin ( (2/3) * \theta), but
 // the standard Kelly error indicator is used to estimate the error.
 // The initial mesh contains three QUAD9 elements which represent the
 // standard quadrants I, II, and III of the domain [-1,1]x[-1,1],
 // i.e.
 // Element 0: [-1,0]x[ 0,1]
 // Element 1: [ 0,1]x[ 0,1]
 // Element 2: [-1,0]x[-1,0]
 // The mesh is provided in the standard libMesh ASCII format file
 // named "lshaped.xda".  In addition, an input file named "ex14.in"
 // is provided which allows the user to set several parameters for
 // the solution so that the problem can be re-run without a
 // re-compile.  The solution technique employed is to have a
 // refinement loop with a linear solve inside followed by a
 // refinement of the grid and projection of the solution to the new grid
 // In the final loop iteration, there is no additional
 // refinement after the solve.  In the input file "ex14.in", the variable
 // "max_r_steps" controls the number of refinement steps,
 // "max_r_level" controls the maximum element refinement level, and
 // "refine_percentage" / "coarsen_percentage" determine the number of
 // elements which will be refined / coarsened at each step.

// LibMesh include files.
#include "mesh.h"
#include "equation_systems.h"
#include "implicit_system.h"
#include "gmv_io.h"
#include "tecplot_io.h"
#include "fe.h"
#include "quadrature_gauss.h"
#include "dense_matrix.h"
#include "dense_vector.h"
#include "sparse_matrix.h"
#include "mesh_refinement.h"
#include "error_vector.h"
#include "kelly_error_estimator.h"
#include "getpot.h"

// Function prototype.  This is the function that will assemble
// the linear system for our Laplace problem.  Note that the
// function will take the \p EquationSystems object and the
// name of the system we are assembling as input.  From the
// \p EquationSystems object we have acess to the \p Mesh and
// other objects we might need.
void assemble_laplace(EquationSystems& es,
                      const std::string& system_name);


// Prototype for calculation of the exact solution.  Useful
// for setting boundary conditions.
Real exact_solution(const Real x,
		    const Real y);
  




int main(int argc, char** argv)
{

  // Initialize libMesh.
  libMesh::init (argc, argv);
  {
    // Set the dimensionality of the mesh = 2
    const unsigned int dim = 2;

    // Create a two-dimensional mesh.
    Mesh mesh (dim);

    // Parse the input file
    GetPot input_file("ex14.in");

    // Read in parameters from the input file
    const unsigned int max_r_steps = input_file("max_r_steps", 3);
    const unsigned int max_r_level = input_file("max_r_level", 3);
    const Real refine_percentage   = input_file("refine_percentage", 0.5);
    const Real coarsen_percentage  = input_file("coarsen_percentage", 0.5);
    
    // Read in the mesh
    mesh.read("lshaped.xda");

    // Mesh Refinement object
    MeshRefinement mesh_refinement(mesh);
    
    // Create an equation systems object.
    EquationSystems equation_systems (mesh);

    // Declare the system and its variables.
    {
      // Creates a system named "Laplace"
      equation_systems.add_system<ImplicitSystem> ("Laplace");

      // Adds the variable "u" to "Laplace".  "u"
      // will be approximated using second-order approximation.
      equation_systems("Laplace").add_variable("u", SECOND);

      // Give the system a pointer to the matrix assembly
      // function.
      equation_systems("Laplace").attach_assemble_function (assemble_laplace);
      
      // Initialize the data structures for the equation system.
      equation_systems.init();

      // Set linear solver max iterations
      equation_systems.set_parameter("linear solver maximum iterations") = 100;

      // Linear solver tolerance.
      equation_systems.set_parameter("linear solver tolerance") = 1.e-12;
      
      // Prints information about the system to the screen.
      equation_systems.print_info();
    }


    // Convenient reference to the system
    ImplicitSystem& system = equation_systems.get_system<ImplicitSystem>("Laplace");

    // A refinement loop.
    for (unsigned int r_step=0; r_step<max_r_steps; r_step++)
      {
	std::cout << "Beginning Solve " << r_step << std::endl;
	
	// Solve the system "Laplace", just like example 2.
	equation_systems("Laplace").solve();

	std::cout << "System has: " << equation_systems.n_dofs()
		  << " degrees of freedom."
		  << std::endl;
	
	std::cout << "Linear solver converged at step: "
		  << system.n_linear_iterations()
		  << ", final residual: "
		  << system.final_linear_residual()
		  << std::endl;
	
	// Possibly refine the mesh
	if (r_step+1 != max_r_steps)
	  {
	    std::cout << "  Refining the mesh..." << std::endl;
		
	    // The \p ErrorVector is a particular \p StatisticsVector
	    // for computing error information on a finite element mesh.
	    ErrorVector error;
		
	    // The \p ErrorEstimator class interrogates a finite element
	    // solution and assigns to each element a positive error value.
	    // This value is used for deciding which elements to refine
	    // and which to coarsen.
	    KellyErrorEstimator error_estimator;
		
	    // Compute the error for each active element using the provided
	    // \p flux_jump indicator.  Note in general you will need to
	    // provide an error estimator specifically designed for your
	    // application.
	    error_estimator.estimate_error (equation_systems,
					    "Laplace",
					    error);
		
	    // This takes the error in \p error and decides which elements
	    // will be coarsened or refined.  Any element within 20% of the
	    // maximum error on any element will be refined, and any
	    // element within 10% of the minimum error on any element might
	    // be coarsened. Note that the elements flagged for refinement
	    // will be refined, but those flagged for coarsening _might_ be
	    // coarsened.
	    mesh_refinement.flag_elements_by_error_fraction (error,
							     refine_percentage,
							     coarsen_percentage,
							     max_r_level);
		
	    // This call actually refines and coarsens the flagged
	    // elements.
	    mesh_refinement.refine_and_coarsen_elements();
		
	    // This call reinitializes the \p EquationSystems object for
	    // the newly refined mesh.  One of the steps in the
	    // reinitialization is projecting the \p solution,
	    // \p old_solution, etc... vectors from the old mesh to
	    // the current one.
	    equation_systems.reinit ();

	    // equation_systems.print_info();
	  }
      }	    




    // Write out the solution
    // After solving the system write the solution
    // to a GMV-formatted plot file.
    GMVIO (mesh).write_equation_systems ("lshaped.gmv",
    					 equation_systems);
  }

  
  // All done.  
  return libMesh::close ();
}




// We now define the exact solution, being careful
// to obtain an angle from atan2 in the correct
// quadrant.
Real exact_solution(const Real x,
		    const Real y)
{
  // The boundary value, given by the exact solution,
  // u_exact = r^(2/3)*sin(2*theta/3).
  Real theta = atan2(y,x);

  // Make sure 0 <= theta <= 2*pi
  if (theta < 0)
    theta += 2. * libMesh::pi;
		  
  return pow(x*x + y*y, 1./3.)*sin(2./3.*theta);
}






// We now define the matrix assembly function for the
// Laplace system.  We need to first compute element
// matrices and right-hand sides, and then take into
// account the boundary conditions, which will be handled
// via a penalty method.
void assemble_laplace(EquationSystems& es,
                      const std::string& system_name)
{
  // It is a good idea to make sure we are assembling
  // the proper system.
  assert (system_name == "Laplace");


  // Declare a performance log.  Give it a descriptive
  // string to identify what part of the code we are
  // logging, since there may be many PerfLogs in an
  // application.
  PerfLog perf_log ("Matrix Assembly",false);
  
    // Get a constant reference to the mesh object.
  const Mesh& mesh = es.get_mesh();

  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to the ImplicitSystem we are solving
  ImplicitSystem& system = es.get_system<ImplicitSystem>("Laplace");
  
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
  
  // A 5th order Gauss quadrature rule for numerical integration.
  QGauss qrule (dim, FIFTH);

  // Tell the finite element object to use our quadrature rule.
  fe->attach_quadrature_rule (&qrule);

  // Declare a special finite element object for
  // boundary integration.
  AutoPtr<FEBase> fe_face (FEBase::build(dim, fe_type));
	      
  // Boundary integration requires one quadraure rule,
  // with dimensionality one less than the dimensionality
  // of the element.
  QGauss qface(dim-1, FIFTH);
  
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
  // const std::vector<Point>& q_point = fe->get_xyz();

  // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector<Real> >& phi = fe->get_phi();

  // The element shape function gradients evaluated at the quadrature
  // points.
  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();

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
  // example 3 for a discussion of the element iterators.  Here we use
  // the \p const_active_local_elem_iterator to indicate we only want
  // to loop over elements that are assigned to the local processor
  // which are "active" in the sense of AMR.  This allows each
  // processor to compute its components of the global matrix for
  // active elements while ignoring parent elements which have been
  // refined.
  const_active_local_elem_iterator           el (mesh.elements_begin());
  const const_active_local_elem_iterator end_el (mesh.elements_end());
  
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
      // summing them.  We use the resize member here because
      // the number of degrees of freedom might have changed from
      // the last element.  Note that this will be the case if the
      // element type is different (i.e. the last element was a
      // triangle, now we are on a quadrilateral).
      Ke.resize (dof_indices.size(),
		 dof_indices.size());

      Fe.resize (dof_indices.size());

      // Stop logging the shape function initialization.
      // If you forget to stop logging an event the PerfLog
      // object will probably catch the error and abort.
      perf_log.stop_event("elem init");      

      // Now we will build the element matrix.  This involves
      // a double loop to integrate the test funcions (i) against
      // the trial functions (j).
      //
      // Now start logging the element matrix computation
      perf_log.start_event ("Ke");

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
	for (unsigned int i=0; i<phi.size(); i++)
	  for (unsigned int j=0; j<phi.size(); j++)
	    Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
	    

      // Stop logging the matrix computation
      perf_log.stop_event ("Ke");


      // At this point the interior element integration has
      // been completed.  However, we have not yet addressed
      // boundary conditions.  For this example we will only
      // consider simple Dirichlet boundary conditions imposed
      // via the penalty method. This is discussed at length in
      // example 3.
      {
	// Start logging the boundary condition computation
	perf_log.start_event ("BCs");

	// The penalty value.  
	const Real penalty = 1.e10;

	// The following loops over the sides of the element.
	// If the element has no neighbor on a side then that
	// side MUST live on a boundary of the domain.
	for (unsigned int s=0; s<elem->n_sides(); s++)
	  if (elem->neighbor(s) == NULL)
	    {
	      // Get a pointer to the side.
	      AutoPtr<Elem> side (elem->build_side(s));
	      
	      // Loop over the nodes on the side with NULL neighbor.
	      for (unsigned int ns=0; ns<side->n_nodes(); ns++)
		{
		  const Real x = side->point(ns)(0);
		  const Real y = side->point(ns)(1);
		  const Real value = exact_solution(x,y);
		  
// 		  std::cout << "(x,y,bc)=("
// 			    << x << ","
// 			    << y << ","
// 			    << value << ")"
// 			    << std::endl;
		  
		  // Find the node on the element matching this node on
		  // the side.  That defines where in the element matrix
		  // the boundary condition will be applied.
		  for (unsigned int n=0; n<elem->n_nodes(); n++)
		    if (elem->node(n) == side->node(ns))
		      {
			// Matrix contribution.
			Ke(n,n) += penalty;
			
			// Right-hand-side contribution.
			Fe(n) += penalty*value;
		      }
		} 
	    } 
	
	// Stop logging the boundary condition computation
	perf_log.stop_event ("BCs");
      } 
      

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

      // Start logging the insertion of the local (element)
      // matrix and vector into the global matrix and vector
      perf_log.stop_event ("matrix insertion");
    }

  // That's it.  We don't need to do anything else to the
  // PerfLog.  When it goes out of scope (at this function return)
  // it will print its log to the screen. Pretty easy, huh?
}
