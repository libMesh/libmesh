/* $Id: ex10.C,v 1.32 2007-08-23 18:07:45 roystgnr Exp $ */

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

  // <h1>Example 10 - Solving a Transient System with Adaptive Mesh Refinement</h1>
  // 
  // This example shows how a simple, linear transient
  // system can be solved in parallel.  The system is simple
  // scalar convection-diffusion with a specified external
  // velocity.  The initial condition is given, and the
  // solution is advanced in time with a standard Crank-Nicolson
  // time-stepping strategy.
  //
  // Also, we use this example to demonstrate writing out and reading
  // in of solutions. We do 25 time steps, then save the solution
  // and do another 25 time steps starting from the saved solution.
 
// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <cmath>

// Basic include file needed for the mesh functionality.
#include "libmesh.h"
#include "mesh.h"
#include "mesh_refinement.h"
#include "gmv_io.h"
#include "equation_systems.h"
#include "fe.h"
#include "quadrature_gauss.h"
#include "dof_map.h"
#include "sparse_matrix.h"
#include "numeric_vector.h"
#include "dense_matrix.h"
#include "dense_vector.h"

#include "getpot.h"

// Some (older) compilers do not offer full stream 
// functionality, \p OStringStream works around this.
// Check example 9 for details.
#include "o_string_stream.h"

// This example will solve a linear transient system,
// so we need to include the \p TransientLinearImplicitSystem definition.
#include "transient_system.h"
#include "linear_implicit_system.h"
#include "vector_value.h"

// To refine the mesh we need an \p ErrorEstimator
// object to figure out which elements to refine.
#include "error_vector.h"
#include "kelly_error_estimator.h"

// The definition of a geometric element
#include "elem.h"

// Function prototype.  This function will assemble the system
// matrix and right-hand-side at each time step.  Note that
// since the system is linear we technically do not need to
// assmeble the matrix at each time step, but we will anyway.
// In subsequent examples we will employ adaptive mesh refinement,
// and with a changing mesh it will be necessary to rebuild the
// system matrix.
void assemble_cd (EquationSystems& es,
		  const std::string& system_name);

// Function prototype.  This function will initialize the system.
// Initialization functions are optional for systems.  They allow
// you to specify the initial values of the solution.  If an
// initialization function is not provided then the default (0)
// solution is provided.
void init_cd (EquationSystems& es,
	      const std::string& system_name);

// Exact solution function prototype.  This gives the exact
// solution as a function of space and time.  In this case the
// initial condition will be taken as the exact solution at time 0,
// as will the Dirichlet boundary conditions at time t.
Real exact_solution (const Real x,
		     const Real y,
		     const Real t);

Number exact_value (const Point& p,
		    const Parameters& parameters,
		    const std::string&,
		    const std::string&)
{
  return exact_solution(p(0), p(1), parameters.get<Real> ("time"));
}



// Begin the main program.  Note that the first
// statement in the program throws an error if
// you are in complex number mode, since this
// example is only intended to work with real
// numbers.
int main (int argc, char** argv)
{
  // Initialize libMesh.
  libMesh::init (argc, argv);

#ifdef ENABLE_PARMESH
  if (libMesh::processor_id() == 0)
    std::cerr << "ERROR: This example uses AMR, which libMesh\n"
              << "does not yet support on parallel meshes!"
              << std::endl;
  return 0;
#endif

#ifndef ENABLE_AMR
  if (libMesh::processor_id() == 0)
    std::cerr << "ERROR: This example requires libMesh to be\n"
              << "compiled with AMR support!"
              << std::endl;
  return 0;
#else

  {    

    // Brief message to the user regarding the program name
    // and command line arguments.

    // Use commandline parameter to specify if we are to
    // read in an initial solution or generate it ourself
    std::cout << "Usage:\n"
      <<"\t " << argv[0] << " -init_timestep 0\n"
      << "OR\n"
      <<"\t " << argv[0] << " -read_solution -init_timestep 26\n"
      << std::endl;

    std::cout << "Running: " << argv[0];

    for (int i=1; i<argc; i++)
      std::cout << " " << argv[i];

    std::cout << std::endl << std::endl;

    // Create a GetPot object to parse the command line
    GetPot command_line (argc, argv);


    // This boolean value is obtained from the command line, it is true
    // if the flag "-read_solution" is present, false otherwise.
    // It indicates whether we are going to read in
    // the mesh and solution files "saved_mesh.xda" and "saved_solution.xda"
    // or whether we are going to start from scratch by just reading
    // "mesh.xda"
    const bool read_solution   = command_line.search("-read_solution");

    // This value is also obtained from the commandline and it specifies the
    // initial value for the t_step looping variable. We must
    // distinguish between the two cases here, whether we read in the 
    // solution or we started from scratch, so that we do not overwrite the
    // gmv output files.
    unsigned int init_timestep = 0;
    
    // Search the command line for the "init_timestep" flag and if it is
    // present, set init_timestep accordingly.
    if(command_line.search("-init_timestep"))
      init_timestep = command_line.next(0);
    else
    {
      if (libMesh::processor_id() == 0)
        std::cerr << "ERROR: Initial timestep not specified\n" << std::endl;

      // This handy function will print the file name, line number,
      // and then abort.  Currrently the library does not use C++
      // exception handling.
      error();
    }



    // This value is also obtained from the command line, and specifies
    // the number of time steps to take.
    unsigned int n_timesteps = 0;

    // Again do a search on the command line for the argument
    if(command_line.search("-n_timesteps"))
      n_timesteps = command_line.next(0);
    else
    {
      std::cout << "ERROR: Number of timesteps not specified\n" << std::endl;
      error();
    }


    // Create a two-dimensional mesh.
    Mesh mesh (2);

    // Create an equation systems object.
    EquationSystems equation_systems (mesh);
    MeshRefinement mesh_refinement (mesh);

    // First we process the case where we do not read in the solution
    if(!read_solution)
    {
      // Read the mesh from file.
      mesh.read ("mesh.xda");

      // Uniformly refine the mesh 5 times
      if(!read_solution)
        mesh_refinement.uniformly_refine (5);

      // Print information about the mesh to the screen.
      mesh.print_info();
      
      
      // Declare the system and its variables.
      // Begin by creating a transient system
      // named "Convection-Diffusion".
      TransientLinearImplicitSystem & system = 
        equation_systems.add_system<TransientLinearImplicitSystem> 
        ("Convection-Diffusion");

      // Adds the variable "u" to "Convection-Diffusion".  "u"
      // will be approximated using first-order approximation.
      system.add_variable ("u", FIRST);

      // Give the system a pointer to the matrix assembly
      // and initialization functions.
      system.attach_assemble_function (assemble_cd);
      system.attach_init_function (init_cd);

      // Initialize the data structures for the equation system.
      equation_systems.init ();
    }

    // Otherwise we read in the solution and mesh
    else 
    {
      // Read in the mesh stored in "saved_mesh.xda"
      mesh.read("saved_mesh.xda");

      // Print information about the mesh to the screen.
      mesh.print_info();

      // Read in the solution stored in "saved_solution.xda"
      equation_systems.read("saved_solution.xda", libMeshEnums::READ);

      // Get a reference to the system so that we can call update() on it
      TransientLinearImplicitSystem & system = 
        equation_systems.get_system<TransientLinearImplicitSystem> 
        ("Convection-Diffusion");

      // We need to call update to put system in a consistent state
      // with the solution that was read in
      system.update();

      // Attach the same matrix assembly function as above. Note, we do not
      // have to attach an init() function since we are initializing the
      // system by reading in "saved_solution.xda"
      system.attach_assemble_function (assemble_cd);
    }

    // Prints information about the system to the screen.
    equation_systems.print_info();

    
    equation_systems.parameters.set<unsigned int>
      ("linear solver maximum iterations") = 250;
    equation_systems.parameters.set<Real>
      ("linear solver tolerance") = TOLERANCE;
      
    if(!read_solution)
      // Write out the initial condition
      GMVIO(mesh).write_equation_systems ("out.gmv.000",
					  equation_systems);
    else
      // Write out the solution that was read in
      GMVIO(mesh).write_equation_systems ("solution_read_in.gmv",
					  equation_systems);

      
    
    // The Convection-Diffusion system requires that we specify
    // the flow velocity.  We will specify it as a RealVectorValue
    // data type and then use the Parameters object to pass it to
    // the assemble function.
    equation_systems.parameters.set<RealVectorValue>("velocity") = 
      RealVectorValue (0.8, 0.8);
    
    // Solve the system "Convection-Diffusion".  This will be done by
    // looping over the specified time interval and calling the
    // \p solve() member at each time step.  This will assemble the
    // system and call the linear solver.
    const Real dt = 0.025;
    Real time     = init_timestep*dt;
    
    // We do 25 timesteps both before and after writing out the
    // intermediate solution
    for(unsigned int t_step=init_timestep; 
                     t_step<(init_timestep+n_timesteps); 
                     t_step++)
      {
	// Increment the time counter, set the time and the
	// time step size as parameters in the EquationSystem.
	time += dt;

	equation_systems.parameters.set<Real> ("time") = time;
	equation_systems.parameters.set<Real> ("dt")   = dt;

	// A pretty update message
	std::cout << " Solving time step ";
	
	// As already seen in example 9, use a work-around
	// for missing stream functionality (of older compilers).
	// Add a set of scope braces to enforce data locality.
	{
	  OStringStream out;

	  OSSInt(out,2,t_step);
	  out << ", time=";
	  OSSRealzeroleft(out,6,3,time);
	  out <<  "...";
	  std::cout << out.str() << std::endl;
	}
	
	// At this point we need to update the old
	// solution vector.  The old solution vector
	// will be the current solution vector from the
	// previous time step.  We will do this by extracting the
	// system from the \p EquationSystems object and using
	// vector assignment.  Since only \p TransientLinearImplicitSystems
	// (and systems derived from them) contain old solutions
	// we need to specify the system type when we ask for it.
	TransientLinearImplicitSystem &  system =
	  equation_systems.get_system<TransientLinearImplicitSystem>("Convection-Diffusion");

	*system.old_local_solution = *system.current_local_solution;
	
	// The number of refinement steps per time step.
	const unsigned int max_r_steps = 2;
	
	// A refinement loop.
	for (unsigned int r_step=0; r_step<max_r_steps; r_step++)
	  {
	    // Assemble & solve the linear system
	    system.solve();
	    
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
		//ErrorEstimator* error_estimator = new KellyErrorEstimator;
		KellyErrorEstimator error_estimator;
		
		// Compute the error for each active element using the provided
		// \p flux_jump indicator.  Note in general you will need to
		// provide an error estimator specifically designed for your
		// application.
		error_estimator.estimate_error (system,
						error);
		
		// This takes the error in \p error and decides which elements
		// will be coarsened or refined.  Any element within 20% of the
		// maximum error on any element will be refined, and any
		// element within 7% of the minimum error on any element might
		// be coarsened. Note that the elements flagged for refinement
		// will be refined, but those flagged for coarsening _might_ be
		// coarsened.
		mesh_refinement.refine_fraction() = 0.80;
		mesh_refinement.coarsen_fraction() = 0.07;
		mesh_refinement.max_h_level() = 5;
		mesh_refinement.flag_elements_by_error_fraction (error);
		
		// This call actually refines and coarsens the flagged
		// elements.
		mesh_refinement.refine_and_coarsen_elements();
		
		// This call reinitializes the \p EquationSystems object for
		// the newly refined mesh.  One of the steps in the
		// reinitialization is projecting the \p solution,
		// \p old_solution, etc... vectors from the old mesh to
		// the current one.
		equation_systems.reinit ();
	      }	    
	  }
	
	// Output evey 10 timesteps to file.
	if ( (t_step+1)%10 == 0)
	  {
	    OStringStream file_name;

	    file_name << "out.gmv.";
	    OSSRealzeroright(file_name,3,0,t_step+1);

	    GMVIO(mesh).write_equation_systems (file_name.str(),
						equation_systems);
	  }
      }

      if(!read_solution)
      {
        mesh.write("saved_mesh.xda");
        equation_systems.write("saved_solution.xda", libMeshEnums::WRITE);
      }
  }
#endif // #ifndef ENABLE_AMR
  
  // All done.  
  return libMesh::close ();
}

// Here we define the initialization routine for the
// Convection-Diffusion system.  This routine is
// responsible for applying the initial conditions to
// the system.
void init_cd (EquationSystems& es,
	      const std::string& system_name)
{
  // It is a good idea to make sure we are initializing
  // the proper system.
  assert (system_name == "Convection-Diffusion");

  // Get a reference to the Convection-Diffusion system object.
  TransientLinearImplicitSystem & system =
    es.get_system<TransientLinearImplicitSystem>("Convection-Diffusion");

  // Project initial conditions at time 0
  es.parameters.set<Real> ("time") = 0;
  
  system.project_solution(exact_value, NULL, es.parameters);
}



// This function defines the assembly routine which
// will be called at each time step.  It is responsible
// for computing the proper matrix entries for the
// element stiffness matrices and right-hand sides.
void assemble_cd (EquationSystems& es,
		  const std::string& system_name)
{
#ifdef ENABLE_AMR
  // It is a good idea to make sure we are assembling
  // the proper system.
  assert (system_name == "Convection-Diffusion");
  
  // Get a constant reference to the mesh object.
  const MeshBase& mesh = es.get_mesh();
  
  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();
  
  // Get a reference to the Convection-Diffusion system object.
  TransientLinearImplicitSystem & system =
    es.get_system<TransientLinearImplicitSystem> ("Convection-Diffusion");
  
  // Get the Finite Element type for the first (and only) 
  // variable in the system.
  FEType fe_type = system.variable_type(0);
  
  // Build a Finite Element object of the specified type.  Since the
  // \p FEBase::build() member dynamically creates memory we will
  // store the object as an \p AutoPtr<FEBase>.  This can be thought
  // of as a pointer that will clean up after itself.
  AutoPtr<FEBase> fe      (FEBase::build(dim, fe_type));
  AutoPtr<FEBase> fe_face (FEBase::build(dim, fe_type));
  
  // A Gauss quadrature rule for numerical integration.
  // Let the \p FEType object decide what order rule is appropriate.
  QGauss qrule (dim,   fe_type.default_quadrature_order());
  QGauss qface (dim-1, fe_type.default_quadrature_order());

  // Tell the finite element object to use our quadrature rule.
  fe->attach_quadrature_rule      (&qrule);
  fe_face->attach_quadrature_rule (&qface);

  // Here we define some references to cell-specific data that
  // will be used to assemble the linear system.  We will start
  // with the element Jacobian * quadrature weight at each integration point.
  const std::vector<Real>& JxW      = fe->get_JxW();
  const std::vector<Real>& JxW_face = fe_face->get_JxW();
  
  // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector<Real> >& phi = fe->get_phi();
  const std::vector<std::vector<Real> >& psi = fe_face->get_phi();

  // The element shape function gradients evaluated at the quadrature
  // points.
  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();

  // The XY locations of the quadrature points used for face integration
  const std::vector<Point>& qface_points = fe_face->get_xyz();
    
  // A reference to the \p DofMap object for this system.  The \p DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.  We will talk more about the \p DofMap
  // in future examples.
  const DofMap& dof_map = system.get_dof_map();
  
  // Define data structures to contain the element matrix
  // and right-hand-side vector contribution.  Following
  // basic finite element terminology we will denote these
  // "Ke" and "Fe".
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;
  
  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<unsigned int> dof_indices;

  // Here we extract the velocity & parameters that we put in the
  // EquationSystems object.
  const RealVectorValue velocity =
    es.parameters.get<RealVectorValue> ("velocity");

  const Real dt = es.parameters.get<Real>   ("dt");
  const Real time = es.parameters.get<Real> ("time");

  // Now we will loop over all the elements in the mesh that
  // live on the local processor. We will compute the element
  // matrix and right-hand-side contribution.  Since the mesh
  // will be refined we want to only consider the ACTIVE elements,
  // hence we use a variant of the \p active_elem_iterator.
  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end(); 
  
  for ( ; el != end_el; ++el)
    {    
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
      
      // Now we will build the element matrix and right-hand-side.
      // Constructing the RHS requires the solution and its
      // gradient from the previous timestep.  This myst be
      // calculated at each quadrature point by summing the
      // solution degree-of-freedom values by the appropriate
      // weight functions.
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
	{
	  // Values to hold the old solution & its gradient.
	  Number   u_old = 0.;
	  Gradient grad_u_old;
	  
	  // Compute the old solution & its gradient.
	  for (unsigned int l=0; l<phi.size(); l++)
	    {
	      u_old      += phi[l][qp]*system.old_solution  (dof_indices[l]);
	      
	      // This will work,
	      // grad_u_old += dphi[l][qp]*system.old_solution (dof_indices[l]);
	      // but we can do it without creating a temporary like this:
	      grad_u_old.add_scaled (dphi[l][qp],system.old_solution (dof_indices[l]));
	    }
	  
	  // Now compute the element matrix and RHS contributions.
	  for (unsigned int i=0; i<phi.size(); i++)
	    {
	      // The RHS contribution
	      Fe(i) += JxW[qp]*(
				// Mass matrix term
				u_old*phi[i][qp] + 
				-.5*dt*(
					// Convection term
					// (grad_u_old may be complex, so the
					// order here is important!)
					(grad_u_old*velocity)*phi[i][qp] +
					
					// Diffusion term
					0.01*(grad_u_old*dphi[i][qp]))     
				);
	      
	      for (unsigned int j=0; j<phi.size(); j++)
		{
		  // The matrix contribution
		  Ke(i,j) += JxW[qp]*(
				      // Mass-matrix
				      phi[i][qp]*phi[j][qp] + 
				      .5*dt*(
					     // Convection term
					     (velocity*dphi[j][qp])*phi[i][qp] +
					     // Diffusion term
					     0.01*(dphi[i][qp]*dphi[j][qp]))      
				      );
		}
	    } 
	} 

      // At this point the interior element integration has
      // been completed.  However, we have not yet addressed
      // boundary conditions.  For this example we will only
      // consider simple Dirichlet boundary conditions imposed
      // via the penalty method. 
      //	
      // The following loops over the sides of the element.
      // If the element has no neighbor on a side then that
      // side MUST live on a boundary of the domain.
      {
	// The penalty value.  
	const Real penalty = 1.e10;

	// The following loops over the sides of the element.
	// If the element has no neighbor on a side then that
	// side MUST live on a boundary of the domain.
	for (unsigned int s=0; s<elem->n_sides(); s++)
	  if (elem->neighbor(s) == NULL)
	    {
	      fe_face->reinit(elem,s);
	      
	      for (unsigned int qp=0; qp<qface.n_points(); qp++)
		{
		  const Number value = exact_solution (qface_points[qp](0),
						       qface_points[qp](1),
						       time);
						       
		  // RHS contribution
		  for (unsigned int i=0; i<psi.size(); i++)
		    Fe(i) += penalty*JxW_face[qp]*value*psi[i][qp];

		  // Matrix contribution
		  for (unsigned int i=0; i<psi.size(); i++)
		    for (unsigned int j=0; j<psi.size(); j++)
		      Ke(i,j) += penalty*JxW_face[qp]*psi[i][qp]*psi[j][qp];
		}
	    } 
      } 

      
      // We have now built the element matrix and RHS vector in terms
      // of the element degrees of freedom.  However, it is possible
      // that some of the element DOFs are constrained to enforce
      // solution continuity, i.e. they are not really "free".  We need
      // to constrain those DOFs in terms of non-constrained DOFs to
      // ensure a continuous solution.  The
      // \p DofMap::constrain_element_matrix_and_vector() method does
      // just that.
      dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);
      
      // The element matrix and right-hand-side are now built
      // for this element.  Add them to the global matrix and
      // right-hand-side vector.  The \p PetscMatrix::add_matrix()
      // and \p PetscVector::add_vector() members do this for us.
      system.matrix->add_matrix (Ke, dof_indices);
      system.rhs->add_vector    (Fe, dof_indices);
      
    }
  // Finished computing the sytem matrix and right-hand side.
#endif // #ifdef ENABLE_AMR
}
