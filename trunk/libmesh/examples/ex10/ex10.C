// $Id: ex10.C,v 1.4 2003-06-04 22:47:44 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2003  Benjamin S. Kirk
  
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




/**
 * C++ include files that we need
 */
#include <iostream>
#include <algorithm>
#include <math.h>

/**
 * Basic include file needed for the mesh functionality.
 */
#include "libmesh.h"
#include "mesh.h"
#include "equation_systems.h"
#include "fe.h"
#include "quadrature_gauss.h"
#include "dof_map.h"
#include "sparse_matrix.h"
#include "numeric_vector.h"
#include "dense_matrix.h"
#include "dense_vector.h"

/**
 * Some (older) compilers do not offer full stream 
 * functionality, \p OStringStream works around this.
 * Check example 9 for details.
 */
#include "o_string_stream.h"

/**
 * This example will solve a linear transient system,
 * so we need to include the \p TransientSystem definition.
 */
#include "transient_system.h"
#include "vector_value.h"

/**
 * To refine the mesh we need an \p ErrorEstimator
 * object to figure out which elements to refine.
 */
#include "error_vector.h"
#include "error_estimator.h"


/**
 * \mainpage Example 10
 *
 * \section Introduction
 *
 * This example shows how a simple, linear transient
 * system can be solved in parallel.  The system is simple
 * scalar convection-diffusion with a specified external
 * velocity.  The initial condition is given, and the
 * solution is advanced in time with a standard Crank-Nicholson
 * time-stepping strategy.
 */


/**
 * Function prototype.  This function will assemble the system
 * matrix and right-hand-side at each time step.  Note that
 * since the system is linear we technically do not need to
 * assmeble the matrix at each time step, but we will anyway.
 * In subsequent examples we will employ adaptive mesh refinement,
 * and with a changing mesh it will be necessary to rebuild the
 * system matrix.
 */
void assemble_cd (EquationSystems& es,
		  const std::string& system_name);


/**
 * Function prototype.  This function will initialize the system.
 * Initialization functions are optional for systems.  They allow
 * you to specify the initial values of the solution.  If an
 * initialization function is not provided then the default (0)
 * solution is provided.
 */
void init_cd (EquationSystems& es,
	      const std::string& system_name);



/**
 * Exact solution function prototype.  This gives the exact
 * solution as a function of space and time.  In this case the
 * initial condition will be taken as the exact solution at time 0,
 * as will the Dirichlet boundary conditions at time t.
 */
Real exact_solution (const Real x,
		     const Real y,
		     const Real t);





int main (int argc, char** argv)
{
#ifndef HAVE_PETSC

  std::cerr << "ERROR: AMR only good with PETSc for some strange reason..."
	    << "       Bug Ben about this."
	    << std::endl;
  here();

  return 0;
  
#endif
  
#ifdef USE_COMPLEX_NUMBERS
  
  std::cerr << "ERROR: Not intended for use with complex numbers."
	    << std::endl;
  here();

  return 0;
  
#else
  
  /**
   * Initialize libMesh.
   */
  libMesh::init (argc, argv);

  {    
    /**
     * Create a two-dimensional mesh.
     */
    Mesh mesh (2);
    
    /**
     * Read the mesh from file.
     */
    mesh.read ("mesh.xda");

    /**
     * Uniformly refine the mesh 5 times.  This is the
     * first time we use the mesh refinement capabilities
     * of the library.
     */
    mesh.mesh_refinement.uniformly_refine (5);
    
    /**
     * Print information about the mesh to the screen.
     */
    mesh.print_info();
    
    /**
     * Create an equation systems object.
     */
    EquationSystems equation_systems (mesh);
    
    /**
     * Declare the system and its variables.
     */
    {
      /**
       * Creates a transient system named "Convection-Diffusion"
       */
      TransientSystem& system = 
	equation_systems.add_system<TransientSystem> ("Convection-Diffusion");
      
      /**
       * Adds the variable "u" to "Convection-Diffusion".  "u"
       * will be approximated using first-order approximation.
       */
      system.add_variable ("u", FIRST);

      /**
       * Give the system a pointer to the matrix assembly
       * and initialization functions.
       */
      system.attach_assemble_function (assemble_cd);
      system.attach_init_function (init_cd);
      
      /**
       * Initialize the data structures for the equation system.
       */
      equation_systems.init ();
      
      /**
       * Prints information about the system to the screen.
       */
      equation_systems.print_info();

      /**
       * Write out the initial conditions.
       */
      mesh.write_gmv ("out_000.gmv",
		      equation_systems);
    }



    /**
     * The Convection-Diffusion system requires that we specify
     * the flow velocity.  We will specify it as a \p RealVectorValue
     * data type.  Check example 9 for details on \p RealVectorValue
     * and \p DataMap.
     */
    RealVectorValue velocity (0.8, 0.8);

    equation_systems.data_map.add_data ("velocity", velocity);
    
    /**
     * Solve the system "Convection-Diffusion".  This will be done by
     * looping over the specified time interval and calling the
     * \p solve() member at each time step.  This will assemble the
     * system and call the linear solver.
     */
    const Real dt = 0.025;
    Real time     = 0.;
    
    for (unsigned int t_step = 0; t_step < 50; t_step++)
      {
	/**
	 * Increment the time counter, set the time and the
	 * time step size as parameters in the EquationSystem.
	 */
	time += dt;

	equation_systems.set_parameter ("time") = time;
	equation_systems.set_parameter ("dt")   = dt;

	// A pretty update message
	std::cout << " Solving time step ";
	/**
	 * As already seen in example 9, use a work-around
	 * for missing stream functionality (of older compilers).
	 */
	{
	  OStringStream out;

	  OSSInt(out,2,t_step);
	  out << ", time=";
	  OSSRealzeroleft(out,6,3,time);
	  out <<  "..." << std::endl;
	  std::cout << out.str();
	}
	
	/**
	 * At this point we need to update the old
	 * solution vector.  The old solution vector
	 * will be the current solution vector from the
	 * previous time step.  We will do this by extracting the
	 * system from the \p EquationSystems object and using
	 * vector assignment.  Since only \p TransientSystems
	 * (and systems derived from them) contain old solutions
	 * we need to specify the system type when we ask for it.
	 */
	TransientSystem&  system =
	  equation_systems.get_system<TransientSystem>("Convection-Diffusion");

	*system.old_local_solution = *system.current_local_solution;


	/**
	 * The number of refinement steps per time step.
	 */
	const unsigned int max_r_steps = 2;

	/**
	 * A refinement loop.
	 */
	for (unsigned int r_step=0; r_step<max_r_steps; r_step++)
	  {
	    /**
	     * Assemble & solve the linear system
	     */
	    system.solve();

	    /**
	     * Possibly refine the mesh
	     */
	    if (r_step+1 != max_r_steps)
	      {
		std::cout << "  Refining the mesh..." << std::endl;

		/**
		 * The \p ErrorVector is a particular \p StatisticsVector
		 * for computing error information on a finite element mesh.
		 */
		ErrorVector error;

		/**
		 * The \p ErrorEstimator class interrogates a finite element
		 * solution and assigns to each element a positive error value.
		 * This value is used for deciding which elements to refine
		 * and which to coarsen.
		 */
		ErrorEstimator error_estimator;

		/**
		 * Compute the error for each active element using the provided
		 * \p flux_jump indicator.  Note in general you will need to
		 * provide an error estimator specifically designed for your
		 * application.
		 */
		error_estimator.flux_jump (equation_systems,
					   "Convection-Diffusion",
					   error);

		/**
		 * This takes the error in \p error and decides which elements
		 * will be coarsened or refined.  Any element within 20% of the
		 * maximum error on any element will be refined, and any
		 * element within 10% of the minimum error on any element might
		 * be coarsened. Note that the elements flagged for refinement
		 * will be refined, but those flagged for coarsening _might_ be
		 * coarsened.
		 */
		mesh.mesh_refinement.flag_elements_by_error_fraction (error,
								      0.80,
								      0.07,
								      5);
		
		/**
		 * This call actually refines and coarsens the flagged
		 * elements.
		 */
		mesh.mesh_refinement.refine_and_coarsen_elements();

		/**
		 * This call reinitializes the \p EquationSystems object for
		 * the newly refined mesh.  One of the steps in the
		 * reinitialization is projecting the \p solution,
		 * \p old_solution, etc... vectors from the old mesh to
		 * the current one.
		 */
		equation_systems.reinit ();
	      }	    
	  }

	/**
	 * Output evey 10 timesteps to file.
	 */
	if ( (t_step+1)%10 == 0)
	  {
	    OStringStream file_name;

	    file_name << "out_";
	    OSSRealzeroright(file_name,3,0,t_step+1);
	    file_name << ".gmv";

// OLD CODE
// 	    std::stringstream file_name;

// 	    file_name << "out_";
// 	    file_name.fill('0');
// 	    file_name.width(3);
// 	    file_name << std::right << t_step+1;
// 	    file_name << ".gmv";

	    mesh.write_gmv (file_name.str(),
			    equation_systems);
	  }
      }
  }


  /**
   * All done.  
   */
  return libMesh::close ();

#endif
}




void init_cd (EquationSystems& es,
	      const std::string& system_name)
{
#ifndef USE_COMPLEX_NUMBERS
  
  /**
   * It is a good idea to make sure we are initializing
   * the proper system.
   */
  assert (system_name == "Convection-Diffusion");
  
  /**
   * Get a constant reference to the mesh object.
   */
  const Mesh& mesh = es.get_mesh();

  /**
   * Get a reference to the Convection-Diffusion system object.
   */
  TransientSystem& system =
    es.get_system<TransientSystem> ("Convection-Diffusion");
  
  /**
   * Get a reference to the \p DofMap for this system.
   */
  const DofMap& dof_map = system.get_dof_map();
  
  /**
   * Get a reference to the solution vector.
   */
  NumericVector<Real>& solution = *system.solution;

  /**
   * A vector to hold the global DOF indices for this element.
   */
  std::vector<unsigned int> dof_indices;

  
  

  /**
   *---------------------------------------------------------------
   * Loop over the local elements and compute the initial value
   * of the solution at the element degrees of freedom.  Assign
   * these initial values to the solution vector.  There is a small
   * catch, however...  We only want to assign the components that
   * live on the local processor, hence there will be an if-test
   * in the loop.
   */

  const_active_local_elem_iterator       elem_it (mesh.elements_begin());
  const const_active_local_elem_iterator elem_end(mesh.elements_end());

  for ( ; elem_it != elem_end; ++elem_it)
    {
      const Elem* elem = *elem_it;

      dof_map.dof_indices (elem, dof_indices);

      /**
       * For these Lagrange-elements the number
       * of degrees of freedom should be <= the number
       * of nodes.
       */
      assert (dof_indices.size() <= elem->n_nodes());

      /**
       * Loop over the element DOFs, compute the initial
       * value if the DOF is local to the processor.
       */
      for (unsigned int i=0; i<dof_indices.size(); i++)
	if ((dof_indices[i] >= solution.first_local_index()) &&
	    (dof_indices[i] <  solution.last_local_index()))
	  {
	    const Point&  p = elem->point (i);
	    const Real    x = p(0);
	    const Real    y = p(1);
	    const Real time = 0.;
	    
	    solution.set (dof_indices[i],
			  exact_solution (x,y,time));	    
	  }	 
    }

  /**
   * The initial solution has now been set for the local
   * solution components.  However, the local matrix assembly
   * will likely depend on solution components that live on
   * other processors.  We need to get those components, and
   * the \p TransientSystem::update() member will do that
   * for us.
   */
  system.update ();

#endif 
}




void assemble_cd (EquationSystems& es,
		  const std::string& system_name)
{
#ifndef USE_COMPLEX_NUMBERS
  
  /**
   * It is a good idea to make sure we are assembling
   * the proper system.
   */
  assert (system_name == "Convection-Diffusion");
  
  /**
   * Get a constant reference to the mesh object.
   */
  const Mesh& mesh = es.get_mesh();

  /**
   * The dimension that we are running
   */
  const unsigned int dim = mesh.mesh_dimension();

  /**
   * Get a reference to the Convection-Diffusion system object.
   */
  TransientSystem& system =
    es.get_system<TransientSystem> ("Convection-Diffusion");
  
  /**
   * Get the Finite Element type for the first (and only) 
   * variable in the system.
   */
  FEType fe_type = system.variable_type(0);

  /**
   * Build a Finite Element object of the specified type.  Since the
   * \p FEBase::build() member dynamically creates memory we will
   * store the object as an \p AutoPtr<FEBase>.  The \p AutoPtr's
   * are covered in more detail in example 5.
   */
  AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));
  
  /**
   * A Gauss quadrature rule for numerical integration.
   * Let the \p FEType object decide what order rule is appropriate.
   */
  QGauss qrule (dim, fe_type.default_quadrature_order());

  /**
   * Tell the finite element object to use our quadrature rule.
   */
  fe->attach_quadrature_rule (&qrule);


  

  /**
   *--------------------------------------------------------------------
   * Here we define some references to cell-specific data that
   * will be used to assemble the linear system.
   */
  
  /**
   * The element Jacobian * quadrature weight at each integration point.   
   */
  const std::vector<Real>& JxW = fe->get_JxW();

  /**
   * The element shape functions evaluated at the quadrature points.
   */
  const std::vector<std::vector<Real> >& phi = fe->get_phi();

  /**
   * The element shape function gradients evaluated at the quadrature
   * points.
   */
  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();
  
  /**
   * A reference to the \p DofMap object for this system.  The \p DofMap
   * object handles the index translation from node and element numbers
   * to degree of freedom numbers.  We will talk more about the \p DofMap
   * in future examples.
   */
  const DofMap& dof_map = system.get_dof_map();

  /**
   * Define data structures to contain the element matrix
   * and right-hand-side vector contribution.  Following
   * basic finite element terminology we will denote these
   * "Ke" and "Fe".
   */
  DenseMatrix<Real> Ke;
  DenseVector<Real> Fe;

  /**
   * This vector will hold the degree of freedom indices for
   * the element.  These define where in the global system
   * the element degrees of freedom get mapped.
   */
  std::vector<unsigned int> dof_indices;



  /**
   *--------------------------------------------------------------------
   * Here we extract the velocity & parameters that we put in the
   * EquationSystems object.
   */
  const RealVectorValue velocity =
    es.data_map.get_data<RealVectorValue> ("velocity");

  const Real dt = es.parameter   ("dt");
  const Real time = es.parameter ("time");



  
  /**
   *--------------------------------------------------------------------
   * Now we will loop over all the elements in the mesh that
   * live on the local processor. We will compute the element
   * matrix and right-hand-side contribution.  Since the mesh
   * will be refined we want to only consider the ACTIVE elements,
   * hence we use a variant of the \p active_elem_iterator.
   */

  const_active_local_elem_iterator           el (mesh.elements_begin());
  const const_active_local_elem_iterator end_el (mesh.elements_end());
  
  for ( ; el != end_el; ++el)
    {    
      /**
       * Store a pointer to the element we are currently
       * working on.  This allows for nicer syntax later.
       */
      const Elem* elem = *el;

      /**
       * Get the degree of freedom indices for the
       * current element.  These define where in the global
       * matrix and right-hand-side this element will
       * contribute to.
       */
      dof_map.dof_indices (elem, dof_indices);

      /**
       * Compute the element-specific data for the current
       * element.  This involves computing the location of the
       * quadrature points (q_point) and the shape functions
       * (phi, dphi) for the current element.
       */
      fe->reinit (elem);

      /**
       * Zero the element matrix and right-hand side before
       * summing them.  We use the resize member here because
       * the number of degrees of freedom might have changed from
       * the last element.  Note that this will be the case if the
       * element type is different (i.e. the last element was a
       * triangle, now we are on a quadrilateral).
       */
      Ke.resize (dof_indices.size(),
		 dof_indices.size());

      Fe.resize (dof_indices.size());

      
      
      /**
       * Now we will build the element matrix and right-hand-side.
       * Constructing the RHS requires the solution and its
       * gradient from the previous timestep.  This myst be
       * calculated at each quadrature point by summing the
       * solution degree-of-freedom values by the appropriate
       * weight functions.
       */
      
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
	{

	  /**
	   * Values to hold the old solution & its gradient.
	   */
	  Real         u_old = 0.;
	  RealGradient grad_u_old;

	  /**
	   * Compute the old solution & its gradient.
	   */
	  for (unsigned int l=0; l<phi.size(); l++)
	    {
	      u_old      += phi[l][qp]*system.old_solution  (dof_indices[l]);

	      /**
	       * This will work,
	       * grad_u_old += dphi[l][qp]*system.old_solution (dof_indices[l]);
	       * but we can do it without creating a temporary like this:
	       */	      
	      grad_u_old.add_scaled (dphi[l][qp],system.old_solution (dof_indices[l]));
	    }

	  /**
	   * Now compute the element matrix and RHS contributions.
	   */
	  for (unsigned int i=0; i<phi.size(); i++)
	    {
	      // The RHS contribution
	      Fe(i) += JxW[qp]*(
				u_old*phi[i][qp] + // Mass matrix term

				-.5*dt*((velocity*grad_u_old)*phi[i][qp] + // Convection
					0.01*(grad_u_old*dphi[i][qp]))     // Diffusion
				);
	      
	      for (unsigned int j=0; j<phi.size(); j++)
		{
		  // The matrix contribution
		  Ke(i,j) += JxW[qp]*(
				      phi[i][qp]*phi[j][qp] + // Mass-matrix

				      .5*dt*((velocity*dphi[j][qp])*phi[i][qp] + // Convection
					     0.01*(dphi[i][qp]*dphi[j][qp]))      // Diffusion
				      );
		  
		} // end of the matrix summation j-loop
	    } // end of the matrix & rhs summation i-loop	  
	} // end of the quadrature point qp-loop


      
      


      
      /**
       *----------------------------------------------------------------
       * At this point the interior element integration has
       * been completed.  However, we have not yet addressed
       * boundary conditions.  For this example we will only
       * consider simple Dirichlet boundary conditions imposed
       * via the penalty method. The penalty method used here
       * is equivalent to lumping the matrix resulting from the
       * L2 projection penalty approach introduced in example 3.
       */      
      {
	
	/**
	 * The following loops over the sides of the element.
	 * If the element has no neighbor on a side then that
	 * side MUST live on a boundary of the domain.
	 */
	for (unsigned int s=0; s<elem->n_sides(); s++)
	  if (elem->neighbor(s) == NULL)
	    {
	      AutoPtr<Elem> side (elem->build_side(s));
	      
	      /**
	       * Loop over the nodes on the side.
	       */
	      for (unsigned int ns=0; ns<side->n_nodes(); ns++)
		{
		  /**
		   * The location on the boundary of the current
		   * node.
		   */
		  const Real xf = side->point(ns)(0);
		  const Real yf = side->point(ns)(1);
		  
		  /**
		   * The penalty value.  \f$ \frac{1}{\epsilon \f$
		   */
		  const Real penalty = 1.e10;
		  
		  /**
		   * The boundary value.
		   */
		  const Real value = exact_solution(xf, yf, time);

		  /**
		   * Find the node on the element matching this node on
		   * the side.  That defined where in the element matrix
		   * the boundary condition will be applied.
		   */
		  for (unsigned int n=0; n<elem->n_nodes(); n++)
		    if (elem->node(n) == side->node(ns))
		      {
			/**
			 * Matrix contribution.
			 */
			Ke(n,n) += penalty;
		  		  
			/**
			 * Right-hand-side contribution.
			 */
			Fe(n) += penalty*value;
		      }
		  
		} // end face node loop	  
	    } // end if (elem->neighbor(side) == NULL)
	
      } // end boundary condition section	  


      

      /**
       * We have now built the element matrix and RHS vector in terms
       * of the element degrees of freedom.  However, it is possible
       * that some of the element DOFs are constrained to enforce
       * solution continuity, i.e. they are not really "free".  We need
       * to constrain those DOFs in terms of non-constrained DOFs to
       * ensure a continuous solution.  The
       * \p DofMap::constrain_element_matrix_and_vector() method does
       * just that.
       */
      dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);
      
      /**
       *----------------------------------------------------------------
       * The element matrix and right-hand-side are now built
       * for this element.  Add them to the global matrix and
       * right-hand-side vector.  The \p PetscMatrix::add_matrix()
       * and \p PetscVector::add_vector() members do this for us.
       */      
      es("Convection-Diffusion").matrix->add_matrix (Ke, dof_indices);
      es("Convection-Diffusion").rhs->add_vector    (Fe, dof_indices);
      
    } // end of element loop



  
  /**
   * That's it.
   */
  return;

#endif
}
