// $Id: ex8.C,v 1.1 2003-04-09 14:08:39 spetersen Exp $
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
#include "newmark_system.h"
#include "equation_systems.h"

/**
 * Define the Finite Element object.
 */
#include "fe.h"

/**
 * Define Gauss quadrature rules.
 */
#include "quadrature_gauss.h"

/**
 * Define useful datatypes for finite element
 * matrix and vector components.
 */
#include "dense_matrix.h"
#include "dense_vector.h"

/**
 * Define matrix and vector data types for the global 
 * equation system.  These are base classes,
 * from which specific implementations, like
 * the PETSc or LASPACK implementations, are derived.
 */
#include "sparse_matrix.h"
#include "numeric_vector.h"

/**
 * Define the DofMap, which handles degree of freedom
 * indexing.
 */
#include "dof_map.h"


/**
 * \mainpage Example 8
 *
 * \section Introduction
 *
 * This is the eighth example program. It builds on
 * the previous example programs. It introduces the
 * NewmarkSystem class. In this example the wave equation
 * is solved using the time integration scheme provided
 * by the NewmarkSystem class.
 *
 * This example comes with a cylindrical mesh given in the
 * universal file pipe-mesh.unv
 * The mesh contains HEX8 and PRISM6 elements.
 */



/**
 * Function prototype.  This is the function that will assemble
 * the linear system for our problem, governed by the linear
 * wave equation
 */
void assemble_wave(EquationSystems<NewmarkSystem>& es,
			const std::string& system_name);



/**
 * Function Prototype. This function will be used to apply the
 * initial conditions.
 */
void apply_initial(EquationSystems<NewmarkSystem>& es,
		   const std::string& system_name);



/**
 * Function Prototype. This function imposes
 * Dirichlet Boundary conditions via the penalty
 * method after the system is assembled.
 */
void fill_dirichlet_bc(EquationSystems<NewmarkSystem>& es,
		       const std::string& system_name,
		       bool do_for_matrix=false);



/**
 *--------------------------------------------------------------------
 * Initialize some (global) handy constants and variables
 *
 * The step size for the time.
 */
Real delta_t = .0000625;

/**
 * The time.
 */
Real t_time;

/**
 * The number of time steps.
 */
unsigned int n_time_steps = 300;

/**
 * The node that should be monitored.
 */
unsigned int result_node = 274;



/**
 * Define the fluid properties (the density rho 
 * and the speed of sound speed).
 */
const Real speed = 1000.;
const Real rho   = 1000.;


/**
 * Define pi.
 */
const Real libmesh_pi = acos(-1.);




int main (int argc, char** argv)
{
  /**
   * Initialize Petsc, like in example 2.
   */
  libMesh::init (argc, argv);

  /**
   * This example is not designed for complex numbers.
   */
#ifdef USE_COMPLEX_NUMBERS

  std::cerr << "ERROR: This example is not intended for " << std::endl
	    << " use with complex numbers." << std::endl;
  here();

  return 0;

#endif


  /**
   * Braces are used to force object scope.
   */   
  {
    /**
     * Check for proper usage.
     */
    if (argc < 2)
      {
	std::cerr << "Usage: " << argv[0] << " [meshfile]"
		  << std::endl;
	
	/**
	 * This handy function will print the file name, line number,
	 * and then abort.  Currrently the library does not use C++
	 * exception handling.
	 */
	error();
      }
    
    /**
     * Tell the user what we are doing.
     */
    else 
      {
	std::cout << "Running " << argv[0];
	
	for (int i=1; i<argc; i++)
	  std::cout << " " << argv[i];
	
	std::cout << std::endl << std::endl;

      };

    /**
     * Get the name of the mesh file
     * from the command line.
     */
    std::string mesh_file = argv[1];
    std::cout << "Mesh file is: " << mesh_file << std::endl;

    /**
     * For now, restrict to dim=3, though this
     * may easily be changed, see example 4
     */
    const unsigned int dim = 3;

    /**
     * Create a dim-dimensional mesh.
     */
    Mesh mesh (dim);

    /**
     * Read the meshfile specified in the command line or
     * use the internal mesh generator to create a uniform
     * grid on square or cube.
     */
    mesh.read(mesh_file);
    // mesh.build_cube (10, 10, 40,
    //		     -1., 1.,
    //		     -1., 1.,
    //	             0., 4.,
    //	             HEX8);

    /**
     * Print information about the mesh to the screen.
     */
    mesh.print_info();
    
    /**
     * Create an equation systems object.
     */
    EquationSystems<NewmarkSystem> equation_systems (mesh);
    
    /**
     * Declare the system and its variables.
     */
    {
      /**
       * Creates a system named "Wave"
       */
      equation_systems.add_system("Wave");

      /**
       * Use a handy reference to this system
       */
      NewmarkSystem & t_system = equation_systems("Wave");
      
      /**
       * Adds the variable "p" to "Wave".  "p"
       * will be approximated using first-order approximation.
       */
      t_system.add_variable("p", FIRST);

      /**
       * Give the system a pointer to the matrix assembly
       * function and the initial condition function defined
       * below.
       */
      t_system.attach_assemble_function (assemble_wave);
      t_system.attach_init_cond_function (apply_initial);
  

      /**
       * Set the Newmark parameters and compute integration
       * constants. Here we simply use the default values
       * of alpha=.25  and delta=.5.
       */
      t_system.set_newmark_parameters(delta_t);
     
      /**
       * Initialize the data structures for the equation system.
       */
      equation_systems.init();

      /**
       * Prints information about the system to the screen.
       */
      equation_systems.print_info();
    };

    /**
     * A file to store the results at certain nodes.
     */
    std::ofstream res_out("result_file.dat");

    /**
     * get the dof_numbers for the nodes that
     * should be monitored.
     */
    unsigned int res_node_no = result_node;
    const Node& res_node = mesh.node(res_node_no-1);
    unsigned int dof_no = res_node.dof_number(0,0,0);


    /**
     * Assemble the time independent system matrices and rhs.
     * This function will also compute the effective system matrix
     * K~=K+a_0*M+a_1*C and apply user specified initial
     * conditions. 
     */
    equation_systems("Wave").assemble();


    /**
     * Now solve for each time step.
     */
    // start with t_time = 0 and write a short header to the
    // nodal result file
    t_time = 0.;
    res_out << "# pressure at node " << res_node_no << "\n"
	    << "# time\tpressure\n"
	    << t_time << "\t" << 0 << std::endl;


    for (unsigned int time_step=0; time_step<n_time_steps; time_step++)
      {

	/**
	 * Update the time.
	 */
	t_time += delta_t;

	/**
	 * Update the rhs.
	 */
	equation_systems("Wave").update_rhs();

	/**
	 * Impose essential boundary conditions.
	 * Not that since the matrix is only assembled once,
	 * the panalty parameter should be added to the matrix
	 * only in the first time step. The applied
	 * boundary conditions may be time-dependent and hence
	 * the rhs vector is considered in each time step. 
	 */
	if (time_step == 0)
	  fill_dirichlet_bc(equation_systems, "Wave", true);
	else
	  fill_dirichlet_bc(equation_systems, "Wave");

	/**
	 * Solve the system "Wave".
	 */
	equation_systems("Wave").solve();

	/**
	 * After solving the system, write the solution
	 * to a GMV-formatted plot file.
	 * Do only for a few time steps.
	 */
	if (time_step == 30 || time_step == 60 ||
	    time_step == 90 || time_step == 120 )
	  {
	    char buf[14];
	    sprintf (buf, "out.%03d.gmv", time_step);
	    mesh.write_gmv (buf, equation_systems);
	  }

	/**
	 * Update the p, v and a.
	 */
	equation_systems("Wave").update_u_v_a();

	/**
	 * Write nodal results to file. The results can then
	 * be viewed with e.g. gnuplot (run gnuplot and type
	 * 'plot "result_file.dat" with lines' in the command line)
	 */
	res_out << t_time << "\t"
		<< equation_systems("Wave").get_vector("displacement")(dof_no)
		<< std::endl;

      }
  };

  
  /**
   * All done.  
   */
  return libMesh::close ();
};




void assemble_wave(EquationSystems<NewmarkSystem>& es,
			const std::string& system_name)
{
  
  /**
   * It is a good idea to make sure we are assembling
   * the proper system.
   */
  assert (system_name == "Wave");

  /**
   * Get a constant reference to the mesh object.
   */
  const Mesh& mesh = es.get_mesh();

  /**
   * The dimension that we are running
   */
  const unsigned int dim = mesh.mesh_dimension();

  /**
   * Get a reference to our system, as before
   */
  NewmarkSystem & t_system = es (system_name);

  /**
   * Get a constant reference to the Finite Element type
   * for the first (and only) variable in the system.
   */
  FEType fe_type = es("Wave").get_dof_map().variable_type(0);

  /**
   * In here, we will add the element matrices to the
   * @e additional matrices "stiffness_mass", "damping",
   * and the additional vector "rhs", not to the members 
   * "matrix" and "rhs".  Therefore, get writable
   * references to them
   */
  SparseMatrix<Number>&   stiffness      = t_system.get_matrix("stiffness");
  SparseMatrix<Number>&   damping        = t_system.get_matrix("damping");
  SparseMatrix<Number>&   mass           = t_system.get_matrix("mass");

  NumericVector<Number>&  force          = t_system.get_vector("force");

  /**
   * Some solver packages (PETSc) are especially picky about
   * allocating sparsity structure and truly assigning values
   * to this structure.  Namely, matrix additions, as performed
   * later, exhibit acceptable performance only for identical
   * sparsity structures.  Therefore, explicitly zero the
   * values in the collective matrix, so that matrix additions
   * encounter identical sparsity structures.
   */
  SparseMatrix<Number>&  matrix           = *t_system.matrix;


  /**
   * Build a Finite Element object of the specified type.  Since the
   * \p FEBase::build() member dynamically creates memory we will
   * store the object as an \p AutoPtr<FEBase>.  This can be thought
   * of as a pointer that will clean up after itself.
   */
  AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));
  
  /**
   * A 2nd order Gauss quadrature rule for numerical integration.
   */
  QGauss qrule (dim, SECOND);

  /**
   * Tell the finite element object to use our quadrature rule.
   */
  fe->attach_quadrature_rule (&qrule);

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
  const std::vector<std::vector<Point> >& dphi = fe->get_dphi();

  /**
   * A reference to the \p DofMap object for this system.  The \p DofMap
   * object handles the index translation from node and element numbers
   * to degree of freedom numbers.  We will talk more about the \p DofMap
   * in future examples.
   */
  const DofMap& dof_map = t_system.get_dof_map();

  /**
   * The mass, damping and stiffness matrices, however, are inherently
   * real.  And since \p DenseMatrix<> offers a method
   * to add a real matrix to a complex matrix, we can safely
   * define element stiffness Ke and mass matrix Me as real. 
   */
  DenseMatrix<Number>   Ae, zero_matrix;
  DenseMatrix<Number>   Ke, Ce, Me;
  DenseVector<Number>   Fe;


  /**
   * This vector will hold the degree of freedom indices for
   * the element.  These define where in the global system
   * the element degrees of freedom get mapped.
   */
  std::vector<unsigned int> dof_indices;


  /**
   *--------------------------------------------------------------------
   * Now we will loop over all the elements in the mesh.
   * We will compute the element matrix and right-hand-side
   * contribution.
   */

  const_elem_iterator           el (mesh.elements_begin());
  const const_elem_iterator end_el (mesh.elements_end());
  
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
       *
       * The \p DenseMatrix::resize() member will automatically
       * zero out the matrix.  Since we are using a \p std::vector
       * for the right-hand-side we will use the \p std::fill algorithm
       * to zero out Fe.
       */
      Ae.resize          (dof_indices.size(),
			  dof_indices.size());
      zero_matrix.resize (dof_indices.size(),
			  dof_indices.size());

      Ke.resize (dof_indices.size(),
		 dof_indices.size());

      Ce.resize (dof_indices.size(),
		 dof_indices.size());

      Me.resize (dof_indices.size(),
		 dof_indices.size());

      Fe.resize (dof_indices.size());


      /**
       *----------------------------------------------------------------
       * Now loop over the quadrature points.  This handles
       * the numeric integration.
       */
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
	{

	  /**
	   * Now we will build the element matrix.  This involves
	   * a double loop to integrate the test funcions (i) against
	   * the trial functions (j).
	   */
	  for (unsigned int i=0; i<phi.size(); i++)
	    for (unsigned int j=0; j<phi.size(); j++)
	      {
		Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
		Me(i,j) += JxW[qp]*phi[i][qp]*phi[j][qp]
		           *1./(speed*speed);
	      }; // end of the matrix summation loop
	  
	}; // end of quadrature point loop


      /**       
       *----------------------------------------------------------------
       * Now compute the contribution to the element matrix and the
       * right-hand-side vector if the current element lies on the
       * boundary. 
       */
      {
	/**
	 * The following loops over the sides of the element.
	 * If the element has no neighbor on a side then that
	 * side MUST live on a boundary of the domain.
	 *
	 * In this example no natural boundary conditions will
	 * be considered. The code is left here so it can easily
	 * be extended.
	 * 
	 */
	// don't do for any side
	for (unsigned int side=0; side<elem->n_sides(); side++)
	  if (!true)
	    // if (elem->neighbor(side) == NULL)
	    {
	      /**
	       * Declare a special finite element object for
	       * boundary integration.
	       */
	      AutoPtr<FEBase> fe_face (FEBase::build(dim, fe_type));
	      
	      /**
	       * Boundary integration requires one quadraure rule,
	       * with dimensionality one less than the dimensionality
	       * of the element.
	       */
	      QGauss qface(dim-1, SECOND);
	      
	      /**
	       * Tell the finte element object to use our
	       * quadrature rule.
	       */
	      fe_face->attach_quadrature_rule (&qface);
	      
	      /**
	       * The value of the shape functions at the quadrature
	       * points.
	       */
	      const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();
	      
	      /**
	       * The Jacobian * Quadrature Weight at the quadrature
	       * points on the face.
	       */
	      const std::vector<Real>& JxW_face = fe_face->get_JxW();
	      
	      /**
	       * Compute the shape function values on the element
	       * face.
	       */
	      fe_face->reinit(elem, side);

	      /**
	       * Here we consider a normal acceleration acc_n=1 applied to
	       * the whole boundary of our mesh.
	       */ 
	      const Real acc_n_value = 1.0;
	      
	      /**
	       * Loop over the face quagrature points for integration.
	       */
	      for (unsigned int qp=0; qp<qface.n_points(); qp++)
		{

		  /**
		   * Right-hand-side contribution due to prescribed
		   * normal acceleration.
		   */
		  for (unsigned int i=0; i<phi_face.size(); i++)
		    {
		      Fe(i) += acc_n_value*rho
			*phi_face[i][qp]*JxW_face[qp];
		    };


		}; // end face quadrature point loop	  
	    }; // end if (elem->neighbor(side) == NULL)

	/**
	 * In this example the Dirichlet boundary conditions will be 
	 * imposed via panalty method after the
	 * system is assembled.
	 */

      }; // end boundary condition section	  



      /**
       * Finally, simply add the contributions to the additional
       * matrices and vector.
       */
      stiffness.add_matrix   (Ke, dof_indices);
      damping.add_matrix     (Ce, dof_indices);
      mass.add_matrix        (Me, dof_indices);

      force.add_vector       (Fe, dof_indices);
    
      /**
       * For the overall matrix, explicitly zero the entries where
       * we added values in the other ones, so that we have 
       * identical sparsity footprints.
       */
      matrix.add_matrix(zero_matrix, dof_indices);
    
    }; // end of element loop
  
  /**
   * All done!
   */
  return;

};

void apply_initial(EquationSystems<NewmarkSystem>& es,
		   const std::string& system_name)
{
  /**
   * Get a reference to our system, as before
   */
  NewmarkSystem & t_system = es (system_name);

  /**
   * Numeric vectors for the pressure, velocity and acceleration
   * values.
   */
  NumericVector<Number>&  pres_vec       = t_system.get_vector("displacement");
  NumericVector<Number>&  vel_vec        = t_system.get_vector("velocity");
  NumericVector<Number>&  acc_vec        = t_system.get_vector("acceleration");

  /**
   * Assume our fluid to be at rest, which would
   * also be the default conditions in class NewmarkSystem.
   */
  pres_vec.zero();
  vel_vec.zero();
  acc_vec.zero();

  // the nodal acceleration values for t=t_0 could be computed
  // in the matrix assembly.
  // acc_vec.add (1., t_system.get_vector("force"));

}




void fill_dirichlet_bc(EquationSystems<NewmarkSystem>& es,
		       const std::string& system_name,
		       bool do_for_matrix)
{
  /**
   * It is a good idea to make sure we are assembling
   * the proper system.
   */
  assert (system_name == "Wave");

  /**
   * Get a reference to our system, as before.
   */
  NewmarkSystem & t_system = es (system_name);

  /**
   * Get writable references to the overall matrix and vector.
   */
  SparseMatrix<Number>&  matrix                = *t_system.matrix;
  NumericVector<Number>& rhs                   = *t_system.rhs;

  /**
   * Get a constant reference to the mesh object.
   */
  const Mesh& mesh = es.get_mesh();

  /**
   * Number of nodes in the mesh.
   */
  unsigned int n_nodes = mesh.n_nodes();

  for (unsigned int n_cnt=0; n_cnt<n_nodes; n_cnt++)
    {
      /**
       * Get a reference to the current node.
       */
      const Node& curr_node = mesh.node(n_cnt);

      /**
       * Check if Dirichlet BCs should be applied to this node.
       * Here a pressure value is applied if the z-coord.
       * is equal to 4.
       */
      const Real coo_tol = 1.e-6;
      const Real z_coo = 4.;

      if (curr_node(2)<z_coo+coo_tol &&
	  curr_node(2)>z_coo-coo_tol)
	{
	  /**
	   * The global number of the respective degree of freedom.
	   */ 
	  unsigned int dn = curr_node.dof_number(0,0,0);

	  /**
	   * The penalty parameter.
	   */
	  const Real penalty = 1.e10;

	  /**
	   * The boundary value that will be applied.
	   */
	  // const Real p_value = 1.;

	  /**
	   * Compute distance from point 0,0,4.
	   */
	  // const Real r_dist = pow(curr_node(0)*curr_node(0)+
          //                        curr_node(1)*curr_node(1),.5);

	  /**
	   * Here we apply sinusoidal pressure values for 0<t<0.002
	   * at one end of the pipe-mesh.
	   */
	  Real p_value;
	  if (t_time < .002 )
	    p_value = sin(2*libmesh_pi*t_time/.002);
	  else
	    p_value = .0;

	  // if (r_dist < .2+coo_tol && t_time < .002 )
	  //   p_value = sin(2*libmesh_pi*t_time/.002);
	  // else
	  //   p_value = .0;

	  /**
	   * Now add the contributions to the matrix and the rhs.
	   */
	  rhs.add(dn, p_value*penalty);

	  if (do_for_matrix)
	    matrix.add(dn, dn, penalty);

	}
    }

}






















