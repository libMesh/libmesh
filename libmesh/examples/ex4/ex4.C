// $Id: ex4.C,v 1.6 2003-02-10 03:55:50 benkirk Exp $

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
#include "mesh_init.h"
#include "mesh.h"
#include "system_data.h"
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
 * Define the DofMap, which handles degree of freedom
 * indexing.
 */
#include "dof_map.h"






/**
 * \mainpage Example 4
 *
 * \section Introduction
 *
 * This is the fourth example program.  It builds on
 * the third example program by showing how to formulate
 * the code in a dimension-independent way.  Very minor
 * changes to the example will allow the problem to be
 * solved in two or three dimensions.
 */


/**
 * Function prototype.  This is the function that will assemble
 * the linear system for our Poisson problem.  Note that the
 * function will take the \p EquationSystems object and the
 * name of the system we are assembling as input.  From the
 * \p EquationSystems object we have acess to the \p Mesh and
 * other objects we might need.
 */
void assemble_poisson(EquationSystems& es,
                      const std::string& system_name);



/**
 * Exact solution function prototype.
 */
Real exact_solution (const Real x,
		     const Real y,
		     const Real z = 0.);





int main (int argc, char** argv)
{
  /**
   * Initialize Petsc, like in example 2.
   */
  libMesh::init (argc, argv);
  
  /**
   * This example is designed for real numbers only.
   */
#ifdef USE_COMPLEX_NUMBERS

  std::cerr << "ERROR: This example is not intended for " << std::endl
	    << " use with complex numbers." << std::endl;
  error();

#endif


  /**
   * Braces are used to force object scope, like in example 2
   */   
  {
    /**
     * Check for proper usage.
     */
    if (argc < 3)
      {
	std::cerr << "Usage: " << argv[0] << " -d 2"
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
     * Get the dimensionality of the mesh from argv[2]
     */
    const unsigned int dim = atoi(argv[2]);     
    
    /**
     * Create a mesh with user-defined dimension.
     */
    Mesh mesh (dim);
    
    /**
     * Use the internal mesh generator to create a uniform
     * grid on the square [-1,1]^D.  We instruct the mesh generator
     * to build a mesh of 5x5 \p Quad9 elements in 2D, or \p Hex27
     * elements in 3D.  Building these higher-order elements allows
     * us to use higher-order approximation, as in example 3.
     */
    mesh.build_cube (5, 5, 5,
		     -1., 1.,
		     -1., 1.,
		     -1., 1.,
		     (dim == 2) ? QUAD9 : HEX27);

    /**
     * Let the elements find their neighbors.
     */
    mesh.find_neighbors();
    
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
       * Creates a system named "Poisson"
       */
      equation_systems.add_system("Poisson");
      
      /**
       * Adds the variable "u" to "Poisson".  "u"
       * will be approximated using second-order approximation.
       */
      equation_systems("Poisson").add_variable("u", SECOND);

      /**
       * Give the system a pointer to the matrix assembly
       * function.
       */
      equation_systems("Poisson").attach_assemble_function (assemble_poisson);
      
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
     * Solve the system "Poisson", just like example 2.
     */
    equation_systems("Poisson").solve();


    /**
     * After solving the system write the solution
     * to a GMV-formatted plot file.
     */
    mesh.write_gmv ((dim == 3) ? "out_3.gmv" : "out_2.gmv",
		    equation_systems);
  };


  libMesh::close ();

  
  /**
   * All done.  
   */
  return 0;
};




void assemble_poisson(EquationSystems& es,
                      const std::string& system_name)
{
  /**
   * It is a good idea to make sure we are assembling
   * the proper system.
   */
  assert (system_name == "Poisson");

  /**
   * Get a constant reference to the mesh object.
   */
  const Mesh& mesh = es.get_mesh();

  /**
   * The dimension that we are running
   */
  const unsigned int dim = mesh.mesh_dimension();

  /**
   * Get a constant reference to the Finite Element type
   * for the first (and only) variable in the system.
   */
  FEType fe_type = es("Poisson").dof_map.component_type(0);

  /**
   * Build a Finite Element object of the specified type.  Since the
   * \p FEBase::build() member dynamically creates memory we will
   * store the object as an \p AutoPtr<FEBase>.  This can be thought
   * of as a pointer that will clean up after itself.
   */
  AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));
  
  /**
   * A 5th order Gauss quadrature rule for numerical integration.
   */
  QGauss qrule (dim, FIFTH);

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
   * The physical XY locations of the quadrature points on the element.
   * These might be useful for evaluating spatially varying material
   * properties at the quadrature points.
   */
  const std::vector<Point>& q_point = fe->get_xyz();

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
  const DofMap& dof_map = es("Poisson").get_dof_map();

  /**
   * Define data structures to contain the element matrix
   * and right-hand-side vector contribution.  Following
   * basic finite element terminology we will denote these
   * "Ke" and "Fe".  Use the complex versions, so that
   * this example compiles successfully, and the error 
   * message in \p main() can catch this irregularity.
   */
  RealDenseMatrix   Ke;
  std::vector<Real> Fe;

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
  for (unsigned int e=0; e<mesh.n_elem(); e++)
    {
      /**
       * Store a pointer to the element we are currently
       * working on.  This allows for nicer syntax later.
       */
      const Elem* elem = mesh.elem(e);

      /**
       * Get the degree of freedom indices for the
       * current element.  These define where in the global
       * matrix and right-hand-side this element will
       * contribute to.
       */
      dof_map.dof_indices (e, dof_indices);

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
      Ke.resize (dof_indices.size(),
		 dof_indices.size());

      Fe.resize (dof_indices.size());

      std::fill (Fe.begin(), Fe.end(), 0.);



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
	      }; // end of the matrix summation loop


	  /**
	   * Now we build the element right-hand-side contribution.
	   * This involves a single loop in which we integrate the
	   * "forcing function" in the PDE against the test functions.
	   */
	  for (unsigned int i=0; i<phi.size(); i++)
	    {
	      const Real x = q_point[qp](0);
	      const Real y = q_point[qp](1);
	      const Real z = q_point[qp](2);
	      const Real eps = 1.e-3;

	      /**
	       * fxy is the forcing function for the Poisson equation.
	       * In this case we set fxy to be a finite difference
	       * Laplacian approximation to the (known) exact solution.
	       *
	       * Note that in 2D the Laplacian of u = u_xx + u_yy,
	       * but in 3D Laplacian of u = u_xx + u_yy + u_zz
	       */
	      const Real uxx = (exact_solution(x-eps,y,z) +
				exact_solution(x+eps,y,z) +
				-2.*exact_solution(x,y,z))/eps/eps;
	      
	      const Real uyy = (exact_solution(x,y-eps,z) +
				exact_solution(x,y+eps,z) +
				-2.*exact_solution(x,y,z))/eps/eps;
	      
	      const Real uzz = (exact_solution(x,y,z-eps) +
				exact_solution(x,y,z+eps) +
				-2.*exact_solution(x,y,z))/eps/eps;

	      const Real fxy = - (uxx + uyy + ((dim==2) ? 0. : uzz));
	      
	      Fe[i] += JxW[qp]*fxy*phi[i][qp];
	    }; // end of the RHS summation loop
	  
	}; // end of quadrature point loop




      
      /**
       *----------------------------------------------------------------
       * At this point the interior element integration has
       * been completed.  However, we have not yet addressed
       * boundary conditions.  For this example we will only
       * consider simple Dirichlet boundary conditions imposed
       * via the penalty method. This is discussed at length in
       * example 3.
       */
      {
	/**
	 * The following loops over the sides of the element.
	 * If the element has no neighbor on a side then that
	 * side MUST live on a boundary of the domain.
	 */
	for (unsigned int side=0; side<elem->n_sides(); side++)
	  if (elem->neighbor(side) == NULL)
	    {
	      /**
	       * Declare a special finite element object for
	       * boundary integration.
	       */
	      AutoPtr<FEBase> fe_face (FEBase::build(dim, fe_type));
	      
	      /**
	       * Boundary integration requires TWO quadraure rules:
	       * one the dimensionality of the element and another
	       * one less than the dimensionality of the element.
	       */
	      QGauss qface0(dim,   FIFTH);
	      QGauss qface1(dim-1, FIFTH);
	      
	      /**
	       * Tell the finte element object to use our
	       * quadrature rule.
	       */
	      fe_face->attach_quadrature_rule (&qface0);
	      
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
	       * The XYZ locations (in physical space) of the
	       * quadrature points on the face.  This is where
	       * we will interpolate the boundary value function.
	       */
	      const std::vector<Point >& qface_point = fe_face->get_xyz();
	      
	      /**
	       * Compute the shape function values on the element
	       * face.
	       */
	      fe_face->reinit(&qface1, elem, side);
	      
	      /**
	       * Loop over the face quagrature points for integration.
	       */
	      for (unsigned int qp=0; qp<qface0.n_points(); qp++)
		{
		  /**
		   * The location on the boundary of the current
		   * face quadrature point.
		   */
		  const Real xf = qface_point[qp](0);
		  const Real yf = qface_point[qp](1);
		  const Real zf = qface_point[qp](2);
		  
		  /**
		   * The penalty value.  \f$ \frac{1}{\epsilon \f$
		   * in the discussion above.
		   */
		  const Real penalty = 1.e10;
		  
		  /**
		   * The boundary value.
		   */
		  const Real value = exact_solution(xf, yf, zf);
		  
		  /**
		   * Matrix contribution of the L2 projection. 
		   */
		  for (unsigned int i=0; i<phi_face.size(); i++)
		    for (unsigned int j=0; j<phi_face.size(); j++)
		      {
			Ke(i,j) += JxW_face[qp]*penalty*phi_face[i][qp]*phi_face[j][qp];
		      };
		  
		  /**
		   * Right-hand-side contribution of the L2
		   * projection.
		   */
		  for (unsigned int i=0; i<phi_face.size(); i++)
		    {
		      Fe[i] += JxW_face[qp]*penalty*value*phi_face[i][qp];
		    };
		  
		}; // end face quadrature point loop	  
	    }; // end if (elem->neighbor(side) == NULL)
      }; // end boundary condition section	  


      

      
      /**
       *----------------------------------------------------------------
       * The element matrix and right-hand-side are now built
       * for this element.  Add them to the global matrix and
       * right-hand-side vector.  The \p PetscMatrix::add_matrix()
       * and \p PetscVector::add_vector() members do this for us.
       * 
       * The preprocessor test for complex numbers is explained
       * in example 5.
       */
#ifndef USE_COMPLEX_NUMBERS

      es("Poisson").matrix.add_matrix (Ke, dof_indices);
      es("Poisson").rhs.add_vector    (Fe, dof_indices);
#endif
      
    }; // end of element loop


  
  /**
   * All done!
   */
  return;
};
