// $Id: ex7.C,v 1.2 2003-02-07 16:19:11 spetersen Exp $

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
 * \mainpage Example 7
 *
 * \section Introduction
 *
 * Warning: This example is still under cunstruction.
 *
 * This is the sixth example program.  It builds on
 * the privious example programs. This example now sows how
 * to deal with complex numbers by solving the Helmholtz
 * equation grad(p)*grad(p)+(omega/c)^2*p=0.
 * For this example the library has to be compiled with
 * complex numbers enabled.
 */


/**
 * Function prototype.  This is the function that will assemble
 * the linear system for our Helmholz problem.
 */
void assemble_helmholtz(EquationSystems& es,
			const std::string& system_name);
/**
 * The frequency for which we are solving the Helmholtz
 * equation.
 */
Real frequency;


int main (int argc, char** argv)
{
  /**
   * Initialize Petsc, like in example 2.
   */
#ifdef HAVE_PETSC
  
  const bool have_petsc = true;
  PetscInitialize (&argc, &argv, NULL, NULL);
  
#else
  
  const bool have_petsc = false;
  
#endif

  /**
   * This example is designed for complex numbers.
   */
#ifndef USE_COMPLEX_NUMBERS

  std::cerr << "ERROR: This example is intended for " << std::endl
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
    if (argc != 5)
      {
	std::cerr << "Usage: " << argv[0] << " -d 2"
		  << " -f [frequency]"
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
     * Get the dimensionality of the mesh from argv[2].
     */
    const unsigned int dim = atoi(argv[2]);     

    /**
     * Get the frequency from argv[4].
     */
    frequency = atoi(argv[4]);

    /**
     * Define the fluid properties. Here (for simplicity) 
     * we define the density rho = 1 and the speed of sound
     * c = 1.
     */
    // const Real   c = 1.0;
    // const Real rho = 1.0;

    /**
     * Create a 2D mesh.
     */
    Mesh mesh (dim);
    
    /**
     * Use the internal mesh generator to create a uniform
     * grid on the square [-1,1]^D.  We instruct the mesh generator
     * to build a mesh of 5x5 \p Quad4 elements in 2D, or \p Hex8
     * elements in 3D.
     */
    mesh.build_cube (5, 5, 5,
		     -1., 1.,
		     -1., 1.,
		     -1., 1.,
		     (dim == 2) ? QUAD4 : HEX8);

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
    EquationSystems equation_systems (mesh, have_petsc);
    
    /**
     * Declare the system and its variables.
     */
    {
      /**
       * Creates a system named "Helmholtz"
       */
      equation_systems.add_system("Helmholtz");
      
      /**
       * Adds the variable "p" to "Helmholtz".  "p"
       * will be approximated using first-order approximation.
       */
      equation_systems("Helmholtz").add_variable("p", FIRST);

      /**
       * Give the system a pointer to the matrix assembly
       * function.
       */
      equation_systems("Helmholtz").attach_assemble_function (assemble_helmholtz);
      
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
     * Solve the system "Helmholtz".
     */
    equation_systems("Helmholtz").solve();


    /**
     * After solving the system write the solution
     * to a GMV-formatted plot file.
     */
    mesh.write_gmv ("out.gmv", equation_systems);
  };


#ifdef HAVE_PETSC

  PetscFinalize();
  
#endif

  
  /**
   * All done.  
   */
  return 0;
};




void assemble_helmholtz(EquationSystems& es,
			const std::string& system_name)
{
  /**
   * It is a good idea to make sure we are assembling
   * the proper system.
   */
  assert (system_name == "Helmholtz");

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
  FEType fe_type = es("Helmholtz").dof_map.component_type(0);

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
  const DofMap& dof_map = es("Helmholtz").get_dof_map();

  /**
   * Define data structures to contain the element matrix
   * [Ae] and right-hand-side vector [Fe] contribution.
   * Here [Ae] is computed using the element stiffness [Ke],
   * mass [Me] and damping [Ce]  matrices, where
   * [Ae] = [Ke]+i*omega*[Ce]-(omega/c)^2*[Me].
   */
  ComplexDenseMatrix   Ae, Ke, Me, Ce;
  std::vector<Complex> Fe;

  /**
   * Calculate the circular frequency omega and define the fluid
   * properties (here: density roh and speed of sound c). 
   */
  const Real omega = frequency*2*3.141592653589793;

  const Real c   = 1.0;
  const Real rho = 1.0;

  /**
   * Define I = 0 + i*1.0
   */ 
  Complex I(0.0, 1.0);

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
      Ae.resize (dof_indices.size(),
		 dof_indices.size());

      Ke.resize (dof_indices.size(),
		 dof_indices.size());

      Ce.resize (dof_indices.size(),
		 dof_indices.size());

      Me.resize (dof_indices.size(),
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
		Me(i,j) += JxW[qp]*phi[i][qp]*phi[j][qp];
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
	      QGauss qface0(dim,   SECOND);
	      QGauss qface1(dim-1, SECOND);
	      
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
	       * Compute the shape function values on the element
	       * face.
	       */
	      fe_face->reinit(&qface1, elem, side);

	      /**
	       * Here we consider a normal velocity vn=1 applied to
	       * the whole boundary of our mesh.
	       */ 
	      const Real vn_value = 1.0;

	      /**
	       * Consider a normal admittance an=1
	       * at some parts of the bounfdary
	       */
	      const Real an_value = 1.0;
	      
	      /**
	       * Loop over the face quagrature points for integration.
	       */
	      for (unsigned int qp=0; qp<qface0.n_points(); qp++)
		{

		  /**
		   * Right-hand-side contribution due to prescribed
		   * normal velocity.
		   */
		  for (unsigned int i=0; i<phi_face.size(); i++)
		    {
		      Fe[i] += -I*vn_value*rho*omega
			*phi_face[i][qp]*JxW_face[qp];
		    };

		  /**
		   * Element matrix contributrion due to precribed
		   * admittance boundary conditions.
		   */
		  for (unsigned int i=0; i<phi_face.size(); i++)
		    for (unsigned int j=0; j<phi_face.size(); j++)
		      {
			Ce(i,j) += rho*an_value
			  *JxW_face[qp]*phi_face[i][qp]*phi_face[j][qp];
		      };

		}; // end face quadrature point loop	  
	    }; // end if (elem->neighbor(side) == NULL)
      }; // end boundary condition section	  

      /**
       * Compute [Ae]
       */
      Ae.add( 1, Ke);
      Ae.add(I*omega, Ce);
      Ae.add(-omega*omega/(c*c), Me);      
      
      /**
       *----------------------------------------------------------------
       * The element matrix and right-hand-side are now built
       * for this element.  Add them to the global matrix and
       * right-hand-side vector.  The \p PetscMatrix::add_matrix()
       * and \p PetscVector::add_vector() members do this for us.
       */
      es("Helmholtz").matrix.add_matrix (Ae, dof_indices);
      es("Helmholtz").rhs.add_vector    (Fe, dof_indices);
      
    }; // end of element loop
  
  /**
   * All done!
   */
  return;
};

