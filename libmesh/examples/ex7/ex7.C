// $Id: ex7.C,v 1.6 2003-02-10 22:03:21 benkirk Exp $

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
 * \mainpage Example 7
 *
 * \section Introduction
 *
 * Warning: This example is still under cunstruction.
 *
 * This is the seventh example program.  It builds on
 * the previous example programs.  It introduces complex
 * numbers and solves a simple Helmholtz equation 
 * grad(p)*grad(p)+(omega/c)^2*p=0.
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
 *--------------------------------------------------------------------
 * Initialize some (global) handy constants and variables
 *
 * The frequency for which we are solving the Helmholtz
 * equation.
 */
Real frequency;


/**
 * Define the fluid properties. Here (for simplicity) 
 * we define the density rho = 1. and the speed of sound
 * speed = 1.
 */
const Real speed = 1.;
const Real rho   = 1.;


/**
 * Define pi.
 */
const Real libmesh_pi = acos(-1.);


/**
 * Define the imaginary unit
 * I = 0. + i*1.
 */ 
#ifdef USE_COMPLEX_NUMBERS
Complex I(0.0, 1.0);
#else
/**
 * Do this for compatibility, so that main() can catch
 * the error of compiling this example without complex support.
 */
Complex I(0.);
#endif








int main (int argc, char** argv)
{
  /**
   * Initialize Petsc, like in example 2.
   */
  libMesh::init (argc, argv);

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
    if (argc != 3)
      {
	std::cerr << "Usage: " << argv[0] << " -f [frequency]"
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
     * For now, restrict to dim=2, though this
     * may easily be changed, see example 4
     */
    const unsigned int dim = 2;

    /**
     * Get the frequency from argv[2] as a @e float
     */
    frequency = atof(argv[2]);

    /**
     * mesh discretization depends on frequency (overestimated)
     *
     * @note  For cool picture (a_p in gmv), try: ./ex7 -f 3
     */
    const unsigned int n_el_per_dim = static_cast<unsigned int>(frequency*40.);

    /**
     * Tell the user the number of elements
     */
    std::cout << " Using " << n_el_per_dim << " x " 
	      << n_el_per_dim << " = " 
	      << n_el_per_dim*n_el_per_dim
	      << " QUAD4 elements"
	      << std::endl << std::endl;


    /**
     * Create a dim-dimensional mesh.
     */
    Mesh mesh (dim);

    /**
     * Use the internal mesh generator to create a uniform
     * grid on the square [-1,1]^2.  We instruct the mesh generator
     * to build a mesh of n x n \p Quad4 elements.
     */
    mesh.build_square (n_el_per_dim, n_el_per_dim,
		       -1., 1.,
		       -1., 1.,
		       QUAD4);

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
     * After solving the system, write the solution
     * to a GMV-formatted plot file.  Now this is 
     * nice ;-) : we have the @e identical interface 
     * to the mesh write method as in the real-only 
     * case, but we output both the real and imaginary 
     * part, the variable "p" prepended with "r_"
     * and "i_", respectively.
     */
    mesh.write_gmv ("out.gmv", equation_systems);
  };


  libMesh::close ();
  
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
   * Now this is slightly different from example 4.
   * The matrix that we finally add to the PETSc matrix
   * must be complex-valued: this is \p Ae.
   * The right-hand-side vector \p Fe is also required to be
   * complex-valued, as we will see later on.
   * The mass, damping and stiffness matrices, however, are inherently
   * real.  And since \p DenseMatrix<> offers a method
   * to add a real matrix to a complex matrix, we can safely
   * define element stiffness Ke and mass matrix Me as real. 
   */
  ComplexDenseMatrix   Ae;
  RealDenseMatrix      Ke, Ce, Me;
  std::vector<Complex> Fe;

  /**
   * Calculate the circular frequency omega and define the fluid
   * properties (here: density roh and speed of sound c). 
   */
  const Real omega = frequency*2*libmesh_pi;


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
       * Compute the total, frequency-dependent element
       * matrix  \f$  Ae = Ke - (\omega / speed)^2 Me \f$.
       * Note that real matrices are added to a complex
       * matrix (see above).  The class \p DenseMatrix<>
       * offers this feature.
       */
      Ae.add( 1., Ke);
      Ae.add(I*omega, Ce);
      Ae.add(-omega*omega/(speed*speed), Me);      
      
      /**
       *----------------------------------------------------------------
       * The element matrix and right-hand-side are now built
       * for this element.  Add them to the global matrix and
       * right-hand-side vector.  The \p PetscMatrix::add_matrix()
       * and \p PetscVector::add_vector() members do this for us.
       */
      es("Helmholtz").matrix->add_matrix (Ae, dof_indices);
      es("Helmholtz").rhs->add_vector    (Fe, dof_indices);
      
    }; // end of element loop
  
  /**
   * All done!
   */
  return;
};

