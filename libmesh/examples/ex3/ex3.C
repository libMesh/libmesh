// $Id: ex3.C,v 1.1 2003-01-31 21:22:11 benkirk Exp $

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
 * \mainpage Example 3
 *
 * \section Introduction
 *
 * This is the third example program.  It builds on
 * the second example program by showing how to solve a simple
 * Poisson system.  Note that we will not comment on things that
 * were already explained in the second example.
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
   * Braces are used to force object scope, like in example 2
   */   
  {
    /**
     * Tell the user what we are doing.
     */
    {
      std::cout << "Running " << argv[0];
      
      for (int i=1; i<argc; i++)
	std::cout << " " << argv[i];
      
      std::cout << std::endl << std::endl;
    };
    
    
    /**
     * Create a 2D mesh.
     */
    Mesh mesh (2);
    
    /**
     * Use the internal mesh generator to create a uniform
     * grid on the square [-1,1]^2.  We instruct the mesh generator
     * to build a mesh of 5x5 \p Quad9 elements.  Building \p Quad9
     * elements instead of the default \p Quad4's we used in example 2
     * allow us to use higher-order approximation.
     */
    mesh.build_cube (5, 5, 0,
		     -1., 1.,
		     -1., 1.,
		      0., 0.,
		     QUAD9);
    
    /**
     * Print information about the mesh to the screen.
     * Note that 5x5 \p Quad9 elements actually has 11x11 nodes,
     * so this mesh is significantly larger than the one in example 2.
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
       * function.  This will be called when needed by the
       * library.
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
     * Solve the system "Poisson"  Note that calling this
     * member will assemble the linear system and invoke
     * the default Petsc solver, however the solver can be
     * controlled from the command line.  For example,
     * you can invoke conjugate gradient with
     *
     * ./ex3 -ksp_type cg
     *
     * and you can get a nice X-window that monitors the solver
     * convergence with
     *
     * ./ex3 -ksp_xmonitor
     */
    equation_systems("Poisson").solve();
  };


#ifdef HAVE_PETSC

  PetscFinalize();
  
#endif

  
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
  const std::vector<real>& JxW = fe->get_JxW();

  /**
   * The physical XY locations of the quadrature points on the element.
   * These might be useful for evaluating spatially varying material
   * properties at the quadrature points.
   */
  const std::vector<Point>& q_point = fe->get_xyz();

  /**
   * The element shape functions evaluated at the quadrature points.
   */
  const std::vector<std::vector<real> >& phi = fe->get_phi();

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
   * "Ke" and "Fe".
   */
  DenseMatrix       Ke;
  std::vector<real> Fe;

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
	      };


	  /**
	   * Now we build the element right-hand-side contribution.
	   * This involves a single loop in which we integrate the
	   * "forcing function" in the PDE against the test functions.
	   */
	  for (unsigned int i=0; i<phi.size(); i++)
	    {
	      const real fxy = 0.;
	      
	      Fe[i] += JxW[qp]*fxy*phi[i][qp];
	    };
	};
      
      /**
       *----------------------------------------------------------------
       * The element matrix and right-hand-side are now built
       * for this element.  Add them to the global matrix and
       * right-hand-side vector.  The \p PetscMatrix::add_matrix()
       * and \p PetscVector::add_vector() members do this for us.
       */
      es("Poisson").matrix.add_matrix (Ke, dof_indices);
      es("Poisson").rhs.add_vector (Fe, dof_indices);
    };


  
  /**
   * All done!
   */
  return;
};
