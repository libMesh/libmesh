// $Id: ex3.C,v 1.15 2003-02-24 14:35:52 benkirk Exp $

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
#include "general_system.h"
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



/**
 * Function prototype for the exact solution.
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
     * to build a mesh of 8x8 \p Quad9 elements.  Building \p Quad9
     * elements instead of the default \p Quad4's we used in example 2
     * allow us to use higher-order approximation.
     */
    mesh.build_square (8, 8,
		       -1., 1.,
		       -1., 1.,
		       QUAD9);

    /**
     * This is the first use of the \p Mesh::find_neighbors() method.
     * When this method is called all the elements are interrogated
     * and neighbors are found.  After this method is called then
     * calling mesh.elem(3)->neighbor(2) will return a pointer to
     * the element in the mesh that borders side 2 of element number
     * 3.  The \p Elem::neighbor(unsigned int s)  member returns
     * \p NULL if the \p s side of the element is on a boundary of
     * the domain.
     */
    mesh.find_neighbors();
    
    /**
     * Print information about the mesh to the screen.
     * Note that 5x5 \p Quad9 elements actually has 11x11 nodes,
     * so this mesh is significantly larger than the one in example 2.
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


    /**
     * After solving the system write the solution
     * to a GMV-formatted plot file.
     */
    mesh.write_gmv ("out.gmv", equation_systems);
  };


  /**
   * All done.  
   */
  return libMesh::close();
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
  FEType fe_type = es("Poisson").get_dof_map().variable_type(0);

  /**
   * Build a Finite Element object of the specified type.  Since the
   * \p FEBase::build() member dynamically creates memory we will
   * store the object as an \p AutoPtr<FEBase>.  This can be thought
   * of as a pointer that will clean up after itself.  Example 4
   * describes some advantages of \p AutoPtr's in the context of
   * quadrature rules.
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
  DenseMatrix<Number> Ke;
  std::vector<Number> Fe;

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
	      const Real eps = 1.e-3;

	      /**
	       * fxy is the forcing function for the Poisson equation.
	       * In this case we set fxy to be a finite difference
	       * Laplacian approximation to the (known) exact solution.
	       *
	       * We will use the second-order accurate FD Laplacian
	       * approximation, which in 2D is
	       *
	       * u_xx + u_yy = (u(i,j-1) + u(i,j+1) +
	       *                u(i-1,j) + u(i+1,j) +
	       *                -4*u(i,j))/h^2
	       */
	      const Real fxy = -(exact_solution(x,y-eps) +
				 exact_solution(x,y+eps) +
				 exact_solution(x-eps,y) +
				 exact_solution(x+eps,y) -
				 4.*exact_solution(x,y))/eps/eps;
	      
	      Fe[i] += JxW[qp]*fxy*phi[i][qp];
	    }; // end of the RHS summation loop
	  
	}; // end of quadrature point loop



      
      /**
       *----------------------------------------------------------------
       * At this point the interior element integration has
       * been completed.  However, we have not yet addressed
       * boundary conditions.  For this example we will only
       * consider simple Dirichlet boundary conditions.
       *
       * There are several ways Dirichlet boundary conditions
       * can be imposed.  A simple approach, which works for
       * interpolary bases like you have with standard Lagrange
       * finite elements, is to assing function values to the
       * degrees of freedom living on the domain boundary. This
       * works well for interpolary bases, but is more difficult
       * when non-interpolary (e.g Legendre or Hierarchic) bases
       * are used.
       *
       * Dirichlet boundary conditions can also be imposed with a
       * "penalty" method.  In this case essentially the L2 projection
       * of the boundary values are added to the matrix. The
       * projection is multiplied by some large factor so that, in
       * floating point arithmetic, the existing (smaller) entries
       * in the matrix and right-hand-side are effectively ignored.
       *
       * This amounts to adding a term of the form
       *
       * \f$ \frac{1}{\epsilon} \int_{\delta \Omega} \phi_i \phi_j = \frac{1}{\epsilon} \int_{\delta \Omega} u \phi_i \f$
       *
       * where \f$ \frac{1}{\epsilon} \f$ is the penalty parameter when \f$ \epsilon \lle 1 \f$
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
	       * Boundary integration requires one quadraure rule,
	       * with dimensionality one less than the dimensionality
	       * of the element.
	       */
	      QGauss qface(dim-1, FIFTH);
	      
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
	       * The XYZ locations (in physical space) of the
	       * quadrature points on the face.  This is where
	       * we will interpolate the boundary value function.
	       */
	      const std::vector<Point >& qface_point = fe_face->get_xyz();
	      
	      /**
	       * Compute the shape function values on the element
	       * face.
	       */
	      fe_face->reinit(elem, side);
	      
	      /**
	       * Loop over the face quagrature points for integration.
	       */
	      for (unsigned int qp=0; qp<qface.n_points(); qp++)
		{
		  /**
		   * The location on the boundary of the current
		   * face quadrature point.
		   */
		  const Real xf = qface_point[qp](0);
		  const Real yf = qface_point[qp](1);
		  
		  /**
		   * The penalty value.  \f$ \frac{1}{\epsilon \f$
		   * in the discussion above.
		   */
		  const Real penalty = 1.e10;
		  
		  /**
		   * The boundary value.
		   */
		  const Real value = exact_solution(xf, yf);
		  
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
		    }		  
		} // end face quadrature point loop	  
	    } // end if (elem->neighbor(side) == NULL)
      } // end boundary condition section	  


      

      
      /**
       *----------------------------------------------------------------
       * The element matrix and right-hand-side are now built
       * for this element.  Add them to the global matrix and
       * right-hand-side vector.  The \p PetscMatrix::add_matrix()
       * and \p PetscVector::add_vector() members do this for us.
       */
      es("Poisson").matrix->add_matrix (Ke, dof_indices);
      es("Poisson").rhs->add_vector (Fe, dof_indices);
      
    }; // end of element loop


  
  /**
   * All done!
   */
  return;
};
