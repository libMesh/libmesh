/* $Id: ex16.C,v 1.9 2006-06-25 00:21:16 benkirk Exp $ */

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



// <h1>Example 16 - Solving an Eigen Problem</h1>
//
// This example introduces the EigenSystem and shows
// how libMesh can be used for eigenvalue analysis.
// 
// For solving eigen problems, libMesh interfaces
// SLEPc (www.grycap.upv.es/slepc/) which again is based on PETSc.
// Hence, this example will only work if the library is compiled
// with SLEPc support enabled.
//
// In this example some eigenvalues for a standard symmetric eigenvalue
// problem A*x=lambda*x are computed, where the matrix A
// is assembled according to a mass matrix.


// libMesh include files.
#include "libmesh.h"
#include "mesh.h"
#include "mesh_generation.h"
#include "gmv_io.h"
#include "eigen_system.h"
#include "equation_systems.h"
#include "fe.h"
#include "quadrature_gauss.h"
#include "dense_matrix.h"
#include "sparse_matrix.h"
#include "numeric_vector.h"
#include "dof_map.h"


// Function prototype.  This is the function that will assemble
// the eigen system. Here, we will simply assemble a mass matrix.
void assemble_mass(EquationSystems& es,
		   const std::string& system_name);



int main (int argc, char** argv)
{
  // Initialize libMesh and the dependent libraries.
  libMesh::init (argc, argv);

  // This example is designed for the SLEPc eigen solver interface.
#ifndef HAVE_SLEPC

  std::cerr << "ERROR: This example requires libMesh to be\n"
	    << "compiled with SLEPc eigen solvers support!"
	    << std::endl;

  return 0;
#else


  // Braces are used to force object scope.  
  {
    // Check for proper usage.
    if (argc < 3)
      {
	std::cerr << "\nUsage: " << argv[0]
		  << " -n <number of eigen values>"
		  << std::endl;
	error();
      }
    
    // Tell the user what we are doing.
    else 
      {
	std::cout << "Running " << argv[0];
	
	for (int i=1; i<argc; i++)
	  std::cout << " " << argv[i];
	
	std::cout << std::endl << std::endl;
      }

    // Set the dimensionality.
    const unsigned int dim = 2;

    // Get the number of eigen values to be computed from argv[2]
    const unsigned int nev = std::atoi(argv[2]);

    // Create a dim-dimensional mesh.
    Mesh mesh (dim);

    // Use the internal mesh generator to create a uniform
    // grid on a square.
    MeshTools::Generation::build_square (mesh, 
					 20, 20,
					 -1., 1.,
					 -1., 1.,
					 QUAD4);

    // Print information about the mesh to the screen.
    mesh.print_info();
    
    // Create an equation systems object.
    EquationSystems equation_systems (mesh);

    // Create a EigenSystem named "Eigensystem" and (for convenience)
    // use a reference to the system we create.
    EigenSystem & eigen_system =
      equation_systems.add_system<EigenSystem> ("Eigensystem");

    // Declare the system variables.
    {
      // Adds the variable "p" to "Eigensystem".   "p"
      // will be approximated using second-order approximation.
      eigen_system.add_variable("p", FIRST);

      // Give the system a pointer to the matrix assembly
      // function defined below.
      eigen_system.attach_assemble_function (assemble_mass);

      // Set necessary parametrs used in EigenSystem::solve(),
      // i.e. the number of requested eigenpairs \p nev and the number
      // of basis vectors \p ncv used in the solution algorithm. Note that
      // ncv >= nev must hold and ncv >= 2*nev is recommended.
      equation_systems.parameters.set<unsigned int>("eigenpairs")    = nev;
      equation_systems.parameters.set<unsigned int>("basis vectors") = nev*3;

      // Set the eigen solver type. SLEPc offers various solvers such as
      // the Arnoldi and subspace method. It
      // also offers interfaces to other solver packages (e.g. ARPACK).
      eigen_system.eigen_solver->set_eigensolver_type(ARNOLDI);

      // Set the solver tolerance and the maximum number of iterations. 
      equation_systems.parameters.set<Real>("linear solver tolerance") = pow(TOLERANCE, 5./3.);
      equation_systems.parameters.set<unsigned int>
	("linear solver maximum iterations") = 1000;

      // Initialize the data structures for the equation system.
      equation_systems.init();

      // Prints information about the system to the screen.
      equation_systems.print_info();

    }
       
    // Solve the system "Eigensystem".
    eigen_system.solve();

    // Get the number of converged eigen pairs.
    unsigned int nconv = eigen_system.get_n_converged();

    std::cout << "Number of converged eigenpairs: " << nconv
	      << "\n" << std::endl;

    // Get the last converged eigenpair
    if (nconv != 0)
      {
	eigen_system.get_eigenpair(nconv-1);
	
	// Write the eigen vector to file.
	char buf[14];
	sprintf (buf, "out.gmv");
	GMVIO (mesh).write_equation_systems (buf, equation_systems);
      }
    else
      {
	std::cout << "WARNING: Solver did not converge!\n" << nconv << std::endl;
      }
  }

#endif // HAVE_SLEPC

  // All done.  
  return libMesh::close ();
}




void assemble_mass(EquationSystems& es,
		   const std::string& system_name)
{
  
  // It is a good idea to make sure we are assembling
  // the proper system.
  assert (system_name == "Eigensystem");

#ifdef HAVE_SLEPC

  // Get a constant reference to the mesh object.
  const Mesh& mesh = es.get_mesh();

  // The dimension that we are running.
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to our system.
  EigenSystem & eigen_system = es.get_system<EigenSystem> (system_name);

  // Get a constant reference to the Finite Element type
  // for the first (and only) variable in the system.
  FEType fe_type = eigen_system.get_dof_map().variable_type(0);

  // A reference to the system matrix
  SparseMatrix<Number>&  matrix_A = *eigen_system.matrix_A;

  // Build a Finite Element object of the specified type.  Since the
  // \p FEBase::build() member dynamically creates memory we will
  // store the object as an \p AutoPtr<FEBase>.  This can be thought
  // of as a pointer that will clean up after itself.
  AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));
  
  // A  Gauss quadrature rule for numerical integration.
  // Use the default quadrature order.
  QGauss qrule (dim, fe_type.default_quadrature_order());

  // Tell the finite element object to use our quadrature rule.
  fe->attach_quadrature_rule (&qrule);

  // The element Jacobian * quadrature weight at each integration point.   
  const std::vector<Real>& JxW = fe->get_JxW();

  // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector<Real> >& phi = fe->get_phi();

  // The element shape function gradients evaluated at the quadrature
  // points.
  // const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();

  // A reference to the \p DofMap object for this system.  The \p DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.
  const DofMap& dof_map = eigen_system.get_dof_map();

  // The element mass matrix.
  DenseMatrix<Number>   Me;

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<unsigned int> dof_indices;


  // Now we will loop over all the elements in the mesh.
  // We will compute the element matrix contribution.

  MeshBase::const_element_iterator       el     = mesh.elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.elements_end();
 
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

      // Zero the element matrices and rhs before
      // summing them.  We use the resize member here because
      // the number of degrees of freedom might have changed from
      // the last element.  Note that this will be the case if the
      // element type is different (i.e. the last element was a
      // triangle, now we are on a quadrilateral).
      Me.resize (dof_indices.size(), dof_indices.size());

      // Now loop over the quadrature points.  This handles
      // the numeric integration.
      // 
      // We will build the element matrix.  This involves
      // a double loop to integrate the test funcions (i) against
      // the trial functions (j).
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
	for (unsigned int i=0; i<phi.size(); i++)
	  for (unsigned int j=0; j<phi.size(); j++)
	      Me(i,j) += JxW[qp]*phi[i][qp]*phi[j][qp];
	  

      // Finally, simply add the element contribution to the
      // overall matrix.
      matrix_A.add_matrix (Me, dof_indices);


    } // end of element loop

#endif // HAVE_SLEPC

  /**
   * All done!
   */
  return;

}
























