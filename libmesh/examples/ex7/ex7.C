/* $Id: ex7.C,v 1.37 2005-01-06 21:54:59 benkirk Exp $ */
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




 // <h1>Example 7 - Introduction to Complex Numbers and the "FrequencySystem"</h1>
 //
 // This is the seventh example program.  It builds on
 // the previous example programs, introduces complex
 // numbers and the FrequencySystem class to solve a 
 // simple Helmholtz equation grad(p)*grad(p)+(omega/c)^2*p=0,
 // for multiple frequencies rather efficiently.
 //
 // The FrequencySystem class offers two solution styles,
 // namely to solve large systems, or to solve
 // moderately-sized systems fast, for multiple frequencies.
 // The latter approach is implemented here.
 //
 // For this example the library has to be compiled with
 // complex numbers enabled. 
 
// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <stdio.h>

// Basic include files needed for overall functionality.
#include "libmesh.h"
#include "libmesh_logging.h"
#include "mesh.h"
#include "mesh_generation.h"
#include "gmv_io.h"
#include "equation_systems.h"

// Include FrequencySystem.  Compared to GeneralSystem,
// this class offers added functionality for the solution of 
// frequency-dependent systems.
#include "frequency_system.h"

// Define the Finite Element object.
#include "fe.h"

// Define Gauss quadrature rules.
#include "quadrature_gauss.h"

// Define useful datatypes for finite element
// matrix and vector components.
#include "dense_matrix.h"
#include "dense_vector.h"

// Define matrix and vector data types for the global 
// equation system.  These are base classes,
// from which specific implementations, like
// the PETSc or LASPACK implementations, are derived.
#include "sparse_matrix.h"
#include "numeric_vector.h"

// Define the DofMap, which handles degree of freedom
// indexing.
#include "dof_map.h"

// Function prototype.  This is the function that will assemble
// the mass, damping and stiffness matrices.  It will <i>not</i>
// form an overall system matrix ready for solution.
void assemble_helmholtz(EquationSystems& es,
			const std::string& system_name);

// Function prototype.  This is the function that will combine
// the previously-assembled mass, damping and stiffness matrices
// to the overall matrix, which then renders ready for solution.
void add_M_C_K_helmholtz(EquationSystems& es,
			 const std::string& system_name);

// Begin the main program.  Note that this example only
// works correctly if complex numbers have been enabled
// in the library.  In order to link against the complex
// PETSc libraries, you must have built PETSc with the same
// C++ compiler that you used to build libMesh.  This is
// so that the name mangling will be the same for the
// routines in both libraries.
int main (int argc, char** argv)
{
  // Initialize Petsc, like in example 2.
  libMesh::init (argc, argv);
  
  // This example is designed for complex numbers.   
#ifndef USE_COMPLEX_NUMBERS

  std::cerr << "ERROR: This example is intended for " << std::endl
	    << " use with complex numbers." << std::endl;
  here();

  return 0;

#else
  
  // Braces are used to force object scope, like in example 2
  {
    // Check for proper usage.
    if (argc < 3)
      {
	std::cerr << "Usage: " << argv[0] << " -f [frequency]"
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
    
    // For now, restrict to dim=2, though this
    // may easily be changed, see example 4
    const unsigned int dim = 2;
    
    // Get the frequency from argv[2] as a <i>float</i>,
    // currently, solve for 1/3rd, 2/3rd and 1/1th of the given frequency
    const Real frequency_in = atof(argv[2]);
    const unsigned int n_frequencies = 3;
    
    // mesh discretization depends on frequency (badly guessed estimate...?)
    const unsigned int n_el_per_dim =
      static_cast<unsigned int>(frequency_in*40.);
    
    // Tell the user the number of elements
    std::cout << " Using " << n_el_per_dim << " x " 
	      << n_el_per_dim << " = " 
	      << n_el_per_dim*n_el_per_dim
	      << " QUAD9 elements"
	      << std::endl << std::endl;
    
    // Create a dim-dimensional mesh.
    Mesh mesh (dim);
    
    // Use the internal mesh generator to create a uniform
    // grid on the square [-1,1]^2.  We instruct the mesh generator
    // to build a mesh of n x n Quad9 elements.
    MeshTools::Generation::build_square (mesh,
					 n_el_per_dim, n_el_per_dim,
					 -1., 1.,
					 -1., 1.,
					 QUAD9);
    
    // Print information about the mesh to the screen.
    mesh.print_info();
    
    // Create an equation systems object, which now handles
    // a frequency system, as opposed to previous examples.
    EquationSystems equation_systems (mesh);
    
    // Create a FrequencySystem named "Helmholtz" & store a
    // reference to it.
    FrequencySystem & f_system =      
      equation_systems.add_system<FrequencySystem> ("Helmholtz");
    
    // Add the variable "p" to "Helmholtz".  "p"
    // will be approximated using second-order approximation.
    f_system.add_variable("p", SECOND);
    
    // Tell the frequency system about the two user-provided
    // functions.  In other circumstances, at least the
    // solve function has to be attached.
    f_system.attach_assemble_function (assemble_helmholtz);
    f_system.attach_solve_function    (add_M_C_K_helmholtz);
    
    // To enable the fast solution scheme, additional
    // <i>global</i> matrices and one global vector, all appropriately sized,
    // have to be added.  The system object takes care of the
    // appropriate size, but the user should better fill explicitly
    // the sparsity structure of the overall matrix, so that the
    // fast matrix addition method can be used, as will be shown later.
    f_system.add_matrix ("stiffness");
    f_system.add_matrix ("damping");
    f_system.add_matrix ("mass");
    f_system.add_vector ("rhs");
    
    // Communicates the frequencies to the system.  Note that
    // the frequency system stores the frequencies as parameters
    // in the equation systems object, so that our assemble and solve
    // functions may directly access them.
    // Will solve for 1/3rd, 2/3rd and 1/1th of the given frequency
    f_system.set_frequencies_by_steps (frequency_in/n_frequencies,
				       frequency_in,
				       n_frequencies);
    
    // Use the parameters of the equation systems object to
    // tell the frequency system about the wave velocity and fluid
    // density.  The frequency system provides default values, but
    // these may be overridden, as shown here.
    equation_systems.parameters.set<Real> ("wave speed") = 1.;
    equation_systems.parameters.set<Real> ("rho")        = 1.;
    
    // Initialize the data structures for the equation system.  <i>Always</i>
    // prior to this, the frequencies have to be communicated to the system.
    equation_systems.init ();
    
    // Prints information about the system to the screen.
    equation_systems.print_info ();

    for (unsigned int n=0; n < n_frequencies; n++)
      {
	// Solve the system "Helmholtz" for the n-th frequency.  
	// Since we attached an assemble() function to the system,
	// the mass, damping and stiffness contributions will only
	// be assembled once.  Then, the system is solved for the
	// given frequencies.  Note that solve() may also solve 
	// the system only for specific frequencies.
	f_system.solve (n,n);
	
	// After solving the system, write the solution
	// to a GMV-formatted plot file, for every frequency.  
	// Now this is nice ;-) : we have the <i>identical</i> 
	// interface to the mesh write method as in the real-only 
	// case, but we output the real and imaginary 
	// part, and the magnitude, where the variable 
	// "p" is prepended with "r_", "i_", and "a_", 
	// respectively.
	char buf[14];
	sprintf (buf, "out%04d.gmv", n);
	GMVIO(mesh).write_equation_systems (buf,
					    equation_systems);
      }
    
    // Alternatively, the whole EquationSystems object can be
    // written to disk.  By default, the additional vectors are also
    // saved.
    equation_systems.write ("eqn_sys.dat", libMeshEnums::WRITE);
  }
  
  // All done.  
  return libMesh::close ();

#endif 
}


// Here we define the matrix assembly routine for
// the Helmholtz system.  This function will be
// called to form the stiffness matrix and right-hand side.
void assemble_helmholtz(EquationSystems& es,
			const std::string& system_name)
{
#ifdef USE_COMPLEX_NUMBERS
    
  // It is a good idea to make sure we are assembling
  // the proper system.
  assert (system_name == "Helmholtz");
  
  // Get a constant reference to the mesh object.
  const Mesh& mesh = es.get_mesh();
  
  // The dimension that we are in
  const unsigned int dim = mesh.mesh_dimension();
  
  // Get a reference to our system, as before
  FrequencySystem & f_system =
    es.get_system<FrequencySystem> (system_name);
  
  // A const reference to the DofMap object for this system.  The DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.
  const DofMap& dof_map = f_system.get_dof_map();
  
  // Get a constant reference to the Finite Element type
  // for the first (and only) variable in the system.
  const FEType fe_type = dof_map.variable_type(0);

  // For the admittance boundary condition,
  // get the fluid density
  const Real rho = es.parameters.get<Real>("rho");
  
  // In here, we will add the element matrices to the
  // <i>additional</i> matrices "stiffness_mass", "damping",
  // and the additional vector "rhs", not to the members 
  // "matrix" and "rhs".  Therefore, get writable
  // references to them
  SparseMatrix<Number>&   stiffness      = f_system.get_matrix("stiffness");
  SparseMatrix<Number>&   damping        = f_system.get_matrix("damping");
  SparseMatrix<Number>&   mass           = f_system.get_matrix("mass");
  NumericVector<Number>&  freq_indep_rhs = f_system.get_vector("rhs");
  
  // Some solver packages (PETSc) are especially picky about
  // allocating sparsity structure and truly assigning values
  // to this structure.  Namely, matrix additions, as performed
  // later, exhibit acceptable performance only for identical
  // sparsity structures.  Therefore, explicitly zero the
  // values in the collective matrix, so that matrix additions
  // encounter identical sparsity structures.
  SparseMatrix<Number>&  matrix           = *f_system.matrix;
  
  // ------------------------------------------------------------------
  // Finite Element related stuff
  //
  // Build a Finite Element object of the specified type.  Since the
  // FEBase::build() member dynamically creates memory we will
  // store the object as an AutoPtr<FEBase>.  This can be thought
  // of as a pointer that will clean up after itself.
  AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));
  
  // A 5th order Gauss quadrature rule for numerical integration.
  QGauss qrule (dim, FIFTH);

  // Tell the finite element object to use our quadrature rule.
  fe->attach_quadrature_rule (&qrule);
  
  // The element Jacobian// quadrature weight at each integration point.   
  const std::vector<Real>& JxW = fe->get_JxW();

  // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector<Real> >& phi = fe->get_phi();
  
  // The element shape function gradients evaluated at the quadrature
  // points.
  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();
  
  // Now this is slightly different from example 4.
  // We will not add directly to the overall (PETSc/LASPACK) matrix,
  // but to the additional matrices "stiffness_mass" and "damping".
  // The same holds for the right-hand-side vector Fe, which we will
  // later on store in the additional vector "rhs". 
  // The zero_matrix is used to explicitly induce the same sparsity
  // structure in the overall matrix.
  // see later on. (At least) the mass, and stiffness matrices, however, 
  // are inherently real.  Therefore, store these as one complex
  // matrix.  This will definitely save memory.
  DenseMatrix<Number> Ke, Ce, Me, zero_matrix;
  DenseVector<Number> Fe;
  
  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<unsigned int> dof_indices;

  // Now we will loop over all the elements in the mesh.
  // We will compute the element matrix and right-hand-side
  // contribution.
  //const_local_elem_iterator           el (mesh.elements_begin());
  //const const_local_elem_iterator end_el (mesh.elements_end());

  MeshBase::const_element_iterator           el = mesh.local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.local_elements_end();
  
  for ( ; el != end_el; ++el)
    {
      // Start logging the element initialization.
      START_LOG("elem init","assemble_helmholtz");
      
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
      
      // Zero & resize the element matrix and right-hand side before
      // summing them, with different element types in the mesh this
      // is quite necessary.
      {
        const unsigned int n_dof_indices = dof_indices.size();

	Ke.resize          (n_dof_indices, n_dof_indices);
	Ce.resize          (n_dof_indices, n_dof_indices);
	Me.resize          (n_dof_indices, n_dof_indices);
	zero_matrix.resize (n_dof_indices, n_dof_indices);
	Fe.resize          (n_dof_indices);
      }
      
      // Stop logging the element initialization.
      STOP_LOG("elem init","assemble_helmholtz");

      // Now loop over the quadrature points.  This handles
      // the numeric integration.
      START_LOG("stiffness & mass","assemble_helmholtz");

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
	{
	  // Now we will build the element matrix.  This involves
	  // a double loop to integrate the test funcions (i) against
	  // the trial functions (j).  Note the braces on the rhs
	  // of Ke(i,j): these are quite necessary to finally compute
	  // Real*(Point*Point) = Real, and not something else...
	  for (unsigned int i=0; i<phi.size(); i++)
	    for (unsigned int j=0; j<phi.size(); j++)
	      {
		Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
		Me(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp]);
	      }	  
	}

      STOP_LOG("stiffness & mass","assemble_helmholtz");

      // Now compute the contribution to the element matrix and the
      // right-hand-side vector if the current element lies on the
      // boundary. 
      //
      // The following loops over the sides of the element.
      // If the element has no neighbor on a side then that
      // side MUST live on a boundary of the domain.
	 
      for (unsigned int side=0; side<elem->n_sides(); side++)
	if (elem->neighbor(side) == NULL)
	  {
	    START_LOG("damping & rhs","assemble_helmholtz");
	      
	    // Declare a special finite element object for
	    // boundary integration.
	    AutoPtr<FEBase> fe_face (FEBase::build(dim, fe_type));
	      
	    // Boundary integration requires one quadraure rule,
	    // with dimensionality one less than the dimensionality
	    // of the element.
	    QGauss qface(dim-1, SECOND);
	      
	    // Tell the finte element object to use our
	    // quadrature rule.
	    fe_face->attach_quadrature_rule (&qface);
	      
	    // The value of the shape functions at the quadrature
	    // points.
	    const std::vector<std::vector<Real> >&  phi_face =
	      fe_face->get_phi();
	      
	    // The Jacobian// Quadrature Weight at the quadrature
	    // points on the face.
	    const std::vector<Real>& JxW_face = fe_face->get_JxW();
	      
	    // Compute the shape function values on the element
	    // face.
	    fe_face->reinit(elem, side);

	    // Here we consider a normal velocity vn=1 applied to
	    // the whole boundary of our mesh.
	    const Real vn_value = 1.;
	      
	    // Consider a normal admittance an=1
	    // at some parts of the bounfdary
	    const Real an_value = 1.;
	      
	    // Loop over the face quadrature points for integration.
	    for (unsigned int qp=0; qp<qface.n_points(); qp++)
	      {
		// Right-hand-side contribution due to prescribed
		// normal velocity.
		for (unsigned int i=0; i<phi_face.size(); i++)
		  Fe(i) += vn_value*phi_face[i][qp]*JxW_face[qp];
		
		// Element matrix contributrion due to precribed
		// admittance boundary conditions.
		for (unsigned int i=0; i<phi_face.size(); i++)
		  for (unsigned int j=0; j<phi_face.size(); j++)
		    Ce(i,j) += rho*an_value*JxW_face[qp]*phi_face[i][qp]*phi_face[j][qp];
	      }

	    STOP_LOG("damping & rhs","assemble_helmholtz");
	  }
      
      // Finally, simply add the contributions to the additional
      // matrices and vector.
      stiffness.add_matrix      (Ke, dof_indices);
      damping.add_matrix        (Ce, dof_indices);
      mass.add_matrix           (Me, dof_indices);
      freq_indep_rhs.add_vector (Fe, dof_indices);
      
      // For the overall matrix, explicitly zero the entries where
      // we added values in the other ones, so that we have 
      // identical sparsity footprints.
      matrix.add_matrix(zero_matrix, dof_indices);
    }
  
  // All done!
#endif
}


// We now define the function which will combine
// the previously-assembled mass, stiffness, and
// damping matrices into a single system matrix.
void add_M_C_K_helmholtz(EquationSystems& es,
			 const std::string& system_name)
{
#ifdef USE_COMPLEX_NUMBERS

  START_LOG("init phase","add_M_C_K_helmholtz");
  
  // It is a good idea to make sure we are assembling
  // the proper system.
  assert (system_name == "Helmholtz");
  
  // Get a reference to our system, as before
  FrequencySystem & f_system =
    es.get_system<FrequencySystem> (system_name);
  
  // Get the frequency, fluid density, and speed of sound
  // for which we should currently solve
  const Real frequency = es.parameters.get<Real> ("current frequency");
  const Real rho       = es.parameters.get<Real> ("rho");
  const Real speed     = es.parameters.get<Real> ("wave speed");
  
  // Compute angular frequency omega and wave number k
  const Real omega = 2.0*libMesh::pi*frequency;
  const Real k     = omega / speed;
  
  // Get writable references to the overall matrix and vector, where the 
  // frequency-dependent system is to be collected
  SparseMatrix<Number>&  matrix          = *f_system.matrix;
  NumericVector<Number>& rhs             = *f_system.rhs;
  
  // Get writable references to the frequency-independent matrices
  // and rhs, though we only need to extract values.  This write access
  // is necessary, since solver packages have to close the data structure 
  // before they can extract values for computation.
  SparseMatrix<Number>&   stiffness      = f_system.get_matrix("stiffness");
  SparseMatrix<Number>&   damping        = f_system.get_matrix("damping");
  SparseMatrix<Number>&   mass           = f_system.get_matrix("mass");
  NumericVector<Number>&  freq_indep_rhs = f_system.get_vector("rhs");
  
  // form the scaling values for the coming matrix and vector axpy's
  const Number scale_stiffness (  1., 0.   );
  const Number scale_damping   (  0., omega);
  const Number scale_mass      (-k*k, 0.   );
  const Number scale_rhs       (  0., -(rho*omega));
  
  // Now simply add the matrices together, store the result
  // in matrix and rhs.  Clear them first.
  matrix.close(); matrix.zero ();
  rhs.close();    rhs.zero    ();
  
  // The matrices from which values are added to another matrix
  // have to be closed.  The add() method does take care of 
  // that, but let us do it explicitly.
  stiffness.close ();
  damping.close   ();
  mass.close      ();

  STOP_LOG("init phase","add_M_C_K_helmholtz");

  START_LOG("global matrix & vector additions","add_M_C_K_helmholtz");
  
  // add the stiffness and mass with the proper frequency to the
  // overall system.  For this to work properly, matrix has
  // to be not only initialized, but filled with the identical
  // sparsity structure as the matrix added to it, otherwise
  // solver packages like PETSc crash.
  //
  // Note that we have to add the mass and stiffness contributions
  // one at a time; otherwise, the real part of matrix would
  // be fine, but the imaginary part cluttered with unwanted products.
  matrix.add (scale_stiffness, stiffness);
  matrix.add (scale_mass,      mass);
  matrix.add (scale_damping,   damping);
  rhs.add    (scale_rhs,       freq_indep_rhs);

  STOP_LOG("global matrix & vector additions","add_M_C_K_helmholtz");
  
  // The "matrix" and "rhs" are now ready for solution   
#endif
}

