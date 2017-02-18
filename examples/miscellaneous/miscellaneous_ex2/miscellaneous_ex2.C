// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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




// <h1>Miscellaneous Example 2 - Complex Numbers and the "FrequencySystem"</h1>
// \author Steffen Petersen
// \date 2003
//
// This example program introduces complex numbers and the
// FrequencySystem class to solve a simple Helmholtz equation,
// Laplacian(p) + (omega/c)^2*p = 0,
// for multiple frequencies rather efficiently.
//
// The FrequencySystem class offers two solution styles:
// 1.) Solve large systems once, and
// 2.) Solve moderately-sized systems for multiple frequencies.
// The latter approach is implemented here.
//
// This example uses an L-shaped domain and nodal boundary data given
// in the files lshape.unv and lshape_data.unv.  For this example, the
// library has to be compiled with complex numbers enabled.

// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <stdio.h>

// Basic include files needed for overall functionality.
#include "libmesh/libmesh.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/unv_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/elem.h"

// Include FrequencySystem.  This class offers added functionality for
// the solution of frequency-dependent systems.
#include "libmesh/frequency_system.h"

// Define the FE object.
#include "libmesh/fe.h"

// Define the QGauss quadrature rule objects.
#include "libmesh/quadrature_gauss.h"

// Useful datatypes for finite element assembly.
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"

// Sparse matrix and vector data types for parallel linear algebra.
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"

// Define the DofMap, which handles degree of freedom indexing.
#include "libmesh/dof_map.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// This problem is only defined on complex-valued fields, for
// which libMesh must be configured with Number == Complex.
#ifdef LIBMESH_USE_COMPLEX_NUMBERS

// Function prototype.  This is the function that will assemble
// the mass, damping and stiffness matrices.  It will not
// form an overall system matrix ready for solution.
void assemble_helmholtz(EquationSystems & es,
                        const std::string & system_name);

// Function prototype.  This is the function that will combine
// the previously-assembled mass, damping and stiffness matrices
// to the overall matrix, which then renders ready for solution.
void add_M_C_K_helmholtz(EquationSystems & es,
                         const std::string & system_name);
#endif

// This example only works correctly if libmesh has been configured
// with --enable-complex.  If you wish to use PETSc, you must also
// build PETSc with complex number support by configuring with
// --with-scalar-type=complex --with-clanguage=cxx, and using the same
// C++ compiler to build both PETSc and libmesh.  This example also
// works with the Eigen sparse linear solvers which are provided in
// libmesh's contrib directory.
int main (int argc, char ** argv)
{
  // The libMeshInit object initializes MPI, PETSc, etc, and must be
  // constructed before all other objects.
  LibMeshInit init (argc, argv);

  // This example is designed for complex numbers.
#ifndef LIBMESH_USE_COMPLEX_NUMBERS
  libmesh_example_requires(false, "--enable-complex");
#else

  // This example requires at least 2D support.
  libmesh_example_requires(2 <= LIBMESH_DIM, "2D support");

  // Check for proper usage.
  if (argc < 3)
    libmesh_error_msg("Usage: " << argv[0] << " -f [frequency]");

  if (init.comm().size() > 1)
    {
      if (init.comm().rank() == 0)
        {
          libMesh::err << "TODO: This example should be able to run in parallel."
                       << std::endl;
        }
      return 0;
    }

  // Tell the user what we are doing.
  else
    {
      libMesh::out << "Running " << argv[0];

      for (int i=1; i<argc; i++)
        libMesh::out << " " << argv[i];

      libMesh::out << std::endl << std::endl;
    }

  // Get the frequency from argv[2] as a float, solve for 1/3rd, 2/3rd
  // and 1/1th of the given frequency.
  const Real frequency_in = atof(argv[2]);
  const unsigned int n_frequencies = 3;

  // Create a mesh, with dimension to be overridden later, distributed
  // across the default MPI communicator.
  Mesh mesh(init.comm());

  // Read the mesh file. Here the file lshape.unv contains
  // an L-shaped domain in .unv format.
  UNVIO unvio(mesh);
  unvio.read("lshape.unv");

  // Manually prepare the mesh for use, this is not done automatically
  // by the UNVIO reader, so we have to do it here. Note that calling
  // this function renumbers the nodes and elements of the Mesh, but
  // the original numbering is still stored in the UNVIO object for
  // correctly mapping dataset values (see below) to the correct nodes.
  mesh.prepare_for_use();

  // Read the dataset accompanying this problem.  The load on the
  // boundary of the domain is stored in the .unv data file
  // lshape_data.unv.  The data are given as complex-valued normal
  // velocities.
  unvio.read_dataset("lshape_data.unv");

  // Print information about the mesh to the screen.
  mesh.print_info();

  // Create an EquationSystems object.
  EquationSystems equation_systems (mesh);

  // Create a FrequencySystem named "Helmholtz" and store a
  // reference to it.
  FrequencySystem & f_system =
    equation_systems.add_system<FrequencySystem> ("Helmholtz");

  // Add the variable "p" to the "Helmholtz" system.  "p"
  // will be approximated using second-order approximation.
  f_system.add_variable("p", SECOND);

  // The FrequencySystem requires two user-provided functions: one for
  // assembling the different operators and another that specifies how
  // to combine them before the solve.
  f_system.attach_assemble_function (assemble_helmholtz);
  f_system.attach_solve_function    (add_M_C_K_helmholtz);

  // To enable the fast solution scheme, additional sparse matrices
  // and one vector have to be added.  The FrequencySystem object
  // takes care of sizing the additional objects.  The user should
  // still set the sparsity structure of the f_system.matrix, so that
  // the fast matrix addition method can be used.  The procedure for
  // this is shown in detail in the assembly function.
  f_system.add_matrix ("stiffness");
  f_system.add_matrix ("damping");
  f_system.add_matrix ("mass");
  f_system.add_vector ("rhs");

  // Communicate the frequencies to the system.  Note that the
  // frequency system stores the frequencies as parameters in the
  // EquationSystems object, so that our assemble and solve functions
  // may directly access them. Must be called before the
  // EquationSystems object is initialized.  Will solve for 1/3,
  // 2/3, and 1 times the given frequency.
  f_system.set_frequencies_by_steps (frequency_in/n_frequencies,
                                     frequency_in,
                                     n_frequencies);

  // Set the wave velocity and fluid density parameter values,
  // otherwise default values will be used.
  equation_systems.parameters.set<Real> ("wave speed") = 1.;
  equation_systems.parameters.set<Real> ("rho")        = 1.;

  // Initialize the EquationSystems object.
  equation_systems.init ();

  // Set values in the "rhs" vector based on the entries stored in the
  // UNVIO object from the dataset we read in.  These values only need
  // to be set once, as they are the same for every frequency. We can
  // only do this once equation_systems.init() has been called...
  {
    NumericVector<Number> & freq_indep_rhs = f_system.get_vector("rhs");

    MeshBase::const_node_iterator       node_it  = mesh.nodes_begin();
    const MeshBase::const_node_iterator node_end = mesh.nodes_end();

    for ( ; node_it != node_end; ++node_it)
      {
        // the current node pointer
        Node * node = *node_it;

        // Get the data read in from the dataset for the current Node, if any.
        const std::vector<Number> * nodal_data = unvio.get_data(node);

        // Set the rhs value based on values read in from the dataset.
        if (nodal_data)
          {
            unsigned int dn = node->dof_number(/*system=*/0,
                                               /*variable=*/0,
                                               /*component=*/0);
            freq_indep_rhs.set(dn, (*nodal_data)[0]);
          }
      }
  }

  // Print information about the system to the screen.
  equation_systems.print_info ();

  for (unsigned int n=0; n < n_frequencies; n++)
    {
      // Solve the "Helmholtz" system for the n-th frequency.
      // Since we attached an assemble() function to the system,
      // the mass, damping, and stiffness contributions will only
      // be assembled once.  Then, the system is solved for the
      // given frequencies.  Note that solve() may also solve
      // the system only for specific frequencies.
      f_system.solve (n, n);

      // After solving the system, write the solution to an ExodusII
      // file for every frequency.
#ifdef LIBMESH_HAVE_EXODUS_API
      std::ostringstream file_name;

      file_name << "out_"
                << std::setw(4)
                << std::setfill('0')
                << std::right
                << n
                << ".e";

      ExodusII_IO(mesh).write_equation_systems (file_name.str(),
                                                equation_systems);
#endif
    }

  // Alternatively, the whole EquationSystems object can be
  // written to disk.  By default, the additional vectors are also
  // saved.
  equation_systems.write ("eqn_sys.dat", WRITE);

  // All done.
  return 0;

#endif
}


#ifdef LIBMESH_USE_COMPLEX_NUMBERS
// Here we define the matrix assembly routine for
// the Helmholtz system.  This function will be
// called to form the stiffness matrix and right-hand side.
void assemble_helmholtz(EquationSystems & es,
                        const std::string & system_name)
{

  // It is a good idea to make sure we are assembling
  // the proper system.
  libmesh_assert_equal_to (system_name, "Helmholtz");

  // Get a constant reference to the mesh object.
  const MeshBase & mesh = es.get_mesh();

  // The maximum dimension of the elements stored in the mesh.
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to our system, as before
  FrequencySystem & f_system =
    es.get_system<FrequencySystem> (system_name);

  // A const reference to the DofMap object for this system.  The DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.
  const DofMap & dof_map = f_system.get_dof_map();

  // Get a constant reference to the finite element type
  // for the first (and only) variable in the system.
  const FEType & fe_type = dof_map.variable_type(0);

  // The fluid density is used by the admittance boundary condition.
  const Real rho = es.parameters.get<Real>("rho");

  // In this example, we assemble the element matrices into the
  // additional matrices "stiffness_mass", "damping", and the element
  // vectors into "rhs". We get writable references to these objects
  // now.
  SparseMatrix<Number> & stiffness = f_system.get_matrix("stiffness");
  SparseMatrix<Number> & damping = f_system.get_matrix("damping");
  SparseMatrix<Number> & mass = f_system.get_matrix("mass");

  // Build a finite element object of the specified type.  Since the
  // FEBase::build() member dynamically creates memory we will
  // store the object as a UniquePtr<FEBase>.  This can be thought
  // of as a pointer that will clean up after itself.
  UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));

  // A 5th-order Gauss quadrature rule for numerical integration.
  QGauss qrule (dim, FIFTH);

  // Tell the finite element object to use our quadrature rule.
  fe->attach_quadrature_rule (&qrule);

  // The element Jacobian times the quadrature weight at each integration point.
  const std::vector<Real> & JxW = fe->get_JxW();

  // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector<Real> > & phi = fe->get_phi();

  // The element shape function gradients evaluated at the quadrature
  // points.
  const std::vector<std::vector<RealGradient> > & dphi = fe->get_dphi();

  // Here we do not assemble directly in the System matrix, but to the
  // additional matrices "stiffness_mass" and "damping".  The same
  // approach is followed for the right-hand-side vector Fe, which we
  // will later on store in the additional vector "rhs".
  // The zero_matrix is used to explicitly induce the same sparsity
  // structure in the overall matrix.  The mass and stiffness matrices
  // are real-valued.  Therefore, it would be possible to store them
  // as a single complex-valued matrix to save on memory, but we do
  // not follow this approach here.
  DenseMatrix<Number> Ke, Ce, Me, zero_matrix;
  DenseVector<Number> Fe;

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;

  // Now we will loop over all the elements in the mesh, and compute
  // the element matrix and right-hand-side contributions.
  MeshBase::const_element_iterator           el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
    {
      // Start logging the element initialization.
      START_LOG("elem init", "assemble_helmholtz");

      // Store a pointer to the element we are currently
      // working on.  This allows for nicer syntax later.
      const Elem * elem = *el;

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
      STOP_LOG("elem init", "assemble_helmholtz");

      // Now loop over the quadrature points.  This handles
      // the numeric integration.
      START_LOG("stiffness & mass", "assemble_helmholtz");

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
          // Now we will build the element matrix.  This involves
          // a double loop to integrate the test funcions (i) against
          // the trial functions (j).
          for (std::size_t i=0; i<phi.size(); i++)
            for (std::size_t j=0; j<phi.size(); j++)
              {
                Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
                Me(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp]);
              }
        }

      STOP_LOG("stiffness & mass", "assemble_helmholtz");

      // Now compute the contribution to the element matrix
      // (due to mixed boundary conditions) if the current
      // element lies on the boundary.
      //
      // The following loops over the sides of the element.
      // If the element has no neighbor on a side then that
      // side MUST live on a boundary of the domain.
      for (unsigned int side=0; side<elem->n_sides(); side++)
        if (elem->neighbor(side) == libmesh_nullptr)
          {
            LOG_SCOPE("damping", "assemble_helmholtz");

            // Declare a special finite element object for
            // boundary integration.
            UniquePtr<FEBase> fe_face (FEBase::build(dim, fe_type));

            // Boundary integration requires one quadraure rule,
            // with dimensionality one less than the dimensionality
            // of the element.
            QGauss qface(dim-1, SECOND);

            // Tell the finte element object to use our
            // quadrature rule.
            fe_face->attach_quadrature_rule (&qface);

            // The value of the shape functions at the quadrature
            // points.
            const std::vector<std::vector<Real> > & phi_face =
              fe_face->get_phi();

            // The Jacobian times the quadrature weight at the quadrature
            // points on the face.
            const std::vector<Real> & JxW_face = fe_face->get_JxW();

            // Compute the shape function values on the element
            // face.
            fe_face->reinit(elem, side);

            // For the Robin BCs, consider a normal admittance an=1
            // at some parts of the bounfdary
            const Real an_value = 1.;

            // Loop over the face quadrature points for integration.
            for (unsigned int qp=0; qp<qface.n_points(); qp++)
              {
                // Element matrix contributrion due to prescribed
                // admittance boundary conditions.
                for (std::size_t i=0; i<phi_face.size(); i++)
                  for (std::size_t j=0; j<phi_face.size(); j++)
                    Ce(i,j) += rho * an_value * JxW_face[qp] * phi_face[i][qp] * phi_face[j][qp];
              }
          }

      // If this assembly program were to be used on an adaptive mesh,
      // we would have to apply any hanging node constraint equations
      // by uncommenting the following lines:
      // std::vector<unsigned int> dof_indicesC = dof_indices;
      // std::vector<unsigned int> dof_indicesM = dof_indices;
      // dof_map.constrain_element_matrix (Ke, dof_indices);
      // dof_map.constrain_element_matrix (Ce, dof_indicesC);
      // dof_map.constrain_element_matrix (Me, dof_indicesM);

      // Finally, simply add the contributions to the additional
      // matrices and vector.
      stiffness.add_matrix (Ke, dof_indices);
      damping.add_matrix (Ce, dof_indices);
      mass.add_matrix (Me, dof_indices);

      // For the overall matrix, explicitly zero the entries where
      // we added values in the other ones, so that we have
      // identical sparsity footprints.
      f_system.matrix->add_matrix(zero_matrix, dof_indices);
    } // end loop over elements
}


// We now define the function which will combine the previously-assembled
// mass, stiffness, and damping matrices into a single system matrix.
void add_M_C_K_helmholtz(EquationSystems & es,
                         const std::string & system_name)
{
  START_LOG("init phase", "add_M_C_K_helmholtz");

  // Verify that we are assembling the system we think we are.
  libmesh_assert_equal_to (system_name, "Helmholtz");

  // Get a reference to the FrequencySystem.
  FrequencySystem & f_system =
    es.get_system<FrequencySystem> (system_name);

  // Get the frequency, fluid density, and sound speed of the current solve.
  const Real frequency = es.parameters.get<Real> ("current frequency");
  const Real rho       = es.parameters.get<Real> ("rho");
  const Real speed     = es.parameters.get<Real> ("wave speed");

  // Compute angular frequency omega and wave number k
  const Real omega = 2.0*libMesh::pi*frequency;
  const Real k     = omega / speed;

  // Get writable references to the system matrix and vector, where the
  // frequency-dependent system is to be collected.
  SparseMatrix<Number> & matrix = *f_system.matrix;
  NumericVector<Number> & rhs = *f_system.rhs;

  // Get writable references to the frequency-independent matrices
  // and rhs, though we only need to extract values.  This write access
  // is necessary, since solver packages have to close the data structure
  // before they can extract values for computation.
  SparseMatrix<Number> & stiffness = f_system.get_matrix("stiffness");
  SparseMatrix<Number> & damping = f_system.get_matrix("damping");
  SparseMatrix<Number> & mass = f_system.get_matrix("mass");
  NumericVector<Number> & freq_indep_rhs = f_system.get_vector("rhs");

  // Compute the scale values for the different operators.
  const Number scale_stiffness (1., 0.);
  const Number scale_damping (0., omega);
  const Number scale_mass (-k*k, 0.);
  const Number scale_rhs (0., -(rho*omega));

  // Now simply add the matrices together, store the result
  // in matrix and rhs.  Clear them first.
  matrix.close ();
  matrix.zero ();
  rhs.close ();
  rhs.zero ();

  // The matrices from which values are added to another matrix
  // have to be closed.  The add() method does take care of
  // that, but let us do it explicitly.
  stiffness.close ();
  damping.close ();
  mass.close ();

  STOP_LOG("init phase", "add_M_C_K_helmholtz");

  START_LOG("global matrix & vector additions", "add_M_C_K_helmholtz");

  // Add the stiffness and mass with the proper frequency to the
  // overall system.  For this to work properly, the matrix has
  // to be not only initialized, but filled with the identical
  // sparsity structure as the matrix added to it, otherwise
  // solver packages like PETSc crash.
  //
  // Note that we have to add the mass and stiffness contributions one
  // at a time, otherwise the real part of matrix would be fine, but
  // the imaginary part would be cluttered with unwanted products.
  matrix.add (scale_stiffness, stiffness);
  matrix.add (scale_mass,      mass);
  matrix.add (scale_damping,   damping);
  rhs.add    (scale_rhs,       freq_indep_rhs);

  STOP_LOG("global matrix & vector additions", "add_M_C_K_helmholtz");

  // The linear system involving "matrix" and "rhs" is now ready to be solved.
}

#endif // LIBMESH_USE_COMPLEX_NUMBERS
