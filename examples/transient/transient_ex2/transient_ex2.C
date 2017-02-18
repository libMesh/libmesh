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



// <h1>Transient Example 2 - The Newmark System and the Wave Equation</h1>
// \author Steffen Petersen
// \date 2003
//
// This is the eighth example program. It builds on
// the previous example programs.  It introduces the
// NewmarkSystem class.  In this example the wave equation
// is solved using the time integration scheme provided
// by the NewmarkSystem class.
//
// This example comes with a cylindrical mesh given in the
// universal file pipe-mesh.unv.
// The mesh contains HEX8 and PRISM6 elements.

// C++ include files that we need
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <stdio.h>
#include <math.h>

// Basic include file needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/replicated_mesh.h"
#include "libmesh/gmv_io.h"
#include "libmesh/vtk_io.h"
#include "libmesh/newmark_system.h"
#include "libmesh/equation_systems.h"

// Define the Finite Element object.
#include "libmesh/fe.h"

// Define Gauss quadrature rules.
#include "libmesh/quadrature_gauss.h"

// Define useful datatypes for finite element
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"

// Define matrix and vector data types for the global
// equation system.  These are base classes,
// from which specific implementations, like
// the PETSc and Eigen implementations, are derived.
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"

// Define the DofMap, which handles degree of freedom
// indexing.
#include "libmesh/dof_map.h"

// The definition of a vertex associated with a Mesh.
#include "libmesh/node.h"

// The definition of a geometric element
#include "libmesh/elem.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Function prototype.  This is the function that will assemble
// the linear system for our problem, governed by the linear
// wave equation.
void assemble_wave(EquationSystems & es,
                   const std::string & system_name);


// Function Prototype.  This function will be used to apply the
// initial conditions.
void apply_initial(EquationSystems & es,
                   const std::string & system_name);

// Function Prototype.  This function imposes
// Dirichlet Boundary conditions via the penalty
// method after the system is assembled.
void fill_dirichlet_bc(EquationSystems & es,
                       const std::string & system_name);

// The main program
int main (int argc, char** argv)
{
  // Initialize libraries, like in example 2.
  LibMeshInit init (argc, argv);

  // Check for proper usage.
  if (argc < 2)
    libmesh_error_msg("Usage: " << argv[0] << " [meshfile]");

  // Tell the user what we are doing.
  else
    {
      libMesh::out << "Running " << argv[0];

      for (int i=1; i<argc; i++)
        libMesh::out << " " << argv[i];

      libMesh::out << std::endl << std::endl;

    }

  // LasPack solvers don't work so well for this example, Trilinos doesn't work at all.
  // PETSc and Eigen both work...
  libmesh_example_requires(libMesh::default_solver_package() == PETSC_SOLVERS || \
                           libMesh::default_solver_package() == EIGEN_SOLVERS, "--enable-petsc");

  // Get the name of the mesh file
  // from the command line.
  std::string mesh_file = argv[1];
  libMesh::out << "Mesh file is: " << mesh_file << std::endl;

  // Skip this 3D example if libMesh was compiled as 1D or 2D-only.
  libmesh_example_requires(3 <= LIBMESH_DIM, "3D support");

  // Create a mesh.
  // This example directly references all mesh nodes and is
  // incompatible with DistributedMesh use.
  //
  // Create a ReplicatedMesh object, with dimension to be overridden
  // later, distributed across the default MPI communicator.
  ReplicatedMesh mesh(init.comm());

  // Read the meshfile specified on the command line.
  mesh.read(mesh_file);

  // Print information about the mesh to the screen.
  mesh.print_info();

  // The node that should be monitored.
  const unsigned int result_node = 274;


  // Time stepping issues
  //
  // Note that the total current time is stored as a parameter
  // in the \pEquationSystems object.
  //
  // the time step size
  const Real delta_t = .0000625;

  // The number of time steps.
  unsigned int n_time_steps = 300;

  // Create an equation systems object.
  EquationSystems equation_systems (mesh);

  // Declare the system and its variables.
  // Create a NewmarkSystem named "Wave"
  equation_systems.add_system<NewmarkSystem> ("Wave");

  // Use a handy reference to this system
  NewmarkSystem & t_system = equation_systems.get_system<NewmarkSystem> ("Wave");

  // Add the variable "p" to "Wave".   "p"
  // will be approximated using first-order approximation.
  t_system.add_variable("p", FIRST);

  // Give the system a pointer to the matrix assembly
  // function and the initial condition function defined
  // below.
  t_system.attach_assemble_function  (assemble_wave);
  t_system.attach_init_function      (apply_initial);

  // Set the time step size, and optionally the
  // Newmark parameters, so that NewmarkSystem can
  // compute integration constants.  Here we simply use
  // pass only the time step and use default values
  // for alpha=.25  and delta=.5.
  t_system.set_newmark_parameters(delta_t);

  // Set the speed of sound and fluid density
  // as EquationSystems parameter,
  // so that assemble_wave() can access it.
  equation_systems.parameters.set<Real>("speed")          = 1000.;
  equation_systems.parameters.set<Real>("fluid density")  = 1000.;

  // Start time integration from t=0
  t_system.time = 0.;

  // Initialize the data structures for the equation system.
  equation_systems.init();

  // Prints information about the system to the screen.
  equation_systems.print_info();

  // A file to store the results at certain nodes.
  std::ofstream res_out("pressure_node.res");

  // get the dof_numbers for the nodes that
  // should be monitored.
  const unsigned int res_node_no = result_node;
  const Node & res_node = mesh.node_ref(res_node_no-1);
  unsigned int dof_no = res_node.dof_number(0, 0, 0);

  // Assemble the time independent system matrices and rhs.
  // This function will also compute the effective system matrix
  // K~=K+a_0*M+a_1*C and apply user specified initial
  // conditions.
  t_system.assemble();

  // Now solve for each time step.
  // For convenience, use a local buffer of the
  // current time.  But once this time is updated,
  // also update the EquationSystems parameter
  // Start with t_time = 0 and write a short header
  // to the nodal result file
  res_out << "# pressure at node " << res_node_no << "\n"
          << "# time\tpressure\n"
          << t_system.time << "\t" << 0 << std::endl;


  for (unsigned int time_step=0; time_step<n_time_steps; time_step++)
    {
      // Update the time.  Both here and in the
      // EquationSystems object
      t_system.time += delta_t;

      // Update the rhs.
      t_system.update_rhs();

      // Impose essential boundary conditions.
      // Not that since the matrix is only assembled once,
      // the penalty parameter should be added to the matrix
      // only in the first time step.  The applied
      // boundary conditions may be time-dependent and hence
      // the rhs vector is considered in each time step.
      if (time_step == 0)
        {
          // The local function fill_dirichlet_bc()
          // may also set Dirichlet boundary conditions for the
          // matrix.  When you set the flag as shown below,
          // the flag will return true.  If you want it to return
          // false, simply do not set it.
          equation_systems.parameters.set<bool>("Newmark set BC for Matrix") = true;

          fill_dirichlet_bc(equation_systems, "Wave");

          // unset the flag, so that it returns false
          equation_systems.parameters.set<bool>("Newmark set BC for Matrix") = false;
        }
      else
        fill_dirichlet_bc(equation_systems, "Wave");

      // Solve the system "Wave".
      t_system.solve();

      // After solving the system, write the solution
      // to a GMV-formatted plot file.
      // Do only for a few time steps.
      if (time_step == 30 || time_step == 60 ||
          time_step == 90 || time_step == 120)
        {
          std::ostringstream file_name;

#ifdef LIBMESH_HAVE_VTK
          file_name << "out_"
                    << std::setw(3)
                    << std::setfill('0')
                    << std::right
                    << time_step
                    << ".pvtu";

          VTKIO(mesh).write_equation_systems (file_name.str(), equation_systems);
#else

          file_name << "out."
                    << std::setw(3)
                    << std::setfill('0')
                    << std::right
                    << time_step
                    << ".gmv";

          GMVIO(mesh).write_equation_systems (file_name.str(), equation_systems);
#endif
        }

      // Update the p, v and a.
      t_system.update_u_v_a();

      // dof_no may not be local in parallel runs, so we may need a
      // global displacement vector
      NumericVector<Number> & displacement = t_system.get_vector("displacement");
      std::vector<Number> global_displacement(displacement.size());
      displacement.localize(global_displacement);

      // Write nodal results to file.  The results can then
      // be viewed with e.g. gnuplot (run gnuplot and type
      // 'plot "pressure_node.res" with lines' in the command line)
      res_out << t_system.time << "\t"
              << global_displacement[dof_no]
              << std::endl;
    }

  // All done.
  return 0;
}

// This function assembles the system matrix and right-hand-side
// for our wave equation.
void assemble_wave(EquationSystems & es,
                   const std::string & system_name)
{
  // It is a good idea to make sure we are assembling
  // the proper system.
  libmesh_assert_equal_to (system_name, "Wave");

  // Get a constant reference to the mesh object.
  const MeshBase & mesh = es.get_mesh();

  // The dimension that we are running.
  const unsigned int dim = mesh.mesh_dimension();

  // Copy the speed of sound and fluid density
  // to a local variable.
  const Real speed = es.parameters.get<Real>("speed");
  const Real rho   = es.parameters.get<Real>("fluid density");

  // Get a reference to our system, as before.
  NewmarkSystem & t_system = es.get_system<NewmarkSystem> (system_name);

  // Get a constant reference to the Finite Element type
  // for the first (and only) variable in the system.
  FEType fe_type = t_system.get_dof_map().variable_type(0);

  // In here, we will add the element matrices to the
  // @e additional matrices "stiffness_mass" and "damping"
  // and the additional vector "force", not to the members
  // "matrix" and "rhs".  Therefore, get writable
  // references to them.
  SparseMatrix<Number> & stiffness = t_system.get_matrix("stiffness");
  SparseMatrix<Number> & damping   = t_system.get_matrix("damping");
  SparseMatrix<Number> & mass      = t_system.get_matrix("mass");
  NumericVector<Number> & force    = t_system.get_vector("force");

  // Some solver packages (PETSc) are especially picky about
  // allocating sparsity structure and truly assigning values
  // to this structure.  Namely, matrix additions, as performed
  // later, exhibit acceptable performance only for identical
  // sparsity structures.  Therefore, explicitly zero the
  // values in the collective matrix, so that matrix additions
  // encounter identical sparsity structures.
  SparseMatrix<Number> & matrix = *t_system.matrix;
  DenseMatrix<Number> zero_matrix;

  // Build a Finite Element object of the specified type.  Since the
  // FEBase::build() member dynamically creates memory we will
  // store the object as a UniquePtr<FEBase>.  This can be thought
  // of as a pointer that will clean up after itself.
  UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));

  // A 2nd order Gauss quadrature rule for numerical integration.
  QGauss qrule (dim, SECOND);

  // Tell the finite element object to use our quadrature rule.
  fe->attach_quadrature_rule (&qrule);

  // The element Jacobian * quadrature weight at each integration point.
  const std::vector<Real> & JxW = fe->get_JxW();

  // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector<Real> > & phi = fe->get_phi();

  // The element shape function gradients evaluated at the quadrature
  // points.
  const std::vector<std::vector<RealGradient> > & dphi = fe->get_dphi();

  // A reference to the DofMap object for this system.  The DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.
  const DofMap & dof_map = t_system.get_dof_map();

  // The element mass, damping and stiffness matrices
  // and the element contribution to the rhs.
  DenseMatrix<Number> Ke, Ce, Me;
  DenseVector<Number> Fe;

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;

  // Now we will loop over all the elements in the mesh.
  // We will compute the element matrix and right-hand-side
  // contribution.
  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
    {
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

      // Zero the element matrices and rhs before
      // summing them.  We use the resize member here because
      // the number of degrees of freedom might have changed from
      // the last element.  Note that this will be the case if the
      // element type is different (i.e. the last element was HEX8
      // and now have a PRISM6).
      {
        const unsigned int n_dof_indices = dof_indices.size();

        Ke.resize          (n_dof_indices, n_dof_indices);
        Ce.resize          (n_dof_indices, n_dof_indices);
        Me.resize          (n_dof_indices, n_dof_indices);
        zero_matrix.resize (n_dof_indices, n_dof_indices);
        Fe.resize          (n_dof_indices);
      }

      // Now loop over the quadrature points.  This handles
      // the numeric integration.
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
          // Now we will build the element matrix.  This involves
          // a double loop to integrate the test funcions (i) against
          // the trial functions (j).
          for (std::size_t i=0; i<phi.size(); i++)
            for (std::size_t j=0; j<phi.size(); j++)
              {
                Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
                Me(i,j) += JxW[qp]*phi[i][qp]*phi[j][qp]
                  *1./(speed*speed);
              } // end of the matrix summation loop
        } // end of quadrature point loop

      // Now compute the contribution to the element matrix and the
      // right-hand-side vector if the current element lies on the
      // boundary.
      {
        // In this example no natural boundary conditions will
        // be considered.  The code is left here so it can easily
        // be extended.
        //
        // don't do this for any side
        for (unsigned int side=0; side<elem->n_sides(); side++)
          if (!true)
            // if (elem->neighbor_ptr(side) == libmesh_nullptr)
            {
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
              const std::vector<std::vector<Real> > &  phi_face = fe_face->get_phi();

              // The Jacobian * Quadrature Weight at the quadrature
              // points on the face.
              const std::vector<Real> & JxW_face = fe_face->get_JxW();

              // Compute the shape function values on the element
              // face.
              fe_face->reinit(elem, side);

              // Here we consider a normal acceleration acc_n=1 applied to
              // the whole boundary of our mesh.
              const Real acc_n_value = 1.0;

              // Loop over the face quadrature points for integration.
              for (unsigned int qp=0; qp<qface.n_points(); qp++)
                {
                  // Right-hand-side contribution due to prescribed
                  // normal acceleration.
                  for (std::size_t i=0; i<phi_face.size(); i++)
                    {
                      Fe(i) += acc_n_value*rho
                        *phi_face[i][qp]*JxW_face[qp];
                    }
                } // end face quadrature point loop
            } // end if (elem->neighbor_ptr(side) == libmesh_nullptr)

        // In this example the Dirichlet boundary conditions will be
        // imposed via panalty method after the
        // system is assembled.

      } // end boundary condition section

      // If this assembly program were to be used on an adaptive mesh,
      // we would have to apply any hanging node constraint equations
      // by uncommenting the following lines:
      // std::vector<unsigned int> dof_indicesC = dof_indices;
      // std::vector<unsigned int> dof_indicesM = dof_indices;
      // dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);
      // dof_map.constrain_element_matrix (Ce, dof_indicesC);
      // dof_map.constrain_element_matrix (Me, dof_indicesM);

      // Finally, simply add the contributions to the additional
      // matrices and vector.
      stiffness.add_matrix (Ke, dof_indices);
      damping.add_matrix   (Ce, dof_indices);
      mass.add_matrix      (Me, dof_indices);

      force.add_vector     (Fe, dof_indices);

      // For the overall matrix, explicitly zero the entries where
      // we added values in the other ones, so that we have
      // identical sparsity footprints.
      matrix.add_matrix(zero_matrix, dof_indices);

    } // end of element loop
}

// This function applies the initial conditions
void apply_initial(EquationSystems & es,
                   const std::string & system_name)
{
  // Get a reference to our system, as before
  NewmarkSystem & t_system = es.get_system<NewmarkSystem> (system_name);

  // Numeric vectors for the pressure, velocity and acceleration
  // values.
  NumericVector<Number> & pres_vec = t_system.get_vector("displacement");
  NumericVector<Number> & vel_vec  = t_system.get_vector("velocity");
  NumericVector<Number> & acc_vec  = t_system.get_vector("acceleration");

  // Assume our fluid to be at rest, which would
  // also be the default conditions in class NewmarkSystem,
  // but let us do it explicetly here.
  pres_vec.zero();
  vel_vec.zero();
  acc_vec.zero();
}

// This function applies the Dirichlet boundary conditions
void fill_dirichlet_bc(EquationSystems & es,
                       const std::string & system_name)
{
  // It is a good idea to make sure we are assembling
  // the proper system.
  libmesh_assert_equal_to (system_name, "Wave");

  // Get a reference to our system, as before.
  NewmarkSystem & t_system = es.get_system<NewmarkSystem> (system_name);

  // Get writable references to the overall matrix and vector.
  SparseMatrix<Number> & matrix = *t_system.matrix;
  NumericVector<Number> & rhs    = *t_system.rhs;

  // Get a constant reference to the mesh object.
  const MeshBase & mesh = es.get_mesh();

  // Get libMesh's pi
  const Real pi = libMesh::pi;

  // Ask the EquationSystems flag whether
  // we should do this also for the matrix
  const bool do_for_matrix =
    es.parameters.get<bool>("Newmark set BC for Matrix");

  // Number of nodes in the mesh.
  unsigned int n_nodes = mesh.n_nodes();

  for (unsigned int n_cnt=0; n_cnt<n_nodes; n_cnt++)
    {
      // Get a reference to the current node.
      const Node & curr_node = mesh.node_ref(n_cnt);

      // Check if Dirichlet BCs should be applied to this node.
      // Use the TOLERANCE from mesh_common.h as tolerance.
      // Here a pressure value is applied if the z-coord.
      // is equal to 4, which corresponds to one end of the
      // pipe-mesh in this directory.
      const Real z_coo = 4.;

      if (std::abs(curr_node(2)-z_coo) < TOLERANCE)
        {
          // The global number of the respective degree of freedom.
          unsigned int dn = curr_node.dof_number(0, 0, 0);

          // The penalty parameter.
          const Real penalty = 1.e10;

          // Here we apply sinusoidal pressure values for 0<t<0.002
          // at one end of the pipe-mesh.
          Real p_value;
          if (t_system.time < .002)
            p_value = sin(2*pi*t_system.time/.002);
          else
            p_value = .0;

          // Now add the contributions to the matrix and the rhs.
          rhs.add(dn, p_value*penalty);

          // Add the panalty parameter to the global matrix
          // if desired.
          if (do_for_matrix)
            matrix.add(dn, dn, penalty);
        }
    } // loop n_cnt
}
