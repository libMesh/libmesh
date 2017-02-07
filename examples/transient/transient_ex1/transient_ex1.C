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



// <h1>Transient Example 1 - Solving a Transient Linear System in Parallel</h1>
// \author Benjamin S. Kirk
// \date 2003
//
// This example shows how a simple, linear transient
// system can be solved in parallel.  The system is simple
// scalar convection-diffusion with a specified external
// velocity.  The initial condition is given, and the
// solution is advanced in time with a standard Crank-Nicolson
// time-stepping strategy.

// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <sstream>
#include <math.h>

// Basic include file needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/gmv_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/exodusII_io.h"

// This example will solve a linear transient system,
// so we need to include the TransientLinearImplicitSystem definition.
#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/vector_value.h"

// The definition of a geometric element
#include "libmesh/elem.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Function prototype.  This function will assemble the system
// matrix and right-hand-side at each time step.  Note that
// since the system is linear we technically do not need to
// assmeble the matrix at each time step, but we will anyway.
// In subsequent examples we will employ adaptive mesh refinement,
// and with a changing mesh it will be necessary to rebuild the
// system matrix.
void assemble_cd (EquationSystems & es,
                  const std::string & system_name);

// Function prototype.  This function will initialize the system.
// Initialization functions are optional for systems.  They allow
// you to specify the initial values of the solution.  If an
// initialization function is not provided then the default (0)
// solution is provided.
void init_cd (EquationSystems & es,
              const std::string & system_name);

// Exact solution function prototype.  This gives the exact
// solution as a function of space and time.  In this case the
// initial condition will be taken as the exact solution at time 0,
// as will the Dirichlet boundary conditions at time t.
Real exact_solution (const Real x,
                     const Real y,
                     const Real t);

Number exact_value (const Point & p,
                    const Parameters & parameters,
                    const std::string &,
                    const std::string &)
{
  return exact_solution(p(0), p(1), parameters.get<Real> ("time"));
}



// We can now begin the main program.  Note that this
// example will fail if you are using complex numbers
// since it was designed to be run only with real numbers.
int main (int argc, char ** argv)
{
  // Initialize libMesh.
  LibMeshInit init (argc, argv);

  // This example requires Adaptive Mesh Refinement support - although
  // it only refines uniformly, the refinement code used is the same
  // underneath
#ifndef LIBMESH_ENABLE_AMR
  libmesh_example_requires(false, "--enable-amr");
#else

  // Skip this 2D example if libMesh was compiled as 1D-only.
  libmesh_example_requires(2 <= LIBMESH_DIM, "2D support");

  // Read the mesh from file.  This is the coarse mesh that will be used
  // in example 10 to demonstrate adaptive mesh refinement.  Here we will
  // simply read it in and uniformly refine it 5 times before we compute
  // with it.
  //
  // Create a mesh object, with dimension to be overridden later,
  // distributed across the default MPI communicator.
  Mesh mesh(init.comm());

  mesh.read ("mesh.xda");

  // Create a MeshRefinement object to handle refinement of our mesh.
  // This class handles all the details of mesh refinement and coarsening.
  MeshRefinement mesh_refinement (mesh);

  // Uniformly refine the mesh 5 times.  This is the
  // first time we use the mesh refinement capabilities
  // of the library.
  mesh_refinement.uniformly_refine (5);

  // Print information about the mesh to the screen.
  mesh.print_info();

  // Create an equation systems object.
  EquationSystems equation_systems (mesh);

  // Add a transient system to the EquationSystems
  // object named "Convection-Diffusion".
  TransientLinearImplicitSystem & system =
    equation_systems.add_system<TransientLinearImplicitSystem> ("Convection-Diffusion");

  // Adds the variable "u" to "Convection-Diffusion".  "u"
  // will be approximated using first-order approximation.
  system.add_variable ("u", FIRST);

  // Give the system a pointer to the matrix assembly
  // and initialization functions.
  system.attach_assemble_function (assemble_cd);
  system.attach_init_function (init_cd);

  // Initialize the data structures for the equation system.
  equation_systems.init ();

  // Prints information about the system to the screen.
  equation_systems.print_info();

  // Write out the initial conditions.
#ifdef LIBMESH_HAVE_EXODUS_API
  // If Exodus is available, we'll write all timesteps to the same file
  // rather than one file per timestep.
  std::string exodus_filename = "transient_ex1.e";
  ExodusII_IO(mesh).write_equation_systems (exodus_filename, equation_systems);
#else
  GMVIO(mesh).write_equation_systems ("out_000.gmv", equation_systems);
#endif

  // The Convection-Diffusion system requires that we specify
  // the flow velocity.  We will specify it as a RealVectorValue
  // data type and then use the Parameters object to pass it to
  // the assemble function.
  equation_systems.parameters.set<RealVectorValue>("velocity") =
    RealVectorValue (0.8, 0.8);

  // Solve the system "Convection-Diffusion".  This will be done by
  // looping over the specified time interval and calling the
  // solve() member at each time step.  This will assemble the
  // system and call the linear solver.
  const Real dt = 0.025;
  system.time = 0.;

  for (unsigned int t_step = 0; t_step < 50; t_step++)
    {
      // Incremenet the time counter, set the time and the
      // time step size as parameters in the EquationSystem.
      system.time += dt;

      equation_systems.parameters.set<Real> ("time") = system.time;
      equation_systems.parameters.set<Real> ("dt")   = dt;

      // A pretty update message
      libMesh::out << " Solving time step ";

      // Do fancy zero-padded formatting of the current time.
      {
        std::ostringstream out;

        out << std::setw(2)
            << std::right
            << t_step
            << ", time="
            << std::fixed
            << std::setw(6)
            << std::setprecision(3)
            << std::setfill('0')
            << std::left
            << system.time
            <<  "...";

        libMesh::out << out.str() << std::endl;
      }

      // At this point we need to update the old
      // solution vector.  The old solution vector
      // will be the current solution vector from the
      // previous time step.  We will do this by extracting the
      // system from the EquationSystems object and using
      // vector assignment.  Since only TransientSystems
      // (and systems derived from them) contain old solutions
      // we need to specify the system type when we ask for it.
      *system.old_local_solution = *system.current_local_solution;

      // Assemble & solve the linear system
      equation_systems.get_system("Convection-Diffusion").solve();

      // Output evey 10 timesteps to file.
      if ((t_step+1)%10 == 0)
        {

#ifdef LIBMESH_HAVE_EXODUS_API
          ExodusII_IO exo(mesh);
          exo.append(true);
          exo.write_timestep (exodus_filename, equation_systems, t_step+1, system.time);
#else
          std::ostringstream file_name;

          file_name << "out_"
                    << std::setw(3)
                    << std::setfill('0')
                    << std::right
                    << t_step+1
                    << ".gmv";


          GMVIO(mesh).write_equation_systems (file_name.str(),
                                              equation_systems);
#endif
        }
    }
#endif // #ifdef LIBMESH_ENABLE_AMR

  // All done.
  return 0;
}

// We now define the function which provides the
// initialization routines for the "Convection-Diffusion"
// system.  This handles things like setting initial
// conditions and boundary conditions.
void init_cd (EquationSystems & es,
              const std::string & libmesh_dbg_var(system_name))
{
  // It is a good idea to make sure we are initializing
  // the proper system.
  libmesh_assert_equal_to (system_name, "Convection-Diffusion");

  // Get a reference to the Convection-Diffusion system object.
  TransientLinearImplicitSystem & system =
    es.get_system<TransientLinearImplicitSystem>("Convection-Diffusion");

  // Project initial conditions at time 0
  es.parameters.set<Real> ("time") = system.time = 0;

  system.project_solution(exact_value, libmesh_nullptr, es.parameters);
}



// Now we define the assemble function which will be used
// by the EquationSystems object at each timestep to assemble
// the linear system for solution.
void assemble_cd (EquationSystems & es,
                  const std::string & libmesh_dbg_var(system_name))
{
#ifdef LIBMESH_ENABLE_AMR
  // It is a good idea to make sure we are assembling
  // the proper system.
  libmesh_assert_equal_to (system_name, "Convection-Diffusion");

  // Get a constant reference to the mesh object.
  const MeshBase & mesh = es.get_mesh();

  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to the Convection-Diffusion system object.
  TransientLinearImplicitSystem & system =
    es.get_system<TransientLinearImplicitSystem> ("Convection-Diffusion");

  // Get a constant reference to the Finite Element type
  // for the first (and only) variable in the system.
  FEType fe_type = system.variable_type(0);

  // Build a Finite Element object of the specified type.  Since the
  // FEBase::build() member dynamically creates memory we will
  // store the object as a UniquePtr<FEBase>.  This can be thought
  // of as a pointer that will clean up after itself.
  UniquePtr<FEBase> fe      (FEBase::build(dim, fe_type));
  UniquePtr<FEBase> fe_face (FEBase::build(dim, fe_type));

  // A Gauss quadrature rule for numerical integration.
  // Let the FEType object decide what order rule is appropriate.
  QGauss qrule (dim,   fe_type.default_quadrature_order());
  QGauss qface (dim-1, fe_type.default_quadrature_order());

  // Tell the finite element object to use our quadrature rule.
  fe->attach_quadrature_rule      (&qrule);
  fe_face->attach_quadrature_rule (&qface);

  // Here we define some references to cell-specific data that
  // will be used to assemble the linear system.  We will start
  // with the element Jacobian * quadrature weight at each integration point.
  const std::vector<Real> & JxW      = fe->get_JxW();
  const std::vector<Real> & JxW_face = fe_face->get_JxW();

  // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector<Real> > & phi = fe->get_phi();
  const std::vector<std::vector<Real> > & psi = fe_face->get_phi();

  // The element shape function gradients evaluated at the quadrature
  // points.
  const std::vector<std::vector<RealGradient> > & dphi = fe->get_dphi();

  // The XY locations of the quadrature points used for face integration
  const std::vector<Point> & qface_points = fe_face->get_xyz();

  // A reference to the DofMap object for this system.  The DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.  We will talk more about the DofMap
  // in future examples.
  const DofMap & dof_map = system.get_dof_map();

  // Define data structures to contain the element matrix
  // and right-hand-side vector contribution.  Following
  // basic finite element terminology we will denote these
  // "Ke" and "Fe".
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;

  // Here we extract the velocity & parameters that we put in the
  // EquationSystems object.
  const RealVectorValue velocity =
    es.parameters.get<RealVectorValue> ("velocity");

  const Real dt = es.parameters.get<Real>   ("dt");

  // Now we will loop over all the elements in the mesh that
  // live on the local processor. We will compute the element
  // matrix and right-hand-side contribution.  Since the mesh
  // will be refined we want to only consider the ACTIVE elements,
  // hence we use a variant of the active_elem_iterator.
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

      // Zero the element matrix and right-hand side before
      // summing them.  We use the resize member here because
      // the number of degrees of freedom might have changed from
      // the last element.  Note that this will be the case if the
      // element type is different (i.e. the last element was a
      // triangle, now we are on a quadrilateral).
      Ke.resize (dof_indices.size(),
                 dof_indices.size());

      Fe.resize (dof_indices.size());

      // Now we will build the element matrix and right-hand-side.
      // Constructing the RHS requires the solution and its
      // gradient from the previous timestep.  This myst be
      // calculated at each quadrature point by summing the
      // solution degree-of-freedom values by the appropriate
      // weight functions.
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
          // Values to hold the old solution & its gradient.
          Number u_old = 0.;
          Gradient grad_u_old;

          // Compute the old solution & its gradient.
          for (std::size_t l=0; l<phi.size(); l++)
            {
              u_old += phi[l][qp]*system.old_solution  (dof_indices[l]);

              // This will work,
              // grad_u_old += dphi[l][qp]*system.old_solution (dof_indices[l]);
              // but we can do it without creating a temporary like this:
              grad_u_old.add_scaled (dphi[l][qp], system.old_solution (dof_indices[l]));
            }

          // Now compute the element matrix and RHS contributions.
          for (std::size_t i=0; i<phi.size(); i++)
            {
              // The RHS contribution
              Fe(i) += JxW[qp]*(
                                // Mass matrix term
                                u_old*phi[i][qp] +
                                -.5*dt*(
                                        // Convection term
                                        // (grad_u_old may be complex, so the
                                        // order here is important!)
                                        (grad_u_old*velocity)*phi[i][qp] +

                                        // Diffusion term
                                        0.01*(grad_u_old*dphi[i][qp]))
                                );

              for (std::size_t j=0; j<phi.size(); j++)
                {
                  // The matrix contribution
                  Ke(i,j) += JxW[qp]*(
                                      // Mass-matrix
                                      phi[i][qp]*phi[j][qp] +

                                      .5*dt*(
                                             // Convection term
                                             (velocity*dphi[j][qp])*phi[i][qp] +

                                             // Diffusion term
                                             0.01*(dphi[i][qp]*dphi[j][qp]))
                                      );
                }
            }
        }

      // At this point the interior element integration has
      // been completed.  However, we have not yet addressed
      // boundary conditions.  For this example we will only
      // consider simple Dirichlet boundary conditions imposed
      // via the penalty method.
      //
      // The following loops over the sides of the element.
      // If the element has no neighbor on a side then that
      // side MUST live on a boundary of the domain.
      {
        // The penalty value.
        const Real penalty = 1.e10;

        // The following loops over the sides of the element.
        // If the element has no neighbor on a side then that
        // side MUST live on a boundary of the domain.
        for (unsigned int s=0; s<elem->n_sides(); s++)
          if (elem->neighbor_ptr(s) == libmesh_nullptr)
            {
              fe_face->reinit(elem, s);

              for (unsigned int qp=0; qp<qface.n_points(); qp++)
                {
                  const Number value = exact_solution (qface_points[qp](0),
                                                       qface_points[qp](1),
                                                       system.time);

                  // RHS contribution
                  for (std::size_t i=0; i<psi.size(); i++)
                    Fe(i) += penalty*JxW_face[qp]*value*psi[i][qp];

                  // Matrix contribution
                  for (std::size_t i=0; i<psi.size(); i++)
                    for (std::size_t j=0; j<psi.size(); j++)
                      Ke(i,j) += penalty*JxW_face[qp]*psi[i][qp]*psi[j][qp];
                }
            }
      }

      // If this assembly program were to be used on an adaptive mesh,
      // we would have to apply any hanging node constraint equations
      dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

      // The element matrix and right-hand-side are now built
      // for this element.  Add them to the global matrix and
      // right-hand-side vector.  The SparseMatrix::add_matrix()
      // and NumericVector::add_vector() members do this for us.
      system.matrix->add_matrix (Ke, dof_indices);
      system.rhs->add_vector    (Fe, dof_indices);
    }

  // That concludes the system matrix assembly routine.
#endif // #ifdef LIBMESH_ENABLE_AMR
}
