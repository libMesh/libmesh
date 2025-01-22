// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// rbOOmit: An implementation of the Certified Reduced Basis method.
// Copyright (C) 2009, 2010 David J. Knezevic
// This file is part of rbOOmit.



// <h1>Reduced Basis Example 4 - Empirical Interpolation Method</h1>
// \author David Knezevic
// \date 2011
//
// In this example problem we develop a reduced basis approximation
// for a parametrized PDE that has "non-affine" parameter
// dependence. This requires the use of the Empirical Interpolation
// Method (EIM).
//
// We first use EIM to construct an affine approximation to the
// non-affine term, which is a parametrized function that is a
// Gaussian with "center" defined by the two parameters (mu_1, mu_2)
// \in [-1,1]^2. We then employ this EIM approximation in order to
// generate a reduced basis approximation for the parametrized PDE:
// -0.05 * Laplacian(u) = f(mu_1, mu_2), with zero Dirichlet boundary
// conditions.

// Basic include file needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/equation_systems.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/getpot.h"
#include "libmesh/rb_data_serialization.h"
#include "libmesh/rb_data_deserialization.h"
#include "libmesh/enum_solver_package.h"

#include "eim_classes.h"
#include "rb_classes.h"

using namespace libMesh;


int main (int argc, char ** argv)
{
  // Initialize libMesh.
  LibMeshInit init (argc, argv);

  // This example requires a linear solver package.
  libmesh_example_requires(libMesh::default_solver_package() != INVALID_SOLVER_PACKAGE,
                           "--enable-petsc, --enable-trilinos, or --enable-eigen");

#if !defined(LIBMESH_HAVE_XDR)
  // We need XDR support to write out reduced bases
  libmesh_example_requires(false, "--enable-xdr");
#elif defined(LIBMESH_DEFAULT_SINGLE_PRECISION)
  // XDR binary support requires double precision
  libmesh_example_requires(false, "double precision");
#elif defined(LIBMESH_DEFAULT_TRIPLE_PRECISION)
  // I have no idea why long double isn't working here... [RHS]
  libmesh_example_requires(false, "double precision");
#endif

  // Skip this 2D example if libMesh was compiled as 1D-only.
  libmesh_example_requires(2 <= LIBMESH_DIM, "2D support");

  // This example requires libmesh to be configured with both
  // DirichletBoundary and Cap'n Proto support.
#if !defined(LIBMESH_ENABLE_DIRICHLET) || !defined(LIBMESH_HAVE_CAPNPROTO)
  libmesh_example_requires(false, "--enable-dirichlet --enable-capnp");
#else

  // Define the names of the input files we will read the problem properties from
  std::string eim_parameters = "eim.in";
  std::string rb_parameters  = "rb.in";
  std::string main_parameters = "reduced_basis_ex4.in";
  GetPot infile(main_parameters);

  unsigned int n_elem = infile("n_elem", 1);       // Determines the number of elements in the "truth" mesh
  const unsigned int dim = 2;                      // The number of spatial dimensions

  // Read the "online_mode" flag from the command line
  const int online_mode =
    libMesh::command_line_next("-online_mode", 0),
            eim_training_function_to_plot =
    libMesh::command_line_next("-eim_training_function_to_plot", 0),
            eim_basis_function_to_plot =
    libMesh::command_line_next("-eim_basis_function_to_plot", 0);

  ReplicatedMesh mesh (init.comm(), dim);
  MeshTools::Generation::build_square (mesh,
                                       n_elem, n_elem,
                                       -1., 1.,
                                       -1., 1.,
                                       QUAD4);

  // Set to true the write binary (XDR) files with EIM data, false to write ASCII files.
  bool eim_binary_io = true;

  if (!online_mode)
    {
      // First run the Offline stage for the EIM approximation.
      {
        libMesh::out << "Perform EIM training" << std::endl << std::endl;

        // Initialize the EquationSystems object for this mesh and attach
        // the EIM and RB Construction objects
        EquationSystems equation_systems (mesh);

        SimpleEIMConstruction & eim_construction =
          equation_systems.add_system<SimpleEIMConstruction> ("EIM");

        // Initialize the data structures for the equation system.
        equation_systems.init ();

        // Print out some information about the "truth" discretization
        mesh.print_info();
        equation_systems.print_info();

        // Initialize the EIM RBEvaluation object
        SimpleEIMEvaluation eim_rb_eval(mesh.comm());

        // Set the rb_eval objects for the RBConstructions
        eim_construction.set_rb_eim_evaluation(eim_rb_eval);

        // Read data from input file and print state
        eim_construction.process_parameters_file(eim_parameters);
        eim_construction.print_info();

        // Perform the EIM Greedy and write out the data
        eim_construction.initialize_eim_construction();
        eim_construction.train_eim_approximation();

        RBDataSerialization::RBEIMEvaluationSerialization rb_eim_eval_writer(eim_rb_eval);
        rb_eim_eval_writer.write_to_file("rb_eim_eval.bin");

        // Write out the EIM basis functions
        eim_rb_eval.write_out_basis_functions("eim_data", eim_binary_io);
      }

      {
        libMesh::out << std::endl << "Perform RB training" << std::endl << std::endl;

        // Initialize the EquationSystems object for this mesh and attach
        // the EIM and RB Construction objects
        EquationSystems equation_systems (mesh);

        SimpleEIMConstruction & eim_construction =
          equation_systems.add_system<SimpleEIMConstruction> ("EIM");
        SimpleRBConstruction & rb_construction =
          equation_systems.add_system<SimpleRBConstruction> ("RB");

        // Initialize the data structures for the equation system.
        equation_systems.init ();

        // Print out some information about the "truth" discretization
        mesh.print_info();
        equation_systems.print_info();

        // Initialize the standard RBEvaluation object
        SimpleRBEvaluation rb_eval(mesh.comm());

        // Initialize the EIM RBEvaluation object
        SimpleEIMEvaluation eim_rb_eval(mesh.comm());

        // Set the rb_eval objects for the RBConstructions
        eim_construction.set_rb_eim_evaluation(eim_rb_eval);
        rb_construction.set_rb_evaluation(rb_eval);

        RBDataDeserialization::RBEIMEvaluationDeserialization rb_eim_eval_reader(eim_rb_eval);
        rb_eim_eval_reader.read_from_file("rb_eim_eval.bin");
        eim_rb_eval.read_in_basis_functions(rb_construction, "eim_data", eim_binary_io);

        // Read data from input file and print state
        rb_construction.process_parameters_file(rb_parameters);

        // attach the EIM theta objects to the RBConstruction and RBEvaluation objects
        eim_rb_eval.initialize_eim_theta_objects();
        rb_eval.get_rb_theta_expansion().attach_multiple_F_theta(eim_rb_eval.get_eim_theta_objects());

        // attach the EIM assembly objects to the RBConstruction object
        eim_construction.initialize_eim_assembly_objects();
        rb_construction.get_rb_assembly_expansion().attach_multiple_F_assembly(eim_construction.get_eim_assembly_objects());

        // Print out the state of rb_construction now that the EIM objects have been attached
        rb_construction.print_info();

        // Need to initialize _after_ EIM greedy so that
        // the system knows how many affine terms there are
        rb_construction.initialize_rb_construction();
        rb_construction.train_reduced_basis();

        RBDataSerialization::RBEvaluationSerialization rb_eval_writer(rb_construction.get_rb_evaluation());
        rb_eval_writer.write_to_file("rb_eval.bin");

        rb_construction.get_rb_evaluation().write_out_basis_functions(rb_construction,
                                                                      "rb_data");
      }
    }
  else // online mode
    {
      SimpleEIMEvaluation eim_rb_eval(mesh.comm());
      SimpleRBEvaluation rb_eval(mesh.comm());

      RBDataDeserialization::RBEIMEvaluationDeserialization rb_eim_eval_reader(eim_rb_eval);
      rb_eim_eval_reader.read_from_file("rb_eim_eval.bin");

      // attach the EIM theta objects to rb_eval objects
      eim_rb_eval.initialize_eim_theta_objects();
      rb_eval.get_rb_theta_expansion().attach_multiple_F_theta(eim_rb_eval.get_eim_theta_objects());

      // Read in the offline data for rb_eval
      RBDataDeserialization::RBEvaluationDeserialization rb_eval_reader(rb_eval);
      rb_eval_reader.read_from_file("rb_eval.bin", /*read_error_bound_data*/ true);

      // Get the parameters at which we will do a reduced basis solve
      Real online_center_x = infile("online_center_x", 0.);
      Real online_center_y = infile("online_center_y", 0.);
      RBParameters online_mu;
      online_mu.push_back_value("center_x", online_center_x);
      online_mu.push_back_value("center_y", online_center_y);

      // Testing: Add secondary center_x and center_y values,
      // corresponding to e.g. a different time step or load step.
      online_mu.push_back_value("center_x", 0.5);
      online_mu.push_back_value("center_y", 0.5);

      // Add 3rd (center_x, center_y) values. For debugging purposes,
      // we want the number of "samples" (3) to be different from the
      // number of parameters (2).
      online_mu.push_back_value("center_x", -0.25);
      online_mu.push_back_value("center_y", -0.25);

      // Here we are going to pre-evaluate the thetas and pass them in
      // to rb_solve() in a loop.  Note that the rb_solve doesn't
      // explicitly need the parameters ("mus"), it just needs the
      // evaluated functions of parameters ("thetas"), but we can
      // perform an rb_solve() with either.
      rb_eval.set_parameters(online_mu);
      rb_eval.print_parameters();

      // When performing an rb_solve() with "mu" values, we only
      // support single-valued RBParameters objects.  In this case,
      // since the RBParameters object stores multiple "samples", we
      // take the approach of pre-evaluating the thetas for each
      // parameter while calling the rb_solve() in a loop.
      const RBThetaExpansion & rb_theta_expansion = rb_eval.get_rb_theta_expansion();

      // Single-entry vector which makes calling the vector-overrides
      // of eval_A_theta(), eval_F_theta(), etc. easier.
      std::vector<RBParameters> mu_vec = {online_mu};

      // 1.) Evaluate and store "A" thetas at all samples:
      std::vector<std::vector<Number>> all_A(rb_theta_expansion.get_n_A_terms());
      for (unsigned int q_a=0; q_a<rb_theta_expansion.get_n_A_terms(); q_a++)
        {
          // Here we call the version of eval_A_theta() taking a
          // vector, so we evaluate theta(mu) at all samples
          // simultaneously.
          all_A[q_a] = rb_theta_expansion.eval_A_theta(q_a, mu_vec);

          // This size of each A_q vector here is:
          // sum_i mu_vec[i].n_samples()
          //
          // Each A_q vector contains the logically 2D ragged array of
          // values (theta[i], sample[j(i)]), arranged in "row-major"
          // format, where the number of samples defined by each mu
          // object can vary, and they are all concatenated together
          // to produce one long list of samples.
          //
          // Example: suppose that there are 2 mu values containing
          // the same sets of parameters, with 5 samples defined in the
          // first and 10 samples defined in the second. Then, the A_q
          // vector will look like:
          //
          // {Theta(mu_vec[0], sample_0), Theta(mu_vec[0], sample_1), ... Theta(mu_vec[0], sample_4),
          //  Theta(mu_vec[1], sample_0), Theta(mu_vec[1], sample_1), ... Theta(mu_vec[1], sample_4), ..., Theta(mu_vec[1], sample_9)}
          //
          // Note: this example is unusual in that it would be
          // conceptually much simpler to have either:
          // (i) A vector of 15 single-sample RBParameters objects, or
          // (ii) A single RBParameters object with 15 samples defined in it.
          // but we can also support a combination of approaches (i)
          // and (ii) by concatenation, if necessary.
        }

      // 2.) Evaluate "F" thetas at all samples:
      std::vector<std::vector<Number>> all_F(rb_theta_expansion.get_n_F_terms());
      for (unsigned int q_f=0; q_f<rb_theta_expansion.get_n_F_terms(); q_f++)
        all_F[q_f] = rb_theta_expansion.eval_F_theta(q_f, mu_vec);

      // The eval_F_theta() calls above use "EIM solves" from eim_rb_eval, hence we
      // can now print out the EIM error indicator values.
      const auto & eim_error_indicators =
        eim_rb_eval.get_rb_eim_error_indicators();

      libMesh::out << std::scientific;
      for (auto idx : index_range(eim_error_indicators))
        libMesh::out << "EIM (error indicator, normalization factor) for step "
                     << idx << ": "
                     << "(" << eim_error_indicators[idx].first
                     << ", " << eim_error_indicators[idx].second
                     << ")" << std::endl;
      libMesh::out << std::endl;

      // 3.) Evaluate "output" thetas at all samples: in this case, the
      // output is particularly simple (does not depend on mu) so this
      // should just be a vector of all 1s.
      std::vector<std::vector<Number>> all_outputs(rb_theta_expansion.get_total_n_output_terms());
      {
        unsigned int output_counter = 0;
        for (unsigned int n=0; n<rb_theta_expansion.get_n_outputs(); n++)
          for (unsigned int q_l=0; q_l<rb_theta_expansion.get_n_output_terms(n); q_l++)
            all_outputs[output_counter++] =
              rb_theta_expansion.eval_output_theta(n, q_l, mu_vec);
      }

      // The total number of thetas is the sum of the "A", "F", and
      // "output" thetas.
      unsigned int n_thetas =
        rb_theta_expansion.get_n_A_terms() +
        rb_theta_expansion.get_n_F_terms() +
        rb_theta_expansion.get_total_n_output_terms();

      // Allocate enough space to store all the evaluated theta values for this RBThetaExpansion
      std::vector<Number> evaluated_thetas(n_thetas);

      // The RBConstruction object is used for plotting reduced-basis solutions (visualization)
      EquationSystems equation_systems (mesh);

      SimpleRBConstruction & rb_construction =
        equation_systems.add_system<SimpleRBConstruction> ("RB");

      equation_systems.init ();

      // Tell the RBConstruction object about the RBEvaluation object
      // and read in the RB basis functions.
      rb_construction.set_rb_evaluation(rb_eval);
      rb_eval.read_in_basis_functions(rb_construction, "rb_data");

      // Loop over each sample, fill the evaluated_thetas array, call rb_solve()
      for (unsigned sample_idx=0; sample_idx<online_mu.n_samples(); ++sample_idx)
        {
          unsigned int counter = 0;

          // Set A Theta values for current sample
          for (unsigned int q_a=0; q_a<rb_theta_expansion.get_n_A_terms(); q_a++)
            evaluated_thetas[counter++] = all_A[q_a][sample_idx];

          // Set F Theta values for current sample
          for (unsigned int q_f=0; q_f<rb_theta_expansion.get_n_F_terms(); q_f++)
            evaluated_thetas[counter++] = all_F[q_f][sample_idx];

          // Set output Theta values for current sample
          {
            unsigned int output_counter = 0;
            for (unsigned int n=0; n<rb_theta_expansion.get_n_outputs(); n++)
              for (unsigned int q_l=0; q_l<rb_theta_expansion.get_n_output_terms(n); q_l++)
                evaluated_thetas[counter++] = all_outputs[output_counter++][sample_idx];
          }

          // Call rb_solve() for the current thetas
          libMesh::out << "Performing solve for step " << sample_idx << std::endl;
          rb_eval.rb_solve(rb_eval.get_n_basis_functions(), &evaluated_thetas);

          // Print the output as well as the corresponding output error bound
          for (auto i : index_range(rb_eval.RB_outputs))
            libMesh::out << "Output value " << i << " = " << rb_eval.RB_outputs[i]
                         << ", error bound " << i << " = " << rb_eval.RB_output_error_bounds[i]
                         << std::endl;

          // Write an exo file to visualize this solution
          rb_construction.load_rb_solution();
#ifdef LIBMESH_HAVE_EXODUS_API
          ExodusII_IO(mesh).write_equation_systems("RB_sol_" + std::to_string(sample_idx) + ".e", equation_systems);
#endif
        } // end for (sample)
    }

#endif // LIBMESH_ENABLE_DIRICHLET
}
