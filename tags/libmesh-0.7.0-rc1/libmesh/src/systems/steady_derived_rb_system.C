// specialization for Base = RBSystem
// which is typedef'd to SteadyDerivedRBSystem

#include "derived_rb_system.h"
#include "libmesh_logging.h"
#include "equation_systems.h"
#include "mesh_base.h"
#include "gmv_io.h"

namespace libMesh
{

template <>
Real DerivedRBSystem<RBSystem>::truth_solve(int plot_solution)
{
  START_LOG("truth_solve()", "DerivedRBSystem");
  
  EquationSystems& es = this->get_equation_systems();
  RBSystem& uber_system = es.get_system<RBSystem>(uber_system_name);

  set_uber_current_parameters();
  
  uber_system.RB_solve(uber_system.get_n_basis_functions());
  
  if(plot_solution > 0)
  {
    uber_system.load_RB_solution();
    *solution = *(uber_system.solution);
    const MeshBase& mesh = get_mesh();
    GMVIO(mesh).write_equation_systems ("unter_uber_truth.gmv",
                                        this->get_equation_systems());
  }

  STOP_LOG("truth_solve()", "DerivedRBSystem");
  
  // Don't bother returning the norm of the uber solution
  return 0.;
}

template <>
void DerivedRBSystem<RBSystem>::enrich_RB_space()
{
  START_LOG("enrich_RB_space()", "DerivedRBSystem");
  
  EquationSystems& es = this->get_equation_systems();
  RBSystem& uber_system = es.get_system<RBSystem>(uber_system_name);
  const unsigned int uber_size = uber_system.get_n_basis_functions();

  DenseVector<Number> new_bf = uber_system.RB_solution;

  // compute Gram-Schmidt orthogonalization
  DenseVector<Number> proj_sum(uber_size);
  for(unsigned int index=0; index<get_n_basis_functions(); index++)
  {
    // orthogonalize using the Identity matrix as the inner product,
    // since the uber basis functions should be orthogonal already
    // (i.e. neglect possible rounding errors in uber orthogonalization)
    Number scalar = new_bf.dot(derived_basis_functions[index]);
    proj_sum.add(scalar, derived_basis_functions[index]);
  }
  new_bf -= proj_sum;
  new_bf.scale(1./new_bf.l2_norm());

  // load the new basis function into the basis_functions vector.
  derived_basis_functions.push_back( new_bf );

  STOP_LOG("enrich_RB_space()", "DerivedRBSystem");
}


template <>
void DerivedRBSystem<RBSystem>::update_RB_system_matrices()
{
  START_LOG("update_RB_system_matrices()", "DerivedRBSystem");
  
  EquationSystems& es = this->get_equation_systems();
  RBSystem& uber_system = es.get_system<RBSystem>(uber_system_name);

  unsigned int derived_RB_size = get_n_basis_functions();
  unsigned int uber_RB_size    = uber_system.get_n_basis_functions();
  
  DenseVector<Number> temp_vector;
  for(unsigned int q_f=0; q_f<get_Q_f(); q_f++)
  {
    for(unsigned int i=(derived_RB_size-delta_N); i<derived_RB_size; i++)
    {
      uber_system.RB_F_q_vector[q_f].get_principal_subvector(uber_RB_size, temp_vector);
      RB_F_q_vector[q_f](i) = temp_vector.dot(derived_basis_functions[i]);
    }
  }

  DenseMatrix<Number> temp_matrix;
  for(unsigned int i=(derived_RB_size-delta_N); i<derived_RB_size; i++)
  {
    for(unsigned int n=0; n<get_n_outputs(); n++)
      for(unsigned int q_l=0; q_l<get_Q_l(n); q_l++)
      {
        uber_system.RB_output_vectors[n][q_l].get_principal_subvector(uber_RB_size, temp_vector);
        RB_output_vectors[n][q_l](i) = temp_vector.dot(derived_basis_functions[i]);
      }

    for(unsigned int j=0; j<derived_RB_size; j++)
    {
      for(unsigned int q_a=0; q_a<get_Q_a(); q_a++)
      {
        // Compute reduced A_q matrix
        uber_system.RB_A_q_vector[q_a].get_principal_submatrix(uber_RB_size, temp_matrix);
        temp_matrix.vector_mult(temp_vector, derived_basis_functions[j]);
        RB_A_q_vector[q_a](i,j) = derived_basis_functions[i].dot(temp_vector);

        if(i!=j)
        {
          temp_vector.zero();
          temp_matrix.vector_mult(temp_vector, derived_basis_functions[i]);
          RB_A_q_vector[q_a](j,i) = (derived_basis_functions[j]).dot(temp_vector);
        }
      }
    }
  }

  STOP_LOG("update_RB_system_matrices()", "DerivedRBSystem");
}


template <>
void DerivedRBSystem<RBSystem>::update_residual_terms(bool compute_inner_products)
{
  START_LOG("update_residual_terms()", "DerivedRBSystem");
  
  EquationSystems& es = this->get_equation_systems();
  RBSystem& uber_system = es.get_system<RBSystem>(uber_system_name);
  
  switch(residual_type_flag)
  {
    case(RESIDUAL_WRT_UBER):
    {
      // First we need to compute the representors of
      // each basis function using inner_product_matrix

      unsigned int derived_RB_size = get_n_basis_functions();
      unsigned int uber_RB_size = uber_system.get_n_basis_functions();
      DenseVector<Number> temp_vector1, temp_vector2;

      // We only need to compute the representors for
      // the right-hand side once.  Note: this will be called at least
      // once even after a restart, in which case F_q_representor[0]
      // will already be set.
      if(!update_residual_terms_called)
      {

        // Assume inner product matrix is the identity, hence don't need to
        // do any solves
        if (compute_inner_products)
        {
	  unsigned int q=0;
	  for(unsigned int q_f1=0; q_f1<get_Q_f(); q_f1++)
	  {
	    for(unsigned int q_f2=q_f1; q_f2<get_Q_f(); q_f2++)
	      {
                uber_system.RB_F_q_vector[q_f2].get_principal_subvector(uber_RB_size, temp_vector1);
                uber_system.RB_F_q_vector[q_f1].get_principal_subvector(uber_RB_size, temp_vector2);
		Fq_representor_norms[q] = temp_vector1.dot( temp_vector2 );
		q++;
	      }
	  }
        } // end if (compute_inner_products)

      } // end if(!update_residual_terms_called)

      // Now compute and store the inner products (if requested)
      if (compute_inner_products)
      {

        DenseMatrix<Number> temp_matrix;
        for(unsigned int q_f=0; q_f<get_Q_f(); q_f++)
	  {
	    for(unsigned int q_a=0; q_a<get_Q_a(); q_a++)
	      {
	        for(unsigned int i=(derived_RB_size-delta_N); i<derived_RB_size; i++)
		  {
                    uber_system.RB_A_q_vector[q_a].get_principal_submatrix(uber_RB_size, temp_matrix);
                    temp_matrix.vector_mult(temp_vector1, derived_basis_functions[i]);
                    uber_system.RB_F_q_vector[q_f].get_principal_subvector(uber_RB_size, temp_vector2);
		    Fq_Aq_representor_norms[q_f][q_a][i] = -temp_vector1.dot(temp_vector2);
		  }
	      }
	  }

        unsigned int q=0;
        for(unsigned int q_a1=0; q_a1<get_Q_a(); q_a1++)
	  {
	    for(unsigned int q_a2=q_a1; q_a2<get_Q_a(); q_a2++)
	      {
	        for(unsigned int i=(derived_RB_size-delta_N); i<derived_RB_size; i++)
		  {
		    for(unsigned int j=0; j<derived_RB_size; j++)
		      {
                        uber_system.RB_A_q_vector[q_a1].get_principal_submatrix(uber_RB_size, temp_matrix);
                        temp_matrix.vector_mult(temp_vector1, derived_basis_functions[i]);
                        uber_system.RB_A_q_vector[q_a2].get_principal_submatrix(uber_RB_size, temp_matrix);
                        temp_matrix.vector_mult(temp_vector2, derived_basis_functions[j]);
		        Aq_Aq_representor_norms[q][i][j] = temp_vector1.dot(temp_vector2);

		        if(i != j)
			  {
                            uber_system.RB_A_q_vector[q_a1].get_principal_submatrix(uber_RB_size, temp_matrix);
                            temp_matrix.vector_mult(temp_vector1, derived_basis_functions[j]);
                            uber_system.RB_A_q_vector[q_a2].get_principal_submatrix(uber_RB_size, temp_matrix);
                            temp_matrix.vector_mult(temp_vector2, derived_basis_functions[i]);
			    Aq_Aq_representor_norms[q][j][i] = temp_vector1.dot(temp_vector2);
			  }
		      }
		  }
	        q++;
	      }
	  }
      } // end if (compute_inner_products)
  
      update_residual_terms_called = true;
      break;
    }
          
    case(RESIDUAL_WRT_TRUTH):
    {
      unsigned int RB_size = get_n_basis_functions();

      if(!update_residual_terms_called)
      {
        // Copy the output terms over from uber_system
        for(unsigned int n=0; n<get_n_outputs(); n++)
        {
          output_dual_norms[n] = uber_system.output_dual_norms[n];
        }

        // Copy the Fq terms over from uber_system
        Fq_representor_norms = uber_system.Fq_representor_norms;
      }

      for(unsigned int q_f=0; q_f<get_Q_f(); q_f++)
      {
        for(unsigned int q_a=0; q_a<get_Q_a(); q_a++)
        {
          for(unsigned int i=(RB_size-delta_N); i<RB_size; i++)
          {
            Fq_Aq_representor_norms[q_f][q_a][i] = 0.;
            for(unsigned int j=0; j<uber_system.get_n_basis_functions(); j++) // Evaluate the dot product
            {
              Fq_Aq_representor_norms[q_f][q_a][i] +=
	        uber_system.Fq_Aq_representor_norms[q_f][q_a][j] * derived_basis_functions[i](j);
            }
          }
        }
      }

      unsigned int q=0;
      for(unsigned int q_a1=0; q_a1<get_Q_a(); q_a1++)
      {
        for(unsigned int q_a2=q_a1; q_a2<get_Q_a(); q_a2++)
        {
          for(unsigned int i=(RB_size-delta_N); i<RB_size; i++)
          {
            for(unsigned int j=0; j<RB_size; j++)
	    {

              Aq_Aq_representor_norms[q][i][j] = 0.;
              if(i != j)
                Aq_Aq_representor_norms[q][j][i] = 0.;

              for(unsigned int k=0; k<uber_system.get_n_basis_functions(); k++)
                for(unsigned int k_prime=0; k_prime<uber_system.get_n_basis_functions(); k_prime++)
                {
                  Aq_Aq_representor_norms[q][i][j] += derived_basis_functions[i](k)*derived_basis_functions[j](k_prime)*
                                                      uber_system.Aq_Aq_representor_norms[q][k][k_prime];
                                                 
                  if(i != j)
                  {
                    Aq_Aq_representor_norms[q][j][i] += derived_basis_functions[j](k)*derived_basis_functions[i](k_prime)*
                                                        uber_system.Aq_Aq_representor_norms[q][k][k_prime];
                  }
                }
            }
          }
          q++;
        }
      }

      update_residual_terms_called = true;

      break;
    }

    default:
    {
      std::cout << "Invalid RESIDUAL_TYPE in update_residual_terms" << std::endl;
      break;
    }
  }

  STOP_LOG("update_residual_terms()", "DerivedRBSystem");
}

template<>
void DerivedRBSystem<RBSystem>::load_RB_solution()
{
  START_LOG("load_RB_solution()", "DerivedRBSystem");

  if(!initialize_calN_dependent_data)
  {
    std::cerr << "Error: We must initialize the calN dependent "
              << "data structures in order to load RB solution."
              << std::endl;
    libmesh_error();
  }

  solution->zero();

  if(RB_solution.size() > get_n_basis_functions())
  {
    std::cerr << "ERROR: System contains " << get_n_basis_functions() << " basis functions."
              << " RB_solution vector constains " << RB_solution.size() << " entries."
              << " RB_solution in RBSystem::load_RB_solution is too long!" << std::endl;
    libmesh_error();
  }

  EquationSystems& es = this->get_equation_systems();
  RBSystem& uber_system = es.get_system<RBSystem>(uber_system_name);

  for(unsigned int i=0; i<RB_solution.size(); i++)
    for(unsigned int j=0; j<uber_system.get_n_basis_functions(); j++)
    {
      solution->add(RB_solution(i)*derived_basis_functions[i](j), uber_system.get_bf(j));
    }

  update();

  STOP_LOG("load_RB_solution()", "DerivedRBSystem");
}

}