// rbOOmit: An implementation of the Certified Reduced Basis method.
// Copyright (C) 2009, 2010 David J. Knezevic

// This file is part of rbOOmit.

// rbOOmit is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
  
// rbOOmit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#ifndef __rb_component_system_h__
#define __rb_component_system_h__

#include "reference_counted_object.h"
#include <vector>
#include "rb_component.h"
#include "dense_matrix.h"
#include "dense_vector.h"

namespace libMesh
{

/**
 * This class is part of the rbOOmit framework.
 *
 * RBComponentSystem provides the functionality to connect a group of
 * RBComponents together to form an assembled system.
 *
 * @author David J. Knezevic, 2011
 */

// ------------------------------------------------------------
// RBComponentSystem class definition

class RBComponentSystem : public ReferenceCountedObject<RBComponentSystem>
{
public:

  /**
   * Constructor.
   */
  RBComponentSystem ();

  /**
   * Destructor.
   */
  virtual ~RBComponentSystem ();

  /**
   * The type of system.
   */
  typedef RBComponentSystem sys_type;

  /**
   * Attach the vector of pointers to the Reference Components
   * used in this component system.
   */
  void attach_reference_components(std::vector<RBReferenceComponent*>& reference_components_in)
  { reference_components = reference_components_in; }
  
  /**
   * Attach the vector of pointers to the source systems used in
   * this component system.
   */
  void attach_source_systems(std::vector<RBSystem*>& source_systems_in)
  { source_systems = source_systems_in; }
  
  /**
   * Add a new type of reference lego to the system.
   * @return a reference to the newly added component.
   */
  RBComponent& add_component(RBReferenceComponent& component,
                             RBSystem* source_system=NULL);

  /**
   * Initialize the reference components and source systems.
   * The boolean \p do_not_assemble determines whether or not
   * we pre-assemble the affine operators in each system.
   * Override in subclasses to perform any "mesh specific"
   * connectivity initialization of physical components.
   */
  virtual void initialize_component_system(bool do_not_assemble);
  
  /**
   * Train the reduced basis spaces for each reference component
   * and source system, and write out all the necessary data.
   */
  void train_and_write_all_components(const std::string& directory_name = "offline_data");

  /**
   * Read in the offline data and initialize all the components accordingly.
   */
  void read_in_all_components(const std::string& directory_name = "offline_data");

  /**
   * @return the number of global degrees of freedom (interface functions).
   */
  unsigned int get_n_global_dofs() const;

  /**
   * @return the number of components in the system.
   */
  unsigned int n_components() const { return components.size(); }
  
  /**
   * @returns a pointer to the specified component.
   */
  RBComponent& get_component(unsigned int id);
  
  /**
   * Loop over the components in the system and assign a unique
   * global ID to each interface function.
   * This function also updates the member n_global_dofs.
   */
  void renumber_interface_functions();

  /**
   * Assemble and solve the truth component system.
   */
  void truth_solve();

  /**
   * Assemble and solve the RB component system.
   * N is the default number of RB functions we use on
   * each component (but we automatically cap N if
   * N > number of basis functions)
   */
  void RB_solve(unsigned int N);

  //----------- PUBLIC DATA MEMBERS -----------//

  typedef std::vector< RBComponent* >::iterator       component_iterator;
  typedef std::vector< RBComponent* >::const_iterator const_component_iterator;

  /**
   * The vector that stores the family of reference legos.
   */
  std::vector< RBComponent* > components;

  /**
   * The vector of reference components used in this system.
   */
  std::vector< RBReferenceComponent* > reference_components;
  
  /**
   * The vector of source systems used in this system.
   */
  std::vector< RBSystem* > source_systems;

  /**
   * The global solution vector that stores the coefficients
   * of the interface functions in the system.
   * Updated by the solve functions.
   */
  DenseVector<Number> global_solution;

  /**
   * The global matrix in which we accumulate the contribution
   * of each component.
   */
  DenseMatrix<Number> global_matrix;
  
  /**
   * The global right-hand side vector in which we accumulate
   * the contribution of each component.
   */
  DenseVector<Number> global_rhs;

protected:

  /**
   * Loop over the components and add the truth local matrices
   * and right-hand side vectors to the global matrix and rhs.
   */
  void assemble_truth_global_matrix_and_rhs();

  /**
   * Loop over the components and add the RB local matrices
   * and right-hand side vectors to the global matrix and rhs.
   */  
  void assemble_RB_global_matrix_and_rhs(unsigned int N_rb);

  /**
   * Add local_matrix and local_rhs to global_matrix and global_rhs.
   * Use \p global_indices to determine the local-to-global mapping.
   */
  void add_local_matrix(const DenseMatrix<Number>& local_matrix,
                        const DenseVector<Number>& local_rhs,
                        std::vector<unsigned int>& global_indices);

private:

  //----------- PRIVATE DATA MEMBERS -----------//
  
  /**
   * The number of global degrees of freedom (interface functions).
   */
  unsigned int n_global_dofs;

};

} // namespace libMesh

#endif
