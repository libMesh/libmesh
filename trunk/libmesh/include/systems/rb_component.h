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

#ifndef __rb_component_h__
#define __rb_component_h__

#include "reference_counted_object.h"
#include <vector>

namespace libMesh
{

class RBReferenceComponent;
class RBSystem;
template <typename T> class DenseMatrix;
template <typename T> class DenseVector;

/**
 * This class is part of the rbOOmit framework.
 *
 * RBComponent provides a simulation component that
 * can be connected to other components via ports.
 * This class points to an RBReferenceComponent, which
 * provides the shape functions etc.
 *
 * @author David J. Knezevic, 2011
 */

// ------------------------------------------------------------
// RBComponent class definition

class RBComponent : public ReferenceCountedObject<RBComponent>
{
public:

  /**
   * Constructor.
   */
  RBComponent (RBReferenceComponent& ref_component_in, RBSystem* source_system=NULL);

  /**
   * Destructor.
   */
  virtual ~RBComponent ();
  
  /**
   * Set a (unique) global ID for this RBComponent.
   */
  void set_global_id(unsigned int id_in);
  
  /**
   * Get a (unique) global ID for this RBComponent.
   */
  unsigned int get_global_id() const;

  /**
   * Get/set the component's parameters.
   */
  const std::vector<Real>& get_component_parameters() const { return component_parameters; };
  void set_component_parameters(const std::vector<Real>& params);

  /**
   * Print the current parameters.
   */
  void print_component_parameters();

  /**
   * Set neighbor component on port \p port.
   */
  void set_neighbor_component(unsigned int port, RBComponent& neighbor);

  /**
   * Set the neighbor components for each port. Here the i^th neighbor
   * is set to the i^th entry of \p neighbors.
   */
  void set_neighbor_components(std::vector< RBComponent* >& neighbors);
  
  /**
   * Get a pointer to the neighbor component on port \p port.
   * NULL indicates that the component has no neighbor on \p port.
   */
  RBComponent* get_neighbor_component(unsigned int port) { return neighbor_components[port]; }

  /**
   * Set the pointer to the RBSystem that provides us with
   * the source bubble function. A NULL source system means
   * that we do not have a source on this component.
   */
  void set_source_system(RBSystem* source_system_in);

  /**
   * Assemble the local contribution to the global matrix
   * and right-hand side vector for this component using the truth
   * approximation.
   */
  void assemble_local_truth_matrix_and_rhs(DenseMatrix<Number>& local_matrix,
                                           DenseVector<Number>& local_rhs);

  /**
   * Assemble the local contribution to the global matrix
   * and right-hand side vector for this component using the RB
   * approximation.
   */
  void assemble_local_RB_matrix_and_rhs(unsigned int N_rb,
                                        DenseMatrix<Number>& local_matrix,
                                        DenseVector<Number>& local_rhs);

  /**
   * Load the solution for this component into ref_component.
   */
  void load_component_solution(DenseVector<Number>& global_solution);

  /**
   * Updates the data source/reference-component extra data
   * that is required to compute the local contribution to
   * the right-hand side vector.
   * If we do not have a source_system, then this function
   * does nothing.
   */
  void update_extra_source_data();

  /**
   * Write out the extra data from the source_system that is required
   * to assemble the global matrix and right-hand side.
   */
  void write_extra_source_data(const std::string& directory_name);
  
  /**
   * Read in the extra data from the source_system that is required
   * to assemble the global matrix and right-hand side.
   */
  void read_extra_source_data(const std::string& directory_name);

  //----------- PUBLIC DATA MEMBERS -----------//

  /**
   * The reference component. This is a reference rather than a pointer
   * to ensure that we already have access to a reference component.
   */
  RBReferenceComponent& ref_component;

  /**
   * The RBSystem that provides us with the source bubble function.
   * If source_system==NULL, then we have no source on this component.
   */
  RBSystem* source_system;

  /**
   * The set of global indices associated with the local interface functions.
   * These determine the local-to-global dof mapping in an RBComponentSystem.
   */
  std::vector< unsigned int > interface_function_global_indices;

  /**
   * The list of neighbor components. The i^th entry
   * of the vector is the neighbor at the i^th port.
   */
  std::vector<RBComponent*> neighbor_components;

  /**
   * Vectors storing the data necessary to assemble the local
   * contribution to the global right-hand side, using the
   * reduced basis approximation provided by source_system.
   */
  std::vector< std::vector<Number> >                source_Fq_g;
  std::vector< std::vector< std::vector<Number> > > source_Aq_bubble_g;

  /**
   * Vector storing the component's parameters.
   */
  std::vector<Real> component_parameters;
  
  /**
   * A unique global ID for this RBComponent.
   */
  unsigned int global_id;

};

} // namespace libMesh

#endif
