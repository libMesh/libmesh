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

#ifndef __rb_reference_component_h__
#define __rb_reference_component_h__

#include "rb_system.h"

namespace libMesh
{

template <typename T> class NumericVector;

/**
 * This class is part of the rbOOmit framework.
 *
 * RBReferenceComponent is an RBSystem that
 * defines a reference component. An RBComponent
 * is analogous to a physical finite element, and
 * and RBReferenceComponent is analogous to a reference
 * element.
 * RBComponents can be connected to other systems via ports
 * and interface functions.
 * RBComponentSystem stores a group of RBComponents,
 * analogous to a mesh of finite elements.
 *
 * @author David J. Knezevic, 2011
 */

// ------------------------------------------------------------
// RBReferenceComponent class definition

class RBReferenceComponent : public RBSystem
{
public:

  /**
   * Constructor.  Optionally initializes required
   * data structures.
   */
  RBReferenceComponent (EquationSystems& es,
                        const std::string& name,
                        const unsigned int number);

  /**
   * Destructor.
   */
  virtual ~RBReferenceComponent ();

  /**
   * The type of the parent.
   */
  typedef RBSystem Parent;

  /**
   * The type of system.
   */
  typedef RBReferenceComponent sys_type;
  
  /**
   * Clear all the data structures associated with
   * the system.
   */
  virtual void clear ();

  /**
   * @returns a clever pointer to the system.
   */
  sys_type & system () { return *this; }
  
  /**
   * Set the number of ports that this components has.
   */
  void set_n_ports(unsigned int n_ports_in);
  
  /**
   * Get the number of ports.
   */
  unsigned int get_n_ports() const
    { return n_ports; }

  /**
   * Get the number of local interface functions.
   */
  unsigned int n_local_interface_functions() const
    { return interface_function_representors.size(); }

  /**
   * Return a reference to the specified interface function representor vector.
   */
  NumericVector<Number>& get_interface_function_representor(unsigned int i)
    { return *interface_function_representors[i]; }

  /**
   * Return the list of local interface functions which have support on port \p i.
   */
  std::vector<unsigned int> get_local_dofs_on_port(unsigned int i) const;

  /**
   * Return the list of ports on which the interface function with local index
   * \p i have support.
   */
  std::vector<unsigned int> get_ports_with_support(unsigned int i) const;

  /**
   * Override attach_A_q to just throw an error. Should use
   * attach_component_affine_operator in RBReferenceComponent
   * and its subclasses.
   */
  virtual void attach_A_q(theta_q_fptr ,
                          affine_assembly_fptr ,
                          affine_assembly_fptr )
  {
    libMesh::out << "Error: Cannot use attach_A_q in RBSystem. "
                 << "Use attach_component_affine_operator instead." << std::endl;
    libmesh_error();
  }

  /**
   * Override attach_F_q to just throw an error. Should use
   * attach_component_affine_operator in RBReferenceComponent
   * and its subclasses.
   */
  virtual void attach_F_q(theta_q_fptr ,
                          affine_assembly_fptr ,
                          affine_assembly_fptr )
  {
    libMesh::out << "Error: Cannot use attach_A_q in RBSystem. "
                 << "Use attach_component_affine_operator instead." << std::endl;
    libmesh_error();
  }

  /**
   * Attach an affine operator for this RBReferenceComponent. This
   * adds the operator to both the lhs and rhs expansions,
   * as required by the bubble function formulation.
   */
  void attach_component_affine_operator(theta_q_fptr theta_q_a,
                                        affine_assembly_fptr A_q_intrr_assembly,
                                        affine_assembly_fptr A_q_bndry_assembly);

  /**
   * Add a new interface function which has support on the port
   * \p port_nums. The interface function will be represented by
   * projecting the function g into the FE space.
   */
  void add_interface_function(unsigned int port_num,
                              Number g(const Point& p,
                                       const Parameters& ,
                                       const std::string& ,
                                       const std::string& ));

  /**
   * Add a new interface function which has support on the ports in
   * \p port_nums. The interface function will be represented by
   * projecting the function g into the FE space.
   */
  void add_interface_function(std::vector<unsigned int>& port_nums,
                              Number g(const Point& p,
                                       const Parameters& ,
                                       const std::string& ,
                                       const std::string& ));

  /**
   * Add a new interface function which has support on the port
   * \p port_num. The interface function will be represented by vec.
   */
  void add_interface_function(unsigned int port_num,
                              NumericVector<Number>& vec);

  /**
   * Add a new interface function which has support on the ports in
   * \p port_nums. The interface function will be represented by vec.
   */
  void add_interface_function(std::vector<unsigned int>& port_nums,
                              NumericVector<Number>& vec);

  /**
   * Overload initialize_RB_system in order to initialize extra data.
   */
  virtual void initialize_RB_system(bool online_mode);

  /**
   * Build a new RBReferenceComponentEvaluation object and add
   * it to the rb_evaluation_objects vector.
   */
  virtual void add_new_rb_evaluation_object();

  /**
   * Assemble the local contribution to the global matrix
   * for this reference component using the truth approximation.
   */
  void assemble_local_truth_matrix(DenseMatrix<Number>& local_matrix);

  /**
   * Assemble the local contribution to the global matrix
   * for this component using the RB approximation.
   */
  void assemble_local_RB_matrix(unsigned int N_rb,
                                DenseMatrix<Number>& local_matrix);

  /**
   * Overload assemble_all_affine_vectors to do nothing,
   * we can only assemble the affine vectors after we
   * select a lifting function.
   */
  virtual void assemble_all_affine_vectors() {}

  /**
   * Load the local portion of the global solution
   * defined by \p global_solution.
   */
  void load_local_solution(DenseVector<Number>& global_solution);

  /**
   * Train a reduced basis space for each interface function, and write out each
   * space to file as we go.
   */
  void train_and_write_all_reduced_bases(const std::string& directory_name = "offline_data");

  /**
   * Read in a reduced basis space for each interface function. These "bubble spaces"
   * are stored in the rb_evaluation_objects vector.
   */
  void read_in_all_reduced_bases(const std::string& directory_name = "offline_data");

  //----------- PUBLIC DATA MEMBERS -----------//

  /**
   * The local indices of the interface functions associated with each port.
   */
  std::vector< std::vector< unsigned int > > local_dofs_on_port;

  /**
   * The ports on which each local interface function has support.
   */
  std::vector< std::vector< unsigned int > > ports_with_support;

  /**
   * The set of interface functions associated with each port.
   */
  std::vector< NumericVector<Number>* > interface_function_representors;
  
  /**
   * We store a set of truth bubble functions in order to
   * compute the truth solution.
   */
  std::vector< NumericVector<Number>* > bubble_functions;

protected:

  /**
   * Helper function that updates the F_q vectors according
   * to the specified interface function.
   * This is necessary because in the RBReferenceComponent formulation
   * the right-hand side operators are a^q(interface_function,v).
   */
  void update_F_q_for_interface_function(unsigned int index);

  /**
   * Helper function that updates the data needed to assemble
   * the local contribution to the matrix.
   */
  void update_data_for_local_matrix();

private:

  //----------- PRIVATE DATA MEMBERS -----------//

  /**
   * The number of ports, i.e. connections, that this system has.
   */
  unsigned int n_ports;

};

} // namespace libMesh

#endif
