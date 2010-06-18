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

#ifndef __rb_param_subdomain_tree_h__
#define __rb_param_subdomain_tree_h__

#include "libmesh_common.h"
#include <vector>

// Forward declaration
class RBSystem;
class RBParamSubdomainNode;

/**
 * This class is part of the rbOOmit framework.
 *
 * RBParamSubdomainTree defines an adaptive partition of parameter
 * space. Each RBParamSubdomainNode in the tree is associated with
 * a reduced basis model that is constructed for the subdomain
 * represented by that node. In this class, the name "hp" refers
 * to a reduced basis approach in which we first subdivide the parameter
 * domain ("h" refinement) and then enrich the reduced basis in each
 * subdomain ("p" refinement).
 *
 * @author David J. Knezevic and Jens L. Eftang, 2009
 */

// ------------------------------------------------------------
// RBParamSubdomainTree class definition

class RBParamSubdomainTree
{
public:

  /**
   * Constructor. Initializes required
   * data structures.
   */
  RBParamSubdomainTree (RBSystem& system, const std::string& parameters_filename);

  /**
   * Destructor.
   */
  virtual ~RBParamSubdomainTree ();

  /**
   * Build an RBParamSubdomainNode. Virtual function so we specify alternative
   * node types in subclasses.
   */
  virtual void build_root_node();

  /**
   * Function that calls hp_greedy on the root node.
   */
  void hp_greedy();

  /**
   * Determine the subdomain containing the parameter new_param
   */
  RBParamSubdomainNode * determine_subdomain(const std::vector<Real>& new_param);

  /**
   * Writes the tree data required for online binary search to a file
   */
  void write_tree_data_to_file(const std::string& directory_name);

  /**
   * Writes data for the current node to specified file and procedes through
   * the tree according to depth first
   */
  void write_tree_data_to_file_recursively(RBParamSubdomainNode * current_node,   std::ofstream& stream,   std::vector<int>& bool_vec);

  /**
   * Reads a tree from file
   */
  void read_tree_data_from_file(const std::string& directory_name);

  /**
   * Reads anchor and model number from file for current node and procedes through
   * the tree according to depth first
   */
  void reconstruct_tree(RBParamSubdomainNode * current_node,  std::ifstream& stream);

  /**
   * Pointer to root RBParamSubdomainNode
   */
  RBParamSubdomainNode * root_node;

  /**
   * Reference to the RBSystem that is used
   * to actually perform the truth and RB solves.
   */
  RBSystem & _rb_system;

  /**
   * Variable that indicates the "h-type" tolerance.
   */
  Real h_tol;

  /**
   * Variable that indicates the "p-type" tolerance.
   */
  Real p_tol;

  /**
   * Conservativity factor to determine POD_tol from h_tol
   */
  Real conserv_factor;

  /**
   * Number of basis functions during "h-type" stage
   */
  unsigned int N_bar;

  /**
   * Flag that indicates whether we use delta_N or the POD-tolerance
   * to govern how many basis functions we add in a time-dependent
   * hp-greedy.
   */
  bool use_delta_N_in_h_stage;

  /**
   * The number of training parameters in a subdomain's
   * training set _after_ the initial call to refine_training_set.
   */
  unsigned int refined_training_set_size;

  /**
   * The minimum number of training parameters in an RBParamSubdomainNode training set.
   */
  unsigned int min_training_set_size;

  /**
   * This value (between 0 and 1) defines an extra margin
   * as a fraction of a subdomain's size in a given dimension
   * to add to the bounding box when enriching the training set.
   */
  Real bbox_margin;

  /**
   * Variable that indicates the (current) number of
   * "leaf nodes"
   */
  int leaf_node_index;

};

#endif
