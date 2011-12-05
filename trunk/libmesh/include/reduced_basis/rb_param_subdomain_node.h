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

#ifndef __rb_param_subdomain_node_h__
#define __rb_param_subdomain_node_h__

#include "libmesh_common.h"
#include <vector>

namespace libMesh
{

// Forward declaration
class RBConstruction;
class RBParamSubdomainTree;

/**
 * This class is part of the rbOOmit framework.
 *
 * RBParamSubdomainNode implements a "node" of an "hp-tree" structure
 * which is used to generate a hierarchical partition of the parameter
 * domain. Each RBParamSubdomainNode is associated with a reduced basis
 * model that is developed on that particular subdomain. In this class,
 * the name "hp" refers to a reduced basis approach in which we first
 * subdivide the parameter domain ("h" refinement) and then enrich the
 * reduced basis in each subdomain ("p" refinement).
 *
 * @author David J. Knezevic and Jens L. Eftang, 2009
 */

// ------------------------------------------------------------
// RBParamSubdomainNode class definition

class RBParamSubdomainNode
{
public:

  /**
   * Define enumeration to define left/right children.
   */
  enum Child { LEFT  = 0,
               RIGHT = 1 };

  /**
   * Constructor. Initializes required data structures.
   */
  RBParamSubdomainNode (RBParamSubdomainTree& tree, const std::vector<Real>& anchor);

  /**
   * Destructor.
   */
  virtual ~RBParamSubdomainNode ();

  /**
   * Test whether *this is the same as the argument.
   */
  bool operator== (const RBParamSubdomainNode& node) const;

  /**
   * Recursive function that performs the "hp" greedy algorithm,
   * i.e. this function recursively constructs the hierarchical
   * subdivision of parameter space and builds the reduced basis
   * spaces on each subdomain.
   * If \p store_basis_functions==true, write out the
   * basis functions associated with each subdomain
   * during the hp-greedy algorithm.
   */
  virtual void hp_greedy(bool store_basis_functions);

  /**
   * Split the current subdomain into two new subdomains.
   * The argument h_stage_split indicates whether or not
   * this is called during the h-stage or the p-stage.
   */
  virtual void split_this_subdomain(bool h_stage_split);

  /**
   * This function performs the "p" stage of the "hp"
   * greedy algorithm.
   */
  virtual Real perform_p_stage(Real greedy_bound);

  /**
   * Write out the offline data for the current subdomain. Overload
   * in subdomain in order to perform extra funcitonality, e.g.
   * a per-subdomain SCM which is not needed during the hp refinement.
   * If \p store_basis_functions==true, write out the basis functions
   * associated with the current subdomain.
   */
  virtual void write_subdomain_data_to_files(bool store_basis_functions);

  /**
   * Returns the number of parameters.
   */
  unsigned int get_n_params() const { return anchor.size(); }

  /**
   * Returns the "weighted Euclidean distance" between two points in the parameter
   * domain. Each component of the each parameter is mapped to the interval
   * (0,1) in order to compute an appropriately scaled distance.
   */
  Real distance(const std::vector<Real>& p1, const std::vector<Real>& p2) const;

  /**
   * Returns the "weighted Euclidean distance" from the anchor point to
   * the input parameter. Each component of the anchor and new_param are mapped to
   * (0,1) before computing the distance.
   */
  Real dist_from_anchor(const std::vector<Real>& new_param) const;

  /**
   * Add a new (LEFT or RIGHT) child RBParamSubdomainNode. Virtual function so
   * subclasses can add children of the same type.
   */
  virtual void add_child(const std::vector<Real>& child_anchor, Child c);

  /**
   * Copy the training set from the associated RBConstruction to this
   * RBParamSubdomainNode.
   */
  void copy_training_set_from_construction();

  /**
   * Refine the training set of the current RBParamSubdomainNode
   */
  void refine_training_set(const unsigned int new_training_set_size);

  /**
   * Initialize the training sets of the child nodes. This splits
   * the current training set based on proximity to the child
   * anchor points.
   */
  void initialize_child_training_sets();

  /**
   * Returns the number of training parameters in training_set
   * local to this processor.
   */
  unsigned int n_local_training_parameters() const;

  /**
   * Returns the number of global training parameters in training_set.
   */
  unsigned int n_global_training_parameters() const;

  /**
   * Returns the i^th training parameter in the parameter set
   * stored on this processor.
   */
  std::vector<Real> get_local_training_parameter(unsigned int i);

  /**
   * Returns the corners of a box containing the current train sample
   * with specified margin. The margin (between 0 and 1) is defined
   * by the bbox_margin member variables and adds a fraction of the
   * parameter range to each side of the box in each dimension.
   */
  std::vector< std::vector<Real> > get_training_bbox();

  /**
   * Returns a new training set which is a random sample of
   * _tree.n_subsampled_training_points from the training set.
   * If the subsample size is greater than the training set size
   * just return the full training set.
   */
  std::vector< std::vector<Number> > get_subsampled_training_set();

  /**
   * Clear the training set associated with this subdomain.
   */
  void clear_training_set();

  /**
   * Pointers to the child RBParamSubdomainNodes.
   */
  RBParamSubdomainNode * left_child;
  RBParamSubdomainNode * right_child;

  /**
   * Reference to the tree containing this RBParamSubdomainNode
   */
  RBParamSubdomainTree& _tree;

  /**
   * Reference to the RBConstruction that is used
   * to actually perform the truth and RB solves.
   */
  RBConstruction& _rb_construction;

  /**
   * The anchor parameter value. The anchor point is used
   * in order to determine which subdomain a parameter
   * belongs to based on proximity. At each step in descending
   * the tree data structure we go "left" or "right"
   * depending on whether a given parameter value is closer
   * to the left child's or right child's anchor point.
   */
  std::vector<Real> anchor;

  /**
   * The distance to the sibling subdomain anchor point. This
   * is used to construct an approximate bounding box around
   * the subdomain when we wish to enrich a subdomain's training
   * set.
   */
  Real distance_between_anchors;

  /**
   * The set of training points for the parameter
   * subdomain corresponding to this node. We copy
   * this training set to _rb_construction in order to
   * generate a reduced basis for this subdomain.
   */
  std::vector< std::vector<Number> > training_set;

  /**
   * An index to identify "leaf nodes" of an hp-tree
   * data structure. The leaf nodes must be distinguished
   * from one another since they will be used for Online
   * reduced basis calculations.
   */
  int model_number;

  /**
   * Boolean flag to indicate whether or not the training
   * set has been initialized.
   */
  bool training_set_initialized;

};

} // namespace libMesh

#endif
