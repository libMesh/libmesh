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

#include "rb_param_subdomain_node.h"
#include "rb_param_subdomain_tree.h"
#include "rb_system.h"
#include "getpot.h"
#include "libmesh_logging.h"

// FOR CREATING A DIRECTORY
#include <sys/types.h>
#include <sys/stat.h>

#include "o_string_stream.h"
#include <fstream>

namespace libMesh
{


RBParamSubdomainTree::RBParamSubdomainTree(RBSystem& rb_system_in, const std::string& parameters_filename)
  : root_node(NULL),
    _rb_system(rb_system_in),
    h_tol(0.),
    p_tol(0.),
    N_bar(1),
    n_subsampled_training_points(0),
    bbox_margin(0.),
    leaf_node_index(0),
    clear_training_sets_during_hp(true)
{
  // Read in properties from parameters_filename
  GetPot infile(parameters_filename);

  // Set hp-Greedy tolerances
  h_tol = infile("h_tol",h_tol);
  p_tol = infile("p_tol",p_tol);

  N_bar = infile("N_bar",N_bar);
  if(N_bar > _rb_system.get_Nmax())
  {
    libMesh::out << "Error: Cannot set N_bar larger than Nmax for an RBParamSubdomainTree."
                 << std::endl;
    libmesh_error();
  }
  
  // by default use all training points
  n_subsampled_training_points = infile( "n_subsampled_training_points", _rb_system.get_n_training_samples() );

  // Set the subdomain bounding-box margin fraction
  bbox_margin = infile("bbox_margin",bbox_margin);

  libMesh::out << std::endl << "RBParamSubdomainTree parameters:" << std::endl;
  libMesh::out << "Tolerance for the h-stage (h_tol): " << h_tol << std::endl;
  libMesh::out << "Tolerance for the p-stage (p_tol): " << p_tol << std::endl;
  libMesh::out << "Nmax during h-stage (N_bar): " << N_bar << std::endl;
  libMesh::out << "n_subsampled_training_points: " << n_subsampled_training_points << std::endl;
  libMesh::out << "Bounding-box margin fraction for subdomain training set enrichment: " << bbox_margin << std::endl;
  libMesh::out << std::endl;
  
  // The rbOOmit code is still in a state of flux
  libmesh_experimental();
}

RBParamSubdomainTree::~RBParamSubdomainTree()
{
  if(root_node)
    {
      delete root_node;
      root_node = NULL;
    }
}

void RBParamSubdomainTree::build_root_node()
{
  libmesh_assert(root_node == NULL);
  root_node = new RBParamSubdomainNode(*this, _rb_system.get_current_parameters());
}

void RBParamSubdomainTree::hp_greedy(bool store_basis_functions)
{
  START_LOG("hp_greedy()", "RBParamSubdomainTree");

  // Build the root node rather than in the ctor
  // so we can use virtual function
  if(!root_node)
    build_root_node();

  root_node->copy_training_set_from_system();

  root_node->hp_greedy(store_basis_functions);

  STOP_LOG("hp_greedy()", "RBParamSubdomainTree");
}

RBParamSubdomainNode * RBParamSubdomainTree::determine_subdomain(const std::vector<Real>& new_param)
{
  START_LOG("determine_subdomain()", "RBParamSubdomainTree");

  RBParamSubdomainNode * current_node = root_node;

  while ( (current_node->left_child != NULL) && (current_node->right_child != NULL) )
    {
      if (current_node->left_child->dist_from_anchor(new_param) <= current_node->right_child->dist_from_anchor(new_param))
	current_node = current_node->left_child;
      else
	current_node = current_node->right_child;
    }

  STOP_LOG("determine_subdomain()", "RBParamSubdomainTree");

  return current_node;
}

void RBParamSubdomainTree::write_tree_data_to_file(const std::string& directory_name)
{
  START_LOG("write_tree_data_to_file()", "RBParamSubdomainTree");

  if(libMesh::processor_id() == 0)
    {

      libMesh::out << std::endl << "Writing tree data to tree_data/tree.dat" << std::endl;

      if( mkdir(directory_name.c_str(), 0777) == -1)
	{
	  libMesh::out << "Directory " << directory_name << " already exists, overwriting contents." << std::endl;
	}
      {
	std::ofstream n_leaves_out;
	OStringStream file_name;
	file_name << directory_name << "/n_leaves.dat";
	n_leaves_out.open( file_name.str().c_str() );
	if ( !n_leaves_out.good() )
	  {
	    libMesh::err << "Error opening n_leaves.dat" << std::endl;
	    libmesh_error();
	  }
	n_leaves_out << leaf_node_index;
	n_leaves_out.close();
      }

      {
	std::ofstream tree_out;
	OStringStream file_name;
	file_name << directory_name << "/tree.dat";
	tree_out.open( file_name.str().c_str() );
	if ( !tree_out.good() )
	  {
	    libMesh::err << "Error opening tree.dat" << std::endl;
	    libmesh_error();
	  }

	std::vector<int> bool_vec_(1);
	bool_vec_[0] = 1;
	leaf_node_index = 0;
	// Actual writing is performed recursively
	write_tree_data_to_file_recursively(root_node, tree_out, bool_vec_);
	tree_out.close();
      }
    }

  STOP_LOG("write_tree_data_to_file()", "RBParamSubdomainTree");
}

void RBParamSubdomainTree::write_tree_data_to_file_recursively(RBParamSubdomainNode * current_node,
                                                               std::ofstream& stream,
                                                               std::vector<int>& bool_vec)
{
  if (current_node != NULL)
    {
      for (unsigned int i = 0; i< bool_vec.size(); i++)
	stream << bool_vec[i];
      stream << " ";
      for (unsigned int i = 0; i< current_node->anchor.size(); i++){
	stream << current_node->anchor[i] << " ";
      }
      if ((current_node->left_child == NULL) && (current_node->right_child == NULL)){
	stream << leaf_node_index;
	leaf_node_index ++;
	}
      else
	{
	  stream << -1;
	}
      stream << std::endl;

      std::vector<int> bool_vec_left(bool_vec.size()+1);
      std::vector<int> bool_vec_right(bool_vec.size()+1);

      for (unsigned int i = 0; i< bool_vec.size(); i++)
	{
	  bool_vec_left[i] = bool_vec[i];
	  bool_vec_right[i] = bool_vec[i];
	}
      bool_vec_left[bool_vec.size()] = 0;
      bool_vec_right[bool_vec.size()] = 1;

      write_tree_data_to_file_recursively(current_node->left_child, stream, bool_vec_left);
      write_tree_data_to_file_recursively(current_node->right_child, stream, bool_vec_right);
    }
}

void RBParamSubdomainTree::read_tree_data_from_file(const std::string& directory_name)
{
  START_LOG("read_tree_data_from_file()", "RBParamSubdomainTree");
  
  libMesh::out << "Reading tree data from tree_data/tree.dat" << std::endl;

  // First, need to build a root node
  if(!root_node)
    build_root_node();

  OStringStream file_name;
  file_name << directory_name << "/tree.dat";
  std::ifstream tree_in(file_name.str().c_str());
  if ( !tree_in.good() )
    {
      libMesh::err << "Error opening tree.dat" << std::endl;
      libmesh_error();
    }
  // Actual reading is performed recursively
  reconstruct_tree(root_node, tree_in);
  tree_in.close();

  STOP_LOG("read_tree_data_from_file()", "RBParamSubdomainTree");
}

void RBParamSubdomainTree::reconstruct_tree(RBParamSubdomainNode * current_node,  std::ifstream& stream)
{
  START_LOG("reconstruct_tree()", "RBParamSubdomainTree");

 std::string bool_vec_string;

 stream >> bool_vec_string;

 for (unsigned i = 0 ; i < _rb_system.get_n_params(); i++)
   stream >> current_node->anchor[i];

 stream >> current_node->model_number;

//  libMesh::out << bool_vec_string << " ";
//  for (unsigned i = 0 ; i < _rb_system.get_n_params(); i++)
//    libMesh::out<< current_node->anchor[i] << " ";
//  libMesh::out << current_node->model_number << std::endl;

 if (current_node->model_number == -1)
   {
    current_node->add_child(current_node->anchor, RBParamSubdomainNode::LEFT);
    current_node->add_child(current_node->anchor, RBParamSubdomainNode::RIGHT);
    reconstruct_tree(current_node->left_child, stream);
    reconstruct_tree(current_node->right_child, stream);
   }

  STOP_LOG("reconstruct_tree()", "RBParamSubdomainTree");
}

} // namespace libMesh
