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

// Configuration data
#include "libmesh_config.h"

#if defined(LIBMESH_HAVE_SLEPC) && (LIBMESH_HAVE_GLPK)

#include "qn_transient_rb_param_subdomain_node.h"
#include "qn_transient_rb_param_subdomain_tree.h"
#include "qn_transient_scm_system.h"
#include "rb_system.h"

namespace libMesh
{

QNTransientRBParamSubdomainNode::QNTransientRBParamSubdomainNode(QNTransientRBParamSubdomainTree& tree,
                                                                 const std::vector<Real>& anchor)
  : Parent(tree, anchor),
    _scm_system(tree._scm_system)
{
}

void QNTransientRBParamSubdomainNode::add_child(const std::vector<Real>& new_anchor, Child c)
{
  // cast the tree reference to a QNTransientRBParamSubdomainTree
  QNTransientRBParamSubdomainTree& qn_tree =
    libmesh_cast_ref<QNTransientRBParamSubdomainTree&>(_tree);

  if(c == LEFT)
  {
    if(left_child != NULL)
    {
      std::cout << "Error: Child already exists!"
                << std::endl;
      libmesh_error();
    }
    else
    {
      // Allocate a new child node
      left_child = new QNTransientRBParamSubdomainNode(qn_tree, new_anchor);
    }
  }

  if(c == RIGHT)
  {
    if(right_child != NULL)
    {
      std::cout << "Error: Child already exists!"
                << std::endl;
      libmesh_error();
    }
    else
    {
      // Allocate a new child node
      right_child = new QNTransientRBParamSubdomainNode(qn_tree, new_anchor);
    }
  }
}

void QNTransientRBParamSubdomainNode::write_subdomain_data_to_files()
{
  RBParamSubdomainNode::write_subdomain_data_to_files();
  
  // Finally, perform the SCM stage on this subdomain
  std::cout << std::endl << "Perform SCM greedy for leaf subdomain " << _tree.leaf_node_index
            << std::endl << std::endl;
  _scm_system.resize_to_new_n_bfs();
  _scm_system.load_training_set(this->training_set);
  _scm_system.perform_SCM_greedy();

  // and write out the SCM data as well.
  std::stringstream dir_name_stream;
  dir_name_stream << "offline_data_hp" << _tree.leaf_node_index;
  const std::string& directory_name = dir_name_stream.str();
  _scm_system.write_offline_data_to_files(directory_name);
}

} // namespace libMesh

#endif // LIBMESH_HAVE_SLEPC && LIBMESH_HAVE_GLPK
