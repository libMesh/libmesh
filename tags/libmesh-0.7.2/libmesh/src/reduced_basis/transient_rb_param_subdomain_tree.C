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

#include "transient_rb_param_subdomain_tree.h"
#include "transient_rb_param_subdomain_node.h"
#include "getpot.h"

namespace libMesh
{

TransientRBParamSubdomainTree::TransientRBParamSubdomainTree
  (TransientRBConstruction& rb_construction_in, const std::string& parameters_filename)
  : Parent(rb_construction_in, parameters_filename),
    conserv_factor(-1.),
    use_delta_N_in_h_stage(false)
{
  // Read in properties from parameters_filename
  GetPot infile(parameters_filename);

  conserv_factor = infile("conserv_factor",conserv_factor);
  libmesh_assert(conserv_factor >= 1.);

  // Set the boolean that indicates whether we use delta_N or the POD_tol
  // in order to determine how many basis functions we add in the h-stage
  // refinement for a time-dependent problem.
  use_delta_N_in_h_stage = infile("use_delta_N_in_h_stage",use_delta_N_in_h_stage);

  libMesh::out << std::endl << "TransientRBParamSubdomainTree parameters:" << std::endl;
  libMesh::out << "Conservativity factor: " << conserv_factor << std::endl;
  libMesh::out << "Using delta_N in h-stage: " << use_delta_N_in_h_stage << std::endl;
  libMesh::out << std::endl;
}

void TransientRBParamSubdomainTree::build_root_node()
{
  libmesh_assert(root_node == NULL);
  root_node = new TransientRBParamSubdomainNode(*this, _rb_construction.get_current_parameters());
}

} // namespace libMesh
