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



// C++ includes


// Local includes
#include "libmesh/system_subset_by_subdomain.h"
#include "libmesh/system.h"
#include "libmesh/dof_map.h"
#include "libmesh/parallel.h"
#include "libmesh/elem.h"

namespace libMesh
{
// ------------------------------------------------------------
// SubdomainSelection implementation
SystemSubsetBySubdomain::SubdomainSelection::SubdomainSelection () = default;
SystemSubsetBySubdomain::SubdomainSelection::~SubdomainSelection () = default;

SystemSubsetBySubdomain::SubdomainSelectionByList::
SubdomainSelectionByList (const std::set<subdomain_id_type> & list):
  _list(list)
{
}

bool
SystemSubsetBySubdomain::SubdomainSelectionByList::
operator()(const subdomain_id_type & subdomain_id)const
{
  return _list.count(subdomain_id); // _list is actually a std::set
}

// ------------------------------------------------------------
// SystemSubsetBySubdomain implementation

SystemSubsetBySubdomain::
SystemSubsetBySubdomain (const System & system,
                         const SubdomainSelection & subdomain_selection,
                         const std::set<unsigned int> * const var_nums):
  SystemSubset(system),
  ParallelObject(system),
  _var_nums(),
  _dof_ids()
{
  this->set_var_nums(var_nums);
  this->init(subdomain_selection);
}

SystemSubsetBySubdomain::
SystemSubsetBySubdomain (const System & system,
                         const std::set<subdomain_id_type> & subdomain_ids,
                         const std::set<unsigned int> * const var_nums):
  SystemSubset(system),
  ParallelObject(system),
  _var_nums(),
  _dof_ids()
{
  this->set_var_nums(var_nums);
  this->init(subdomain_ids);
}

SystemSubsetBySubdomain::~SystemSubsetBySubdomain () = default;

const std::vector<unsigned int> &
SystemSubsetBySubdomain::dof_ids() const
{
  return _dof_ids;
}

void
SystemSubsetBySubdomain::
set_var_nums (const std::set<unsigned int> * const var_nums)
{
  _var_nums.clear();

  if (var_nums != nullptr)
    _var_nums = *var_nums;

  else
    for (auto i : make_range(_system.n_vars()))
      _var_nums.insert(i);
}

void
SystemSubsetBySubdomain::
init (const SubdomainSelection & subdomain_selection)
{
  _dof_ids.clear();

  std::vector<std::vector<dof_id_type>> dof_ids_per_processor(this->n_processors());

  const DofMap & dof_map = _system.get_dof_map();
  std::vector<dof_id_type> dof_indices;

  const MeshBase & mesh = _system.get_mesh();
  for (const auto & elem : mesh.active_local_element_ptr_range())
    if (subdomain_selection(elem->subdomain_id()))
      {
        for (const auto & var_num : _var_nums)
          {
            dof_map.dof_indices (elem, dof_indices, var_num);
            for (const auto & dof : dof_indices)
              for (auto proc : make_range(this->n_processors()))
                if ((dof>=dof_map.first_dof(proc)) && (dof<dof_map.end_dof(proc)))
                  dof_ids_per_processor[proc].push_back(dof);
          }
      }

  /* Distribute information among processors.  */
  std::vector<Parallel::Request> request_per_processor(this->n_processors());
  for (auto proc : make_range(this->n_processors()))
    if (proc!=this->processor_id())
      this->comm().send(proc,dof_ids_per_processor[proc],request_per_processor[proc]);
  for (auto proc : make_range(this->n_processors()))
    {
      std::vector<dof_id_type> received_dofs;
      if (proc==this->processor_id())
        received_dofs = dof_ids_per_processor[proc];
      else
        this->comm().receive(proc,received_dofs);
      _dof_ids.insert(_dof_ids.end(), received_dofs.begin(), received_dofs.end());
    }

  /* Sort and unique the vector (using the same mechanism as in \p
     DofMap::prepare_send_list()).  */
  std::sort(_dof_ids.begin(), _dof_ids.end());
  std::vector<unsigned int>::iterator new_end = std::unique (_dof_ids.begin(), _dof_ids.end());
  std::vector<unsigned int> (_dof_ids.begin(), new_end).swap (_dof_ids);

  /* Wait for sends to be complete.  */
  for (auto proc : make_range(this->n_processors()))
    {
      if (proc!=this->processor_id())
        {
          request_per_processor[proc].wait();
        }
    }
}

void
SystemSubsetBySubdomain::
init (const std::set<subdomain_id_type> & subdomain_ids)
{
  SubdomainSelectionByList selection(subdomain_ids);
  this->init(selection);
}

} // namespace libMesh
