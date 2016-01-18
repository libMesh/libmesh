// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef DTKSOLUTIONTRANSFER_H
#define DTKSOLUTIONTRANSFER_H

#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_TRILINOS_HAVE_DTK

#include "libmesh/solution_transfer.h"
#include "libmesh/dtk_adapter.h"

// Trilinos
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_DefaultComm.hpp>

// DTK
#include <DTK_SharedDomainMap.hpp>

#include <string>

namespace libMesh
{

/**
 * Implementation of a SolutionTransfer object that uses the
 * DataTransferKit (https://github.com/CNERG/DataTransferKit) to
 * transfer variables back and forth between systems.
 *
 * \author Derek Gaston
 * \date 2013
 */
class DTKSolutionTransfer : public SolutionTransfer
{
public:
  DTKSolutionTransfer(const libMesh::Parallel::Communicator & comm
                      LIBMESH_CAN_DEFAULT_TO_COMMWORLD);
  virtual ~DTKSolutionTransfer();

  /**
   * Transfer the values of a variable to another.
   *
   * This is meant for transferring values from one EquationSystems to
   * another even in the case of having different meshes.
   *
   * Note that the first time this function is called for one
   * combination of EquationSystems a lot of setup and caching is
   * done.  Subsequent transfers between the same EquationSystems will
   * be _much_ faster.
   */
  virtual void transfer(const Variable & from_var, const Variable & to_var) libmesh_override;

protected:
  typedef DataTransferKit::SharedDomainMap<DTKAdapter::MeshContainerType,DTKAdapter::MeshContainerType> shared_domain_map_type;

  /// COMM_WORLD for now
  Teuchos::RCP<const Teuchos::Comm<int> > comm_default;

  /// The DTKAdapter associated with each EquationSystems
  std::map<EquationSystems *, DTKAdapter *> adapters;

  /// The dtk shared domain maps for pairs of EquationSystems (from, to)
  std::map<std::pair<EquationSystems *, EquationSystems *>, shared_domain_map_type * > dtk_maps;
};

} // namespace libMesh

#endif // #ifdef LIBMESH_TRILINOS_HAVE_DTK

#endif // #define DTKSOLUTIONTRANSFER_H
