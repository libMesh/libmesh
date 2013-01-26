#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_HAVE_DTK

#ifndef DTKSOLUTIONTRANSFER_H
#define DTKSOLUTIONTRANSFER_H

#include "libmesh/solution_transfer.h"
#include "libmesh/dtk_adapter.h"

// Trilinos
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_DefaultComm.hpp>

// DTK
#include <DTK_SharedDomainMap.hpp>

#include <string>

/**
 * Implementation of a SolutionTransfer object that uses the DataTransferKit (https://github.com/CNERG/DataTransferKit) to transfer variables back and forth between systems.
 */
class DTKSolutionTransfer : public SolutionTransfer
{
public:
  DTKSolutionTransfer();
  virtual ~DTKSolutionTransfer();
  
  /**
   * Transfer the values of a variable to another.
   * 
   * This is meant for transferring values from one EquationSystems to another
   * even in the case of having different meshes.
   *
   * Note that the first time this function is called for one combination of EquationSystems
   * a lot of setup and caching is done.  Subsequent transfers between the same EquationSystems
   * will be _much_ faster.
   */
  virtual void transfer(const Variable & from_var, const Variable & to_var);

protected:
  typedef DataTransferKit::SharedDomainMap<DTKAdapter::MeshContainerType,DTKAdapter::MeshContainerType> shared_domain_map_type;

  /// COMM_WORLD for now
  Teuchos::RCP<const Teuchos::Comm<int> > comm_default;
  
  /// The DTKAdapter associated with each EquationSystems
  std::map<EquationSystems *, DTKAdapter *> adapters;

  /// The dtk shared domain maps for pairs of EquationSystems (from, to)
  std::map<std::pair<EquationSystems *, EquationSystems *>, shared_domain_map_type * > dtk_maps;
};

#endif

#endif
