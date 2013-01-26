#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_HAVE_DTK

#ifndef DTKADAPTER_H
#define DTKADAPTER_H

#include "libmesh/dtk_evaluator.h"

#include <DTK_MeshManager.hpp>
#include <DTK_MeshContainer.hpp>
#include <DTK_MeshTraits.hpp>
#include <DTK_MeshTraitsFieldAdapter.hpp>
#include <DTK_FieldManager.hpp>
#include <DTK_FieldContainer.hpp>
#include <DTK_FieldEvaluator.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>

class DTKAdapter
{
public:
  DTKAdapter(Teuchos::RCP<const Teuchos::Comm<int> > in_comm, EquationSystems & in_es);

  typedef DataTransferKit::MeshContainer<int>                                  MeshContainerType;
  typedef DataTransferKit::FieldContainer<double>                              FieldContainerType;
  typedef DataTransferKit::MeshTraits<MeshContainerType>::global_ordinal_type  GlobalOrdinal;
  typedef DataTransferKit::FieldEvaluator<GlobalOrdinal,FieldContainerType>    EvaluatorType;
  typedef Teuchos::RCP<EvaluatorType>                                          RCP_Evaluator;


  Teuchos::RCP<DataTransferKit::MeshManager<MeshContainerType> > get_mesh_manager() { return mesh_manager; }
  RCP_Evaluator get_variable_evaluator(std::string var_name);
  Teuchos::RCP<DataTransferKit::FieldManager<MeshContainerType> > get_target_coords() { return target_coords; }
  Teuchos::RCP<DataTransferKit::FieldManager<FieldContainerType> > get_values_to_fill(std::string var_name);

  /**
   * After computing values for a variable in this EquationSystems
   * we need to take those values and put them in the actual solution vector.
   */
  void update_variable_values(std::string var_name);

protected:
  /**
   * Small helper function for finding the system containing the variable.
   *
   * Note that this implies that variable names are unique across all systems!
   */
  System * find_sys(std::string var_name);

  /**
   * Helper that returns the DTK ElementTopology for a given Elem
   */
  DataTransferKit::DTK_ElementTopology get_element_topology(const Elem * elem);

  /**
   * Helper function that fills the std::set with all of the node numbers of
   * nodes connected to local elements.
   */
  void get_semi_local_nodes(std::set<unsigned int> & semi_local_nodes);
  
  Teuchos::RCP<const Teuchos::Comm<int> > comm;
  EquationSystems & es;
  const MeshBase & mesh;
  unsigned int dim;

  unsigned int num_local_nodes;
  Teuchos::ArrayRCP<int> vertices;
  
  Teuchos::RCP<DataTransferKit::MeshManager<MeshContainerType> > mesh_manager;
  RCP_Evaluator field_evaluator;
  Teuchos::RCP<DataTransferKit::FieldManager<MeshContainerType> > target_coords;

  /// Map of variable names to arrays to be filled by a transfer
  std::map<std::string, Teuchos::RCP<DataTransferKit::FieldManager<FieldContainerType> > > values_to_fill;

  /// Map of variable names to RCP_Evaluator objects
  std::map<std::string, RCP_Evaluator> evaluators;
};

#endif

#endif
