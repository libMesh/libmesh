#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_HAVE_DTK

#ifndef DTKEVALUATOR_H
#define DTKEVALUATOR_H

#include "libmesh/equation_systems.h"
#include "libmesh/mesh.h"
#include "libmesh/system.h"

#include <DTK_MeshContainer.hpp>
#include <DTK_FieldEvaluator.hpp>
#include <DTK_FieldContainer.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>

#include <string>

namespace libMesh {

class DTKEvaluator : public DataTransferKit::FieldEvaluator<int,DataTransferKit::FieldContainer<double> >
{
public:
  typedef DataTransferKit::MeshContainer<int>	      MeshContainerType;
  typedef DataTransferKit::FieldContainer<Number>     FieldContainerType;

  DTKEvaluator(System & in_sys, std::string var_name);

  FieldContainerType evaluate( const Teuchos::ArrayRCP<int>& elements,
                               const Teuchos::ArrayRCP<double>& coords );

protected:
  System & sys;
  NumericVector<Number> & current_local_solution;
  EquationSystems & es;
  MeshBase & mesh;
  unsigned int dim;
  DofMap & dof_map;
  unsigned int var_num;
  const FEType& fe_type;
};

} // namespace libMesh

#endif

#endif
