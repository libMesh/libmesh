// Function Definitions for NoTimeHistory

// Local includes
#include "libmesh/no_solution_history.h"

namespace libMesh
{

void NoSolutionHistory::store()
{
  // Do nothing
}

void NoSolutionHistory::retrieve()
{
  // Nothing was stored, so nothing can be retrieved
  libmesh_not_implemented();
}

}
