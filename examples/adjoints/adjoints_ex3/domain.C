
#include "libmesh/mesh_base.h"

#include "domain.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Mesh construction
void build_domain (MeshBase & mesh,
                   FEMParameters & param)
{
  mesh.read(param.domainfile);

  libMesh::out << "Making elements 2nd order" << std::endl;

  // Right now we are setting approximation orders in the code, rather than reading them in
  // That needs to be fixed and the second ordering should be done only if one of the
  // approximation orders is greater than 1
  mesh.all_second_order();
}
