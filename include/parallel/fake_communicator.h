
#ifndef LIBMESH_FAKE_COMMUNICATOR_H
#define LIBMESH_FAKE_COMMUNICATOR_H

#include "parallel.h"

namespace libMesh
{

class FakeCommunicator
{
  operator const Communicator& () {
    libmesh_error();
    return libMesh::CommWorld;
  }

  operator Communicator& () {
    libmesh_error();
    return libMesh::CommWorld;
  }

};

} // namespace libMesh

#endif // LIBMESH_FAKE_COMMUNICATOR_H
