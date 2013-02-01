// C++ include files that we need
#include "libmesh/parameters.h"
#include "libmesh/point.h"
#include "libmesh/vector_value.h"

using namespace libMesh;

void adjoint_read_initial_parameters();
void adjoint_finish_initialization();

Number adjoint_initial_value(const Point& p,
                     const Parameters&,
                     const std::string&,
                     const std::string&);

Gradient adjoint_initial_grad(const Point& p,
                      const Parameters&,
                      const std::string&,
                      const std::string&);
