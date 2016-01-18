#include "libmesh/parameters.h"
#include "libmesh/point.h"
#include "libmesh/vector_value.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

void read_initial_parameters();
void finish_initialization();

Number initial_value(const Point & p,
                     const Parameters &,
                     const std::string &,
                     const std::string &);

Gradient initial_grad(const Point & p,
                      const Parameters &,
                      const std::string &,
                      const std::string &);
