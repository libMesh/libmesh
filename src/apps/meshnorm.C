// Open the mesh and solution file named on standard input, find all
// variables therein, and output their norms and seminorms.

#include "libmesh/libmesh.h"

#include "libmesh/mesh.h"
#include "libmesh/equation_systems.h"

using namespace libMesh;

void output_norms(const System & sys, const NumericVector<Number> & vec, const std::string & vecname)
{
  for (unsigned int k = 0; k != sys.n_vars(); ++k)
    {
      libMesh::out << "Norms in system " << sys.name() <<
        ", vector " << vecname <<
        ", variable " << sys.variable_name(k) << ":" << std::endl;
      Real l1_vecnorm = sys.calculate_norm(vec, k, DISCRETE_L1);
      libMesh::out << "l1     norm: " << l1_vecnorm << std::endl;
      if (l1_vecnorm)
        {
          libMesh::out << "l2     norm: " << sys.calculate_norm(vec, k, DISCRETE_L2) << std::endl;
          libMesh::out << "linf   norm: " << sys.calculate_norm(vec, k, DISCRETE_L_INF) << std::endl;
          libMesh::out << "H1     norm: " << sys.calculate_norm(vec, k, H1) << std::endl;
          libMesh::out << "L1     norm: " << sys.calculate_norm(vec, k, L1) << std::endl;
          libMesh::out << "L2     norm: " << sys.calculate_norm(vec, k, L2) << std::endl;
          libMesh::out << "Linf   norm: " << sys.calculate_norm(vec, k, L_INF) << std::endl;
          libMesh::out << "H1 seminorm: " << sys.calculate_norm(vec, k, H1_SEMINORM) << std::endl;
          libMesh::out << "H1     norm: " << sys.calculate_norm(vec, k, H1) << std::endl;
        }
    }
}

int main(int argc, char ** argv)
{
  LibMeshInit init (argc, argv);

  Mesh mesh(init.comm());
  EquationSystems es(mesh);

  libMesh::out << "Usage: " << argv[0]
               << " mesh solution" << std::endl;

  libMesh::out << "Loading..." << std::endl;

  mesh.read(argv[1]);
  libMesh::out << "Loaded mesh " << argv[1] << std::endl;
  mesh.print_info();

  es.read(argv[2], EquationSystems::READ_HEADER |
          EquationSystems::READ_DATA |
          EquationSystems::READ_ADDITIONAL_DATA |
          EquationSystems::READ_BASIC_ONLY);
  libMesh::out << "Loaded solution " << argv[2] << std::endl;
  es.print_info();

  libMesh::out.precision(16);

  for (unsigned int i = 0; i != es.n_systems(); ++i)
    {
      System & sys = es.get_system(i);

      output_norms(sys, *sys.solution, std::string("solution"));
      for (unsigned int j = 0; j != sys.n_vectors(); ++j)
        output_norms(sys, sys.get_vector(j), sys.vector_name(j));
    }
}
