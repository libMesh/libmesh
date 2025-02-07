// local includes
#include "advection_system.h"

// LibMesh includes
#include "libmesh/equation_systems.h" // EquationSystems::comm()
#include "libmesh/getpot.h" // GetPot input file parsing
#include "libmesh/libmesh_logging.h" // LOG_SCOPE
#include "libmesh/numeric_vector.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/string_to_enum.h"

namespace libMesh
{

AdvectionSystem::AdvectionSystem (EquationSystems& es,
		                  const std::string& name,
		                  const unsigned int number)
  : Parent(es, name, number),
    _q1_var(0),
    _fe_order(CONSTANT),
    _fe_family(MONOMIAL)
{
  // Allocate NumericVectors in _Fh. I could not figure out
  // how to do this in the initialization list, I don't think it's
  // possible.
  // TODO: the number of problem dimensions is hard-coded here, we should
  // make this depend on the input file parameters instead.
  for (unsigned int i=0; i<2; ++i)
    _Fh.push_back(NumericVector<Number>::build(es.comm()));
}


AdvectionSystem::~AdvectionSystem() = default;

void AdvectionSystem::assemble_claw_rhs (NumericVector<Number> & q)
{
  LOG_SCOPE("assemble_claw_rhs()", "AdvectionSystem");

  // The input to this function is the solution vector, "q", from
  // either the initial condition or the previous timestep.
  this->update_Fh(q);

  // Allocate storage for temporary vector used in the computations below
  auto temp = NumericVector<Number>::build(this->comm());
  temp->init(this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);

  rhs->zero();

  // rhs += (Ai - AVG_i)*Fi
  for (auto i : index_range(_Fh))
  {
    this->get_advection_matrix(i).vector_mult(*temp, *_Fh[i]); // temp = Ai*Fi
    rhs->add(+1., *temp);
    this->get_avg_matrix(i).vector_mult(*temp, *_Fh[i]); // temp = AVG_i*Fi
    rhs->add(-1., *temp);
  }

  // Jump term
  // rhs -= LxF * (J*q)
  this->get_jump_matrix().vector_mult(*temp, q);
  rhs->add(-get_LxF_constant(), *temp);

  // Apply boundary conditions. This is also expressed as a matrix-vector product
  // rhs -= BC_i * Fh[i]
  for (auto i : index_range(_Fh))
  {
    this->get_boundary_condition_matrix(i).vector_mult(*temp, *_Fh[i]); // temp = BC_i * Fi
    rhs->add(-1., *temp);
  }
}

void AdvectionSystem::update_Fh (NumericVector<Number> & q)
{
  LOG_SCOPE("update_Fh()", "AdvectionSystem");

  // Fi = u(i)*q
  for (auto i : index_range(_Fh))
  {
    _Fh[i]->zero();
    _Fh[i]->add(_u(i), q);
    _Fh[i]->close();
  }
}

void AdvectionSystem::init_data ()
{
  // Print info about the type of discretization being used.  Note: we
  // can use the same exact assembly code while changing only the
  // approximation order to test both "DG" and "FV" discretizations.
  // In the FV case, the "interior" integral contributions are zero
  // since they depend on "dphi".
  libMesh::out << "Adding q1 variable using ("
               << Utility::enum_to_string(_fe_order)
               << ", "
               << Utility::enum_to_string(_fe_family)
               << ") approximation to the system."
               << std::endl;

  _q1_var = this->add_variable ("q1", _fe_order, _fe_family);

  Parent::init_data();

  for (auto & vec : _Fh)
    vec->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);
}

void AdvectionSystem::process_parameters_file (const std::string& parameters_filename)
{
  Parent::process_parameters_file(parameters_filename);

  // First read in data from parameters_filename
  GetPot infile(parameters_filename);
  const Real u1_in = infile("u1", 1.);
  const Real u2_in = infile("u2", 1.);
  _u = Point(u1_in, u2_in, 0.);

  std::string fe_order_str = infile("fe_order", std::string("CONSTANT"));
  std::string fe_family_str = infile("fe_family", std::string("MONOMIAL"));

  // Convert input strings to enumerations
  _fe_order = Utility::string_to_enum<Order>(fe_order_str);
  _fe_family = Utility::string_to_enum<FEFamily>(fe_family_str);
}

void AdvectionSystem::print_info ()
{
  Parent::print_info();

  libMesh::out << std::endl << "==== AdvectionSystem ====" << std::endl;
  libMesh::out << "u1 = " << _u(0) << ", u2 = " << _u(1) << std::endl;
  libMesh::out << std::endl;
}

} // namespace libMesh



