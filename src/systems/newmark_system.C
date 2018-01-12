// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



// C++ includes

// Local includes
#include "libmesh/newmark_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/numeric_vector.h"

namespace libMesh
{




// ------------------------------------------------------------
// NewmarkSystem static members
const Real NewmarkSystem::_default_alpha    = .25;
const Real NewmarkSystem::_default_delta    = .5;
const Real NewmarkSystem::_default_timestep = 1.;



// ------------------------------------------------------------
// NewmarkSystem implementation
NewmarkSystem::NewmarkSystem (EquationSystems & es,
                              const std::string & name_in,
                              const unsigned int number_in) :
  LinearImplicitSystem (es, name_in, number_in),
  _a_0                 (1./(_default_alpha*_default_timestep*_default_timestep)),
  _a_1                 (_default_delta/(_default_alpha*_default_timestep)),
  _a_2                 (1./(_default_alpha*_default_timestep)),
  _a_3                 (1./(2.*_default_alpha)-1.),
  _a_4                 (_default_delta/_default_alpha -1.),
  _a_5                 (_default_timestep/2.*(_default_delta/_default_alpha-2.)),
  _a_6                 (_default_timestep*(1.-_default_delta)),
  _a_7                 (_default_delta*_default_timestep),
  _finished_assemble   (false)

{
  // default values of the newmark parameters
  es.parameters.set<Real>("Newmark alpha") = _default_alpha;
  es.parameters.set<Real>("Newmark delta") = _default_delta;

  // time step size.
  // should be handled at a later stage through EquationSystems?
  es.parameters.set<Real>("Newmark time step") = _default_timestep;

  // add additional matrices and vectors that will be used in the
  // newmark algorithm to the data structure
  // functions LinearImplicitSystem::add_matrix and LinearImplicitSystem::add_vector
  // are used so we do not have to bother about initialization and
  // dof mapping

  // system matrices
  this->add_matrix ("stiffness");
  this->add_matrix ("damping");
  this->add_matrix ("mass");

  // load vector
  this->add_vector ("force");

  // the displacement and the time derivatives
  this->add_vector ("displacement");
  this->add_vector ("velocity");
  this->add_vector ("acceleration");

  // contributions to the rhs
  this->add_vector ("rhs_m");
  this->add_vector ("rhs_c");

  // results from the previous time step
  this->add_vector ("old_solution");
  this->add_vector ("old_acceleration");
}



NewmarkSystem::~NewmarkSystem ()
{
  this->clear();
}



void NewmarkSystem::clear ()
{
  // use parent clear this will also clear the
  // matrices and vectors added in the constructor
  LinearImplicitSystem::clear();

  // Get a reference to the EquationSystems
  EquationSystems & es =
    this->get_equation_systems();

  // default values of the newmark parameters
  es.parameters.set<Real>("Newmark alpha") = _default_alpha;
  es.parameters.set<Real>("Newmark delta") = _default_delta;

  // time step size.  should be handled at a later stage through EquationSystems?
  es.parameters.set<Real>("Newmark time step") = _default_timestep;

  // set bool to false
  _finished_assemble = false;
}



void NewmarkSystem::reinit ()
{
  libmesh_not_implemented();
}



void NewmarkSystem::assemble ()
{
  if (!_finished_assemble)
    {
      // prepare matrix with the help of the _dof_map,
      // fill with sparsity pattern
      LinearImplicitSystem::assemble();

      // compute the effective system matrix
      this->compute_matrix();

      // apply initial conditions
      this->initial_conditions();

      _finished_assemble = true;
    }
}


void NewmarkSystem::initial_conditions ()
{
  // libmesh_assert(init_cond_fptr);

  // Log how long the user's matrix assembly code takes
  LOG_SCOPE("initial_conditions ()", "NewmarkSystem");

  // Set all values to 0, then
  // call the user-specified function for initial conditions.
  this->get_vector("displacement").zero();
  this->get_vector("velocity").zero();
  this->get_vector("acceleration").zero();
  this->user_initialization();
}



void NewmarkSystem::compute_matrix ()
{
  // close the component matrices
  this->get_matrix ("stiffness").close();
  this->get_matrix ("mass"     ).close();
  this->get_matrix ("damping"  ).close();

  // close & zero the system matrix
  this->matrix->close (); this->matrix->zero();

  // add up the matrices
  this->matrix->add (1.,   this->get_matrix ("stiffness"));
  this->matrix->add (_a_0, this->get_matrix ("mass"));
  this->matrix->add (_a_1, this->get_matrix ("damping"));

}



void NewmarkSystem::update_rhs ()
{
  LOG_SCOPE("update_rhs ()", "NewmarkSystem");

  // zero the rhs-vector
  NumericVector<Number> & the_rhs = *this->rhs;
  the_rhs.zero();

  // get writable references to some vectors
  NumericVector<Number> & rhs_m = this->get_vector("rhs_m");
  NumericVector<Number> & rhs_c = this->get_vector("rhs_c");

  // zero the vectors for matrix-vector product
  rhs_m.zero();
  rhs_c.zero();

  // compute auxiliary vectors rhs_m and rhs_c
  rhs_m.add(_a_0, this->get_vector("displacement"));
  rhs_m.add(_a_2, this->get_vector("velocity"));
  rhs_m.add(_a_3, this->get_vector("acceleration"));

  rhs_c.add(_a_1, this->get_vector("displacement"));
  rhs_c.add(_a_4, this->get_vector("velocity"));
  rhs_c.add(_a_5, this->get_vector("acceleration"));

  // compute rhs
  the_rhs.add(this->get_vector("force"));
  the_rhs.add_vector(rhs_m, this->get_matrix("mass"));
  the_rhs.add_vector(rhs_c, this->get_matrix("damping"));
}



void NewmarkSystem::update_u_v_a ()
{
  LOG_SCOPE("update_u_v_a ()", "NewmarkSystem");

  // get some references for convenience
  const NumericVector<Number> &  solu = *this->solution;

  NumericVector<Number> &  disp_vec   = this->get_vector("displacement");
  NumericVector<Number> &  vel_vec    = this->get_vector("velocity");
  NumericVector<Number> &  acc_vec    = this->get_vector("acceleration");
  NumericVector<Number> &  old_acc    = this->get_vector("old_acceleration");
  NumericVector<Number> &  old_solu   = this->get_vector("old_solution");

  // copy data
  old_solu = disp_vec;
  disp_vec = solu;
  old_acc  = acc_vec;

  // compute the new acceleration vector
  acc_vec.scale(-_a_3);
  acc_vec.add(_a_0, disp_vec);
  acc_vec.add(-_a_0, old_solu);
  acc_vec.add(-_a_2,vel_vec);

  // compute the new velocity vector
  vel_vec.add(_a_6,old_acc);
  vel_vec.add(_a_7,acc_vec);
}



void NewmarkSystem::set_newmark_parameters (const Real delta_T,
                                            const Real alpha,
                                            const Real delta)
{
  libmesh_assert_not_equal_to (delta_T, 0.);

  // Get a reference to the EquationSystems
  EquationSystems & es =
    this->get_equation_systems();

  // the newmark parameters
  es.parameters.set<Real>("Newmark alpha") = alpha;
  es.parameters.set<Real>("Newmark delta") = delta;

  // time step size.
  // should be handled at a later stage through EquationSystems?
  es.parameters.set<Real>("Newmark time step") = delta_T;

  // the constants for time integration
  _a_0 = 1./(alpha*delta_T*delta_T);
  _a_1 = delta/(alpha*delta_T);
  _a_2 = 1./(alpha*delta_T);
  _a_3 = 1./(2.*alpha)-1.;
  _a_4 = delta/alpha -1.;
  _a_5 = delta_T/2.*(delta/alpha-2.);
  _a_6 = delta_T*(1.-delta);
  _a_7 = delta*delta_T;
}

} // namespace libMesh
