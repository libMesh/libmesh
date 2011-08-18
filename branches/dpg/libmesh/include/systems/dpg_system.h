
// $Id: fem_system.h 4278 2011-03-21 15:23:30Z roystgnr $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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



#ifndef __dpg_system_h__
#define __dpg_system_h__

// C++ includes

// Local Includes
#include "fem_system.h"

namespace libMesh
{

// Forward Declarations
class DiffContext;
class FEMContext;
class DPGContext;


/**
 * This class provides a specific system class.  It aims
 * at nonlinear implicit systems, requiring only a
 * cell residual calculation from the user.  Note
 * that still additional vectors/matrices may be added,
 * as offered in the class \p ExplicitSystem.
 *
 * This class is part of the new DifferentiableSystem framework,
 * which is still experimental.  Users of this framework should
 * beware of bugs and future API changes.
 *
 * @author Truman E. Ellis
 */

// ------------------------------------------------------------
// DPGSystem class definition

class DPGSystem : public FEMSystem
{
public:

  /**
   * Constructor.  Optionally initializes required
   * data structures.
   */
  DPGSystem (EquationSystems& es,
	         const std::string& name,
	         const unsigned int number);

  /**
   * Destructor.
   */
  virtual ~DPGSystem ();

  /**
   * The type of system.
   */
  typedef DPGSystem sys_type;

  /**
   * The type of the parent.
   */
  typedef FEMSystem Parent;
  
  /**
   * Clear all the data structures associated with
   * the system. 
   */
  virtual void clear ();

  /**
   * Prepares \p matrix or \p rhs for matrix assembly.
   * Users may reimplement this to add pre- or post-assembly
   * code before or after calling DPGSystem::assembly()
   */
  virtual void assembly (bool get_residual, bool get_jacobian);

  /**
   * Invokes the solver associated with the system.  For steady state
   * solvers, this will find a root x where F(x) = 0.  For transient
   * solvers, this will integrate dx/dt = F(x).
   *
   * For moving mesh systems, this also translates the mesh to the
   * solution position.
   */
  virtual void solve ();

  /**
   * Builds a FEMContext object with enough information to do
   * evaluations on each element.
   *
   * For most problems, the default DPGSystem implementation is correct; users
   * who subclass FEMContext will need to also reimplement this method to build
   * it.
   */
  virtual AutoPtr<DiffContext> build_context();

  /*
   * Prepares the result of a build_context() call for use.
   * 
   * Most DPGSystem-based problems will need to reimplement this in order to
   * call FE::get_*() as their particular physics requires.
   */
  virtual void init_context(DiffContext &);

  /**
   * @returns the number of test variables in the system
   */
  unsigned int n_test_vars() const;

  /**
   * Adds the test variable \p var to the list of test variables
   * for this system.  Returns the index number for the new test variable.
   */
  unsigned int add_test_variable (const std::string& var,
      const FEType& type,
      const std::set<subdomain_id_type> * const active_subdomains = NULL);

  /**
   * Adds the test variable \p var to the list of test variables
   * for this system.  Same as before, but assumes \p LAGRANGE
   * as default value for \p FEType.family.
   */
  unsigned int add_test_variable (const std::string& var,
      const Order order = FIRST,
      const FEFamily = LAGRANGE,
      const std::set<subdomain_id_type> * const active_subdomains = NULL);

  /** 
   * Return a constant reference to test \p Variable \p var.
   */
  const Variable & test_variable (unsigned int var) const;

  /**
   * @returns true if a test variable named \p var exists in this DPGSystem
   */
  bool has_test_variable(const std::string& var) const;
  
  /**
   * @returns the name of test variable \p i.
   */
  const std::string & test_variable_name(const unsigned int i) const;
  
  /**
   * @returns the variable number assoicated with
   * the user-specified test variable named \p var.
   */
  unsigned short int test_variable_number (const std::string& var) const;

  /**
   * @returns the finite element type test variable number \p i.
   */
  const FEType & test_variable_type (const unsigned int i) const;

  /**
   * @returns the finite element type for test variable \p var.
   */
  const FEType & test_variable_type (const std::string& var) const;
 
  /**
   * Fills the std::set with the degrees of freedom on the local
   * processor corresponding the the variable number passed in.
   */
  // void local_test_dof_indices (const unsigned int var,
  //                         std::set<unsigned int> & var_indices) const;
 
protected:
  /**
   * Initializes the member data fields associated with
   * the system, so that, e.g., \p assemble() may be used.
   */
  virtual void init_data ();

private:
  /**
   * The test \p Variables in this \p DPGSystem.
   */
  std::vector<Variable> _test_variables;

  /**
   * The test variable numbers corresponding to user-specified
   * names, useful for name-based lookups.
   */
  std::map<std::string, unsigned short int> _test_variable_numbers;
};



// ------------------------------------------------------------
// DPGSystem inline methods



inline
unsigned int DPGSystem::n_test_vars() const
{
  return _test_variables.size();
}



inline
const Variable & DPGSystem::test_variable (const unsigned int i) const
{
  libmesh_assert (i < _test_variables.size());

  return _test_variables[i];
}



inline
const std::string & DPGSystem::test_variable_name (const unsigned int i) const
{
  libmesh_assert (i < _test_variables.size());

  return _test_variables[i].name();
}



inline
const FEType & DPGSystem::test_variable_type (const unsigned int i) const
{
  libmesh_assert (i < _test_variables.size());
  
  return _test_variables[i].type();
}



inline
const FEType & DPGSystem::test_variable_type (const std::string& var) const
{
  return _test_variables[this->test_variable_number(var)].type();
}

} // namespace libMesh

#endif
