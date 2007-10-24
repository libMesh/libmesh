// $Id: system_base.h,v 1.1 2003-11-05 22:26:42 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002-2003  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __system_base_h__
#define __system_base_h__

// C++ includes

// Local Includes
#include "libmesh_common.h"
#include "fe_type.h"
#include "dof_map.h"
#include "auto_ptr.h"
#include "reference_counted_object.h"


// Forward Declarations
class SystemBase;
class EquationSystems;
class Mesh;
class Xdr;
template <typename T> class SparseMatrix;
template <typename T> class NumericVector;
template <typename T> class LinearSolverInterface;

/**
 * This is the base class for classes which contain 
 * information related to any physical process that might be simulated.  
 * Such information may range from the actual solution values to 
 * algorithmic flags that may be used to control the numerical methods
 * employed.  In general, use an \p EqnSystems<T_sys> object to handle
 * one or more of the children of this class.
 * Note that templating \p EqnSystems relaxes the use of virtual members.
 *
 * @author Benjamin S. Kirk, 2002-2003.
 */

// ------------------------------------------------------------
// SystemBase class definition
class SystemBase : public ReferenceCountedObject<SystemBase>
{
protected:

  /**
   * Constructor.  Optionally initializes required
   * data structures.  Protected so that this base class
   * cannot be explicitly instantiated.
   */
  SystemBase (EquationSystems& es,
	      const std::string& name,
	      const unsigned int number);
  
public:

  
  /**
   * Destructor.
   */
  virtual ~SystemBase ();

  /**
   * Clear all the data structures associated with
   * the system.  Protected so that only children can
   * use it.
   */
  virtual void clear ();
  
  /**
   * Initializes degrees of freedom and other 
   * required data on the current mesh.  Note that the matrix
   * is not initialized at this time since it may not be required
   * for all applications. @e Should be overloaded in derived classes.
   * Protected so that only children can use it.
   */
  void init ();
  
  /**
   * Reinitializes degrees of freedom and other 
   * required data on the current mesh.  Note that the matrix
   * is not initialized at this time since it may not be required
   * for all applications. @e Should be overloaded in derived classes.
   * Protected so that only children can use it.
   */
  virtual void reinit ();
  
  /**
   * Update the local values to reflect the solution
   * on neighboring processors.
   */
  virtual void update () = 0;
  
  /**
   * Prepares \p matrix and \p _dof_map for matrix assembly.
   * Does not actually assemble anything.  For matrix assembly,
   * use the \p assemble() in derived classes.
   * @e Should be overloaded in derived classes.
   * Protected so that only children can use it.
   */
  virtual void assemble ();

  /**
   * Assembles & solves the linear system Ax=b. 
   */
  virtual std::pair<unsigned int, Real> solve ();
  
  /**
   * @returns \p true when the other system contains
   * identical data, up to the given threshold.  Outputs
   * some diagnostic info when \p verbose is set.
   */
  virtual bool compare (const SystemBase& other_system, 
			const Real threshold,
			const bool verbose) const;

  /**
   * @returns the system name.
   */
  const std::string & name () const;

  /**
   * @returns the type of system, helpful in identifying
   * which system type to use when reading equation system
   * data from file.  Has be overloaded in derived classes.
   */
  virtual std::string system_type () const = 0;

  /**
   * Projects the vector defined on the old mesh onto the
   * new mesh. 
   */
  void project_vector (NumericVector<Number>&) const;

  /**
   * Projects the vector defined on the old mesh onto the
   * new mesh. The original vector is unchanged and the new vector
   * is passed through the second argument.
   */
  void project_vector (const NumericVector<Number>&,
		       NumericVector<Number>&) const;
  
  /**
   * @returns the system number.   
   */
  unsigned int number () const;
  
  /**
   * Fill the input vector \p global_soln so that it contains
   * the global solution on all processors.
   * Requires communication with all other processors.
   */
  void update_global_solution (std::vector<Number>& global_soln) const;

  /**
   * Fill the input vector \p global_soln so that it contains
   * the global solution on processor \p dest_proc.
   * Requires communication with all other processors.
   */
  void update_global_solution (std::vector<Number>& global_soln,
			       const unsigned int dest_proc) const;

  /**
   * @returns a constant reference to this systems's \p _mesh.
   */
  const Mesh & get_mesh() const;
  
  /**
   * @returns a constant reference to this system's \p _dof_map.
   */
  const DofMap & get_dof_map() const;
  
  /**
   * @returns a writeable reference to this system's \p _dof_map.
   */
  DofMap & get_dof_map();

  /**
   * @returns \p true if the system is active, false otherwise.
   * An active system will be solved.
   */
  bool active () const;

  /**
   * Activates the system.  Only active systems are solved.
   */
  void activate ();

  /**
   * Deactivates the system.  Only active systems are solved.
   */
  void deactivate ();
    
  /**
   * Adds the additional matrix \p mat_name to this system.  Only
   * allowed @e prior to \p assemble().  All additional matrices
   * have the same sparsity pattern as the matrix used during
   * solution.  When not \p SystemBase but the @e user wants to
   * initialize the mayor matrix, then all the additional matrices,
   * if existent, have to be initialized by the user, too.
   */
  void add_matrix (const std::string& mat_name);

  /**
   * @returns a const reference to this system's @e additional matrix
   * named \p mat_name.  @e None of these matrices is involved in the 
   * solution process.  Access is only granted when the matrix is already
   * properly initialized.
   */
  const SparseMatrix<Number> & get_matrix(const std::string& mat_name) const;

  /**
   * @returns a writeable reference to this system's @e additional matrix
   * named \p mat_name.  @e None of these matrices is involved in the 
   * solution process.  Access is only granted when the matrix is already
   * properly initialized.
   */
  SparseMatrix<Number> & get_matrix(const std::string& mat_name);

  /**
   * @returns the number of additional vectors handled by this system
   */
  unsigned int n_additional_matrices () const;

  /**
   * Adds the additional vector \p vec_name to this system.  Only
   * allowed @e prior to \p init().  All the additional vectors
   * are similarly distributed, like \p rhs and \p solution,
   * and inititialized to zero.
   */
  void add_vector (const std::string& vec_name);

  /**
   * @returns a const reference to this system's @e additional vector
   * named \p vec_name.  Access is only granted when the vector is already
   * properly initialized.
   */
  const NumericVector<Number> & get_vector(const std::string& vec_name) const;

  /**
   * @returns a writeable reference to this system's @e additional vector
   * named \p vec_name.  Access is only granted when the vector is already
   * properly initialized.
   */
  NumericVector<Number> & get_vector(const std::string& vec_name);

  /**
   * @returns the number of additional vectors handled by this system
   */
  unsigned int n_additional_vectors () const;

  /**
   * @returns the number of variables in the system
   */
  unsigned int n_vars() const;

  /**
   * @returns the number of degrees of freedom in the system
   */
  unsigned int n_dofs() const;

  /**
   * @returns the number of constrained degrees of freedom
   * in the system.
   */
  unsigned int n_constrained_dofs() const;
  
  /**
   * @returns the number of degrees of freedom local
   * to this processor
   */
  unsigned int n_local_dofs() const;
  
  /**
   * Adds the variable \p var to the list of variables
   * for this system.
   */
  void add_variable (const std::string& var,
		     const FEType& type);

  /**
   * Adds the variable \p var to the list of variables
   * for this system.  Same as before, but assumes \p LAGRANGE
   * as default value for \p FEType.family.
   */
  void add_variable (const std::string& var,
		     const Order order);

  /**
   * @returns the name of variable \p i.
   */
  const std::string & variable_name(const unsigned int i) const;
  
  /**
   * @returns the variable number assoicated with
   * the user-specified variable named \p var.
   */
  unsigned short int variable_number (const std::string& var) const;

  /**
   * @returns the finite element type variable number \p i.
   */
  const FEType & variable_type (const unsigned int i) const;

  /**
   * @returns the finite element type for variable \p var.
   */
  const FEType & variable_type (const std::string& var) const;

  /**
   * Reads the basic data for this System.
   */
  void read (Xdr& io, 
	     const bool read_header=true,
	     const bool read_additional_data=true);

  /**
   * Reads additional data, namely vectors, for this System.
   */
  void read_data (Xdr& io,
		  const bool read_additional_data=true);

  /**
   * Writes the basic data for this System.
   */
  void write (Xdr& io) const;

  /**
   * Writes additional data, namely vectors, for this System.
   */
  void write_data (Xdr& io,
		   const bool write_additional_data = true) const;

  /**
   * @returns a string containing information about the
   * system.
   */
  std::string get_info () const;

  /**
   * Register a user function to use in initializing the system.
   */
  void attach_init_function (void fptr(EquationSystems& es,
				       const std::string& name));
  
  /**
   * Register a user function to use in assembling the system
   * matrix and RHS.
   */
  void attach_assemble_function (void fptr(EquationSystems& es,
					   const std::string& name));
  
  /**
   * Function that initializes the system.
   */
  void (* init_system) (EquationSystems& es,
			const std::string& name);
  
  /**
   * Function that assembles the system.
   */
  void (* assemble_system) (EquationSystems& es,
			    const std::string& name);



  //--------------------------------------------------
  // The data & solver interface for solving the linear
  // system Ax=b
  
  /**
   * Data structure to hold solution values.
   */
  AutoPtr<NumericVector<Number> > solution;

  /**
   * Data structure to hold the system right-hand-side.
   */
  AutoPtr<NumericVector<Number> > rhs;
  
  /**
   * Data structure to hold the system matrix.
   */
  AutoPtr<SparseMatrix<Number> > matrix;

  /**
   * Interface to the Petsc solvers.
   */
  AutoPtr<LinearSolverInterface<Number> > linear_solver_interface;


protected:

  
  /**
   * Initializes the data for the system.  Note that this is called
   * before any user-supplied intitialization function so that all
   * required storage will be available.
   */
  virtual void init_data ();

  /**
   * Re-update the local values when the mesh has changed.
   * This method takes the data updated by \p update() and
   * makes it up-to-date on the current mesh.
   */
  virtual void re_update () = 0;
  
  /**
   * A name associated with this system.
   */
  const std::string _sys_name;

  /**
   * The number associated with this system
   */
  const unsigned int _sys_number;

  /**
   * The names of the variables associated with this
   * system.
   */
  std::vector<std::string> _var_names;
  
  /**
   * The finite element type for the variables associated
   * with this system.
   */
  std::map<std::string, FEType> _var_type;
  
  /**
   * The variable numbers corresponding to user-specified
   * names.
   */
  std::map<std::string, unsigned short int> _var_num;

  /**
   * Flag stating if the system is active or not.
   */
  bool _active;
  
  /**
   * Some systems need multiple matrices.
   */
  std::map<std::string, SparseMatrix<Number>* > _other_matrices;

  /**
   * \p true when additional matrices may still be added, \p false otherwise.
   */
  bool _can_add_matrices;
  
  /**
   * Some systems need multiple vectors.  
   */
  std::map<std::string, NumericVector<Number>* > _other_vectors;

  /**
   * \p true when additional vectors may still be added, \p false otherwise.
   */
  bool _can_add_vectors;



  // -------------------------------------------------
  // Necessary classes
  //
  
  /**
   * Data structure describing the relationship between
   * nodes, variables, etc... and degrees of freedom.
   */
  DofMap _dof_map;

  /**
   * Constant reference to the \p EquationSystems object
   * used for the simulation.
   */
  EquationSystems& _equation_systems;
  
  /**
   * Constant reference to the \p mesh data structure used
   * for the simulation.   
   */
  const Mesh& _mesh;
};



// ------------------------------------------------------------
// SystemBase inline methods
inline
const std::string & SystemBase::name() const
{
  return _sys_name;
}



inline
unsigned int SystemBase::number() const
{
  return _sys_number;
}



inline
const Mesh& SystemBase::get_mesh() const
{
  return _mesh;
}



inline
const DofMap& SystemBase::get_dof_map() const
{
  return _dof_map;
}



inline
DofMap& SystemBase::get_dof_map()
{
  return _dof_map;
}



inline
bool SystemBase::active() const
{
  return _active;  
}



inline
void SystemBase::activate ()
{
  _active = true;
}



inline
void SystemBase::deactivate ()
{
  _active = false;
}



inline
unsigned int SystemBase::n_vars() const
{
//  assert (!_var_names.empty());

  return _var_names.size();
}



inline
const std::string & SystemBase::variable_name (const unsigned int i) const
{
  assert (i < n_vars());

  return _var_names[i];
}



inline
const FEType & SystemBase::variable_type (const unsigned int i) const
{
  return variable_type(_var_names[i]);
}



inline
const FEType & SystemBase::variable_type (const std::string& var) const
{
  std::map<std::string, FEType>::const_iterator
    pos = _var_type.find(var);

  assert (pos != _var_type.end());
  
  return pos->second;
}



inline
unsigned int SystemBase::n_dofs() const
{
  return _dof_map.n_dofs();
}



inline
unsigned int SystemBase::n_constrained_dofs() const
{
#ifdef ENABLE_AMR

  return _dof_map.n_constrained_dofs();

#else

  std::cerr << "This feature requires AMR!" << std::endl;
  error();
  return 0;

#endif
}



inline
unsigned int SystemBase::n_local_dofs() const
{
  return _dof_map.n_dofs_on_processor(libMesh::processor_id());
}



inline
unsigned int SystemBase::n_additional_matrices () const
{
 return _other_matrices.size(); 
}



inline
unsigned int SystemBase::n_additional_vectors () const
{
 return _other_vectors.size(); 
}


#endif
