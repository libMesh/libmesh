// $Id$

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



#ifndef __system_h__
#define __system_h__

// C++ includes
#include <vector>

// Local Includes
#include "auto_ptr.h"
#include "enum_xdr_mode.h"
#include "fe_type.h"
#include "libmesh_common.h"
#include "reference_counted_object.h"
#include "system_norm.h"
#include "elem_range.h"

// Forward Declarations
class System;
class EquationSystems;
class MeshBase;
class Xdr;
class DofMap;
class Parameters;
class Point;
template <typename T> class NumericVector;
template <typename T> class VectorValue;
typedef VectorValue<Number> NumberVectorValue;
typedef NumberVectorValue Gradient;

/**
 * This is the base class for classes which contain 
 * information related to any physical process that might be simulated.  
 * Such information may range from the actual solution values to 
 * algorithmic flags that may be used to control the numerical methods
 * employed.  In general, use an \p EqnSystems<T_sys> object to handle
 * one or more of the children of this class.
 * Note that templating \p EqnSystems relaxes the use of virtual members.
 *
 * @author Benjamin S. Kirk, 2003-2004.
 */

// ------------------------------------------------------------
// System class definition
class System : public ReferenceCountedObject<System>
{
protected:

  /**
   * Constructor.  Optionally initializes required
   * data structures.  Protected so that this base class
   * cannot be explicitly instantiated.
   */
  System (EquationSystems& es,
	  const std::string& name,
	  const unsigned int number);
  
public:

  
  /**
   * Destructor.
   */
  virtual ~System ();

  /**
   * The type of system.
   */
  typedef System sys_type;

  /**
   * @returns a clever pointer to the system.
   */
  sys_type & system () { return *this; }
  
  /**
   * Clear all the data structures associated with
   * the system.
   */
  virtual void clear ();
  
  /** 
   * Initializes degrees of freedom on the current mesh.
   * Sets the 
   */
  void init ();
  
  /**
   * Reinitializes degrees of freedom and other 
   * required data on the current mesh.  Note that the matrix
   * is not initialized at this time since it may not be required
   * for all applications. @e Should be overloaded in derived classes.
   */
  virtual void reinit ();
  
  /**
   * Update the local values to reflect the solution
   * on neighboring processors.
   */
  virtual void update ();
  
  /**
   * Prepares \p matrix and \p _dof_map for matrix assembly.
   * Does not actually assemble anything.  For matrix assembly,
   * use the \p assemble() in derived classes.
   * @e Should be overloaded in derived classes.
   */
  virtual void assemble ();

  /**
   * Solves the system.  Must be overloaded in derived systems.
   */
  virtual void solve () = 0;
  
  /**
   * @returns \p true when the other system contains
   * identical data, up to the given threshold.  Outputs
   * some diagnostic info when \p verbose is set.
   */
  virtual bool compare (const System& other_system, 
			const Real threshold,
			const bool verbose) const;

  /**
   * @returns the system name.
   */
  const std::string & name () const;

  /**
   * @returns the type of system, helpful in identifying
   * which system type to use when reading equation system
   * data from file.  Should be overloaded in derived classes.
   */
  virtual std::string system_type () const { return "BasicSystem"; }

  /**
   * Projects the continuous functions onto the current solution. 
   */
  void project_solution (Number fptr(const Point& p,
				     const Parameters& parameters,
                                     const std::string& sys_name,
				     const std::string& unknown_name),
                         Gradient gptr(const Point& p,
				       const Parameters& parameters,
                                       const std::string& sys_name,
				       const std::string& unknown_name),
			 Parameters& parameters) const;

  /**
   * Projects the continuous functions onto the current mesh.
   */
  void project_vector (Number fptr(const Point& p,
				   const Parameters& parameters,
                                   const std::string& sys_name,
				   const std::string& unknown_name),
                       Gradient gptr(const Point& p,
				     const Parameters& parameters,
                                     const std::string& sys_name,
				     const std::string& unknown_name),
		       Parameters& parameters,
		       NumericVector<Number>& new_vector) const;
  
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
  const MeshBase & get_mesh() const;
  
  /**
   * @returns a reference to this systems's \p _mesh.
   */
  MeshBase & get_mesh();
  
  /**
   * @returns a constant reference to this system's \p _dof_map.
   */
  const DofMap & get_dof_map() const;
  
  /**
   * @returns a writeable reference to this system's \p _dof_map.
   */
  DofMap & get_dof_map();

  /**
   * @returns a constant reference to this system's parent EquationSystems object.
   */
  const EquationSystems & get_equation_systems() const { return _equation_systems; }

  /**
   * @returns a reference to this system's parent EquationSystems object.
   */
  EquationSystems & get_equation_systems() { return _equation_systems; }

  /**
   * @returns \p true if the system is active, \p false otherwise.
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
   * Vector iterator typedefs.
   */
  typedef std::map<std::string, NumericVector<Number>* >::iterator       vectors_iterator;
  typedef std::map<std::string, NumericVector<Number>* >::const_iterator const_vectors_iterator;

  /**
   * Beginning of vectors container
   */
  vectors_iterator vectors_begin ();

  /**
   * Beginning of vectors container
   */
  const_vectors_iterator vectors_begin () const;

  /**
   * End of vectors container
   */
  vectors_iterator vectors_end ();

  /**
   * End of vectors container
   */
  const_vectors_iterator vectors_end () const;

  /**
   * Adds the additional vector \p vec_name to this system.  Only
   * allowed @e prior to \p init().  All the additional vectors
   * are similarly distributed, like the \p solution,
   * and inititialized to zero.
   * 
   * By default vectors added by add_vector are projected to changed grids by
   * reinit().  To zero them instead (more efficient), pass "false" as the
   * second argument
   */
  NumericVector<Number> & add_vector (const std::string& vec_name,
				      const bool projections=true);

  /**
   * Tells the System whether or not to project the solution vector onto new
   * grids when the system is reinitialized.  The solution will be projected
   * unless project_solution_on_reinit() = false is called.
   */
  bool& project_solution_on_reinit (void)
    { return _solution_projection; }
  
  /**
   * @returns \p true if this \p System has a vector associated with the
   * given name, \p false otherwise.
   */
  bool have_vector (const std::string& vec_name) const;
  
  /**
   * @returns a const reference to this system's @e additional vector
   * named \p vec_name.  Access is only granted when the vector is already
   * properly initialized.
   */
  const NumericVector<Number> & get_vector (const std::string& vec_name) const;

  /**
   * @returns a writeable reference to this system's @e additional vector
   * named \p vec_name.  Access is only granted when the vector is already
   * properly initialized.
   */
  NumericVector<Number> & get_vector (const std::string& vec_name);

  /**
   * @returns the number of vectors (in addition to the solution)
   * handled by this system
   * This is the size of the \p _vectors map
   */
  unsigned int n_vectors () const;

  /**
   * @returns the number of variables in the system
   */
  unsigned int n_vars() const;

  /**
   * @returns the number of degrees of freedom in the system
   */
  unsigned int n_dofs() const;

  /**
   * Returns the number of active degrees of freedom
   * for this System.
   */
  unsigned int n_active_dofs() const;

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
   * for this system.  Returns the index number for the new variable.
   */
  unsigned int add_variable (const std::string& var,
		             const FEType& type);

  /**
   * Adds the variable \p var to the list of variables
   * for this system.  Same as before, but assumes \p LAGRANGE
   * as default value for \p FEType.family.
   */
  unsigned int add_variable (const std::string& var,
		             const Order order = FIRST,
		             const FEFamily = LAGRANGE);

  /**
   * @returns true if a variable named \p var exists in this System
   */
  bool has_variable(const std::string& var) const;
  
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
   * @returns a norm of variable \p var in the vector \p v, in the specified
   * norm (e.g. L2, H0, H1)
   */
  Real calculate_norm(NumericVector<Number>& v,
		      unsigned int var = 0,
		      FEMNormType norm_type = L2) const;

  /**
   * @returns a norm of the vector \p v, using \p component_norm and \p
   * component_scale to choose and weight the norms of each variable.
   */
  Real calculate_norm(NumericVector<Number>& v,
		      const SystemNorm &norm) const;

  /**
   * Reads the basic data header for this System.
   */
  void read_header (Xdr& io, 
		    const bool read_header=true,
		    const bool read_additional_data=true,
		    const bool read_legacy_format=false);

  /**
   * Reads additional data, namely vectors, for this System.
   */
  void read_legacy_data (Xdr& io,
			 const bool read_additional_data=true);

  /**
   * Reads additional data, namely vectors, for this System.
   * This method may safely be called on a distributed-memory mesh.
   */
  void read_serialized_data (Xdr& io,
			     const bool read_additional_data=true);

  /**
   * Reads additional data, namely vectors, for this System.
   * This method may safely be called on a distributed-memory mesh.
   * This method will read an individual file for each processor in the simulation
   * where the local solution components for that processor are stored.
   */
  void read_parallel_data (Xdr &io,
			   const bool read_additional_data);
  /**
   * Writes the basic data header for this System.
   */
  void write_header (Xdr& io,
		     const bool write_additional_data) const;

  /**
   * Writes additional data, namely vectors, for this System.
   * This method may safely be called on a distributed-memory mesh.
   */
  void write_serialized_data (Xdr& io,
			      const bool write_additional_data = true) const;

  /**
   * Writes additional data, namely vectors, for this System.
   * This method may safely be called on a distributed-memory mesh.
   * This method will create an individual file for each processor in the simulation
   * where the local solution components for that processor will be stored.
   */
  void write_parallel_data (Xdr &io,
			    const bool write_additional_data) const;
  
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
   * Register a user function for imposing constraints.
   */
  void attach_constraint_function (void fptr(EquationSystems& es,
					     const std::string& name));
  
  /**
   * User-provided initialization function.  Can be overloaded
   * in derived classes.
   */
  virtual void user_initialization ();
  
  /**
   * User-provided  assembly function.  Can be overloaded
   * in derived classes.
   */
  virtual void user_assembly ();
  
  /**
   * User provided constraint function.
   */
  virtual void user_constrain ();

  /**
   * Re-update the local values when the mesh has changed.
   * This method takes the data updated by \p update() and
   * makes it up-to-date on the current mesh.
   */
  virtual void re_update ();

  /**
   * Restrict vectors after the mesh has coarsened
   */
  virtual void restrict_vectors ();

  /**
   * Prolong vectors after the mesh has refined
   */
  virtual void prolong_vectors ();

  /**
   * Flag which tells the system to whether or not to
   * call the user assembly function during each call to solve().
   * By default, every call to solve() begins with a call to the
   * user assemble, so this flag is true.  (For explicit systems,
   * "solving" the system occurs during the assembly step, so this
   * flag is always true for explicit systems.)  
   * 
   * You will only want to set this to false if you need direct
   * control over when the system is assembled, and are willing to
   * track the state of its assembly yourself.  An example of such a
   * case is an implicit system with multiple right hand sides.  In
   * this instance, a single assembly would likely be followed with
   * multiple calls to solve.
   *
   * The frequency system and Newmark system have their own versions
   * of this flag, called _finished_assemble, which might be able to
   * be replaced with this more general concept.
   */
  bool assemble_before_solve;


  //--------------------------------------------------
  // The solution and solution access members
  
  /**
   * @returns the current solution for the specified global
   * DOF.
   */
  Number current_solution (const unsigned int global_dof_number) const;

  /**
   * Data structure to hold solution values.
   */
  AutoPtr<NumericVector<Number> > solution;

  /**
   * All the values I need to compute my contribution
   * to the simulation at hand.  Think of this as the
   * current solution with any ghost values needed from
   * other processors.  This vector is necessarily larger
   * than the \p solution vector in the case of a parallel
   * simulation.  The \p update() member is used to synchronize
   * the contents of the \p solution and \p current_local_solution
   * vectors.
   */
  AutoPtr<NumericVector<Number> > current_local_solution;


protected:

  /**
   * Initializes the data for the system.  Note that this is called
   * before any user-supplied intitialization function so that all
   * required storage will be available.
   */
  virtual void init_data ();

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

  
private:

  

  /**
   * Reads an input vector from the stream \p io and assigns
   * the values to a set of \p DofObjects.  This method uses
   * blocked input and is safe to call on a distributed memory-mesh.
   */   
  template <typename iterator_type>
  unsigned int read_serialized_blocked_dof_objects (const unsigned int var,
						    const unsigned int n_objects,
						    const iterator_type begin,
						    const iterator_type end,
						    Xdr &io,
						    NumericVector<Number> &vec) const;

  /**
   * Reads a vector for this System.
   * This method may safely be called on a distributed-memory mesh.
   */
  void read_serialized_vector (Xdr& io,
			       NumericVector<Number> &vec);

  /**
   * Writes an output vector to the stream \p io for a set of \p DofObjects.
   * This method uses blocked output and is safe to call on a distributed memory-mesh.
   */   
  template <typename iterator_type>
  unsigned int write_serialized_blocked_dof_objects (const NumericVector<Number> &vec,
						     const unsigned int var,
						     const unsigned int n_objects,
						     const iterator_type begin,
						     const iterator_type end,
						     Xdr &io) const;

  /**
   * Writes a vector for this System.
   * This method may safely be called on a distributed-memory mesh.
   */
  void write_serialized_vector (Xdr& io,
				const NumericVector<Number> &vec) const;

  /**
   * Function that initializes the system.
   */
  void (* _init_system) (EquationSystems& es,
			 const std::string& name);
  
  /**
   * Function that assembles the system.
   */
  void (* _assemble_system) (EquationSystems& es,
			     const std::string& name);

  /**
   * Function to impose constraints.
   */
  void (* _constrain_system) (EquationSystems& es, 
			      const std::string& name
			      );

  /**
   * Data structure describing the relationship between
   * nodes, variables, etc... and degrees of freedom.
   */
  AutoPtr<DofMap> _dof_map;

  /**
   * Constant reference to the \p EquationSystems object
   * used for the simulation.
   */
  EquationSystems& _equation_systems;
  
  /**
   * Constant reference to the \p mesh data structure used
   * for the simulation.   
   */
  MeshBase& _mesh;
  
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
   * Some systems need an arbitrary number of vectors.
   * This map allows names to be associated with arbitrary
   * vectors.  All the vectors in this map will be distributed
   * in the same way as the solution vector.
   */
  std::map<std::string, NumericVector<Number>* > _vectors;

  /**
   * Holds true if a vector by that name should be projected
   * onto a changed grid, false if it should be zeroed.
   */
  std::map<std::string, bool> _vector_projections;

  /**
   * Holds true if the solution vector should be projected
   * onto a changed grid, false if it should be zeroed.
   * This is true by default.
   */
  bool _solution_projection;

  /**
   * \p true when additional vectors may still be added, \p false otherwise.
   */
  bool _can_add_vectors;

  /**
   * This flag is used only when *reading* in a system from file.
   * Based on the system header, it keeps track of whether or not
   * additional vectors were actually written for this file.
   */
  bool _additional_data_written;


  /**
   * This class implements projecting a vector from 
   * an old mesh to the newly refined mesh.  This
   * may be executed in parallel on multiple threads.
   */
  class ProjectVector 
  {
  private:
    const System                &system;
    const NumericVector<Number> &old_vector;
    NumericVector<Number>       &new_vector;

  public:
    ProjectVector (const System &system_in,
		   const NumericVector<Number> &old_v_in,
		   NumericVector<Number> &new_v_in) :
    system(system_in),
    old_vector(old_v_in),
    new_vector(new_v_in)
    {}
    
    void operator()(const ConstElemRange &range) const;
  };


  /**
   * This class builds the send_list of old dof indices
   * whose coefficients are needed to perform a projection.
   * This may be executed in parallel on multiple threads.
   * The end result is a \p send_list vector which is 
   * unsorted and may contain duplicate elements.
   * The \p unique() method can be used to sort and
   * create a unique list.
   */
  class BuildProjectionList
  {
  private:
    const System              &system;

  public:
    BuildProjectionList (const System &system_in) :
    system(system_in)
    {}

    BuildProjectionList (BuildProjectionList &other, Threads::split) :
     system(other.system)
    {}
    
    void unique();
    void operator()(const ConstElemRange &range);
    void join (const BuildProjectionList &other);
    std::vector<unsigned int> send_list;
  };



  /**
   * This class implements projecting a vector from 
   * an old mesh to the newly refined mesh.  This
   * may be exectued in parallel on multiple threads.
   */
  class ProjectSolution
  {
  private:
    const System                &system;

    Number (* fptr)(const Point& p,
		    const Parameters& parameters,
		    const std::string& sys_name,
		    const std::string& unknown_name);

    Gradient (* gptr)(const Point& p,
		      const Parameters& parameters,
		      const std::string& sys_name,
		      const std::string& unknown_name);

    Parameters &parameters;
    NumericVector<Number>       &new_vector;

  public:
    ProjectSolution (const System &system_in,
		     Number fptr_in(const Point& p,
				    const Parameters& parameters,
				    const std::string& sys_name,
				    const std::string& unknown_name),
		     Gradient gptr_in(const Point& p,
				      const Parameters& parameters,
				      const std::string& sys_name,
				      const std::string& unknown_name),
		     Parameters &parameters_in,
		     NumericVector<Number> &new_v_in) :
    system(system_in),
    fptr(fptr_in),
    gptr(gptr_in),
    parameters(parameters_in),
    new_vector(new_v_in)
    {}
    
    void operator()(const ConstElemRange &range) const;
  };
};



// ------------------------------------------------------------
// System inline methods
inline
const std::string & System::name() const
{
  return _sys_name;
}



inline
unsigned int System::number() const
{
  return _sys_number;
}



inline
const MeshBase & System::get_mesh() const
{
  return _mesh;
}



inline
MeshBase & System::get_mesh()
{
  return _mesh;
}



inline
const DofMap & System::get_dof_map() const
{
  return *_dof_map;
}



inline
DofMap & System::get_dof_map()
{
  return *_dof_map;
}



inline
bool System::active() const
{
  return _active;  
}



inline
void System::activate ()
{
  _active = true;
}



inline
void System::deactivate ()
{
  _active = false;
}



inline
unsigned int System::n_vars() const
{
  return _var_names.size();
}



inline
const std::string & System::variable_name (const unsigned int i) const
{
  libmesh_assert (i < n_vars());

  return _var_names[i];
}



inline
const FEType & System::variable_type (const unsigned int i) const
{
  return variable_type(_var_names[i]);
}



inline
const FEType & System::variable_type (const std::string& var) const
{
  std::map<std::string, FEType>::const_iterator
    pos = _var_type.find(var);

  libmesh_assert (pos != _var_type.end());
  
  return pos->second;
}



inline
unsigned int System::n_active_dofs() const
{
  return this->n_dofs() - this->n_constrained_dofs();
}




inline
bool System::have_vector (const std::string& vec_name) const
{
  return (_vectors.count(vec_name));
}



inline
unsigned int System::n_vectors () const
{
  return _vectors.size(); 
}

inline
System::vectors_iterator System::vectors_begin ()
{
  return _vectors.begin();
}

inline
System::const_vectors_iterator System::vectors_begin () const
{
  return _vectors.begin();
}

inline
System::vectors_iterator System::vectors_end ()
{
  return _vectors.end();
}

inline
System::const_vectors_iterator System::vectors_end () const
{
  return _vectors.end();
}



#endif
