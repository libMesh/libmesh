// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_EQUATION_SYSTEMS_H
#define LIBMESH_EQUATION_SYSTEMS_H

// Local Includes
#include "libmesh/libmesh_common.h"
#include "libmesh/parameters.h"
#include "libmesh/system.h"
#include "libmesh/parallel_object.h"

#ifdef LIBMESH_FORWARD_DECLARE_ENUMS
namespace libMesh
{
enum XdrMODE : int;
}
#else
#include "libmesh/enum_xdr_mode.h"
#endif

// HP aCC needs these for some reason
#ifdef __HP_aCC
# include "libmesh/frequency_system.h"
# include "libmesh/transient_system.h"
# include "libmesh/newmark_system.h"
# include "libmesh/steady_system.h"
#endif

// C++ includes
#include <cstddef>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace libMesh
{

// Forward Declarations
template <typename> class ElemTempl;
typedef ElemTempl<Real> Elem;
template <typename> class MeshBaseTempl;
typedef MeshBaseTempl<Real> MeshBase;

/**
 * This is the \p EquationSystems class.  It is in charge
 * of handling all the various equation systems defined
 * for a \p MeshBase.  It may have multiple systems, which may
 * be active or inactive, so that at different solution
 * stages only a sub-set may be solved for.  Also, through
 * the templated access, different types of systems
 * may be handled.  Also other features, like flags,
 * parameters, I/O etc are provided.
 *
 * \author Benjamin S. Kirk
 * \date 2002-2007
 * \brief Manages multiples systems of equations.
 */
class EquationSystems : public ReferenceCountedObject<EquationSystems>,
                        public ParallelObject

{
public:

  /**
   * Define enumeration to set properties in EquationSystems::read()
   */
  enum ReadFlags { READ_HEADER           = 1,
                   READ_DATA             = 2,
                   READ_ADDITIONAL_DATA  = 4,
                   READ_LEGACY_FORMAT    = 8,
                   TRY_READ_IFEMS        = 16,
                   READ_BASIC_ONLY       = 32 };

  /**
   * Define enumeration to set properties in EquationSystems::write()
   */
  enum WriteFlags { WRITE_DATA             = 1,
                    WRITE_ADDITIONAL_DATA  = 2,
                    WRITE_PARALLEL_FILES   = 4,
                    WRITE_SERIAL_FILES     = 8 };

  /**
   * Constructor.
   */
  EquationSystems (MeshBase & mesh);

  /**
   * Destructor.  Should be virtual, since the user may want to derive
   * subclasses of EquationSystems.
   */
  virtual ~EquationSystems ();

  /**
   * Restores the data structure to a pristine state.
   */
  virtual void clear ();

  /**
   * Initialize all the systems
   */
  virtual void init ();

  /**
   * Handle any mesh changes and reinitialize all the systems on the
   * updated mesh
   */
  virtual void reinit ();

  /**
   * Enable or disable default ghosting functors on the Mesh and on
   * all Systems.  Standard ghosting is enabled by default.  If
   * disabled, default ghosting will also be disabled on any later
   * added systems.
   *
   * Unless other equivalent ghosting functors have been added,
   * removing the default coupling functor is only safe for explicit
   * solves, and removing the default algebraic ghosting functor is
   * only safe for codes where no evaluations on neighbor cells (e.g.
   * no jump error estimators) are done.
   */
  virtual void enable_default_ghosting (bool enable);

  /**
   * Updates local values for all the systems
   */
  void update ();

  /**
   * \returns The number of equation systems.
   */
  unsigned int n_systems() const;

  /**
   * \returns \p true if the system named \p name exists within
   * this EquationSystems object.
   */
  bool has_system (const std::string & name) const;

  /**
   * \returns A constant reference to the system named \p name.
   * The template argument defines the return type.  For example,
   * const SteadySystem & sys = eq.get_system<SteadySystem> ("sys");
   * is an example of how the method might be used
   */
  template <typename T_sys>
  const T_sys & get_system (const std::string & name) const;

  /**
   * \returns A writable reference to the system named \p name.
   * The template argument defines the return type.  For example,
   * const SteadySystem & sys = eq.get_system<SteadySystem> ("sys");
   * is an example of how the method might be used
   */
  template <typename T_sys>
  T_sys & get_system (const std::string & name);

  /**
   * \returns A constant reference to system number \p num.
   * The template argument defines the return type.  For example,
   * const SteadySystem & sys = eq.get_system<SteadySystem> (0);
   * is an example of how the method might be used
   */
  template <typename T_sys>
  const T_sys & get_system (const unsigned int num) const;

  /**
   * \returns A writable reference to the system number \p num.
   * The template argument defines the return type.  For example,
   * const SteadySystem & sys = eq.get_system<SteadySystem> (0);
   * is an example of how the method might be used
   */
  template <typename T_sys>
  T_sys & get_system (const unsigned int num);

  /**
   * \returns A constant reference to the system named \p name.
   */
  const System & get_system (const std::string & name) const;

  /**
   * \returns A writable reference to the system named \p name.
   */
  System & get_system (const std::string & name);

  /**
   * \returns A constant reference to system number \p num.
   */
  const System & get_system (const unsigned int num) const;

  /**
   * \returns A writable reference to the system number \p num.
   */
  System & get_system (const unsigned int num);

  /**
   * Add the system of type \p system_type named \p name to the
   * systems array.
   */
  virtual System & add_system (const std::string & system_type,
                               const std::string & name);

  /**
   * Add the system named \p name to the systems array.
   */
  template <typename T_sys>
  T_sys & add_system (const std::string & name);

  /**
   * Remove the system named \p name from the systems array.
   *
   * \deprecated This function may not work as intended and has not
   * been actively tested over the years. If you need the ability to
   * delete a System from an EquationSystems object, it could probably
   * be added.
   */
#ifdef LIBMESH_ENABLE_DEPRECATED
  void delete_system (const std::string & name);
#endif

  /**
   * \returns The total number of variables in all
   * systems.
   */
  unsigned int n_vars () const;

  /**
   * \returns The total number of degrees of freedom
   * in all systems.
   */
  std::size_t n_dofs () const;

  /**
   * \returns The number of active degrees of freedom
   * for the EquationSystems object.
   */
  std::size_t n_active_dofs() const;

  /**
   * Call \p solve on all the individual equation systems.
   *
   * By default this function solves each equation system once,
   * in the order they were added.  For more sophisticated decoupled
   * problems the user may with to override this behavior in a derived
   * class.
   */
  virtual void solve ();

  /**
   * Call \p adjoint_solve on all the individual equation systems.
   *
   * By default this function solves each system's adjoint once,
   * in the reverse order from that in which they were added.  For
   * more sophisticated decoupled problems the user may with to
   * override this behavior in a derived class.
   */
  virtual void adjoint_solve (const QoISet & qoi_indices = QoISet());

  /**
   * Call \p sensitivity_solve on all the individual equation systems.
   *
   * By default this function solves each sensitivity system once,
   * in the order in which in which they were added.  For
   * more sophisticated decoupled problems the user may with to
   * override this behavior in a derived class.
   */
  virtual void sensitivity_solve (const ParameterVector & parameters);

  /**
   * Fill the input vector \p var_names with the names
   * of the variables for each system. If \p type is passed,
   * only variables of the specified type will be populated.
   * If systems_names!=nullptr, only include names from the
   * specified systems.
   */
  void build_variable_names (std::vector<std::string> & var_names,
                             const FEType * type=nullptr,
                             const std::set<std::string> * system_names=nullptr) const;

  /**
   * Fill the input vector \p soln with the solution values for the
   * system named \p name.
   *
   * \note The input vector \p soln will only be assembled on
   * processor 0, so this method is only applicable to outputting plot
   * files from processor 0.
   */
  void build_solution_vector (std::vector<Number> & soln,
                              const std::string & system_name,
                              const std::string & variable_name = "all_vars") const;

  /**
   * Fill the input vector \p soln with solution values.  The
   * entries will be in variable-major format (corresponding to
   * the names from \p build_variable_names()).
   * If systems_names!=nullptr, only include data from the
   * specified systems.
   */
  void build_solution_vector (std::vector<Number> & soln,
                              const std::set<std::string> * system_names=nullptr) const;

  /**
   * A version of build_solution_vector which is appropriate for
   * "parallel" output formats like Nemesis.
   *
   * \returns A std::unique_ptr to a node-major NumericVector of total
   * length n_nodes*n_vars that various I/O classes can then use to
   * get the local values they need to write on each processor.
   */
  std::unique_ptr<NumericVector<Number>>
  build_parallel_solution_vector(const std::set<std::string> * system_names=nullptr) const;

  /**
   * Retrieve \p vars_active_subdomains, which indicates the active
   * subdomains for each variable in \p names.
   */
  void get_vars_active_subdomains(const std::vector<std::string> & names,
                                  std::vector<std::set<subdomain_id_type>> & vars_active_subdomains) const;

  /**
   * Retrieve the solution data for CONSTANT MONOMIALs.  If \p names
   * is populated, only the variables corresponding to those names will
   * be retrieved.  This can be used to filter which variables are retrieved.
   *
   * \deprecated Call the more appropriately-named build_elemental_solution_vector()
   * instead.
   */
  void get_solution (std::vector<Number> & soln,
                     std::vector<std::string> & names) const;

  /**
   * Retrieve the solution data for CONSTANT MONOMIALs.  If \p names
   * is populated, only the variables corresponding to those names will
   * be retrieved.  This can be used to filter which variables are retrieved.
   *
   * This is the more appropriately-named replacement for the get_solution()
   * function defined above.
   */
  void build_elemental_solution_vector (std::vector<Number> & soln,
                                        std::vector<std::string> & names) const;

  /**
   * Finds system and variable numbers for any variables of \p type
   * corresponding to the entries in the input 'names' vector.
   */
  std::vector<std::pair<unsigned int, unsigned int>>
  find_variable_numbers (std::vector<std::string> & names,
                         const FEType * type=nullptr) const;

  /**
   * Builds a parallel vector of CONSTANT MONOMIAL solution values
   * corresponding to the entries in the input 'names' vector.  This
   * vector is approximately uniformly distributed across all of the
   * available processors.
   *
   * The related function build_elemental_solution_vector() is
   * implemented by calling this function and then calling
   * localize_to_one() on the resulting vector.
   *
   * \returns A nullptr (if no CONSTANT, MONOMIAL variables exist on
   * the system) or a std::unique_ptr to a var-major numeric vector of
   * total length n_elem * n_vars ordered according to:
   * [u0, u1, ... uN, v0, v1, ... vN, w0, w1, ... wN]
   * for constant monomial variables (u, v, w) on a mesh with N elements.
   */
  std::unique_ptr<NumericVector<Number>>
  build_parallel_elemental_solution_vector (std::vector<std::string> & names) const;

  /**
   * Fill the input vector \p soln with solution values.  The
   * entries will be in variable-major format (corresponding to
   * the names from \p build_variable_names()).
   * If systems_names!=nullptr, only include data from the
   * specified systems.
   * If vertices_only == true, then for higher-order elements only
   * the solution at the vertices is computed. This can be useful
   * when plotting discontinous solutions on higher-order elements
   * and only a lower-order representation is required.
   */
  void build_discontinuous_solution_vector
  (std::vector<Number> & soln,
   const std::set<std::string> * system_names = nullptr,
   const std::vector<std::string> * var_names = nullptr,
   bool vertices_only = false) const;

  /**
   * Read & initialize the systems from disk using the XDR data format.
   * This format allows for machine-independent binary output.
   *
   * Set which sections of the file to read by bitwise OR'ing the
   * EquationSystems::ReadFlags enumeration together. For example, to
   * read all sections of the file, set read_flags to:
   * (READ_HEADER | READ_DATA | READ_ADDITIONAL_DATA)
   *
   * \note The equation system can be defined without initializing
   * the data vectors to any solution values.  This can be done
   * by omitting READ_DATA in the read_flags parameter.
   *
   * If XdrMODE is omitted, it will be inferred as READ for filenames
   * containing .xda or as DECODE for filenames containing .xdr
   *
   * \param name Name of the file to be read.
   * \param read_flags Single flag created by bitwise-OR'ing several flags together.
   * \param mode Controls whether reading is done in binary or ascii mode.
   * \param partition_agnostic If true then the mesh and degrees of freedom
   * will be temporarily renumbered in a partition agnostic way so that
   * files written using "n" mpi processes can be re-read on "m" mpi
   * processes.  This renumbering is not compatible with meshes
   * that have two nodes in exactly the same position!
   */
  template <typename InValType>
  void read (const std::string & name,
             const XdrMODE,
             const unsigned int read_flags=(READ_HEADER | READ_DATA),
             bool partition_agnostic = true);

  void read (const std::string & name,
             const XdrMODE mode,
             const unsigned int read_flags=(READ_HEADER | READ_DATA),
             bool partition_agnostic = true)
  { read<Number>(name, mode, read_flags, partition_agnostic); }

  template <typename InValType>
  void read (const std::string & name,
             const unsigned int read_flags=(READ_HEADER | READ_DATA),
             bool partition_agnostic = true);

  void read (const std::string & name,
             const unsigned int read_flags=(READ_HEADER | READ_DATA),
             bool partition_agnostic = true)
  { read<Number>(name, read_flags, partition_agnostic); }


  /**
   * Write the systems to disk using the XDR data format.
   * This format allows for machine-independent binary output.
   *
   * Set the writing properties using the EquationSystems::WriteFlags
   * enumeration. Set which sections to write out by bitwise OR'ing
   * the enumeration values. Write everything by setting write_flags to:
   * (WRITE_DATA | WRITE_ADDITIONAL_DATA)
   *
   * \note The solution data can be omitted by calling
   * this routine with WRITE_DATA omitted in the write_flags argument.
   *
   * If XdrMODE is omitted, it will be inferred as WRITE for filenames
   * containing .xda or as ENCODE for filenames containing .xdr
   *
   * \param name Name of the file to be read.
   * \param write_flags Single flag created by bitwise-OR'ing several flags together.
   * \param mode Controls whether reading is done in binary or ascii mode.
   * \param partition_agnostic If true then the mesh and degrees of freedom
   * will be temporarily renumbered in a partition agnostic way so that
   * files written using "n" mpi processes can be re-read on "m" mpi
   * processes.  This renumbering is not compatible with meshes
   * that have two nodes in exactly the same position!
   */
  void write (const std::string & name,
              const XdrMODE,
              const unsigned int write_flags=(WRITE_DATA),
              bool partition_agnostic = true) const;

  void write (const std::string & name,
              const unsigned int write_flags=(WRITE_DATA),
              bool partition_agnostic = true) const;

  /**
   * \returns \p true when this equation system contains
   * identical data, up to the given threshold.  Delegates
   * most of the comparisons to perform to the responsible
   * systems
   */
  virtual bool compare (const EquationSystems & other_es,
                        const Real threshold,
                        const bool verbose) const;

  /**
   * \returns A string containing information about the
   * systems, flags, and parameters.
   */
  virtual std::string get_info() const;

  /**
   * Prints information about the equation systems, by default to
   * libMesh::out.
   */
  void print_info (std::ostream & os=libMesh::out) const;

  /**
   * Same as above, but allows you to also use stream syntax.
   */
  friend std::ostream & operator << (std::ostream & os,
                                     const EquationSystems & es);

  /**
   * \returns A constant reference to the mesh
   */
  const MeshBase & get_mesh() const;

  /**
   * \returns A reference to the mesh
   */
  MeshBase & get_mesh();

  /**
   * Serializes a distributed mesh and its associated
   * degree of freedom numbering for all systems
   **/
  void allgather ();

  /**
   * Calls to reinit() will also do two-step coarsen-then-refine
   **/
  void enable_refine_in_reinit() { this->_refine_in_reinit = true; }

  /**
   * Calls to reinit() will not try to coarsen or refine the mesh
   **/
  void disable_refine_in_reinit() { this->_refine_in_reinit = false; }

  /**
   * \returns Whether or not calls to reinit() will try to coarsen/refine the mesh
   **/
  bool refine_in_reinit_flag() { return this->_refine_in_reinit; }

  /**
   * Handle any mesh changes and project any solutions onto the
   * updated mesh.
   *
   * \returns Whether or not the mesh may have changed.
   */
  bool reinit_solutions ();

  /**
   * Reinitialize all systems on the current mesh.
   */
  virtual void reinit_systems ();

  /**
   * Data structure holding arbitrary parameters.
   */
  Parameters parameters;


protected:


  /**
   * The mesh data structure
   */
  MeshBase & _mesh;

  /**
   * Data structure holding the systems.
   */
  std::map<std::string, System *> _systems;

  /**
   * Typedef for system iterators
   */
  typedef std::map<std::string, System *>::iterator       system_iterator;

  /**
   * Typedef for constant system iterators
   */
  typedef std::map<std::string, System *>::const_iterator const_system_iterator;

  /**
   * Flag for whether to call coarsen/refine in reinit().
   * Default value: true
   */
  bool _refine_in_reinit;

  /**
   * Flag for whether to enable default ghosting on newly added Systems.
   * Default value: true
   */
  bool _enable_default_ghosting;

private:

  /**
   * Actual read implementation.  This can be called repeatedly
   * inside a try-catch block in an attempt to read broken files.
   *
   * \param name Name of the file to be read.
   * \param read_flags Single flag created by bitwise-OR'ing several flags together.
   * \param partition_agnostic If true then the mesh and degrees of freedom
   * will be temporarily renumbered in a partition agnostic way so that
   * files written using "n" mpi processes can be re-read on "m" mpi
   * processes.
   *
   * \note This renumbering is not compatible with meshes that have
   * two nodes in exactly the same position!
   */
  template <typename InValType>
  void _read_impl (const std::string & name,
                   const XdrMODE,
                   const unsigned int read_flags,
                   bool partition_agnostic = true);

  /**
   * This function is used in the implementation of add_system,
   * it loops over the nodes and elements of the Mesh, adding the
   * system to each one.  The main reason to separate this part
   * is to avoid coupling this header file to mesh.h, and elem.h.
   */
  void _add_system_to_nodes_and_elems();

  /**
   * This just calls DofMap::remove_default_ghosting() but using a
   * shim lets us forward-declare DofMap.
   */
  void _remove_default_ghosting(unsigned int sys_num);
};



// ------------------------------------------------------------
// EquationSystems inline methods
inline
const MeshBase & EquationSystems::get_mesh () const
{
  return _mesh;
}



inline
MeshBase & EquationSystems::get_mesh ()
{
  return _mesh;
}


inline
unsigned int EquationSystems::n_systems () const
{
  return cast_int<unsigned int>(_systems.size());
}




template <typename T_sys>
inline
T_sys & EquationSystems::add_system (const std::string & name)
{
  T_sys * ptr = nullptr;

  if (!_systems.count(name))
    {
      const unsigned int sys_num = this->n_systems();
      ptr = new T_sys(*this, name, sys_num);

      _systems.insert (std::make_pair(name, ptr));

      if (!_enable_default_ghosting)
        this->_remove_default_ghosting(sys_num);

      // Tell all the \p DofObject entities to add a system.
      this->_add_system_to_nodes_and_elems();
    }
  else
    {
      // We now allow redundant add_system calls, to make it
      // easier to load data from files for user-derived system
      // subclasses
      ptr = &(this->get_system<T_sys>(name));
    }

  // Return a dynamically casted reference to the newly added System.
  return *ptr;
}



inline
bool EquationSystems::has_system (const std::string & name) const
{
  if (_systems.find(name) == _systems.end())
    return false;
  return true;
}




template <typename T_sys>
inline
const T_sys & EquationSystems::get_system (const unsigned int num) const
{
  libmesh_assert_less (num, this->n_systems());


  const_system_iterator       pos = _systems.begin();
  const const_system_iterator end = _systems.end();

  for (; pos != end; ++pos)
    if (pos->second->number() == num)
      break;

  // Check for errors
  if (pos == end)
    libmesh_error_msg("ERROR: no system number " << num << " found!");

  // Attempt dynamic cast
  return *cast_ptr<T_sys *>(pos->second);
}




template <typename T_sys>
inline
T_sys & EquationSystems::get_system (const unsigned int num)
{
  libmesh_assert_less (num, this->n_systems());

  const_system_iterator       pos = _systems.begin();
  const const_system_iterator end = _systems.end();

  for (; pos != end; ++pos)
    if (pos->second->number() == num)
      break;

  // Check for errors
  if (pos == end)
    libmesh_error_msg("ERROR: no system number " << num << " found!");

  // Attempt dynamic cast
  return *cast_ptr<T_sys *>(pos->second);
}






template <typename T_sys>
inline
const T_sys & EquationSystems::get_system (const std::string & name) const
{
  const_system_iterator pos = _systems.find(name);

  // Check for errors
  if (pos == _systems.end())
    libmesh_error_msg("ERROR: no system named \"" << name << "\" found!");

  // Attempt dynamic cast
  return *cast_ptr<T_sys *>(pos->second);
}






template <typename T_sys>
inline
T_sys & EquationSystems::get_system (const std::string & name)
{
  system_iterator pos = _systems.find(name);

  // Check for errors
  if (pos == _systems.end())
    libmesh_error_msg("ERROR: no system named " << name << " found!");

  // Attempt dynamic cast
  return *cast_ptr<T_sys *>(pos->second);
}







inline
const System & EquationSystems::get_system (const std::string & name) const
{
  return this->get_system<System>(name);
}



inline
System & EquationSystems::get_system (const std::string & name)
{
  return this->get_system<System>(name);
}



inline
const System & EquationSystems::get_system (const unsigned int num) const
{
  return this->get_system<System>(num);
}



inline
System & EquationSystems::get_system (const unsigned int num)
{
  return this->get_system<System>(num);
}


} // namespace libMesh


#endif // LIBMESH_EQUATION_SYSTEMS_H
