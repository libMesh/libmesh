// $Id: dof_map.h,v 1.1.1.1 2003-01-10 16:17:48 libmesh Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __dof_map_h__
#define __dof_map_h__

// C++ Includes   -----------------------------------
#include <vector>
#include <map>


class MeshBase;

// Local Includes -----------------------------------
#include "mesh_config.h"
#include "mesh_common.h"
#include "elem.h"
#include "order.h"
#include "fe_type.h"
#include "perf_log.h"
#include "dense_matrix.h"


/**
 * This class handles the numbering of degrees of freedom on a mesh.
 * For systems of equations the class supports a fixed number of components.
 * The degrees of freedom are numbered such that sequential, contiguous blocks
 * correspond to distinct subdomains.  This is so that the resulting data
 * structures will work well with parallel linear algebra packages.
 *
 * @author Benjamin S. Kirk, 2002
 */

// ------------------------------------------------------------
// Dof Map class definition

class DofMap
{
 public:

  /**
   * Constructor.  Requires a reference to the mesh.
   */
  DofMap(const MeshBase& mesh);

  /**
   * Destructor
   */
  ~DofMap();



#ifdef ENABLE_AMR

  /**
   * A row of the Dof constraint matrix.  Store as a float
   * to save space since the factors are simple rational
   * fractions.
   */
  typedef std::map<unsigned int, float> DofConstraintRow;

  /** 
   * The constraint matrix storage format. 
   */
  typedef std::map <unsigned int, DofConstraintRow> DofConstraints;

#endif

  /**
   * Distrubute dofs on the current mesh.  Also builds the send list for
   * processor \p proc_id, which defaults to 0 for ease of use in serial
   * applications. 
   */
  void distribute_dofs(const unsigned int proc_id);

  /**
   * Computes the sparsity pattern for the matrix corresponding
   * to \p proc_id.  Produces data that can be fed to Petsc for
   * preallocation of sparse matrices.
   */
  void compute_sparsity (const unsigned int proc_id);

  /**
   * Returns a constant reference to the \p send_list for this processor.  The
   * \p send_list contains the global indices of all the components in the
   * global solution vector that influence the current processor.  This
   * information can be used for gathers at each solution step to retrieve
   * solution values needed for computation.
   */
  const std::vector<unsigned int>& get_send_list() const { return send_list; };
  
  /**
   * Returns a constant reference to the \p n_nz list for this processor.
   * The vector contains the bandwidth of the on-processor coupling for each
   * row of the global matrix that the current processor owns.  This
   * information can be used to preallocate space for a parallel sparse matrix.
   */
  const std::vector<unsigned int>& get_n_nz() const { return n_nz; };
  
  /**
   * Returns a constant reference to the \p n_oz list for this processor.
   * The vector contains the bandwidth of the off-processor coupling for each
   * row of the global matrix that the current processor owns.  This
   * information can be used to preallocate space for a parallel sparse matrix.
   */
  const std::vector<unsigned int>& get_n_oz() const { return n_oz; };

  /**
   * Add an unknown of order \p order and finite element type
   * \p type to the system of equations.
   */
  void add_component (const FEType& type)
  {
    component_types.push_back (type);
    return;
  };
  
  /**
   * @returns the approximation order for component \p c.
   */
  Order component_order (const unsigned int c) const
  { return component_types[c].order; };
  
  /**
   * @returns the finite element type for component \p c.
   */
  FEType component_type (const unsigned int c) const
  { return component_types[c]; };
  
  /**
   * Returns the number of components in the global solution vector. Defaults
   * to 1, should be 1 for a scalar equation, 3 for 2D incompressible Navier
   * Stokes (u,v,p), etc...
   */  
  unsigned int n_components() const
  {
    return component_types.size();
  };

  /**
   * Return the global degree of freedom number for \p component located
   * at \p node.
   */
  unsigned int node_dof_number(const unsigned int node,
			       const unsigned int component=0,
			       const unsigned int index=0) const;

  /**
   * Return the global degree of freedom number for \p component located
   * at \p elem.
   */
  unsigned int elem_dof_number(const unsigned int elem,
			       const unsigned int component=0,
			       const unsigned int index=0) const;

  /**
   * @returns the total number of degrees of freedom in the problem.
   */
  unsigned int n_dofs() const { return n_dfs; };
  
  /**
   * Returns the number of degrees of freedom on subdomain \p proc.
   */
  unsigned int n_dofs_on_processor(const unsigned int proc) const
  { assert(proc < first_df.size()); return (last_df[proc] - first_df[proc]+1); };

  /**
   * Returns the first dof index that is local to subdomain \p proc.
   */
  unsigned int first_dof(unsigned int proc) const
  { assert(proc < first_df.size()); return first_df[proc]; };
  
  /**
   * Returns the last dof index that is local to subdomain \p proc.
   */
  unsigned int last_dof(unsigned int proc) const
  { assert(proc < last_df.size()); return last_df[proc]; };

  /**
   * @returns the total number of degrees of freedom at node n.
   * If no component is specified then the sum over all components
   * is returned.
   */
  unsigned int n_dofs_at_node(const unsigned int n,
  			      const unsigned int c=invalid_number) const;
  
  /**
   * @returns the total number of degrees of freedom internl to element e.
   * If no component is specified then the sum over all components
   * is returned.
   */
  unsigned int n_dofs_on_elem(const unsigned int e,
  			      const unsigned int c=invalid_number) const;
  


#ifdef ENABLE_AMR

  /**
   * @returns the total number of constrained degrees of freedom
   * in the problem.
   */
  unsigned int n_constrained_dofs() const { return dof_constraints.size(); };

#endif

  
  /**
   * Fills the vector di with the global degree of freedom indices
   * for the element. If no component number is specified then all
   * components are returned.
   */
  void dof_indices (const unsigned int e,
		    std::vector<unsigned int>& di,
		    const unsigned int cn=invalid_number) const;



#ifdef ENABLE_AMR

  /**
   * Rebuilds the degree of freedom constraints.
   */ 
  void create_dof_constraints ();

  /**
   * Adds the user-defined row to the constraint matrix.
   * Issues a warning if the DOF was already constrained.
   */
  void add_constraint_row (const unsigned int dof_number,
			   const DofConstraintRow& constraint_row);
  
  /**
   * @returns true if the degree of freedom dof is constrained,
   * false otherwise.
   */
  bool is_constrained_dof (const unsigned int dof) const;
  
  /**
   * Constrains the element matrix.  This method requires the
   * element matrix to be square, in which case the elem_dofs
   * correspond to the global DOF indices of both the rows and
   * columns of the element matrix.  For this case the rows
   * and columns of the matrix necessarily correspond to variables
   * of the same approximation order.
   */
  void constrain_element_matrix (DenseMatrix& matrix,
				 std::vector<unsigned int>& elem_dofs) const;
  
  /**
   * Constrains the element matrix.  This method allows the
   * element matrix to be non-square, in which case the row_dofs
   * and col_dofs may be of different size and correspond to
   * variables approximated in different spaces.
   */
  void constrain_element_matrix (DenseMatrix& matrix,
				 std::vector<unsigned int>& row_dofs,
				 std::vector<unsigned int>& col_dofs) const;
  
  /**
   * Constrains the element vector.
   */
  void constrain_element_vector (std::vector<real>&         rhs,
				 std::vector<unsigned int>& dofs) const;
  
  /**
   * Constrains the element matrix and vector.  This method requires
   * the element matrix to be square, in which case the elem_dofs
   * correspond to the global DOF indices of both the rows and
   * columns of the element matrix.  For this case the rows
   * and columns of the matrix necessarily correspond to variables
   * of the same approximation order.
   */
  void constrain_element_matrix_and_vector (DenseMatrix& matrix,
					    std::vector<real>& rhs,
					    std::vector<unsigned int>& elem_dofs) const;
  
#endif



  /**
   * Reinitialize the underlying data strucures conformal to the current mesh.
   */
  void reinit ();
  
  /**
   * Free all memory associated with the object, but keep the mesh pointer.
   */
  void clear ();
  


#ifdef ENABLE_AMR

  /**
   * Prints the dof_constraints data structure.
   */
  void print_dof_constraints() const;

#endif
  
  /**
   * An invalid DOF index used to initialize data structures.  ONLY refer to
   * \p invalid_entry, and the numerical value may change.
   */
  static const unsigned int invalid_number;

  /**
   * Degree of freedom coupling.  If left empty each DOF
   * couples to all others.  Can be used to reduce memory
   * requirements for sparse matrices.  DOF 0 might only
   * couple to itself, in which case \p dof_coupling(0,0)
   * should be 1 and \p dof_coupling(0,j) = 0 for j not equal
   * to 0.
   */
  CouplingMatrix dof_coupling;

  
 private:


  /**
   * Return the global degree of freedom number for \p component located
   * at \p node as a writeable reference.
   */
  unsigned int & set_node_dof_number(const unsigned int node,
				     const unsigned int component,
				     const unsigned int index=0);
  
  /**
   * Return the global degree of freedom number for \p component located
   * at \p elem as a writeable reference.
   */
  unsigned int & set_elem_dof_number(const unsigned int elem,
				     const unsigned int component,
				     const unsigned int index=0);
  

#ifdef ENABLE_AMR

  /**
   * Build the constraint matrix C associated with the element
   * degree of freedom indices elem_dofs.
   */
  void build_constraint_matrix (DenseMatrix& C,
				std::vector<unsigned int>& elem_dofs) const;

#endif


  /**
   * Finds all the DOFS associated with the element DOFs elem_dofs.
   * This will account for off-element couplings via hanging nodes.
   */
  void find_connected_dofs (std::vector<unsigned int>& elem_dofs) const;
  
#ifdef ENABLE_EXPENSIVE_DATA_STRUCTURES

  /**
   * Data structure defining the nodal degrees of freedom.  This structure
   * allows for multiple DOFS per node for each component.
   */
  std::vector<std::vector<std::vector<unsigned int> > > node_dofs;
  
  /**
   * Data structure defining the element degrees of freedom.  This structure
   * allows for multiple DOFS per element for each component.
   */
  std::vector<std::vector<std::vector<unsigned int> > > elem_dofs;
  
#else
  
  /**
   * Data structure that holds DOF numbers for entities in the mesh
   * that live at the nodes.
   */
  std::vector<unsigned int> node_id_map;
  
  /**
   * Data structure that holds DOF numbers for entities in the mesh
   * that live on the element.
   */
  std::vector<unsigned int> elem_id_map;
  
#endif
  
  /**
   * The finite element type for each component.
   */
  std::vector<FEType> component_types;
  
  /**
   * Constant reference to the computational mesh.
   */
  const MeshBase& mesh;

  /**
   * The dimensionality of the mesh.
   */
  const unsigned int dim;
  
  /**
   * The number of nodes that _were_ in the mesh
   * when this DofMap was initialized.
   */
  unsigned int n_nodes;

  /**
   * The number of elements that _were_ in the mesh
   * when this DofMap was initialized.
   */
  unsigned int n_elem;

  /**
   * Total number of degrees of freedom.
   */
  unsigned int n_dfs;

  /**
   * First DOF index on processor \p p.
   */
  std::vector<unsigned int> first_df;

  /**
   * Last DOF index (plus 1) on processor \p p.
   */
  std::vector<unsigned int> last_df;

  /**
   * A list containing all the global DOF indicies that affect the
   * solution on my subdomain.
   */
  std::vector<unsigned int> send_list;

  /**
   * The number of on-processor nonzeros in my portion of the
   * global matrix.
   */
  std::vector<unsigned int> n_nz;
  
  /**
   * The number of off-processor nonzeros in my portion of the
   * global matrix.
   */
  std::vector<unsigned int> n_oz;


#ifdef ENABLE_AMR

  /**
   * Data structure containing DOF constraints.  The ith
   * entry is the constraint matrix row for DOF i.
   */
  DofConstraints dof_constraints;

#endif

  /**
   * Object to use for performance logging.
   */
  PerfLog perf_log;

  // make the mesh base class (that will contain us) a friend
 friend class MeshBase;
 friend class Mesh;
 friend class BoundaryMesh;
};



// ------------------------------------------------------------
// Dof Map member functions
inline
DofMap::DofMap(const MeshBase& m) :
  mesh(m),
  dim(mesh.mesh_dimension()),
  n_nodes(0),
  n_elem(0),
  n_dfs(0),
#ifdef ENABLE_PERFORMANCE_LOGGING
  perf_log("DofMap", true)
#else
  perf_log("DofMap", false)
#endif
  
{
};




inline
DofMap::~DofMap()
{
  clear();
};



inline
unsigned int DofMap::node_dof_number(const unsigned int node,
				     const unsigned int component,
				     const unsigned int index) const
{
#ifdef ENABLE_EXPENSIVE_DATA_STRUCTURES

  assert (node < n_nodes);
  assert (component < n_components());
  

  if (node_dofs[component].empty())
    {
      return invalid_number;
    }
  else
    {
      assert (index < node_dofs[component][node].size());
      return node_dofs[component][node][index];
    };
  
#else
  
  // If this assert fails then you are trying to use
  // objects with mulitple DOFs per node and the library
  // has not been compiled to support this!
  assert (index == 0);
  assert (node < n_nodes);
  assert (component < n_components());

  return node_id_map[(node) + (component)*(n_nodes)];
  
#endif
};



inline
unsigned int & DofMap::set_node_dof_number(const unsigned int node,
					   const unsigned int component,
					   const unsigned int index)
{
#ifdef ENABLE_EXPENSIVE_DATA_STRUCTURES

  assert (node < n_nodes);
  assert (component < n_components());
  

  if (node_dofs[component].empty())
    {
      error();
      return node_dofs[component][node][index];
    }
  else
    {
      assert (index < node_dofs[component][node].size());
      return node_dofs[component][node][index];
    };
  
#else
  
  // If this assert fails then you are trying to use
  // objects with mulitple DOFs per node and the library
  // has not been compiled to support this!
  assert (index == 0);
  assert (node < n_nodes);
  assert (component < n_components());

  return node_id_map[(node) + (component)*(n_nodes)];

#endif
};



inline
unsigned int DofMap::elem_dof_number(const unsigned int elem,
				     const unsigned int component,
				     const unsigned int index) const
{
#ifdef ENABLE_EXPENSIVE_DATA_STRUCTURES

  assert (elem < n_elem);
  assert (component < n_components());
 
  if (elem_dofs[component].empty())
    {
      return invalid_number;
    }
  else
    {
      assert (index < elem_dofs[component][elem].size()); 
      return elem_dofs[component][elem][index];
    }
  
#else
  
  // If this assert fails then you are trying to use
  // objects with mulitple DOFs per node and the library
  // has not been compiled to support this!
  assert (index == 0);
  assert (elem < n_elem);
  assert (component < n_components());

  return elem_id_map[(elem) + (component)*(n_elem)];
  
#endif
};



inline
unsigned int & DofMap::set_elem_dof_number(const unsigned int elem,
					   const unsigned int component,
					   const unsigned int index)
{
#ifdef ENABLE_EXPENSIVE_DATA_STRUCTURES

  assert (elem < n_elem);
  assert (component < n_components());
 
  if (elem_dofs[component].empty())
    {
      error();
      return elem_dofs[component][elem][index];
    }
  else
    {
      assert (index < elem_dofs[component][elem].size()); 
      return elem_dofs[component][elem][index];
    }
  
#else
  
  // If this assert fails then you are trying to use
  // objects with mulitple DOFs per node and the library
  // has not been compiled to support this!
  assert (index == 0);
  assert (elem < n_elem);
  assert (component < n_components());

  return elem_id_map[(elem) + (component)*(n_elem)];

#endif
};



inline
unsigned int DofMap::n_dofs_at_node (const unsigned int n,
				     const unsigned int cn) const
{
#ifdef ENABLE_EXPENSIVE_DATA_STRUCTURES

  assert (!node_dofs.empty());

  unsigned int sum=0;

  for (unsigned int c=0; c<n_components(); c++)
    if (!node_dofs[c].empty())
      if ((c == cn) || (cn == invalid_number))
	{
	  assert (n < node_dofs[c].size());
	  
	  sum += node_dofs[c][n].size();
	};
  

  return sum;
  
#else

  assert (n < n_nodes);
  assert (n*n_components() < node_id_map.size());

  unsigned int sum = 0;

  for (unsigned int c=0; c<n_components(); c++)
    if ((c == cn) || (cn == invalid_number))
      if (node_dof_number(n,c) != invalid_number)
	sum++;

  return sum;
  
#endif  
};



inline
unsigned int DofMap::n_dofs_on_elem (const unsigned int e,
				     const unsigned int cn) const
{
#ifdef ENABLE_EXPENSIVE_DATA_STRUCTURES

  assert (!elem_dofs.empty());

  unsigned int sum=0;

  for (unsigned int c=0; c<n_components(); c++)
    if (!elem_dofs[c].empty())
      if ((c == cn) || (cn == invalid_number))
	{
	  assert (e < elem_dofs[c].size());
	  
	  sum += elem_dofs[c][e].size();
	};

  return sum;
  
#else

  assert (e < n_elem);
  assert (e*n_components() < elem_id_map.size());

  unsigned int sum = 0;

  for (unsigned int c=0; c<n_components(); c++)
    if ((c == cn) || (cn == invalid_number))
      if (elem_dof_number(e,c) != invalid_number)
	sum++;

  return sum;
  
#endif  
};



#ifdef ENABLE_AMR

inline
bool DofMap::is_constrained_dof (const unsigned int dof) const
{
  if (dof_constraints.count(dof) != 0)
    return true;

  return false;
};

#endif


#endif
