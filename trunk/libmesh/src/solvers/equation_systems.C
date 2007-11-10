// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson
  
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


// System includes
#include <sstream>

// Local Includes
#include "explicit_system.h"
#include "fe_interface.h"
#include "frequency_system.h"
#include "linear_implicit_system.h"
#include "mesh_refinement.h"
#include "newmark_system.h"
#include "nonlinear_implicit_system.h"
#include "parallel.h"
#include "transient_system.h"
#include "dof_map.h"
#include "mesh_base.h"
#include "elem.h"

// Include the systems before this one to avoid
// overlapping forward declarations.
#include "equation_systems.h"

// Forward Declarations




// ------------------------------------------------------------
// EquationSystems class implementation
EquationSystems::EquationSystems (MeshBase& m, MeshData* mesh_data) :
  _mesh      (m),
  _mesh_data (mesh_data)
{
  // Set default parameters
  this->parameters.set<Real>        ("linear solver tolerance") = TOLERANCE * TOLERANCE;
  this->parameters.set<unsigned int>("linear solver maximum iterations") = 5000;
}



EquationSystems::~EquationSystems ()
{
  this->clear ();
}



void EquationSystems::clear ()
{
  // Clear any additional parameters
  parameters.clear ();

  // clear the systems.  We must delete them
  // since we newed them!
  while (!_systems.empty())
    {
      system_iterator pos = _systems.begin();
      
      System *sys = pos->second;
      delete sys;
      sys = NULL;
      
      _systems.erase (pos);
    }
}



void EquationSystems::init ()
{
  const unsigned int n_sys = this->n_systems();

  assert (n_sys != 0);

  // Distribute the mesh if possible
  if (libMesh::n_processors() > 1)
    _mesh.delete_remote_elements();
  
  // Tell all the \p DofObject entities how many systems
  // there are.
  {
    MeshBase::node_iterator       node_it  = _mesh.nodes_begin();
    const MeshBase::node_iterator node_end = _mesh.nodes_end();

    for ( ; node_it != node_end; ++node_it)
      (*node_it)->set_n_systems(n_sys);
    
    MeshBase::element_iterator       elem_it  = _mesh.elements_begin();
    const MeshBase::element_iterator elem_end = _mesh.elements_end();
    
    for ( ; elem_it != elem_end; ++elem_it)
      (*elem_it)->set_n_systems(n_sys);
  }

  for (unsigned int i=0; i != this->n_systems(); ++i)
    this->get_system(i).init();
}



void EquationSystems::reinit ()
{
  const unsigned int n_sys = this->n_systems();

  assert (n_sys != 0);
  
#ifdef DEBUG
  // Make sure all the \p DofObject entities know how many systems
  // there are.
  {
    // All the nodes
    MeshBase::node_iterator       node_it  = _mesh.nodes_begin();
    const MeshBase::node_iterator node_end = _mesh.nodes_end();

    for ( ; node_it != node_end; ++node_it)
      assert((*node_it)->n_systems() == n_sys);
    
    // All the elements
    MeshBase::element_iterator       elem_it  = _mesh.elements_begin();
    const MeshBase::element_iterator elem_end = _mesh.elements_end();
    
    for ( ; elem_it != elem_end; ++elem_it)
      assert((*elem_it)->n_systems() == n_sys);
  }
#endif

  // Localize each system's vectors
  for (unsigned int i=0; i != this->n_systems(); ++i)
    this->get_system(i).re_update();

#ifdef ENABLE_AMR

  bool dof_constraints_created = false;
  bool mesh_changed = false;

  // FIXME: For backwards compatibility, assume
  // refine_and_coarsen_elements or refine_uniformly have already
  // been called
  {
    for (unsigned int i=0; i != this->n_systems(); ++i)
      {
	System &sys = this->get_system(i);
        sys.get_dof_map().distribute_dofs(_mesh);
        sys.get_dof_map().create_dof_constraints(_mesh);
        sys.prolong_vectors();
      }
    mesh_changed = true;
    dof_constraints_created = true;
  }
  
  // FIXME: Where should the user set maintain_level_one now??
  // Don't override previous settings, for now

  MeshRefinement mesh_refine(_mesh);

  mesh_refine.face_level_mismatch_limit() = false;

  // Try to coarsen the mesh, then restrict each system's vectors
  // if necessary
  if (mesh_refine.coarsen_elements())
    {
      for (unsigned int i=0; i != this->n_systems(); ++i)
        {
	  System &sys = this->get_system(i);
          if (!dof_constraints_created)
            {
              sys.get_dof_map().distribute_dofs(_mesh);
              sys.get_dof_map().create_dof_constraints(_mesh);
            }
          sys.restrict_vectors();
        }
      mesh_changed = true;
      dof_constraints_created = true;
    }

  // Once vectors are all restricted, we can delete
  // children of coarsened elements
  if (mesh_changed)
    this->get_mesh().contract();
  
  // Try to refine the mesh, then prolong each system's vectors
  // if necessary
  if (mesh_refine.refine_elements())
    {
      for (unsigned int i=0; i != this->n_systems(); ++i)
        {
	  System &sys = this->get_system(i);
          if (!dof_constraints_created)
            {
              sys.get_dof_map().distribute_dofs(_mesh);
              sys.get_dof_map().create_dof_constraints(_mesh);
            }
          sys.prolong_vectors();
        }
      mesh_changed = true;
      dof_constraints_created = true;
    }

  // If the mesh has changed, systems will need to create new dof
  // constraints and update their global solution vectors
  if (mesh_changed)
    {
      for (unsigned int i=0; i != this->n_systems(); ++i)
        this->get_system(i).reinit();
    }
#endif // #ifdef ENABLE_AMR
}



void EquationSystems::allgather ()
{
  // A serial mesh means nothing needs to be done
  if (_mesh.is_serial())
    return;

  const unsigned int n_sys = this->n_systems();

  assert (n_sys != 0);

  // Gather the mesh
  _mesh.allgather();
  
  // Tell all the \p DofObject entities how many systems
  // there are.
  {
    MeshBase::node_iterator       node_it  = _mesh.nodes_begin();
    const MeshBase::node_iterator node_end = _mesh.nodes_end();

    for ( ; node_it != node_end; ++node_it)
      (*node_it)->set_n_systems(n_sys);
    
    MeshBase::element_iterator       elem_it  = _mesh.elements_begin();
    const MeshBase::element_iterator elem_end = _mesh.elements_end();
    
    for ( ; elem_it != elem_end; ++elem_it)
      (*elem_it)->set_n_systems(n_sys);
  }

  // And distribute each system's dofs
  for (unsigned int i=0; i != this->n_systems(); ++i)
    this->get_system(i).get_dof_map().distribute_dofs(_mesh);
}




void EquationSystems::update ()
{
  // Localize each system's vectors
  for (unsigned int i=0; i != this->n_systems(); ++i)
    this->get_system(i).update();
}



System & EquationSystems::add_system (const std::string& sys_type,
				      const std::string& name)
{
  // Build a Newmark system
  if      (sys_type == "Newmark")
    this->add_system<NewmarkSystem> (name);

  // Build an Explicit system
  else if ((sys_type == "Explicit"))
    this->add_system<ExplicitSystem> (name);

  // Build an Implicit system
  else if ((sys_type == "Implicit") ||
	   (sys_type == "Steady"  ))
    this->add_system<ImplicitSystem> (name);

  // build a transient implicit linear system
  else if ((sys_type == "Transient") ||
	   (sys_type == "TransientImplicit") ||
	   (sys_type == "TransientLinearImplicit"))
    this->add_system<TransientLinearImplicitSystem> (name);

  // build a transient implicit nonlinear system
  else if (sys_type == "TransientNonlinearImplicit")
    this->add_system<TransientNonlinearImplicitSystem> (name);

  // build a transient explicit system
  else if (sys_type == "TransientExplicit")
    this->add_system<TransientExplicitSystem> (name);

  // build a linear implicit sytsem
  else if (sys_type == "LinearImplicit")
    this->add_system<LinearImplicitSystem> (name);

#if defined(USE_COMPLEX_NUMBERS)
  // build a frequency system
  else if (sys_type == "Frequency")
    this->add_system<FrequencySystem> (name);
#endif

  else
    {
      std::cerr << "ERROR: Unknown system type: "
		<< sys_type
		<< std::endl;
      error();
    }

  // Return a reference to the new system
  //return (*this)(name);
  return this->get_system(name);
}






void EquationSystems::delete_system (const std::string& name)
{
  deprecated();

  if (!_systems.count(name))
    {
      std::cerr << "ERROR: no system named "
                << name  << std::endl;

      error();
    }
  
  delete _systems[name];
  
  _systems.erase (name);
}



void EquationSystems::solve ()
{
  assert (this->n_systems());

  for (unsigned int i=0; i != this->n_systems(); ++i)
    this->get_system(i).solve();
}
 
 
void EquationSystems::build_variable_names (std::vector<std::string>& var_names) const
{
  assert (this->n_systems());
  
  var_names.resize (this->n_vars());

  unsigned int var_num=0;
  
  const_system_iterator       pos = _systems.begin();
  const const_system_iterator end = _systems.end();

  for (; pos != end; ++pos)
    for (unsigned int vn=0; vn<pos->second->n_vars(); vn++)
      var_names[var_num++] = pos->second->variable_name(vn);       
}



void EquationSystems::build_solution_vector (std::vector<Number>&,
					     const std::string&,
					     const std::string&) const
{
  //TODO:[BSK] re-implement this from the method below
  error();

//   // Get a reference to the named system
//   const System& system = this->get_system(system_name);

//   // Get the number associated with the variable_name we are passed
//   const unsigned short int variable_num = system.variable_number(variable_name);

//   // Get the dimension of the current mesh
//   const unsigned int dim = _mesh.mesh_dimension();

//   // If we're on processor 0, allocate enough memory to hold the solution.
//   // Since we're only looking at one variable, there will be one solution value
//   // for each node in the mesh.
//   if (_mesh.processor_id() == 0)
//     soln.resize(_mesh.n_nodes());

//   // Vector to hold the global solution from all processors
//   std::vector<Number> sys_soln;
  
//   // Update the global solution from all processors
//   system.update_global_solution (sys_soln, 0);
  
//   // Temporary vector to store the solution on an individual element.
//   std::vector<Number>       elem_soln;   

//   // The FE solution interpolated to the nodes
//   std::vector<Number>       nodal_soln;  

//   // The DOF indices for the element
//   std::vector<unsigned int> dof_indices; 

//   // Determine the finite/infinite element type used in this system 
//   const FEType& fe_type    = system.variable_type(variable_num);

//   // Define iterators to iterate over all the elements of the mesh
//   const_active_elem_iterator       it (_mesh.elements_begin());
//   const const_active_elem_iterator end(_mesh.elements_end());

//   // Loop over elements
//   for ( ; it != end; ++it)
//     {
//       // Convenient shortcut to the element pointer
//       const Elem* elem = *it;

//       // Fill the dof_indices vector for this variable
//       system.get_dof_map().dof_indices(elem,
// 				       dof_indices,
// 				       variable_num);

//       // Resize the element solution vector to fit the
//       // dof_indices for this element.
//       elem_soln.resize(dof_indices.size());

//       // Transfer the system solution to the element
//       // solution by mapping it through the dof_indices vector.
//       for (unsigned int i=0; i<dof_indices.size(); i++)
// 	elem_soln[i] = sys_soln[dof_indices[i]];

//       // Using the FE interface, compute the nodal_soln
//       // for the current elemnt type given the elem_soln
//       FEInterface::nodal_soln (dim,
// 			       fe_type,
// 			       elem,
// 			       elem_soln,
// 			       nodal_soln);

//       // Sanity check -- make sure that there are the same number
//       // of entries in the nodal_soln as there are nodes in the
//       // element!
//       assert (nodal_soln.size() == elem->n_nodes());

//       // Copy the nodal solution over into the correct place in
//       // the global soln vector which will be returned to the user.
//       for (unsigned int n=0; n<elem->n_nodes(); n++)
// 	soln[elem->node(n)] = nodal_soln[n];
//     }
}




void EquationSystems::build_solution_vector (std::vector<Number>& soln) const
{
  // This function must be run on all processors at once
  parallel_only();

  assert (this->n_systems());

  const unsigned int dim = _mesh.mesh_dimension();
  const unsigned int nn  = _mesh.n_nodes();
  const unsigned int nv  = this->n_vars();

  // We'd better have a contiguous node numbering
  assert (nn == _mesh.max_node_id());

  // allocate storage to hold
  // (number_of_nodes)*(number_of_variables) entries.
  soln.resize(nn*nv);

  std::vector<Number>             sys_soln;

  // (Note that we use an unsigned short int here even though an
  // unsigned char would be more that sufficient.  The MPI 1.1
  // standard does not require that MPI_SUM, MPI_PROD etc... be
  // implemented for char data types. 12/23/2003 - BSK)  
  std::vector<unsigned short int> node_conn(nn);

  
  // Zero out the soln vector
  std::fill (soln.begin(),       soln.end(),       libMesh::zero);

  
  // Get the number of elements that share each node.  We will
  // compute the average value at each node.  This is particularly
  // useful for plotting discontinuous data.
  MeshBase::element_iterator       e_it  = _mesh.active_local_elements_begin();
  const MeshBase::element_iterator e_end = _mesh.active_local_elements_end(); 

  for ( ; e_it != e_end; ++e_it)
    for (unsigned int n=0; n<(*e_it)->n_nodes(); n++)
      node_conn[(*e_it)->node(n)]++;

  // Gather the distributed node_conn arrays in the case of
  // multiple processors
  Parallel::sum(node_conn);

  unsigned int var_num=0;

  // For each system in this EquationSystems object,
  // update the global solution and if we are on processor 0,
  // loop over the elements and build the nodal solution
  // from the element solution.  Then insert this nodal solution
  // into the vector passed to build_solution_vector.
  const_system_iterator       pos = _systems.begin();
  const const_system_iterator end = _systems.end();

  for (; pos != end; ++pos)
    {  
      const System& system  = *(pos->second);
      const unsigned int nv_sys = system.n_vars();
      
      system.update_global_solution (sys_soln);
      
      std::vector<Number>       elem_soln;   // The finite element solution
      std::vector<Number>       nodal_soln;  // The FE solution interpolated to the nodes
      std::vector<unsigned int> dof_indices; // The DOF indices for the finite element 
      
      for (unsigned int var=0; var<nv_sys; var++)
	{
	  const FEType& fe_type    = system.variable_type(var);
	  
	  MeshBase::element_iterator       it  = _mesh.active_local_elements_begin();
	  const MeshBase::element_iterator end = _mesh.active_local_elements_end(); 

	  for ( ; it != end; ++it)
	    {
	      const Elem* elem = *it;
	      system.get_dof_map().dof_indices (elem, dof_indices, var);
	      
	      elem_soln.resize(dof_indices.size());
	      
	      for (unsigned int i=0; i<dof_indices.size(); i++)
		elem_soln[i] = sys_soln[dof_indices[i]];
		  
	      FEInterface::nodal_soln (dim,
				       fe_type,
				       elem,
				       elem_soln,
				       nodal_soln);

#ifdef ENABLE_INFINITE_ELEMENTS
	      // infinite elements should be skipped...
	      if (!elem->infinite())
#endif
		{ 
		  assert (nodal_soln.size() == elem->n_nodes());
		  
		  for (unsigned int n=0; n<elem->n_nodes(); n++)
		    {
		      assert (node_conn[elem->node(n)] != 0);
		      soln[nv*(elem->node(n)) + (var + var_num)] +=
			nodal_soln[n]/static_cast<Real>(node_conn[elem->node(n)]);
		    }
		}
	    }	 
	}

      var_num += nv_sys;
    }

  // Now each processor has computed contriburions to the
  // soln vector.  Gather them all up.
  Parallel::sum(soln);
}




void EquationSystems::build_discontinuous_solution_vector (std::vector<Number>& soln) const
{
  assert (this->n_systems());

  const unsigned int dim = _mesh.mesh_dimension();
  const unsigned int nv  = this->n_vars();
  unsigned int tw=0;

  // get the total weight
  {
    MeshBase::element_iterator       it  = _mesh.active_elements_begin();
    const MeshBase::element_iterator end = _mesh.active_elements_end(); 

    for ( ; it != end; ++it)
      tw += (*it)->n_nodes();
  }
  

  // Only if we are on processor zero, allocate the storage
  // to hold (number_of_nodes)*(number_of_variables) entries.
  if (_mesh.processor_id() == 0)
    soln.resize(tw*nv);

  std::vector<Number> sys_soln; 

  
  unsigned int var_num=0;

  // For each system in this EquationSystems object,
  // update the global solution and if we are on processor 0,
  // loop over the elements and build the nodal solution
  // from the element solution.  Then insert this nodal solution
  // into the vector passed to build_solution_vector.
  const_system_iterator       pos = _systems.begin();
  const const_system_iterator end = _systems.end();

  for (; pos != end; ++pos)
    {  
      const System& system  = *(pos->second);
      const unsigned int nv_sys = system.n_vars();
      
      system.update_global_solution (sys_soln, 0);
      
      if (_mesh.processor_id() == 0)
	{
	  std::vector<Number>       elem_soln;   // The finite element solution
	  std::vector<Number>       nodal_soln;  // The FE solution interpolated to the nodes
	  std::vector<unsigned int> dof_indices; // The DOF indices for the finite element 
	      
	  for (unsigned int var=0; var<nv_sys; var++)
	    {
	      const FEType& fe_type    = system.variable_type(var);

	      MeshBase::element_iterator       it  = _mesh.active_elements_begin();
	      const MeshBase::element_iterator end = _mesh.active_elements_end(); 

	      unsigned int nn=0;
	      
	      for ( ; it != end; ++it)
		{
		  const Elem* elem = *it;
		  system.get_dof_map().dof_indices (elem, dof_indices, var);
		  
		  elem_soln.resize(dof_indices.size());
		  
		  for (unsigned int i=0; i<dof_indices.size(); i++)
		    elem_soln[i] = sys_soln[dof_indices[i]];
		  
		  FEInterface::nodal_soln (dim,
					   fe_type,
					   elem,
					   elem_soln,
					   nodal_soln);

#ifdef ENABLE_INFINITE_ELEMENTS
		  // infinite elements should be skipped...
		  if (!elem->infinite())
#endif
		    { 
		      assert (nodal_soln.size() == elem->n_nodes());
		  
		      for (unsigned int n=0; n<elem->n_nodes(); n++)
			{
			  soln[nv*(nn++) + (var + var_num)] +=
			    nodal_soln[n];
			}
		    }
		}
	    }	 
	}

      var_num += nv_sys;
    }
}



bool EquationSystems::compare (const EquationSystems& other_es, 
			       const Real threshold,
			       const bool verbose) const
{
  // safety check, whether we handle at least the same number
  // of systems
  std::vector<bool> os_result;

  if (this->n_systems() != other_es.n_systems())
    {
      if (verbose)
        {
	  std::cout << "  Fatal difference. This system handles " 
		    << this->n_systems() << " systems," << std::endl
		    << "  while the other system handles "
		    << other_es.n_systems() 
		    << " systems." << std::endl
		    << "  Aborting comparison." << std::endl;
	}
      return false;
    }
  else
    {
      // start comparing each system      
      const_system_iterator       pos = _systems.begin();
      const const_system_iterator end = _systems.end();

      for (; pos != end; ++pos)
        {
	  const std::string& sys_name = pos->first;
	  const System&  system        = *(pos->second);
      
	  // get the other system
	  const System& other_system   = other_es.get_system (sys_name);

	  os_result.push_back (system.compare (other_system, threshold, verbose));

	}

    }
  

  // sum up the results
  if (os_result.size()==0)
    return true;
  else
    {
      bool os_identical;
      unsigned int n = 0;
	  do
	    {
	      os_identical = os_result[n];
	      n++;
	    }
	  while (os_identical && n<os_result.size());
	  return os_identical;
	}
}



std::string EquationSystems::get_info () const
{
  std::ostringstream out;
  
  out << " EquationSystems\n"
      << "  n_systems()=" << this->n_systems() << '\n';

  // Print the info for the individual systems
  const_system_iterator       pos = _systems.begin();
  const const_system_iterator end = _systems.end();

  for (; pos != end; ++pos)
    out << pos->second->get_info();

  
//   // Possibly print the parameters  
//   if (!this->parameters.empty())
//     {  
//       out << "  n_parameters()=" << this->n_parameters() << '\n';
//       out << "   Parameters:\n";
      
//       for (std::map<std::string, Real>::const_iterator
// 	     param = _parameters.begin(); param != _parameters.end();
// 	   ++param)
// 	out << "    "
// 	    << "\""
// 	    << param->first
// 	    << "\""
// 	    << "="
// 	    << param->second
// 	    << '\n';
//     }
  
  return out.str();
}



void EquationSystems::print_info (std::ostream& os) const
{
  os << this->get_info()
     << std::endl;
}



std::ostream& operator << (std::ostream& os, const EquationSystems& es)
{
  es.print_info(os);
  return os;
}



unsigned int EquationSystems::n_vars () const
{
  unsigned int tot=0;
  
  const_system_iterator       pos = _systems.begin();
  const const_system_iterator end = _systems.end();

  for (; pos != end; ++pos)
    tot += pos->second->n_vars();

  return tot;
}



unsigned int EquationSystems::n_dofs () const
{
  unsigned int tot=0;

  const_system_iterator       pos = _systems.begin();
  const const_system_iterator end = _systems.end();

  for (; pos != end; ++pos)
    tot += pos->second->n_dofs();

  return tot;      
}




unsigned int EquationSystems::n_active_dofs () const
{
  unsigned int tot=0;

  const_system_iterator       pos = _systems.begin();
  const const_system_iterator end = _systems.end();

  for (; pos != end; ++pos)
    tot += pos->second->n_active_dofs();

  return tot;      
}


void EquationSystems::_add_system_to_nodes_and_elems()
{
  // All the nodes
  MeshBase::node_iterator       node_it  = _mesh.nodes_begin();
  const MeshBase::node_iterator node_end = _mesh.nodes_end();
 
  for ( ; node_it != node_end; ++node_it)
    (*node_it)->add_system();
 
  // All the elements
  MeshBase::element_iterator       elem_it  = _mesh.elements_begin();
  const MeshBase::element_iterator elem_end = _mesh.elements_end();
 
  for ( ; elem_it != elem_end; ++elem_it)
    (*elem_it)->add_system();
}
