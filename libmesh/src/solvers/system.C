// $Id: system.C,v 1.11 2005-02-22 22:17:43 jwpeterson Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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
#include <sstream>   // for std::ostringstream


// Local includes
#include "system.h"
#include "equation_systems.h"
#include "libmesh.h"
#include "sparse_matrix.h"
#include "libmesh_logging.h"
#include "utility.h"



// ------------------------------------------------------------
// System implementation
System::System (EquationSystems& es,
		const std::string& name,
		const unsigned int number) :
  
  init_system              (NULL),
  assemble_system          (NULL),
  solution                 (NumericVector<Number>::build()),
  current_local_solution   (NumericVector<Number>::build()),
  _sys_name                (name),
  _sys_number              (number),
  _active                  (true),
  _can_add_vectors         (true),
  _additional_data_written (false),
  _dof_map                 (number),
  _equation_systems        (es),
  _mesh                    (es.get_mesh())
{
}



System::~System ()
{
  // Null-out the function pointers.  Since this
  // class is getting destructed it is pointless,
  // but a good habit.
  init_system = assemble_system = NULL;

  // Clear data
  this->clear ();

  assert (!libMesh::closed());
}



void System::clear ()
{
  _var_names.clear ();

  _var_type.clear ();
  
  _var_num.clear ();
  
  _dof_map.clear ();
  
  solution->clear ();

  current_local_solution->clear ();
  
  // clear any user-added vectors
  {
    for (vectors_iterator pos = _vectors.begin(); pos != _vectors.end(); ++pos)
      {
	pos->second->clear ();
	delete pos->second;
	pos->second = NULL;
      }
    
    _vectors.clear();
    _can_add_vectors = true;
  }
}



void System::init ()
{
  // First initialize any required data
  this->init_data();

  // Then call the user-provided intialization function,
  // if it was provided
  if (init_system != NULL)
    this->init_system (_equation_systems, this->name()); 
}



void System::init_data ()
{
  // Distribute the degrees of freedom on the mesh
  _dof_map.distribute_dofs (this->get_mesh());

#ifdef ENABLE_AMR

  // Recreate any hanging node constraints
  _dof_map.create_dof_constraints(this->get_mesh());

#endif

  // Resize the solution conformal to the current mesh
  solution->init (this->n_dofs(), this->n_local_dofs());

  // Resize the current_local_solution for the current mesh
  current_local_solution->init (this->n_dofs());

  // from now on, no chance to add additional vectors
  _can_add_vectors = false;

  // initialize & zero other vectors, if necessary
  for (vectors_iterator pos = _vectors.begin(); pos != _vectors.end(); ++pos)
      pos->second->init (this->n_dofs(), this->n_local_dofs());
}



void System::reinit ()
{
#ifdef ENABLE_AMR

//   // Augment the send_list for the old mesh to handle
//   // additional data dependencies introduced by refinement
//   _dof_map.augment_send_list_for_projection (this->get_mesh());

  // Re-update the data in this system using the
  // augmented send_list
  this->re_update ();
  
  // Distribute the degrees of freedom on the mesh
  _dof_map.distribute_dofs (this->get_mesh());
  
  // Recreate any hanging node constraints
  _dof_map.create_dof_constraints(this->get_mesh());

  // Clear the _vectors and reinitialize them
  for (vectors_iterator pos = _vectors.begin(); pos != _vectors.end(); ++pos)
    {
      pos->second->clear ();
      pos->second->init (this->n_dofs(), this->n_local_dofs());
    }
  
  // Project the solution to the new mesh
  {
    // Project the current local solution
    this->project_vector (*current_local_solution);

    // Update the solution based on the projected
    // current_local_solution.
    {
      solution->init (this->n_dofs(), this->n_local_dofs());
	
      assert (solution->size() == current_local_solution->size());
      assert (solution->size() == current_local_solution->local_size());

      const unsigned int first_local_dof = solution->first_local_index();
      const unsigned int local_size      = solution->local_size();
      
      for (unsigned int i=0; i<local_size; i++)
	solution->set(i+first_local_dof,
		      (*current_local_solution)(i+first_local_dof));
    } 
  }
  
#endif
}



void System::update ()
{
  const std::vector<unsigned int>& send_list = _dof_map.get_send_list ();

  // Check sizes
  assert (current_local_solution->local_size() == solution->size());
  assert (!send_list.empty());
  assert (send_list.size() <= solution->size());

  // Create current_local_solution from solution.  This will
  // put a local copy of solution into current_local_solution.
  // Only the necessary values (specified by the send_list)
  // are copied to minimize communication
  solution->localize (*current_local_solution, send_list); 
}



void System::re_update ()
{
  //const std::vector<unsigned int>& send_list = _dof_map.get_send_list ();

  // Explicitly build a send_list
  std::vector<unsigned int> send_list(solution->size());
  Utility::iota (send_list.begin(), send_list.end(), 0);
  
  // Check sizes
  assert (current_local_solution->local_size() == solution->size());
  assert (current_local_solution->size()       == solution->size());
  assert (!send_list.empty());
  assert (send_list.size() <= solution->size());

  // Create current_local_solution from solution.  This will
  // put a local copy of solution into current_local_solution.
  solution->localize (*current_local_solution, send_list);
}



void System::assemble ()
{
  // Log how long the user's matrix assembly code takes
  START_LOG("assemble()", "System");
  
  // Call the user-specified matrix assembly function
  if (this->assemble_system != NULL)
    this->assemble_system (_equation_systems, this->name());

  // Stop logging the user code
  STOP_LOG("assemble()", "System");
}



bool System::compare (const System& other_system, 
			  const Real threshold,
			  const bool verbose) const
{
  // we do not care for matrices, but for vectors
  assert (!_can_add_vectors);
  assert (!other_system._can_add_vectors);

  if (verbose)
    {
      std::cout << "  Systems \"" << _sys_name << "\"" << std::endl;
      std::cout << "   comparing matrices not supported." << std::endl;
      std::cout << "   comparing names...";
    }

  // compare the name: 0 means identical
  const int name_result = _sys_name.compare(other_system.name());
  if (verbose)
    {
      if (name_result == 0)
	std::cout << " identical." << std::endl;
      else
	std::cout << "  names not identical." << std::endl;
      std::cout << "   comparing solution vector...";
    }


  // compare the solution: -1 means identical
  const int solu_result = solution->compare (*other_system.solution.get(),
					     threshold);

  if (verbose)
    {
      if (solu_result == -1)
	std::cout << " identical up to threshold." << std::endl;
      else
	std::cout << "  first difference occured at index = " 
		  << solu_result << "." << std::endl;
    }


  // safety check, whether we handle at least the same number
  // of vectors
  std::vector<int> ov_result;

  if (this->n_vectors() != other_system.n_vectors())
    {
      if (verbose)
        {
	  std::cout << "   Fatal difference. This system handles " 
		    << this->n_vectors() << " add'l vectors," << std::endl
		    << "   while the other system handles "
		    << other_system.n_vectors() 
		    << " add'l vectors." << std::endl
		    << "   Aborting comparison." << std::endl;
	}
      return false;
    }
  else if (this->n_vectors() == 0)
    {
      // there are no additional vectors...
      ov_result.clear ();
    }
  else
    {
      // compare other vectors
      for (const_vectors_iterator pos = _vectors.begin();
	   pos != _vectors.end(); ++pos)
        {
	  if (verbose)
	      std::cout << "   comparing vector \""
			<< pos->first << "\" ...";

	  // assume they have the same name
	  const NumericVector<Number>& other_system_vector = 
	      other_system.get_vector(pos->first);

	  ov_result.push_back(pos->second->compare (other_system_vector,
						    threshold));

	  if (verbose)
	    {
	      if (ov_result[ov_result.size()-1] == -1)
		std::cout << " identical up to threshold." << std::endl;
	      else
		std::cout << " first difference occured at" << std::endl
			  << "   index = " << ov_result[ov_result.size()-1] << "." << std::endl;
	    }

	}

    } // finished comparing additional vectors


  bool overall_result;
       
  // sum up the results
  if ((name_result==0) && (solu_result==-1))
    {
      if (ov_result.size()==0)
	overall_result = true;
      else
        {
	  bool ov_identical;
	  unsigned int n    = 0;
	  do
	    {
	      ov_identical = (ov_result[n]==-1);
	      n++;
	    }
	  while (ov_identical && n<ov_result.size());
	  overall_result = ov_identical;
	}
    }
  else
    overall_result = false;

  if (verbose)
    {
      std::cout << "   finished comparisons, ";
      if (overall_result)
	std::cout << "found no differences." << std::endl << std::endl;
      else 
	std::cout << "found differences." << std::endl << std::endl;
    }
	  
  return overall_result;
}



void System::update_global_solution (std::vector<Number>& global_soln) const
{
  global_soln.resize (solution->size());

  solution->localize (global_soln);
}



void System::update_global_solution (std::vector<Number>& global_soln,
				     const unsigned int   dest_proc) const
{
  global_soln.resize        (solution->size());

  solution->localize_to_one (global_soln, dest_proc);
}



NumericVector<Number> & System::add_vector (const std::string& vec_name)
{
  // only add vectors before initializing...
  if (!_can_add_vectors)
    {
      std::cerr << "ERROR: Too late.  Cannot add vectors to the system after initialization"
		<< std::endl
		<< " any more.  You should have done this earlier."
		<< std::endl;
      error();
    }

  // Return the vector if it is already there.
  if (this->have_vector(vec_name))
    {
      return *(_vectors[vec_name]);
    }

  // Otherwise build the vector and return it.
  NumericVector<Number>* buf = NumericVector<Number>::build().release();
  _vectors.insert (std::make_pair (vec_name, buf));

  return *buf;
}



const NumericVector<Number> & System::get_vector (const std::string& vec_name) const
{
  // Make sure the vector exists
  const_vectors_iterator pos = _vectors.find(vec_name);
  
  if (pos == _vectors.end())
    {
      std::cerr << "ERROR: vector "
		<< vec_name
		<< " does not exist in this system!"
		<< std::endl;      
      error();
    }
  
  return *(pos->second);
}



NumericVector<Number> & System::get_vector (const std::string& vec_name)
{
  // Make sure the vector exists
  vectors_iterator pos = _vectors.find(vec_name);
  
  if (pos == _vectors.end())
    {
      std::cerr << "ERROR: vector "
		<< vec_name
		<< " does not exist in this system!"
		<< std::endl;      
      error();
    }
  
  return *(pos->second);
}



void System::add_variable (const std::string& var,
			   const FEType& type)
{  
  // Make sure the variable isn't there already
  if (_var_num.count(var))
    {
      std::cerr << "ERROR: variable "
		<< var
		<< " has already been added for this system!"
		<< std::endl;
      error();
    }
  
  // Add the variable to the list
  _var_names.push_back (var);
  _var_type[var]  = type;
  _var_num[var]   = (this->n_vars()-1);

  // Add the variable to the _dof_map
  _dof_map.add_variable (type);
}



void System::add_variable (const std::string& var,
			       const Order order,
			       const FEFamily family)
{
  FEType fe_type(order, family);
  
  this->add_variable(var, fe_type);
}



unsigned short int System::variable_number (const std::string& var) const
{
  // Make sure the variable exists
  std::map<std::string, unsigned short int>::const_iterator
    pos = _var_num.find(var);
  
  if (pos == _var_num.end())
    {
      std::cerr << "ERROR: variable "
		<< var
		<< " does not exist in this system!"
		<< std::endl;      
      error();
    }
  
  return pos->second;
}



std::string System::get_info() const
{
  std::ostringstream out;

  
  const std::string& sys_name = this->name();
      
  out << "   System \"" << sys_name << "\"\n"
      << "    Type \""  << this->system_type() << "\"\n"
      << "    Variables=";
  
  for (unsigned int vn=0; vn<this->n_vars(); vn++)
      out << "\"" << this->variable_name(vn) << "\" ";
     
  out << '\n';

#ifndef ENABLE_INFINITE_ELEMENTS
  out << "    Finite Element Types=";
  for (unsigned int vn=0; vn<this->n_vars(); vn++)
  {
      out << "\"" << this->get_dof_map().variable_type(vn).family << "\" ";
  }
#else
  out << "    Finite Element Types=";
  for (unsigned int vn=0; vn<this->n_vars(); vn++)
  {
      out << "\"" << this->get_dof_map().variable_type(vn).family << "\", ";
      out << "\"" << this->get_dof_map().variable_type(vn).radial_family << "\" ";
  }

  out << '\n' << "    Infinite Element Mapping=";
  for (unsigned int vn=0; vn<this->n_vars(); vn++)
  {
      out << "\"" << this->get_dof_map().variable_type(vn).inf_map << "\" ";
  }
#endif      

  out << '\n';
      
  out << "    Approximation Orders=";
  for (unsigned int vn=0; vn<this->n_vars(); vn++)
  {
#ifndef ENABLE_INFINITE_ELEMENTS
      out << "\"" << this->get_dof_map().variable_type(vn).order << "\" ";
#else
      out << "\"" << this->get_dof_map().variable_type(vn).order << "\", ";
      out << "\"" << this->get_dof_map().variable_type(vn).radial_order << "\" ";
#endif
  }

  out << '\n';
      
  out << "    n_dofs()="             << this->n_dofs()             << '\n';
  out << "    n_local_dofs()="       << this->n_local_dofs()       << '\n';
#ifdef ENABLE_AMR
  out << "    n_constrained_dofs()=" << this->n_constrained_dofs() << '\n';
#endif

  out << "    " << "n_vectors()="  << this->n_vectors()  << '\n';
//   out << "    " << "n_additional_matrices()=" << this->n_additional_matrices() << '\n';
  
  return out.str();
}



void System::attach_init_function (void fptr(EquationSystems& es,
					     const std::string& name))
{
  assert (fptr != NULL);
  
  init_system = fptr;
}



void System::attach_assemble_function (void fptr(EquationSystems& es,
						 const std::string& name))
{
  assert (fptr != NULL);
  
  assemble_system = fptr;  
}
