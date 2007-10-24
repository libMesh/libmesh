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



// C++ includes

// Local includes
#include "dof_object.h"




// ------------------------------------------------------------
// DofObject class static member initialization
const unsigned int       DofObject::invalid_id           = libMesh::invalid_uint;
const unsigned short int DofObject::invalid_processor_id = static_cast<unsigned short int>(-1);



// ------------------------------------------------------------
// DofObject class members
// Copy Constructor
DofObject::DofObject (const DofObject& dof_obj) :
  ReferenceCountedObject<DofObject>(),
#ifdef ENABLE_AMR
  old_dof_object (NULL),
#endif
  _id            (dof_obj._id),
  _processor_id  (dof_obj._processor_id),
  _n_systems     (dof_obj._n_systems),
  _n_vars        (NULL),
  _n_comp        (NULL),
  _dof_ids       (NULL)
{

  // Allocate storage for the dof numbers and copy
  // the values. 
  // IT IS UNDEFINED BEHAVIOR TO ALLOCATE AN ARRAY WITH ZERO ENTRIES,
  // IF n_systems==0, leave _n_vars, _n_comp, and _dof_ids NULL.
  if (this->n_systems() > 0)
  {
    _n_vars  = new unsigned char  [this->n_systems()];
    _n_comp  = new unsigned char* [this->n_systems()];
    _dof_ids = new unsigned int*  [this->n_systems()];
  }

  // If n_systems==0, we don't enter this for loop.
  for (unsigned int s=0; s<this->n_systems(); s++)
    {      
      _n_vars[s]  = dof_obj.n_vars(s);

      // In case you have a system with no variables, it is undefined
      // behavior (UB) to allocate a zero-length array here.
      if (this->n_vars(s) > 0)
	{
	  _n_comp[s]  = new unsigned char [this->n_vars(s)]; 
	  _dof_ids[s] = new unsigned int  [this->n_vars(s)];
	}
      
      for (unsigned int v=0; v<this->n_vars(s); v++)
	{
	  _n_comp[s][v]  = dof_obj.n_comp(s,v);

	  if (this->n_comp(s,v) > 0)
	    _dof_ids[s][v] = dof_obj.dof_number(s,v,0);
          else
	    _dof_ids[s][v] = invalid_id;
	}
    }

  // Check that everything worked
#ifdef DEBUG

  assert (this->n_systems() == dof_obj.n_systems());

  for (unsigned int s=0; s<this->n_systems(); s++)
    {
      assert (this->n_vars(s) == dof_obj.n_vars(s));

      for (unsigned int v=0; v<this->n_vars(s); v++)
	{
	  assert (this->n_comp(s,v) == dof_obj.n_comp(s,v));

	  for (unsigned int c=0; c<this->n_comp(s,v); c++)
	    assert (this->dof_number(s,v,c) == dof_obj.dof_number(s,v,c));
	}
    }
  
#endif
}




#ifdef ENABLE_AMR

void  DofObject::clear_old_dof_object ()
{
  // If we have been called before...
  // prevent a memory leak
  if (old_dof_object != NULL)
    {
      delete this->old_dof_object;
      this->old_dof_object = NULL;
    }  
}



void DofObject::set_old_dof_object ()
{
  this->clear_old_dof_object();

  assert (this->old_dof_object == NULL);
  
  // Make a new DofObject, assign a copy of \p this.
  // Make sure the copy ctor for DofObject works!!
  this->old_dof_object = new DofObject(*this);
}

#endif



void DofObject::set_n_systems (const unsigned int ns)
{
  // Check for trivial return
  if (ns == this->n_systems()) return;
 
#ifdef DEBUG

  if (ns != static_cast<unsigned int>(static_cast<unsigned char>(ns)))
    {
      std::cerr << "Unsigned char not big enough to hold ns!" << std::endl
		<< "Recompile with _n_systems set to a bigger type!"
		<< std::endl;
      
      error();
    }
					
#endif


  // Clear any existing data.  This is safe to call
  // even if we don't have any data.
  this->clear_dofs();

  // Set the new number of systems
  _n_systems = static_cast<unsigned char>(ns);
  
  // Allocate storage for the systems
  _n_vars    = new unsigned char  [this->n_systems()];
  _n_comp    = new unsigned char* [this->n_systems()];
  _dof_ids   = new unsigned int*  [this->n_systems()];

  // No variables have been declared yet.
  for (unsigned int s=0; s<this->n_systems(); s++)
    {
      _n_vars[s]  = 0;
      _n_comp[s]  = NULL;
      _dof_ids[s] = NULL;
    }
}



void DofObject::add_system()
{
  if (this->n_systems() > 0)
    {
      // Copy the old systems to temporary storage
      unsigned char  *old_n_vars  = new unsigned char  [this->n_systems()];
      unsigned char **old_n_comp  = new unsigned char* [this->n_systems()];
      unsigned int  **old_dof_ids = new unsigned int*  [this->n_systems()]; 
      
      for (unsigned int s=0; s<this->n_systems(); s++)
	{
	  old_n_vars[s]  = _n_vars[s];
	  old_n_comp[s]  = _n_comp[s];
	  old_dof_ids[s] = _dof_ids[s];
	}
  
      // Delete old storage
      assert (_n_vars  != NULL); delete [] _n_vars;  _n_vars  = NULL;
      assert (_n_comp  != NULL); delete [] _n_comp;  _n_comp  = NULL;
      assert (_dof_ids != NULL); delete [] _dof_ids; _dof_ids = NULL;
  
      // Allocate space for new system
      _n_vars  = new unsigned char  [this->n_systems()+1];
      _n_comp  = new unsigned char* [this->n_systems()+1];
      _dof_ids = new unsigned int*  [this->n_systems()+1];
      
      // Copy the other systems
      for (unsigned int s=0; s<this->n_systems(); s++)
	{
	  _n_vars[s]  = old_n_vars[s];
	  _n_comp[s]  = old_n_comp[s];
	  _dof_ids[s] = old_dof_ids[s];
	}
      
      // Delete temporary storage
      assert (old_n_vars  != NULL); delete [] old_n_vars;  old_n_vars  = NULL;
      assert (old_n_comp  != NULL); delete [] old_n_comp;  old_n_comp  = NULL;
      assert (old_dof_ids != NULL); delete [] old_dof_ids; old_dof_ids = NULL;
    }
  else
    {
      assert (_n_vars  == NULL);
      assert (_n_comp  == NULL);
      assert (_dof_ids == NULL);
      
      // Allocate space for new system
      _n_vars  = new unsigned char  [this->n_systems()+1];
      _n_comp  = new unsigned char* [this->n_systems()+1];
      _dof_ids = new unsigned int*  [this->n_systems()+1];      
    }
  
  // Initialize the new system
  _n_vars[this->n_systems()]  = 0;
  _n_comp[this->n_systems()]  = NULL;
  _dof_ids[this->n_systems()] = NULL;
  
  // Done. Don't forget to increment the number of systems!
  _n_systems++;
}



void DofObject::set_n_vars(const unsigned int s,
			   const unsigned int nvars)
{
  assert (s < this->n_systems());
  assert (_n_vars != NULL);

#ifdef DEBUG

  if (nvars != static_cast<unsigned int>(static_cast<unsigned char>(nvars)))
    {
      std::cerr << "Unsigned char not big enough to hold nvar!" << std::endl
		<< "Recompile with _n_vars set to a bigger type!"
		<< std::endl;
      
      error();
    }
					
#endif

  
  
  // If we already have memory allocated clear it.
  if (this->n_vars(s) != 0)
    {
      assert (_n_comp[s]  != NULL); delete [] _n_comp[s];  _n_comp[s]  = NULL;
      assert (_dof_ids    != NULL);
      assert (_dof_ids[s] != NULL); delete [] _dof_ids[s]; _dof_ids[s] = NULL;
    }

  // Reset the number of variables in the system  
  _n_vars[s] = static_cast<unsigned char>(nvars);

  if (this->n_vars(s) > 0)
    {
      _n_comp[s]  = new unsigned char [this->n_vars(s)];
      _dof_ids[s] = new unsigned int  [this->n_vars(s)];
      
      for (unsigned int v=0; v<this->n_vars(s); v++)
	{
	  _n_comp[s][v]  = 0;
	  _dof_ids[s][v] = invalid_id - 1;
	}
    }
}



void DofObject::set_n_comp(const unsigned int s,
			   const unsigned int var,
			   const unsigned int ncomp)
{
  assert (s < this->n_systems());
  assert (var < this->n_vars(s));
  assert (_dof_ids != NULL);
  assert (_dof_ids[s] != NULL);
  
  // Check for trivial return
  if (ncomp == this->n_comp(s,var)) return;

#ifdef DEBUG

  if (ncomp != static_cast<unsigned int>(static_cast<unsigned char>(ncomp)))
    {
      std::cerr << "Unsigned char not big enough to hold ncomp!" << std::endl
		<< "Recompile with _n_comp set to a bigger type!"
		<< std::endl;
      
      error();
    }
  
#endif

  // We use (invalid_id - 1) to signify no
  // components for this object
  if (ncomp == 0)
    {
      _dof_ids[s][var] = (invalid_id - 1);
    }
  
  assert (_n_comp != NULL);
  assert (_n_comp[s] != NULL);
    
  _n_comp[s][var]  = static_cast<unsigned char>(ncomp);
}



void DofObject::set_dof_number(const unsigned int s,
			       const unsigned int var,
			       const unsigned int comp,
			       const unsigned int dn)
{
  assert (s < this->n_systems());
  assert (var  < this->n_vars(s));
  assert (_dof_ids != NULL);
  assert (_dof_ids[s] != NULL);
  assert (comp < this->n_comp(s,var));
  
  //We intend to change all dof numbers together or not at all
  if (comp)
    assert ((dn == invalid_id && _dof_ids[s][var] == invalid_id) || 
            (dn == _dof_ids[s][var] + comp));
  else
    _dof_ids[s][var] = dn;


  assert(this->dof_number(s, var, comp) == dn);
}
