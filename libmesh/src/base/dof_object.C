// $Id: dof_object.C,v 1.9 2004-01-11 15:56:46 benkirk Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2004  Benjamin S. Kirk, John W. Peterson
  
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




// Local includes
#include "dof_object.h"




// ------------------------------------------------------------
// DofObject class static member initialization
const unsigned int       DofObject::invalid_id           = libMesh::invalid_uint;
const unsigned short int DofObject::invalid_processor_id = static_cast<unsigned short int>(-1);



// ------------------------------------------------------------
// DofObject class members
#ifdef ENABLE_AMR

void  DofObject::clear_old_dof_object ()
{
  // If we have been called before...
  // prevent a memory leak
  if (old_dof_object != NULL)
    {
      delete old_dof_object;
      old_dof_object = NULL;
    }  
}



void DofObject::set_old_dof_object ()
{
  this->clear_old_dof_object();

  assert (old_dof_object == NULL);
  
  // Make a new DofObject, assign a copy of \p this
  old_dof_object = new DofObject(*this);
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


  
#ifdef ENABLE_EXPENSIVE_DATA_STRUCTURES

  // Allocate storage for the systems
  _n_vars    = new unsigned char  [this->n_systems()];
  _n_comp    = new unsigned char* [this->n_systems()];
  _dof_ids   = new unsigned int** [this->n_systems()];

  // No variables have been declared yet.
  for (unsigned int s=0; s<this->n_systems(); s++)
    {
      _n_vars[s]  = 0;
      _n_comp[s]  = NULL;
      _dof_ids[s] = NULL;
    }
  
#else
  
  // Allocate storage for the systems
  _n_vars    = new unsigned char  [this->n_systems()];
  _dof_ids   = new unsigned int*  [this->n_systems()];
 
  // No variables have been declared yet.
  for (unsigned int s=0; s<this->n_systems(); s++)
    {
      _n_vars[s]  = 0;
      _dof_ids[s] = NULL;
    }
    
#endif
}



void DofObject::add_system()
{
#ifdef ENABLE_EXPENSIVE_DATA_STRUCTURES

  if (this->n_systems() > 0)
    {
      // Copy the old systems to temporary storage
      unsigned char  *old_n_vars  = new unsigned char  [this->n_systems()];
      unsigned char **old_n_comp  = new unsigned char* [this->n_systems()];
      unsigned int ***old_dof_ids = new unsigned int** [this->n_systems()]; 
      
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
      _dof_ids = new unsigned int** [this->n_systems()+1];
      
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
      _dof_ids = new unsigned int** [this->n_systems()+1];      
    }
  
  // Initialize the new system
  _n_vars[this->n_systems()]  = 0;
  _n_comp[this->n_systems()]  = NULL;
  _dof_ids[this->n_systems()] = NULL;

#else

  if (this->n_systems() > 0)
    {
      // Copy the old systems to temporary storage
      unsigned char *old_n_vars  = new unsigned char [this->n_systems()];
      unsigned int **old_dof_ids = new unsigned int* [this->n_systems()]; 
      
      for (unsigned int s=0; s<this->n_systems(); s++)
	{
	  old_n_vars[s]  = _n_vars[s];
	  old_dof_ids[s] = _dof_ids[s];
	}
      
      // Delete old storage
      assert (_n_vars  != NULL); delete [] _n_vars;  _n_vars  = NULL;
      assert (_dof_ids != NULL); delete [] _dof_ids; _dof_ids = NULL;
      
      // Allocate space for new system
      _n_vars  = new unsigned char [this->n_systems()+1];
      _dof_ids = new unsigned int* [this->n_systems()+1];
      
      // Copy the other systems
      for (unsigned int s=0; s<this->n_systems(); s++)
	{
	  _n_vars[s]  = old_n_vars[s];
	  _dof_ids[s] = old_dof_ids[s];
	}

      // Delete temporary storage
      assert (old_n_vars  != NULL); delete [] old_n_vars;  old_n_vars  = NULL;
      assert (old_dof_ids != NULL); delete [] old_dof_ids; old_dof_ids = NULL;
    }
  else
    {
      assert (_n_vars  == NULL);
      assert (_dof_ids == NULL);
      
      // Allocate space for new system
      _n_vars  = new unsigned char [this->n_systems()+1];
      _dof_ids = new unsigned int* [this->n_systems()+1];
    }
  
  // Initialize the new system
  _n_vars[this->n_systems()]  = 0;
  _dof_ids[this->n_systems()] = NULL;

#endif

  
  // Done. Don't forget to increment the number of systems!
  _n_systems++;
}



void DofObject::set_n_vars(const unsigned int s,
			   const unsigned int nvars)
{
  assert (s < this->n_systems());
  assert (_n_vars != NULL);

  // Check for trivial return
  if (nvars == this->n_vars(s)) return;

#ifdef DEBUG

  if (nvars != static_cast<unsigned int>(static_cast<unsigned char>(nvars)))
    {
      std::cerr << "Unsigned char not big enough to hold nvar!" << std::endl
		<< "Recompile with _n_vars set to a bigger type!"
		<< std::endl;
      
      error();
    }
					
#endif

  
  
#ifdef ENABLE_EXPENSIVE_DATA_STRUCTURES
  
  // If we already have memory allocated clear it.
  if (this->n_vars(s) != 0)
    {
      assert (_dof_ids    != NULL);
      assert (_dof_ids[s] != NULL);
      assert (_n_comp     != NULL);
      assert (_n_comp[s]  != NULL);

      for (unsigned int v=0; v<this->n_vars(s); v++)
	if (this->n_comp(s,v) != 0)
	  {
	    assert (_dof_ids[s][v] != NULL); delete [] _dof_ids[s][v]; _dof_ids[s][v] = NULL;
	    
	    _n_comp[s][v] = 0;
	  }
      
      assert (_n_comp[s]  != NULL); delete [] _n_comp[s];  _n_comp[s]  = NULL;
      assert (_dof_ids[s] != NULL); delete [] _dof_ids[s]; _dof_ids[s] = NULL;
    }

  // Reset the number of variables in the system  
  _n_vars[s] = static_cast<unsigned char>(nvars);

  if (this->n_vars(s) > 0)
    {
      _n_comp[s]  = new unsigned char [this->n_vars(s)];
      _dof_ids[s] = new unsigned int* [this->n_vars(s)];
      
      for (unsigned int v=0; v<this->n_vars(s); v++)
	{
	  _n_comp[s][v]  = 0;
	  _dof_ids[s][v] = NULL;
	}
    }


#else

  // If we already have memory allocated clear it.
  if (this->n_vars(s) > 0)
    {
      assert (_dof_ids[s] != NULL); delete [] _dof_ids[s]; _dof_ids[s] = NULL;
    }

  // Reset the number of variables in the system    
  _n_vars[s] = static_cast<unsigned char>(nvars);

  if (this->n_vars(s) > 0)
    {
      _dof_ids[s] = new unsigned int [this->n_vars(s)];
      
      // We use (invalid_id - 1) to signify n_comp(s,var) = 0.
      // This eliminates an additional array that would determine
      // when a variable has no components or one component active
      // on this object. Note that without expensive data structures
      // these are the only possiblities.
      for (unsigned int v=0; v<this->n_vars(s); v++)
	_dof_ids[s][v] = invalid_id - 1;
    }

#endif
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


  
#ifdef ENABLE_EXPENSIVE_DATA_STRUCTURES

  assert (_n_comp != NULL);
  assert (_n_comp[s] != NULL);
    
  // If we already have memory allocated clear it.
  if (this->n_comp(s,var))
    {
      assert (_dof_ids[s][var] != NULL); delete [] _dof_ids[s][var]; _dof_ids[s][var] = NULL;
    }
  
  _n_comp[s][var]  = static_cast<unsigned char>(ncomp);

  if (this->n_comp(s,var) > 0)
    {
      _dof_ids[s][var] = new unsigned int [ncomp];
      
      for (unsigned int c=0; c<ncomp; c++)
	_dof_ids[s][var][c] = invalid_id;
    }

#else

  // Remeber...  use (invalid_id - 1) to signify no
  // components for this object
  if (ncomp == 0)
    {
      _dof_ids[s][var] = (invalid_id - 1);
    }
  else if (ncomp == 1)
    {      
      _dof_ids[s][var] = invalid_id;
    }
  else
    {
      std::cerr << "ERROR: You must compile with --enable-expensive to have" << std::endl
		<< "multiple components in a variable!" << std::endl;
      
      error();
    }
  
#endif
}



void DofObject::set_dof_number(const unsigned int s,
			       const unsigned int var,
			       const unsigned int comp,
			       const unsigned int dn)
{
  //assert (dn != invalid_id);
  assert (s < this->n_systems());
  assert (var  < this->n_vars(s));
  assert (comp < this->n_comp(s,var));
  assert (_dof_ids != NULL);
  assert (_dof_ids[s] != NULL);
    
#ifdef ENABLE_EXPENSIVE_DATA_STRUCTURES
  
  assert (_dof_ids[s][var] != NULL);
  
  _dof_ids[s][var][comp] = dn;

#else

  // Reserve (invalid_id-1) to signify n_com(var) = 0.
  assert (dn != (invalid_id-1));
  
  _dof_ids[s][var] = dn;

#endif
}
