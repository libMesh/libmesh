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
  _n_v_comp      (NULL),
  _dof_ids       (NULL)
{

  // Allocate storage for the dof numbers and copy
  // the values. 
  // IT IS UNDEFINED BEHAVIOR TO ALLOCATE AN ARRAY WITH ZERO ENTRIES,
  // IF n_systems==0, leave _n_v_comp, and _dof_ids NULL.
  if (this->n_systems() > 0)
    {
      _n_v_comp = new unsigned char* [this->n_systems()];
      _dof_ids  = new unsigned int*  [this->n_systems()];

      // gotta specifically NULL these - we rely later that
      // _n_v_comp[s] == NULL is synonymous with no variables in the system.
      for (unsigned int s=0; s<this->n_systems(); s++)
	{
	  _n_v_comp[s] = NULL;
	  _dof_ids[s]  = NULL;
	}      
    }

  // If n_systems==0, we don't enter this for loop.
  for (unsigned int s=0; s<this->n_systems(); s++)
    {
      // In case you have a system with no variables, it is undefined
      // behavior (UB) to allocate a zero-length array here.
      if (dof_obj.n_vars(s) > 0)
	{
	  _n_v_comp[s] = new unsigned char [dof_obj.n_vars(s)+1]; 
	  _dof_ids[s]  = new unsigned int  [dof_obj.n_vars(s)];
	  
	  _n_v_comp[s][0] = dof_obj.n_vars(s);

	}
      for (unsigned int v=0; v<this->n_vars(s); v++)
	{
	  _n_v_comp[s][v+1]  = dof_obj.n_comp(s,v);

	  if (this->n_comp(s,v) > 0)
	    _dof_ids[s][v] = dof_obj.dof_number(s,v,0);
          else
	    _dof_ids[s][v] = invalid_id;
	}
    }

  // Check that everything worked
#ifdef DEBUG

  libmesh_assert (this->n_systems() == dof_obj.n_systems());

  for (unsigned int s=0; s<this->n_systems(); s++)
    {
      libmesh_assert (this->n_vars(s) == dof_obj.n_vars(s));

      for (unsigned int v=0; v<this->n_vars(s); v++)
	{
	  libmesh_assert (this->n_comp(s,v) == dof_obj.n_comp(s,v));

	  for (unsigned int c=0; c<this->n_comp(s,v); c++)
	    libmesh_assert (this->dof_number(s,v,c) == dof_obj.dof_number(s,v,c));
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

  libmesh_assert (this->old_dof_object == NULL);
  
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
      
      libmesh_error();
    }
					
#endif


  // Clear any existing data.  This is safe to call
  // even if we don't have any data.
  this->clear_dofs();

  // Set the new number of systems
  _n_systems = static_cast<unsigned char>(ns);
  
  // Allocate storage for the systems
  _n_v_comp = new unsigned char* [this->n_systems()];
  _dof_ids  = new unsigned int*  [this->n_systems()];

  // No variables have been declared yet.
  for (unsigned int s=0; s<this->n_systems(); s++)
    {
      _n_v_comp[s] = NULL;
      _dof_ids[s]  = NULL;
    }
}



void DofObject::add_system()
{
  if (this->n_systems() > 0)
    {
      // Copy the old systems to temporary storage
      unsigned char **old_n_v_comp = new unsigned char* [this->n_systems()];
      unsigned int  **old_dof_ids  = new unsigned int*  [this->n_systems()]; 
      
      for (unsigned int s=0; s<this->n_systems(); s++)
	{
	  old_n_v_comp[s] = _n_v_comp[s];
	  old_dof_ids[s]  = _dof_ids[s];
	}
      
      // Delete old storage
      libmesh_assert (_n_v_comp != NULL); delete [] _n_v_comp; _n_v_comp = NULL;
      libmesh_assert (_dof_ids  != NULL); delete [] _dof_ids;  _dof_ids  = NULL;
  
      // Allocate space for new system
      _n_v_comp= new unsigned char* [this->n_systems()+1];
      _dof_ids = new unsigned int*  [this->n_systems()+1];
      
      // Copy the other systems
      for (unsigned int s=0; s<this->n_systems(); s++)
	{
	  _n_v_comp[s] = old_n_v_comp[s];
	  _dof_ids[s]  = old_dof_ids[s];
	}
	       
      // Delete temporary storage
      libmesh_assert (old_n_v_comp != NULL); delete [] old_n_v_comp; old_n_v_comp = NULL;
      libmesh_assert (old_dof_ids  != NULL); delete [] old_dof_ids;  old_dof_ids  = NULL;
    }
  else
    {
      libmesh_assert (_n_v_comp == NULL);
      libmesh_assert (_dof_ids  == NULL);
      
      // Allocate space for new system
      _n_v_comp = new unsigned char* [this->n_systems()+1];
      _dof_ids  = new unsigned int*  [this->n_systems()+1];      
    }
  
  // Initialize the new system
  _n_v_comp[this->n_systems()] = NULL;
  _dof_ids[this->n_systems()]  = NULL;
  
  // Done. Don't forget to increment the number of systems!
  _n_systems++;
}



void DofObject::set_n_vars(const unsigned int s,
			   const unsigned int nvars)
{
  libmesh_assert (s < this->n_systems());
  libmesh_assert (_n_v_comp != NULL);
  libmesh_assert (_dof_ids  != NULL);

#ifdef DEBUG

  if (nvars != static_cast<unsigned int>(static_cast<unsigned char>(nvars)))
    {
      std::cerr << "Unsigned char not big enough to hold nvar!" << std::endl
		<< "Recompile with _n_vars set to a bigger type!"
		<< std::endl;
      
      libmesh_error();
    }
					
#endif

  
  
  // If we already have memory allocated clear it.
  if (this->n_vars(s) != 0)
    {
      libmesh_assert (_n_v_comp[s] != NULL); delete [] _n_v_comp[s]; _n_v_comp[s] = NULL;
      libmesh_assert (_dof_ids[s]  != NULL); delete [] _dof_ids[s];  _dof_ids[s]  = NULL;
    }

  // Reset the number of variables in the system  
  if (nvars > 0)
    {
      libmesh_assert (_n_v_comp[s] == NULL);
      libmesh_assert (_dof_ids[s]  == NULL);
      
      _n_v_comp[s] = new unsigned char [nvars+1];
      _dof_ids[s]  = new unsigned int  [nvars];
      
      _n_v_comp[s][0] = static_cast<unsigned char>(nvars);

      libmesh_assert (nvars == this->n_vars(s));
      
      for (unsigned int v=0; v<this->n_vars(s); v++)
	{
	  _n_v_comp[s][v+1] = 0;
	  _dof_ids[s][v]    = invalid_id - 1;
	}
    }
  else // (nvars == 0)
    {
      libmesh_assert (_n_v_comp[s] == NULL);
      libmesh_assert (_dof_ids[s]  == NULL);
    }
}



void DofObject::set_n_comp(const unsigned int s,
			   const unsigned int var,
			   const unsigned int ncomp)
{
  libmesh_assert (s < this->n_systems());
  libmesh_assert (var < this->n_vars(s));
  libmesh_assert (_dof_ids != NULL);
  libmesh_assert (_dof_ids[s] != NULL);
  
  // Check for trivial return
  if (ncomp == this->n_comp(s,var)) return;

#ifdef DEBUG

  if (ncomp != static_cast<unsigned int>(static_cast<unsigned char>(ncomp)))
    {
      std::cerr << "Unsigned char not big enough to hold ncomp!" << std::endl
		<< "Recompile with _n_v_comp set to a bigger type!"
		<< std::endl;
      
      libmesh_error();
    }
  
#endif

  // We use (invalid_id - 1) to signify no
  // components for this object
  if (ncomp == 0)
    {
      _dof_ids[s][var] = (invalid_id - 1);
    }
  
  libmesh_assert (_n_v_comp    != NULL);
  libmesh_assert (_n_v_comp[s] != NULL);
    
  _n_v_comp[s][var+1]  = static_cast<unsigned char>(ncomp);
}



void DofObject::set_dof_number(const unsigned int s,
			       const unsigned int var,
			       const unsigned int comp,
			       const unsigned int dn)
{
  libmesh_assert (s < this->n_systems());
  libmesh_assert (var  < this->n_vars(s));
  libmesh_assert (_dof_ids != NULL);
  libmesh_assert (_dof_ids[s] != NULL);
  libmesh_assert (comp < this->n_comp(s,var));
  
  //We intend to change all dof numbers together or not at all
  if (comp)
    libmesh_assert ((dn == invalid_id && _dof_ids[s][var] == invalid_id) || 
		    (dn == _dof_ids[s][var] + comp));
  else
    _dof_ids[s][var] = dn;


  libmesh_assert(this->dof_number(s, var, comp) == dn);
}
