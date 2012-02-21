// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


namespace libMesh
{



// ------------------------------------------------------------
// DofObject class static member -now initialized in header
const unsigned int DofObject::invalid_id;
const processor_id_type DofObject::invalid_processor_id;



// ------------------------------------------------------------
// DofObject class members
// Copy Constructor
DofObject::DofObject (const DofObject& dof_obj) :
  ReferenceCountedObject<DofObject>(),
#ifdef LIBMESH_ENABLE_AMR
  old_dof_object (NULL),
#endif
  _id            (dof_obj._id),
  _processor_id  (dof_obj._processor_id),
  _idx_buf       (dof_obj._idx_buf)
{

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


// Deep-copying assignment operator
DofObject& DofObject::operator= (const DofObject& dof_obj)
{
#ifdef LIBMESH_ENABLE_AMR
  this->clear_old_dof_object();

  this->old_dof_object = new DofObject(*(dof_obj.old_dof_object));
#endif

  _id           = dof_obj._id;
  _processor_id = dof_obj._processor_id;
  _idx_buf      = dof_obj._idx_buf;


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

  return *this;
}





#ifdef LIBMESH_ENABLE_AMR

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
  if (ns == this->n_systems())
    return;

  // Clear any existing data.  This is safe to call
  // even if we don't have any data.
  this->clear_dofs();

  // Set the new number of systems
  _idx_buf.resize(ns, ns);
  _idx_buf[0] = ns;


#ifdef DEBUG

  // check that all systems now exist and that they have 0 size
  libmesh_assert (ns == this->n_systems());
  for (unsigned int s=0; s<this->n_systems(); s++)
    libmesh_assert (this->n_vars(s) == 0);

#endif
}



void DofObject::add_system()
{
  // quick return?
  if (this->n_systems() == 0)
    {
      this->set_n_systems(1);
      return;
    }

  DofObject::index_buffer_t::iterator it = _idx_buf.begin();

  std::advance(it, this->n_systems());

  // this inserts the current vector size at the position for the new system - creating the
  // entry we need for the new system indicating there are 0 variables.
  _idx_buf.insert(it, _idx_buf.size());

  // cache this value before we screw it up!
  const unsigned int ns_orig = this->n_systems();

  // incriment the number of systems and the offsets for each of
  // the systems including the new one we just added.
  for (unsigned int i=0; i<ns_orig+1; i++)
    _idx_buf[i]++;

  libmesh_assert (this->n_systems() == (ns_orig+1));
  libmesh_assert (this->n_vars(ns_orig) == 0);
}



void DofObject::set_n_vars(const unsigned int s,
			   const unsigned int nvars)
{
  libmesh_assert (s < this->n_systems());

  // BSK - note that for compatibility with the previous implementation
  // calling this method when (nvars == this->n_vars()) requires that
  // we invalidate the DOF indices and set the number of components to 0.
  // Note this was a bit of a suprise to me - there was no quick return in
  // the old method, which caused removal and readdition of the DOF indices
  // even in the case of (nvars == this->n_vars()), resulting in n_comp(s,v)
  // implicitly becoming 0 regardless of any previous value.
  // quick return?
  if (nvars == this->n_vars(s))
    {
      for (unsigned int v=0; v<nvars; v++)
	this->set_n_comp(s,v,0);
      return;
    }

  // since there is ample opportunity to screw up other systems, let us
  // cache their current sizes and later assert that they are unchanged.
#ifdef DEBUG
  DofObject::index_buffer_t old_system_sizes;
  old_system_sizes.reserve(this->n_systems());

  for (unsigned int s_ctr=0; s_ctr<this->n_systems(); s_ctr++)
    old_system_sizes.push_back(this->n_vars(s_ctr));
#endif

  // remove current indices if we have some
  if (this->n_vars(s) != 0)
    {
      const unsigned int old_nvars_s = this->n_vars(s);

      DofObject::index_buffer_t::iterator
	it  = _idx_buf.begin(),
	end = _idx_buf.begin();

      std::advance(it,  this->start_idx(s));
      std::advance(end, this->end_idx(s));
      _idx_buf.erase(it,end);

      for (unsigned int ctr=(s+1); ctr<this->n_systems(); ctr++)
	_idx_buf[ctr] -= 2*old_nvars_s;
    }

  // better not have any now!
  libmesh_assert (this->n_vars(s) == 0);

  // had better not screwed up any of our sizes!
#ifdef DEBUG
  for (unsigned int s_ctr=0; s_ctr<this->n_systems(); s_ctr++)
    if (s_ctr != s)
      libmesh_assert(this->n_vars(s_ctr) == old_system_sizes[s_ctr]);
#endif

  // OK, if the user requested 0 that is what we have
  if (nvars == 0)
    return;

  {
    // array to hold new indices
    DofObject::index_buffer_t var_idxs(2*nvars);
    for (unsigned int v=0; v<nvars; v++)
      {
	var_idxs[2*v    ] = 0;
	var_idxs[2*v + 1] = invalid_id - 1;
      }

    DofObject::index_buffer_t::iterator it = _idx_buf.begin();
    std::advance(it, this->end_idx(s));
    _idx_buf.insert(it, var_idxs.begin(), var_idxs.end());

    for (unsigned int ctr=(s+1); ctr<this->n_systems(); ctr++)
      _idx_buf[ctr] += 2*nvars;

    // resize _idx_buf to fit so no memory is wasted.
    DofObject::index_buffer_t(_idx_buf).swap(_idx_buf);
  }

  // that better had worked.  Assert stuff.
  libmesh_assert (nvars == this->n_vars(s));

#ifdef DEBUG
  for (unsigned int v=0; v<this->n_vars(s); v++)
    libmesh_assert (this->n_comp(s,v) == 0);
  // again, all other system sizes shoudl be unchanged!
  for (unsigned int s_ctr=0; s_ctr<this->n_systems(); s_ctr++)
    if (s_ctr != s)
      libmesh_assert(this->n_vars(s_ctr) == old_system_sizes[s_ctr]);

//   std::cout << " [ ";
//   for (unsigned int i=0; i<_idx_buf.size(); i++)
//     std::cout << _idx_buf[i] << " ";
//   std::cout << "]\n";

#endif
}



void DofObject::set_n_comp(const unsigned int s,
			   const unsigned int var,
			   const unsigned int ncomp)
{
  libmesh_assert (s   < this->n_systems());
  libmesh_assert (var < this->n_vars(s));

  // Check for trivial return
  if (ncomp == this->n_comp(s,var)) return;

  const unsigned int
    start_idx_sys = this->start_idx(s),
    base_offset  = start_idx_sys + 2*var;

  libmesh_assert ((base_offset + 1) < _idx_buf.size());

  // set the number of components
  _idx_buf[base_offset] = ncomp;

  // We use (invalid_id - 1) to signify no
  // components for this object
  _idx_buf[base_offset + 1] = (ncomp == 0) ? invalid_id - 1 : invalid_id;

  libmesh_assert (ncomp == this->n_comp(s,var));
}



void DofObject::set_dof_number(const unsigned int s,
			       const unsigned int var,
			       const unsigned int comp,
			       const unsigned int dn)
{
  libmesh_assert (s < this->n_systems());
  libmesh_assert (var  < this->n_vars(s));
  libmesh_assert (comp < this->n_comp(s,var));

  const unsigned int
    start_idx_sys = this->start_idx(s);

  libmesh_assert ((start_idx_sys + 2*var + 1) < _idx_buf.size());

  unsigned int
    &base_idx(_idx_buf[start_idx_sys + 2*var + 1]);

  //We intend to change all dof numbers together or not at all
  if (comp)
    libmesh_assert ((dn == invalid_id && base_idx == invalid_id) ||
		    (dn == base_idx + comp));
  else
    base_idx = dn;

// #ifdef DEBUG
//   std::cout << " [ ";
//   for (unsigned int i=0; i<_idx_buf.size(); i++)
//     std::cout << _idx_buf[i] << " ";
//   std::cout << "]\n";

// #endif

  libmesh_assert(this->dof_number(s, var, comp) == dn);
}

} // namespace libMesh
