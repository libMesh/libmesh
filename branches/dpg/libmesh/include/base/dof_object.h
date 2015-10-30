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



#ifndef __dof_object_h__
#define __dof_object_h__

// C++ includes
#include <vector>

// Local includes
#include "libmesh_config.h"
#include "libmesh_common.h"
#include "libmesh.h" // libMesh::invalid_uint
#include "reference_counted_object.h"

namespace libMesh
{

// Forward declarations
class DofObject;


// Define processor id storage type.  We default to short to save
// space, but expanding to support more than 2^16-2 procs should work
// too.
typedef unsigned short int processor_id_type;
//typedef unsigned int processor_id_type;

/**
 * The \p DofObject defines an abstract base class for objects that
 * have degrees of freedom associated with them.  Examples of such
 * objects are the \p Node and \p Elem classes.  This class can
 * not be instantiated, only derived from.
 *
 * \author Benjamin S. Kirk
 * \date 2003, 2011
 * \version $Revision$
 */

class DofObject : public ReferenceCountedObject<DofObject>
{
  
protected:
  
  /**
   * Constructor. Protected so that you can't instantiate one of these.
   */
  DofObject ();

public:
  
  /**
   * Copy-constructor.
   */
  DofObject (const DofObject&);

  /**
   * Destructor.
   */ 
  virtual ~DofObject ();
  
  /**
   * Deep-copying assignment operator
   */ 
  DofObject& operator= (const DofObject& dof_obj);

#ifdef LIBMESH_ENABLE_AMR

  /**
   * This object on the last mesh.  Useful for projecting
   * solutions from one mesh to another.
   */
  DofObject* old_dof_object;

  /**
   * Sets the \p old_dof_object to NULL
   */
  void clear_old_dof_object ();

  /**
   * Sets the \p old_dof_object to a copy of \p this
   */
  void set_old_dof_object ();

#endif
  
  /**
   * Clear the \p DofMap data structures and return to
   * a pristine state.
   */
  void clear_dofs ();

  /**
   * Sets all degree of freedom numbers to \p invalid_id
   */
  void invalidate_dofs (const unsigned int sys_num = libMesh::invalid_uint);

  /**
   * Sets the id to \p invalid_id
   */
  void invalidate_id ();

  /**
   * Sets the processor id to \p invalid_processor_id
   */
  void invalidate_processor_id ();

  /**
   * Invalidates all the indices for this \p DofObject
   */
  void invalidate ();
  
  /**
   * @returns the number of degrees of freedom associated with
   * system \p s for this object. Optionally only counts degrees
   * of freedom for variable number \p var
   */
  unsigned int n_dofs (const unsigned int s, 
		       const unsigned int var =
		       libMesh::invalid_uint) const;
  
  /**
   * \returns the \p id for this \p DofObject
   */
  unsigned int id () const;

  /**
   * \returns the \p id for this \p DofObject as a writeable reference.
   */
  unsigned int & set_id ();

  /**
   * Sets the \p id for this \p DofObject
   */
  void set_id (const unsigned int id)
  { this->set_id() = id; }

  /**
   * @returns \p true if this \p DofObject has a valid \p id set,
   * \p false otherwise.
   */
  bool valid_id () const;
  
  /**
   * @returns the processor that this element belongs to.
   */
  processor_id_type processor_id () const;
  
  /**
   * @returns the processor that this element belongs to as a
   * writeable reference.
   */
  processor_id_type& processor_id ();

  /**
   * Sets the \p processor_id for this \p DofObject.
   */  
  void processor_id (const processor_id_type id);

  /**
   * @returns \p true if this \p DofObject has a valid \p id set,
   * \p false otherwise.
   */
  bool valid_processor_id () const;
  
  /**
   * @returns the number of systems associated with this
   * \p DofObject
   */
  unsigned int n_systems() const;

  /**
   *  Sets the number of systems for this \p DofObject
   */
  void set_n_systems (const unsigned int s);

  /**
   * Adds an additional system to the \p DofObject
   */
  void add_system ();
  
  /**
   * @returns the number of variables associated with system \p s
   * for this \p DofObject
   */
  unsigned int n_vars(const unsigned int s) const;
  
  /**
   * Sets number of variables associated with system \p s for this
   * \p DofObject.  Has the effect of setting the number of components
   * to 0 even when called even with (nvars == this->n_vars(s)).
   */
  void set_n_vars(const unsigned int s,
		  const unsigned int nvars);
  
  /**
   * @returns the number of components for variable \p var
   * of system \p s associated with this \p DofObject. 
   * For example, the \p HIERARCHIC shape functions may
   * have @e multiple dof's associated with @e one node.  Another
   * example is the \p MONOMIALs, where only the elements
   * hold the dof's, but for the different spatial directions,
   * and orders, see \p FE.
   */
  unsigned int n_comp(const unsigned int s,
		      const unsigned int var) const;
  
  /**
   * Sets the number of components for variable \p var
   * of system \p s associated with this \p DofObject
   */
  void set_n_comp(const unsigned int s,
		  const unsigned int var,
		  const unsigned int ncomp);
  
  /**
   * @returns the global degree of freedom number variable \p var,
   * component \p comp for system \p s associated with this \p DofObject
   */
  unsigned int dof_number(const unsigned int s,
			  const unsigned int var,
			  const unsigned int comp) const;
  
  /**
   * Sets the global degree of freedom number variable \p var,
   * component \p comp for system \p s associated with this \p DofObject
   */
  void set_dof_number(const unsigned int s,
		      const unsigned int var,
		      const unsigned int comp,
		      const unsigned int dn);

  /**
   * @returns true if any system has variables which have been assigned,
   * false otherwise
   */
  bool has_dofs(const unsigned int s=libMesh::invalid_uint) const;
  
  /**
   * Implemented in Elem and Node.
   */
//  virtual bool operator==(const DofObject& ) const
//  { libmesh_error(); return false; }

  
  /**
   * An invaild \p id to distinguish an uninitialized \p DofObject
   */
  static const unsigned int invalid_id = libMesh::invalid_uint;

  /**
   * An invalid \p processor_id to distinguish DOFs that have
   * not been assigned to a processor.
   */
  static const processor_id_type invalid_processor_id = static_cast<processor_id_type>(-1);

  
private:

  
  /**
   * The \p id of the \p DofObject
   */
  unsigned int _id;

  /**
   * The \p processor_id of the \p DofObject.
   * Degrees of freedom are wholly owned by processors,
   * however they may be duplicated on other processors.
   *
   * This is stored as an unsigned short int since we cannot
   * expect to be solving on 65000+ processors any time soon,
   * can we??
   */
  processor_id_type _processor_id;

  /**
   * DOF index information.  This is packed into a contiguous buffer of the following format:
   *
   * [ns end_0 end_1 ... end_{ns-1} (nc_0 idx_0 nc_1 idx_1 ... nc_nv idx_nv)_0 (nc_0 idx_0 nc_1 idx_1 ... nc_nv idx_nv)_1 ... (nc_0 idx_0 nc_1 idx_1 ... nc_nv idx_nv)_ns ]
   *
   * where 'end_s' is the index past the end of the variable storage for system \p s.
   * Note that we specifically do not store the end for the last system - this always _idx_buf.size().
   *
   * Specifically, consider the case of 4 systems, with 3, 0, 1, 2 DOFs respectively.  The _idx_buf then looks like
   *
   * [4 10 10 12 () (nc_0 idx_0 nc_1 idx_1 nc_2 idx_2) () (nc_0 idx_0) (nc_0 idx_0 nc_1 idx_1)]
   * [0  1  2  3        4     5    6     7    8     9        10    11     12    13   14    15]
   *
   * The ending index is then given by
   *
   * end_s = _idx_buf.size(), s == (ns-1),
   *       = _idx_buf[s+1]    otherwise.
   *
   * The starting indices are not specifically stored, but rather inferred as follows:
   *
   * start_s = _idx_buf[s];
   */
  typedef std::vector<unsigned int> index_buffer_t;
  index_buffer_t _idx_buf;

  /**
   * The starting index for system \p s.
   */
  unsigned int start_idx(const unsigned int s) const;

  /**
   * The ending index for system \p s.
   */
  unsigned int end_idx(const unsigned int s) const;
};



//------------------------------------------------------
// Inline functions
inline
DofObject::DofObject () :
#ifdef LIBMESH_ENABLE_AMR
  old_dof_object(NULL),
#endif
  _id (invalid_id),
  _processor_id (invalid_processor_id)
{
  this->invalidate();
}





inline
DofObject::~DofObject ()
{
  // Free all memory.
#ifdef LIBMESH_ENABLE_AMR
  this->clear_old_dof_object ();
#endif
  this->clear_dofs ();
}



inline
void DofObject::invalidate_dofs (const unsigned int sys_num)
{
  // If the user does not specify the system number...
  if (sys_num >= this->n_systems()) 
    {
      for (unsigned int s=0; s<this->n_systems(); s++)
        for (unsigned int v=0; v<this->n_vars(s); v++)
	  if (this->n_comp(s,v))
	    this->set_dof_number(s,v,0,invalid_id);
    }
  // ...otherwise invalidate the dofs for all systems
  else
    for (unsigned int v=0; v<this->n_vars(sys_num); v++)
      if (this->n_comp(sys_num,v))
        this->set_dof_number(sys_num,v,0,invalid_id);
}



inline
void DofObject::invalidate_id ()
{
  this->set_id (invalid_id);
}



inline
void DofObject::invalidate_processor_id ()
{
  this->processor_id (invalid_processor_id);
}



inline
void DofObject::invalidate ()
{
  this->invalidate_dofs ();
  this->invalidate_id ();
  this->invalidate_processor_id ();
}



inline
void DofObject::clear_dofs ()
{
  // vector swap trick to force deallocation
  index_buffer_t().swap(_idx_buf);

  libmesh_assert (this->n_systems() == 0);
  libmesh_assert (_idx_buf.empty());
}
  


inline
unsigned int DofObject::n_dofs (const unsigned int s,
				const unsigned int var) const
{
  libmesh_assert (s < this->n_systems());
  
  unsigned int num = 0;

  // Count all variables
  if (var == libMesh::invalid_uint)
    for (unsigned int v=0; v<this->n_vars(s); v++)
      num += this->n_comp(s,v);
  
  // Only count specified variable
  else
    num = this->n_comp(s,var);

  return num;
}



inline
unsigned int DofObject::id () const
{
  libmesh_assert (this->valid_id());
  return _id;
}



inline
unsigned int & DofObject::set_id ()
{
  return _id;
}



inline
bool DofObject::valid_id () const
{
  return (DofObject::invalid_id != _id);
}


inline
processor_id_type DofObject::processor_id () const
{
  return _processor_id;
}



inline
processor_id_type& DofObject::processor_id ()
{
  return _processor_id;
}



inline
void DofObject::processor_id (const processor_id_type id)
{
  this->processor_id() = id;
}



inline
bool DofObject::valid_processor_id () const
{
  return (DofObject::invalid_processor_id != _processor_id);
}



inline
unsigned int DofObject::n_systems () const
{
  return _idx_buf.empty() ? 0 : _idx_buf[0];
}



inline
unsigned int DofObject::n_vars(const unsigned int s) const
{
  libmesh_assert (s < this->n_systems());

  return (this->end_idx(s) - this->start_idx(s)) / 2;
}




inline
unsigned int DofObject::n_comp(const unsigned int s,
			       const unsigned int var) const
{
  libmesh_assert (s   < this->n_systems());
  libmesh_assert (var < this->n_vars(s));

# ifdef DEBUG
  // Does this ever happen?  I doubt it... 3/7/2003 (BSK)
  if (var >= this->n_vars(s))
    {
      libMesh::err << "s=" << s << ", var=" << var << std::endl
		    << "this->n_vars(s)=" << this->n_vars(s) << std::endl
		    << "this->n_systems()=" << this->n_systems() << std::endl;
      libmesh_error();
    }
# endif

  const unsigned int
    start_idx_sys = this->start_idx(s);

  libmesh_assert ((start_idx_sys + 2*var) < _idx_buf.size());
  
  return _idx_buf[start_idx_sys + 2*var];
}



inline
unsigned int DofObject::dof_number(const unsigned int s,
				   const unsigned int var,
				   const unsigned int comp) const
{
  libmesh_assert (s    < this->n_systems());
  libmesh_assert (var  < this->n_vars(s));
  libmesh_assert (comp < this->n_comp(s,var));

  const unsigned int
    start_idx_sys = this->start_idx(s);

  libmesh_assert ((start_idx_sys + 2*var + 1) < _idx_buf.size());

  const unsigned int
    base_idx = _idx_buf[start_idx_sys + 2*var + 1];

  // if the first component is invalid, they
  // are all invalid
  if (base_idx == invalid_id)
    return invalid_id;

  // otherwise the index is the first component
  // index augemented by the component number
  else
    return (base_idx + comp);
}



inline
bool DofObject::has_dofs (const unsigned int sys) const
{
  if (sys == libMesh::invalid_uint)
    {
      for (unsigned int s=0; s<this->n_systems(); s++)
	if (this->n_vars(s))
	  return true;
    }
  
  else
    {
      libmesh_assert (sys < this->n_systems());

      if (this->n_vars(sys))
	return true;
    }
  
  return false;
}


  
inline
unsigned int DofObject::start_idx (const unsigned int s) const
{
  libmesh_assert (s < this->n_systems());
  libmesh_assert (s < _idx_buf.size());

  return _idx_buf[s];
}

  

inline
unsigned int DofObject::end_idx (const unsigned int s) const
{
  libmesh_assert (s < this->n_systems());
  libmesh_assert (s < _idx_buf.size());

  return ((s+1) == this->n_systems()) ? _idx_buf.size() : _idx_buf[s+1];
}


} // namespace libMesh


#endif // #ifndef __dof_object_h__
