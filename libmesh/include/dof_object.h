// $Id: dof_object.h,v 1.1 2003-02-13 01:49:48 benkirk Exp $

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



#ifndef __dof_object_h__
#define __dof_object_h__

// C++ includes

// Local includes
#include "mesh_common.h"




/**
 * The \p DofObject defines an abstract base class for objects that
 * have degrees of freedom associated with them.  Examples of such
 * objects are the \p Node and \p Elem classes.  This class can
 * not be instantiated, only derived from.
 *
 * This class is intended to be extremely lightweight.  To help acheive
 * this goal no \p std::vector<> or anything else that might be heavy
 * is used to store the degree of freedom indices.  If the library is
 * configured with \p --enable-expensive then there can be multiple
 * components in any variable, if the library is configured with
 * \p --disable-expensive then there can be at MOST one component
 * per variable, i.e n_comp(var) <= 1.
 *
 * \author Benjamin S. Kirk
 * \date 2003
 * \version $Revision: 1.1 $
 */

class DofObject
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
   * Clear the \p DofMap data structures and return to
   * a pristine state.
   */
  void clear_dofs ();

  /**
   * Sets all degree of freedom numbers to \p invalid_id
   */
  void invalidate_dofs ();

  /**
   * @returns the number of degrees of freedom associated with
   * this object.  Optionally only counts degrees of freedom for
   * variable number \p var
   */
  unsigned int n_dofs (const unsigned int var =
		       static_cast<unsigned int>(-1)) const;
  
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
  { set_id() = id; };

  /**
   * @returns the number of variables associated with this
   * \p DofObject
   */
  unsigned int n_vars() const;
  
  /**
   * @returns the number of variables associated with this
   * \p DofObject as a writeable reference.
   */
  void set_n_vars(const unsigned int nvars);
  
  /**
   * @returns the number of components for variable \p var
   * associated with this \p DofObject
   */
  unsigned int n_comp(const unsigned int var) const;
  
  /**
   * @returns the number of components for variable \p var
   * associated with this \p DofObject as a writeable reference
   */
  void set_n_comp(const unsigned int var,
		  const unsigned int ncomp);
  
  /**
   * @returns the global degree of freedom number variable \p var,
   * component \p comp associated with this \p DofObject
   */
  unsigned int dof_number(const unsigned int var,
			  const unsigned int comp) const;
  
  /**
   * @returns the global degree of freedom number variable \p var,
   * component \p comp associated with this \p DofObject as a
   * writeable reference.
   */
  void set_dof_number(const unsigned int var,
		      const unsigned int comp,
		      const unsigned int dn);
  
  /**
   * An invaild \p id to distinguish an uninitialized \p Node
   */
  static const unsigned int invalid_id;

  
private:

  
  /**
   * The \p id of the \p Node
   */
  unsigned int _id;

  /**
   * The number of variables associated with this
   * \p DofObject.  This is stored as an unsigned char
   * for storage efficiency.
   */
  unsigned char _n_vars;

#ifdef ENABLE_EXPENSIVE_DATA_STRUCTURES
  
  /**
   * The number of components for each variable associated
   * with this \p DofObject.  This is stored as an unsigned char
   * for storage efficiency.
   */
  unsigned char *_n_comp;

  /**
   * The global degree of freedom numbers for each component
   * of each variable.
   */
  unsigned int **_dof_ids;

#else

  /**
   * The global degree of freedom numbers for each variable.
   * (Only one component allowed.)
   */
  unsigned int *_dof_ids;

#endif

};



//------------------------------------------------------
// Inline functions
inline
DofObject::DofObject () :
  _id (invalid_id),
  _n_vars (0),
#ifdef ENABLE_EXPENSIVE_DATA_STRUCTURES
  _n_comp (NULL),
#endif
  _dof_ids (NULL)
{
};



inline
DofObject::DofObject (const DofObject& dof_obj) :
  _id (dof_obj._id),
  _n_vars (dof_obj._n_vars),
#ifdef ENABLE_EXPENSIVE_DATA_STRUCTURES
  _n_comp (NULL),
#endif
  _dof_ids (NULL)
{
  // Allocate storage for the dof numbers and copy
  // the values.
  
#ifdef ENABLE_EXPENSIVE_DATA_STRUCTURES
  
  _n_comp  = new unsigned char [n_vars()];
  _dof_ids = new unsigned int* [n_vars()];
  
  for (unsigned int v=0; v<n_vars(); v++)
    {
      _n_comp[v]  = dof_obj.n_comp(v);

      _dof_ids[v] = new unsigned int [dof_obj.n_comp(v)];

      for (unsigned int c=0; c<n_comp(v); c++)
	set_dof_number(v,c,dof_obj.dof_number(v,c));
    };

#else
  
  _dof_ids = new unsigned int [n_vars()];
  
  for (unsigned int v=0; v<n_vars(); v++)
    _dof_ids[v] = dof_obj._dof_ids[v];
  
#endif
};



inline
DofObject::~DofObject ()
{
  clear_dofs ();
};



inline
void DofObject::clear_dofs ()
{
#ifdef ENABLE_EXPENSIVE_DATA_STRUCTURES
  
  for (unsigned int v=0; v<n_vars(); v++)
    {
      delete [] _dof_ids[v]; _dof_ids[v] = NULL;
    };

  delete [] _dof_ids; _dof_ids = NULL;
  delete [] _n_comp;  _n_comp  = NULL;

#else

  delete [] _dof_ids; _dof_ids = NULL;
  
#endif

  _n_vars = 0;
  _id = invalid_id;
 };



inline
void DofObject::invalidate_dofs ()
{
  for (unsigned int v=0; v<n_vars(); v++)
    for (unsigned int c=0; c<n_comp(v); c++)
      set_dof_number(v,c,invalid_id);
};



inline
unsigned int DofObject::n_dofs (const unsigned int var) const
{
  unsigned int num = 0;

  // Count all variables
  if (var == static_cast<unsigned int>(-1))
    for (unsigned int v=0; v<n_vars(); v++)
      for (unsigned int c=0; c<n_comp(v); c++)
	num++;
  
  // Only count specified variable
  else
    {
      assert (var < n_vars());

      for (unsigned int c=0; c<n_comp(var); c++)
	num++;	
    };

  return num;
};



inline
unsigned int DofObject::id () const
{
  assert (_id != invalid_id);

  return _id;
};



inline
unsigned int & DofObject::set_id ()
{
  return _id;
};



inline
unsigned int DofObject::n_vars() const
{
  return static_cast<unsigned int>(_n_vars);
};



inline
void DofObject::set_n_vars(const unsigned int nvars)
{

//#if !defined(NDEBUG)

  if (nvars != static_cast<unsigned int>(static_cast<unsigned char>(nvars)))
    {
      std::cerr << "Unsigned char not big enough to hold nvar!" << std::endl
		<< "Recompile with _n_vars set to a bigger type!"
		<< std::endl;
      
      error();
    };
					
//#endif

  
#ifdef ENABLE_EXPENSIVE_DATA_STRUCTURES
  
  // If we already have memory allocated clear it.
  if (n_vars())
    {
      for (unsigned int v=0; v<n_vars(); v++)
	if (n_comp(v) != 0)
	  {
	    delete [] _dof_ids[v]; _dof_ids[v] = NULL;
	    
	    _n_comp[v] = 0;
	  };
      
      delete [] _n_comp;  _n_comp  = NULL;
      delete [] _dof_ids; _dof_ids = NULL;
    };
  
  _n_vars = static_cast<unsigned char>(nvars);

  _n_comp  = new unsigned char [n_vars()];
  _dof_ids = new unsigned int* [n_vars()];

  for (unsigned int v=0; v<n_vars(); v++)
    {
      _n_comp[v]  = 0;
      _dof_ids[v] = NULL;
    };

#else

  // If we already have memory allocated clear it.
  if (n_vars())
    {
      delete [] _dof_ids; _dof_ids = NULL;
    };
  
  _n_vars = static_cast<unsigned char>(nvars);

  _dof_ids = new unsigned int [n_vars()];

  // We use (invalid_id - 1) to signify n_comp(var) = 0.
  for (unsigned int v=0; v<n_vars(); v++)
    _dof_ids[v] = invalid_id - 1;
  
#endif
};



inline
unsigned int DofObject::n_comp(const unsigned int var) const
{
#ifdef ENABLE_EXPENSIVE_DATA_STRUCTURES
  
  if (var >= n_vars())
    return 0;
  
  return static_cast<unsigned int>(_n_comp[var]);

#else

  // If the dof isn't numbered yet there are effectively
  // zero dofs for this variable
  if (_dof_ids[var] == (invalid_id-1))
    return 0;
  
  return 1;
  
#endif
};



inline
void DofObject::set_n_comp(const unsigned int var,
			   const unsigned int ncomp)
{

//#if !defined(NDEBUG)

  if (ncomp != static_cast<unsigned int>(static_cast<unsigned char>(ncomp)))
    {
      std::cerr << "Unsigned char not big enough to hold ncomp!" << std::endl
		<< "Recompile with _n_comp set to a bigger type!"
		<< std::endl;
      
      error();
    };
					
//#endif

  
#ifdef ENABLE_EXPENSIVE_DATA_STRUCTURES
  
  assert (var < n_vars());

  // If we already have memory allocated clear it.
  if (n_comp(var))
    {
      delete [] _dof_ids[var]; _dof_ids[var] = NULL;
    };
  
  _n_comp[var]  = static_cast<unsigned char>(ncomp);
  _dof_ids[var] = new unsigned int [ncomp];

  for (unsigned int c=0; c<ncomp; c++)
    _dof_ids[var][c] = invalid_id;

#else

  if (ncomp == 0)
    _dof_ids[var] = (invalid_id - 1);

  else if (ncomp != 1)
    {
      std::cerr << "ERROR: You must compile with --enable-expensive to have" << std::endl
		<< "multiple components in a variable!" << std::endl;
      
      error();
    }
  
  assert (ncomp == 1);

  _dof_ids[var] = invalid_id;
  
#endif
};



inline
unsigned int DofObject::dof_number(const unsigned int var,
				   const unsigned int comp) const
{
  assert (var  < n_vars());
  assert (comp < n_comp(var));

#ifdef ENABLE_EXPENSIVE_DATA_STRUCTURES
  
  return _dof_ids[var][comp];

#else

  return _dof_ids[var];

#endif
};



inline
void DofObject::set_dof_number(const unsigned int var,
			       const unsigned int comp,
			       const unsigned int dn)
{
  assert (dn != invalid_id);
  assert (var  < n_vars());
  assert (comp < n_comp(var));
    
#ifdef ENABLE_EXPENSIVE_DATA_STRUCTURES
  
  
  _dof_ids[var][comp] = dn;

#else

  // Reserve (invalid_id-1) to signify n_com(var) = 0.
  assert (dn != (invalid_id-1));
  
  _dof_ids[var] = dn;

#endif
};
  



#endif
