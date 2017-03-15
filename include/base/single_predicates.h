// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_SINGLE_PREDICATES_H
#define LIBMESH_SINGLE_PREDICATES_H

// Local includes
#include <cstddef>         // for NULL with gcc 4.6.2 - I'm serious!
#include "libmesh/libmesh_common.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/id_types.h"

// C++ includes
#include <vector>
#include <set>

namespace libMesh
{

// Forward declarations
class BoundaryInfo;
class DofMap;


/**
 * This file declares several predicates in the Predicates namespace.  They
 * are called "single predicates" since the purpose of each one is to act
 * as a single functor which returns true or false depending on the result
 * of the operator() function.  The single predicates are used together
 * as building blocks to create the "multi predicates" which can be found
 * in the multi_predicates.h header file.
 *
 * \author John W. Peterson
 * \date 2004
 */
namespace Predicates
{
// Forward declaration
template <typename T> struct abstract_multi_predicate;

// abstract single predicate.  Derived classes must implement the clone()
// function.  Be careful using it since it allocates memory!  The clone()
// function is necessary since the predicate class has pure virtual
// functions.
template <typename T>
struct predicate
{
  virtual ~predicate() {}
  virtual bool operator()(const T & it) const = 0;

protected:
  friend struct abstract_multi_predicate<T>;
  virtual predicate * clone() const = 0;
};


// The is_null predicate returns true if the underlying pointer is NULL.
template <typename T>
struct is_null : predicate<T>
{
  virtual ~is_null() {}
  virtual bool operator()(const T & it) const libmesh_override { return *it == libmesh_nullptr; }

protected:
  virtual predicate<T> * clone() const libmesh_override { return new is_null<T>(*this); }
};

// The not_null predicate simply returns true if the pointer is not NULL.
template <typename T>
struct not_null : is_null<T>
{
  virtual bool operator()(const T & it) const libmesh_override { return !is_null<T>::operator()(it); }

protected:
  virtual predicate<T> * clone() const libmesh_override { return new not_null<T>(*this); }
};


// The active predicate returns true if the pointer is active.
template <typename T>
struct active : predicate<T>
{
  virtual ~active() {}
  virtual bool operator()(const T & it) const libmesh_override { return (*it)->active(); }

protected:
  virtual predicate<T> * clone() const libmesh_override { return new active<T>(*this); }
};

// The not active predicate returns true when the pointer is inactive
template <typename T>
struct not_active : active<T>
{
  virtual bool operator()(const T & it) const libmesh_override { return !active<T>::operator()(it); }

protected:
  virtual predicate<T> * clone() const libmesh_override { return new not_active<T>(*this); }
};


// The ancestor predicate returns true if the pointer is ancestor.
template <typename T>
struct ancestor : predicate<T>
{
  virtual ~ancestor() {}
  virtual bool operator()(const T & it) const libmesh_override { return (*it)->ancestor(); }

protected:
  virtual predicate<T> * clone() const libmesh_override { return new ancestor<T>(*this); }
};

// The not_ancestor predicate returns true when the pointer is not ancestor
template <typename T>
struct not_ancestor : ancestor<T>
{
  virtual bool operator()(const T & it) const libmesh_override { return !ancestor<T>::operator()(it); }

protected:
  virtual predicate<T> * clone() const libmesh_override { return new not_ancestor<T>(*this); }
};


// The subactive predicate returns true if the pointer is subactive.
template <typename T>
struct subactive : predicate<T>
{
  virtual ~subactive() {}
  virtual bool operator()(const T & it) const libmesh_override { return (*it)->subactive(); }

protected:
  virtual predicate<T> * clone() const libmesh_override { return new subactive<T>(*this); }
};

// The not_subactive predicate returns true when the pointer is not subactive
template <typename T>
struct not_subactive : subactive<T>
{
  virtual bool operator()(const T & it) const libmesh_override { return !subactive<T>::operator()(it); }

protected:
  virtual predicate<T> * clone() const libmesh_override { return new not_subactive<T>(*this); }
};



// The pid predicate returns true if the pointers
// processor id matches a given processor id.
template <typename T>
struct pid : predicate<T>
{
  // Constructor
  pid(processor_id_type p) : _pid(p) {}
  virtual ~pid() {}

  // op()
  virtual bool operator()(const T & it) const libmesh_override { return (*it)->processor_id() == _pid; }

protected:
  virtual predicate<T> * clone() const libmesh_override { return new pid<T>(*this); }
  const processor_id_type _pid;
};



// The bid predicate returns true if has_boundary_id(node, id) returns true.
template <typename T>
struct bid : predicate<T>
{
  // Constructor
  bid(boundary_id_type bid,
      const BoundaryInfo & bndry_info) :
    _bid(bid),
    _bndry_info(bndry_info)
  {}
  virtual ~bid() {}

  // op()
  virtual bool operator()(const T & it) const libmesh_override;

protected:
  virtual predicate<T> * clone() const libmesh_override { return new bid<T>(*this); }
  const boundary_id_type _bid;
  const BoundaryInfo & _bndry_info;
};



// The bnd predicate returns true if n_boundary_ids(node) > 0.
template <typename T>
struct bnd : predicate<T>
{
  // Constructor
  bnd(const BoundaryInfo & bndry_info) :
    _bndry_info(bndry_info)
  {}
  virtual ~bnd() {}

  // op()
  virtual bool operator()(const T & it) const libmesh_override;

protected:
  virtual predicate<T> * clone() const libmesh_override { return new bnd<T>(*this); }
  const BoundaryInfo & _bndry_info;
};



// The semilocal_pid predicate returns true if the element
// pointed to is semilocal to (has nodes shared with an element of) a
// given processor id.
template <typename T>
struct semilocal_pid : predicate<T>
{
  // Constructor
  semilocal_pid(processor_id_type p) : _pid(p) {}
  virtual ~semilocal_pid() {}

  // op()
  virtual bool operator()(const T & it) const libmesh_override { return (*it)->is_semilocal(_pid); }

protected:
  virtual predicate<T> * clone() const libmesh_override { return new semilocal_pid<T>(*this); }
  const processor_id_type _pid;
};



// The facelocal_pid predicate returns true if the element
// pointed to is face-local to (is on or has a neighbor on the
// partition of) a given processor id.
template <typename T>
struct facelocal_pid : predicate<T>
{
  // Constructor
  facelocal_pid(processor_id_type p) : _pid(p) {}
  virtual ~facelocal_pid() {}

  // op()
  virtual bool operator()(const T & it) const libmesh_override
  {
    if ((*it)->processor_id() == _pid)
      return true;
    for (unsigned int n = 0; n != (*it)->n_neighbors(); ++n)
      if ((*it)->neighbor_ptr(n) &&
          (*it)->neighbor_ptr(n)->processor_id() == _pid)
        return true;
    return false;
  }

protected:
  virtual predicate<T> * clone() const libmesh_override { return new facelocal_pid<T>(*this); }
  const processor_id_type _pid;
};



// The not_pid predicate returns ture if the pointers
// processor id does _not_ match p.
template <typename T>
struct not_pid : pid<T>
{
  not_pid(processor_id_type p) : pid<T>(p) {}

  virtual bool operator()(const T & it) const libmesh_override { return !pid<T>::operator()(it); }

protected:
  virtual predicate<T> * clone() const libmesh_override { return new not_pid<T>(*this); }
};


// The elem_type predicate returns true if the pointers type matches
// the given type.  Of course, this one can only be instantiated for
// objects which return Elem pointers when dereferenced.
template <typename T>
struct elem_type : predicate<T>
{
  // Constructor
  elem_type (ElemType t) : _elem_type(t) {}
  virtual ~elem_type() {}

  virtual bool operator()(const T & it) const libmesh_override { return (*it)->type() == _elem_type; }

protected:
  virtual predicate<T> * clone() const libmesh_override { return new elem_type<T>(*this); }
  const ElemType _elem_type;
};



#ifdef LIBMESH_ENABLE_AMR
// The flagged predicate returns true if the pointers refinement flag
// matches the given rflag.  Of course, this one can only be
// instantiated for objects which return Elem pointers when
// dereferenced.
template <typename T>
struct flagged : predicate<T>
{
  // Constructor
  flagged (unsigned char rflag) : _rflag(rflag) {}
  virtual ~flagged() {}

  virtual bool operator()(const T & it) const libmesh_override { return (*it)->refinement_flag() == _rflag; }

protected:
  virtual predicate<T> * clone() const libmesh_override { return new flagged<T>(*this); }
  const unsigned char _rflag;
};
#endif // LIBMESH_ENABLE_AMR






// The level predicate returns true if the pointers level
// matches the given level.
template <typename T>
struct level : predicate<T>
{
  // Constructor
  level (unsigned int l) : _level(l) {}
  virtual ~level() {}

  virtual bool operator()(const T & it) const libmesh_override { return (*it)->level() == _level; }

protected:
  virtual predicate<T> * clone() const libmesh_override { return new level<T>(*this); }
  const unsigned int _level;
};



// The not_level predicate returns true if the pointers level
// _does not_ match the given level.
template <typename T>
struct not_level : level<T>
{
  // Constructor
  not_level(unsigned int l) : level<T>(l) {}

  virtual bool operator()(const T & it) const libmesh_override { return !level<T>::operator()(it); }

protected:
  virtual predicate<T> * clone() const libmesh_override { return new not_level<T>(*this); }
};




// The null_neighbor predicate returns true if the pointer has any
// NULL neigbors.
template <typename T>
struct null_neighbor : predicate<T>
{
  virtual ~null_neighbor() {}
  virtual bool operator()(const T & it) const libmesh_override
  {
    return (*it)->on_boundary();
  }

protected:
  virtual predicate<T> * clone() const libmesh_override { return new null_neighbor<T>(*this); }
};



// This predicate simply forwards the work of determining whether
// a particular side is on the boundary to the iterator itself, which
// has more information.
template <typename T>
struct boundary_side : predicate<T>
{
  virtual ~boundary_side() {}
  virtual bool operator()(const T & it) const libmesh_override
  {
    return it.side_on_boundary();
  }

protected:
  virtual predicate<T> * clone() const libmesh_override { return new boundary_side<T>(*this); }
};

// The subdomain predicate returns true if the pointers
// subdimain id matches a given subdomain id.
template <typename T>
struct subdomain : predicate<T>
{
  // Constructor
  subdomain(subdomain_id_type sid) : _subdomain(sid) {}
  virtual ~subdomain() {}

  // op()
  virtual bool operator()(const T & it) const libmesh_override { return (*it)->subdomain_id() == _subdomain; }

protected:
  virtual predicate<T> * clone() const libmesh_override { return new subdomain<T>(*this); }
  const subdomain_id_type _subdomain;
};


// The subdomain_set predicate returns true if the pointer's
// subdomain id is in the provided std::set<subdomain_id_type>.
template <typename T>
struct subdomain_set : predicate<T>
{
  // Constructor
  subdomain_set(std::set<subdomain_id_type> sset) : _subdomain_set(sset) {}
  virtual ~subdomain_set() {}

  // op()
  virtual bool operator()(const T & it) const libmesh_override { return _subdomain_set.count((*it)->subdomain_id()); }

protected:
  virtual predicate<T> * clone() const libmesh_override { return new subdomain_set<T>(*this); }
  const std::set<subdomain_id_type> _subdomain_set;
};


// The evaluable predicate returns true if the pointer (which must be
// a local element) has degrees of freedom which can be evaluated for
// the specified DofMap and variable.
template <typename T>
struct evaluable : predicate<T>
{
  // Constructor
  evaluable(const DofMap & dof_map,
            unsigned int var_num) :
    _dof_map(dof_map), _var_num(var_num) {}
  virtual ~evaluable() {}

  // op()
  virtual bool operator()(const T & it) const libmesh_override;

protected:
  virtual predicate<T> * clone() const libmesh_override { return new evaluable<T>(*this); }
  const DofMap & _dof_map;
  unsigned int _var_num;
};


} // namespace Predicates


} // namespace libMesh

#endif // LIBMESH_SINGLE_PREDICATES_H
