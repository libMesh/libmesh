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

#ifndef LIBMESH_MULTI_PREDICATES_H
#define LIBMESH_MULTI_PREDICATES_H

// Local includes
#include "libmesh/single_predicates.h"

// C++ includes
#include <vector>

namespace libMesh {
class Elem;
}

namespace libMesh
{

// Forward declarations
class BoundaryInfo;

/**
 * This namespace defines several multi_predicates which are used by
 * the element and node iterators.  These classes are not in general
 * used by the user, although they could be.
 *
 * \author John W. Peterson
 * \date 2004
 */
namespace Predicates
{

// Empty place-holder base class for multi_predicates
struct multi_predicate {};


// This class represents a generic combination of more than one predicate.
// It is meant to be derived from to actually be used.
template <typename T>
struct abstract_multi_predicate : multi_predicate
{
  // virtual destructor.
  virtual ~abstract_multi_predicate()
  {
    // Clean-up vector
    for (std::size_t i=0; i<_predicates.size(); ++i)
      delete _predicates[i];
  }

  // operator= (perform deep copy of entries in _predicates vector
  abstract_multi_predicate & operator=(const abstract_multi_predicate & rhs)
  {
    // First clear out the predicates vector
    for (std::size_t i=0; i<_predicates.size(); ++i)
      delete _predicates[i];

    // Now copy over the information from the rhs.
    this->deep_copy(rhs);

    return *this;
  }

  // operator() checks all the predicates in the vector.
  virtual bool operator()(const T & it) const
  {
    for (std::size_t i=0; i<_predicates.size(); ++i)
      {
        const predicate<T> * pred = _predicates[i];

        libmesh_assert (pred);

        if ( ! (*pred)(it) )
          return false;
      }

    return true;
  }

protected:
  // Do not instantiate the base class.
  abstract_multi_predicate() {}

  // Copy constructor.
  abstract_multi_predicate(const abstract_multi_predicate & rhs)
  {
    this->deep_copy(rhs);
  }

  // The deep_copy function is used by both the op= and
  // copy constructors.  This function uses the default (empty)
  // copy constructor for the predicate class.
  void deep_copy(const abstract_multi_predicate & rhs)
  {
    for (std::size_t i=0; i<rhs._predicates.size(); ++i)
      _predicates.push_back(rhs._predicates[i]->clone());
  }

  // Predicates to be evaluated.
  std::vector<predicate<T> *> _predicates;
};



/**
 * Used to iterate over NULL entries in a container.
 */
template <typename T>
struct IsNull : abstract_multi_predicate<T>
{
  // Constructor, pushes back a single predicate
  IsNull()
  {
    this->_predicates.push_back(new is_null<T>);
  }
};



/**
 * Used to iterate over non-NULL entries in a container.
 */
template <typename T>
struct NotNull : abstract_multi_predicate<T>
{
  // Constructor, pushes back a single predicate
  NotNull()
  {
    this->_predicates.push_back(new not_null<T>);
  }
};



/**
 * Used to iterate over non-NULL, active entries in a container.
 */
template <typename T>
struct Active : abstract_multi_predicate<T>
{
  // Constructor, pushes back two single predicates
  Active()
  {
    this->_predicates.push_back(new not_null<T>);
    this->_predicates.push_back(new active<T>);
  }
};



/**
 * Used to iterate over non-NULL, inactive entries in a container.
 */
template <typename T>
struct NotActive : abstract_multi_predicate<T>
{
  // Constructor, pushes back two single predicates
  NotActive()
  {
    this->_predicates.push_back(new not_null<T>);
    this->_predicates.push_back(new not_active<T>);
  }
};




/**
 * Used to iterate over non-NULL, entries that have children (i.e. are
 * ancestors) in a container.
 */
template <typename T>
struct Ancestor : abstract_multi_predicate<T>
{
  // Constructor, pushes back two single predicates
  Ancestor()
  {
    this->_predicates.push_back(new not_null<T>);
    this->_predicates.push_back(new ancestor<T>);
  }
};




/**
 * Used to iterate over non-NULL, entries that have no children (i.e. are not
 * ancestors) in a container.
 */
template <typename T>
struct NotAncestor : abstract_multi_predicate<T>
{
  // Constructor, pushes back two single predicates
  NotAncestor()
  {
    this->_predicates.push_back(new not_null<T>);
    this->_predicates.push_back(new not_ancestor<T>);
  }
};




/**
 * Used to iterate over non-NULL, subactive entries (i.e. has no
 * active children) in a container.
 */
template <typename T>
struct SubActive : abstract_multi_predicate<T>
{
  // Constructor, pushes back two single predicates
  SubActive()
  {
    this->_predicates.push_back(new not_null<T>);
    this->_predicates.push_back(new subactive<T>);
  }
};




/**
 * Used to iterate over non-NULL, non-subactive entries (i.e. has one
 * or more active children) in a container.
 */
template <typename T>
struct NotSubActive : abstract_multi_predicate<T>
{
  // Constructor, pushes back two single predicates
  NotSubActive()
  {
    this->_predicates.push_back(new not_null<T>);
    this->_predicates.push_back(new not_subactive<T>);
  }
};



/**
 * Used to iterate over non-NULL, local entries (i.e. owned by the
 * current processor) in a container.
 */
template <typename T>
struct Local : abstract_multi_predicate<T>
{
  // Constructor, pushes back two single predicates
  Local(processor_id_type my_pid)
  {
    this->_predicates.push_back(new not_null<T>);
    this->_predicates.push_back(new pid<T>(my_pid));
  }
};


/**
 * Used to iterate over non-NULL, semi-local entries (i.e. are not
 * subactive and have are owned by an attached processor) in a
 * container.
 *
 * FIXME: This is not currently safe to use on adaptively-refined
 * grids, it should be added back when Elem::is_semilocal() has been
 * patched to not require the Elem to be active.
 */
// template <typename T>
// struct SemiLocal : abstract_multi_predicate<T>
// {
//   // Constructor, pushes back two single predicates
//   SemiLocal(processor_id_type my_pid)
//   {
//     this->_predicates.push_back(new not_null<T>);
//     this->_predicates.push_back(new not_subactive<T>);
//     this->_predicates.push_back(new semilocal_pid<T>(my_pid));
//   }
// };


/**
 * Used to iterate over non-NULL, active, non sub-active, semi-local
 * elements in a container.
 */
template <typename T>
struct ActiveSemiLocal : abstract_multi_predicate<T>
{
  // Constructor, pushes back two single predicates
  ActiveSemiLocal(processor_id_type my_pid)
  {
    this->_predicates.push_back(new not_null<T>);
    this->_predicates.push_back(new active<T>);
    this->_predicates.push_back(new not_subactive<T>);
    this->_predicates.push_back(new semilocal_pid<T>(my_pid));
  }
};


/**
 * Used to iterate over non-NULL, face-local entries (i.e. are not
 * subactive and are on or have a neighbor on processor my_pid) in a
 * container.
 */
template <typename T>
struct FaceLocal : abstract_multi_predicate<T>
{
  // Constructor, pushes back two single predicates
  FaceLocal(processor_id_type my_pid)
  {
    this->_predicates.push_back(new not_null<T>);
    this->_predicates.push_back(new not_subactive<T>);
    this->_predicates.push_back(new facelocal_pid<T>(my_pid));
  }
};



/**
 * Used to iterate over non-NULL, non-local entries in a
 * container.
 */
template <typename T>
struct NotLocal : abstract_multi_predicate<T>
{
  // Constructor, pushes back two single predicates
  NotLocal(processor_id_type my_pid)
  {
    this->_predicates.push_back(new not_null<T>);
    this->_predicates.push_back(new not_pid<T>(my_pid));
  }
};


/**
 * Used to iterate over non-NULL, active, non-local entries in a
 * container.
 */
template <typename T>
struct ActiveNotLocal : abstract_multi_predicate<T>
{
  // Constructor, pushes back two single predicates
  ActiveNotLocal(processor_id_type my_pid)
  {
    this->_predicates.push_back(new not_null<T>);
    this->_predicates.push_back(new active<T>);
    this->_predicates.push_back(new not_pid<T>(my_pid));
  }
};


/**
 * Used to iterate over non-NULL, elements of a given geometric type.
 */
template <typename T>
struct Type : abstract_multi_predicate<T>
{
  Type(ElemType type)
  {
    this->_predicates.push_back(new not_null<T>);
    this->_predicates.push_back(new elem_type<T>(type));
  }
};



/**
 * Used to iterate over non-NULL, active elements of a given geometric type.
 */
template <typename T>
struct ActiveType : abstract_multi_predicate<T>
{
  ActiveType(ElemType type)
  {
    this->_predicates.push_back(new not_null<T>);
    this->_predicates.push_back(new active<T>);
    this->_predicates.push_back(new elem_type<T>(type));
  }
};



#ifdef LIBMESH_ENABLE_AMR
/**
 * Used to iterate over non-NULL, elements with a given refinement
 * flag.
 */
template <typename T>
struct Flagged : abstract_multi_predicate<T>
{
  Flagged(unsigned char rflag)
  {
    this->_predicates.push_back(new not_null<T>);
    this->_predicates.push_back(new flagged<T>(rflag));
  }
};



/**
 * Used to iterate over non-NULL, elements with a given refinement
 * flag belonging to a given processor.
 */
template <typename T>
struct FlaggedPID : abstract_multi_predicate<T>
{
  FlaggedPID(unsigned char rflag, processor_id_type proc_id)
  {
    this->_predicates.push_back(new not_null<T>);
    this->_predicates.push_back(new flagged<T>(rflag));
    this->_predicates.push_back(new pid<T>(proc_id));
  }
};

#endif // LIBMESH_ENABLE_AMR




/**
 * Used to iterate over non-NULL, active elements owned by a given
 * processor.
 */
template <typename T>
struct ActivePID : abstract_multi_predicate<T>
{
  ActivePID(processor_id_type proc_id)
  {
    this->_predicates.push_back(new not_null<T>);
    this->_predicates.push_back(new active<T>);
    this->_predicates.push_back(new pid<T>(proc_id));
  }
};





/**
 * Used to iterate over non-NULL, active, local elements owned by a
 * given processor.
 */
template <typename T>
struct ActiveLocal : abstract_multi_predicate<T>
{
  ActiveLocal(processor_id_type my_pid)
  {
    this->_predicates.push_back(new not_null<T>);
    this->_predicates.push_back(new active<T>);
    this->_predicates.push_back(new pid<T>(my_pid));
  }
};





/**
 * Used to iterate over non-NULL elements owned by a given processor.
 */
template <typename T>
struct PID : abstract_multi_predicate<T>
{
  PID(processor_id_type proc_id)
  {
    this->_predicates.push_back(new not_null<T>);
    this->_predicates.push_back(new pid<T>(proc_id));
  }
};



/**
 * Used to iterate over non-NULL elements on the boundary with a given
 * ID.
 */
template <typename T>
struct BID : abstract_multi_predicate<T>
{
  BID(boundary_id_type bndry_id, const BoundaryInfo & bndry_info)
  {
    this->_predicates.push_back(new not_null<T>);
    this->_predicates.push_back(new bid<T>(bndry_id, bndry_info));
  }
};



/**
 * Used to iterate over non-NULL elements on the boundary.
 */
template <typename T>
struct BND : abstract_multi_predicate<T>
{
  BND(const BoundaryInfo & bndry_info)
  {
    this->_predicates.push_back(new not_null<T>);
    this->_predicates.push_back(new bnd<T>(bndry_info));
  }
};



/**
 * Used to iterate over non-NULL elements *not* owned by a given
 * processor.
 */
template <typename T>
struct NotPID : abstract_multi_predicate<T>
{
  NotPID(processor_id_type proc_id)
  {
    this->_predicates.push_back(new not_null<T>);
    this->_predicates.push_back(new not_pid<T>(proc_id));
  }
};



/**
 * Used to iterate over non-NULL elements of a specified (refinement) level.
 */
template <typename T>
struct Level : abstract_multi_predicate<T>
{
  Level(unsigned int l)
  {
    this->_predicates.push_back(new not_null<T>);
    this->_predicates.push_back(new level<T>(l));
  }
};



/**
 * Used to iterate over non-NULL elements *not* of a specified
 * (refinement) level.
 */
template <typename T>
struct NotLevel : abstract_multi_predicate<T>
{
  NotLevel(unsigned int l)
  {
    this->_predicates.push_back(new not_null<T>);
    this->_predicates.push_back(new not_level<T>(l));
  }
};



/**
 * Used to iterate over non-NULL local elements with a specified
 * (refinement) level.
 */
template <typename T>
struct LocalLevel : abstract_multi_predicate<T>
{
  LocalLevel(processor_id_type my_pid,
             unsigned int l)
  {
    this->_predicates.push_back(new not_null<T>);
    this->_predicates.push_back(new pid<T>(my_pid));
    this->_predicates.push_back(new level<T>(l));
  }
};



/**
 * Used to iterate over non-NULL local elements *not* of a specified
 * (refinement) level.
 */
template <typename T>
struct LocalNotLevel : abstract_multi_predicate<T>
{
  LocalNotLevel(processor_id_type my_pid,
                unsigned int l)
  {
    this->_predicates.push_back(new not_null<T>);
    this->_predicates.push_back(new pid<T>(my_pid));
    this->_predicates.push_back(new not_level<T>(l));
  }
};



/**
 * Used to iterate over non-NULL, active elements which are on the
 * boundary.
 */
template <typename T>
struct ActiveOnBoundary : abstract_multi_predicate<T>
{
  ActiveOnBoundary()
  {
    this->_predicates.push_back(new not_null<T>);
    this->_predicates.push_back(new active<T>);
    this->_predicates.push_back(new null_neighbor<T>);
  }
};



/**
 * Used to iterate over the sides of an element which are on the
 * boundary of the Mesh.
 */
template <typename T>
struct BoundarySide : abstract_multi_predicate<T>
{
  BoundarySide()
  {
    this->_predicates.push_back(new boundary_side<T>);
  }
};



/**
 * Used to iterate over non-NULL, active elements with a given PID on
 * a given subdomain.
 */
template <typename T>
struct ActiveLocalSubdomain : abstract_multi_predicate<T>
{
  ActiveLocalSubdomain(processor_id_type my_pid,
                       subdomain_id_type subdomain_id)
  {
    this->_predicates.push_back(new not_null<T>);
    this->_predicates.push_back(new active<T>);
    this->_predicates.push_back(new pid<T>(my_pid));
    this->_predicates.push_back(new subdomain<T>(subdomain_id));
  }
};



/**
 * Used to iterate over non-NULL, active elements on a given
 * subdomain.
 */
template <typename T>
struct ActiveSubdomain : abstract_multi_predicate<T>
{
  ActiveSubdomain(subdomain_id_type subdomain_id)
  {
    this->_predicates.push_back(new not_null<T>);
    this->_predicates.push_back(new active<T>);
    this->_predicates.push_back(new subdomain<T>(subdomain_id));
  }
};



/**
 * Used to iterate over non-NULL, active elements whose
 * subdomains are in a user-specified set.
 */
template <typename T>
struct ActiveSubdomainSet : abstract_multi_predicate<T>
{
  ActiveSubdomainSet(std::set<subdomain_id_type> sset)
  {
    this->_predicates.push_back(new not_null<T>);
    this->_predicates.push_back(new active<T>);
    this->_predicates.push_back(new subdomain_set<T>(sset));
  }
};



/**
 * Used to iterate over non-NULL elements not owned by a given
 * processor but semi-local to that processor, i.e. ghost elements.
 */
template <typename T>
struct Ghost : abstract_multi_predicate<T>
{
  Ghost(processor_id_type my_pid)
  {
    this->_predicates.push_back(new not_null<T>);
    this->_predicates.push_back(new active<T>);
    this->_predicates.push_back(new not_pid<T>(my_pid));
    this->_predicates.push_back(new semilocal_pid<T>(my_pid));
  }
};



/**
 * Used to iterate over elements where solutions indexed by a given
 * DofMap are evaluable for a given variable var_num.
 */
template <typename T>
struct Evaluable: abstract_multi_predicate<T>
{
  Evaluable(processor_id_type my_pid,
            const DofMap & dof_map,
            unsigned int var_num = libMesh::invalid_uint)
  {
    this->_predicates.push_back(new not_null<T>);
    this->_predicates.push_back(new active<T>);
    this->_predicates.push_back(new pid<T>(my_pid));
    this->_predicates.push_back(new evaluable<T>(dof_map, var_num));
  }
};

}


} // namespace libMesh

#endif // LIBMESH_MULTI_PREDICATES_H
