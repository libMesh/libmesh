// $Id: mesh_iterators.C,v 1.8 2006-11-15 23:51:00 jwpeterson Exp $

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

// Local includes
// #include "mesh_base.h"
#include "mesh.h"
#include "elem.h"

// This file contains the implementation of all the different iterator
// functions for the mesh class.  They were put here to save space in the
// header files.


// default begin() accessor
Mesh::element_iterator
Mesh::elements_begin ()
{
  Predicates::NotNull<elem_iterator_imp> p;
  return element_iterator(_elements.begin(), _elements.end(), p);
}


// active elements begin() accessor
Mesh::element_iterator
Mesh::active_elements_begin ()
{
  Predicates::Active<elem_iterator_imp> p;
  return element_iterator(_elements.begin(), _elements.end(), p);
}


// not active elements begin() accessor
Mesh::element_iterator
Mesh::not_active_elements_begin ()
{
  Predicates::NotActive<elem_iterator_imp> p;
  return element_iterator(_elements.begin(), _elements.end(), p);
}


// subactive elements begin() accessor
Mesh::element_iterator
Mesh::subactive_elements_begin ()
{
  Predicates::SubActive<elem_iterator_imp> p;
  return element_iterator(_elements.begin(), _elements.end(), p);
}


// not subactive elements begin() accessor
Mesh::element_iterator
Mesh::not_subactive_elements_begin ()
{
  Predicates::NotSubActive<elem_iterator_imp> p;
  return element_iterator(_elements.begin(), _elements.end(), p);
}


// local elements begin() accessor
Mesh::element_iterator
Mesh::local_elements_begin ()
{
  Predicates::Local<elem_iterator_imp> p;
  return element_iterator(_elements.begin(), _elements.end(), p);
}


// active local elements begin() accessor
Mesh::element_iterator
Mesh::active_local_elements_begin ()
{
  Predicates::ActiveLocal<elem_iterator_imp> p;
  return element_iterator(_elements.begin(), _elements.end(), p);
}


// level elements begin() accessor
Mesh::element_iterator
Mesh::level_elements_begin (const unsigned int level)
{
  Predicates::Level<elem_iterator_imp> p(level);
  return element_iterator(_elements.begin(), _elements.end(), p);
}


// not level elements begin() accessor
Mesh::element_iterator
Mesh::not_level_elements_begin (const unsigned int level)
{
  Predicates::NotLevel<elem_iterator_imp> p(level);
  return element_iterator(_elements.begin(), _elements.end(), p);
}



// pid elements begin() accessor
Mesh::element_iterator
Mesh::pid_elements_begin (const unsigned int proc_id)
{
  Predicates::PID<elem_iterator_imp> p(proc_id);
  return element_iterator(_elements.begin(), _elements.end(), p);
}


// type elements begin() accessor
Mesh::element_iterator
Mesh::type_elements_begin (const ElemType type)
{
  Predicates::Type<elem_iterator_imp> p(type);
  return element_iterator(_elements.begin(), _elements.end(), p);
}



// active type elements begin() accessor
Mesh::element_iterator
Mesh::active_type_elements_begin (const ElemType type)
{
  Predicates::ActiveType<elem_iterator_imp> p(type);
  return element_iterator(_elements.begin(), _elements.end(), p);
}



// active pid elements begin() accessor
Mesh::element_iterator
Mesh::active_pid_elements_begin (const unsigned int proc_id)
{
  Predicates::ActivePID<elem_iterator_imp> p(proc_id);
  return element_iterator(_elements.begin(), _elements.end(), p);
}










// default const begin() accessor
Mesh::const_element_iterator
Mesh::elements_begin () const
{
  Predicates::NotNull<const_elem_iterator_imp> p;
  return const_element_iterator(_elements.begin(), _elements.end(), p);
}



// const active begin() accessor
Mesh::const_element_iterator
Mesh::active_elements_begin () const
{
  Predicates::Active<const_elem_iterator_imp> p;
  return const_element_iterator(_elements.begin(), _elements.end(), p);
}


// const not active begin() accessor
Mesh::const_element_iterator
Mesh::not_active_elements_begin () const
{
  Predicates::NotActive<const_elem_iterator_imp> p;
  return const_element_iterator(_elements.begin(), _elements.end(), p);
}


// const subactive begin() accessor
Mesh::const_element_iterator
Mesh::subactive_elements_begin () const
{
  Predicates::SubActive<const_elem_iterator_imp> p;
  return const_element_iterator(_elements.begin(), _elements.end(), p);
}


// const not subactive begin() accessor
Mesh::const_element_iterator
Mesh::not_subactive_elements_begin () const
{
  Predicates::NotSubActive<const_elem_iterator_imp> p;
  return const_element_iterator(_elements.begin(), _elements.end(), p);
}


// const local begin() accessor

Mesh::const_element_iterator
Mesh::local_elements_begin () const
{
  Predicates::Local<const_elem_iterator_imp> p;
  return const_element_iterator(_elements.begin(), _elements.end(), p);
}


// const active local begin() accessor
Mesh::const_element_iterator
Mesh::active_local_elements_begin () const
{
  Predicates::ActiveLocal<const_elem_iterator_imp> p;
  return const_element_iterator(_elements.begin(), _elements.end(), p);
}



// const level begin() accessor
Mesh::const_element_iterator
Mesh::level_elements_begin (const unsigned int level) const
{
  Predicates::Level<const_elem_iterator_imp> p(level);
  return const_element_iterator(_elements.begin(), _elements.end(), p);
}



// const not level begin() accessor
Mesh::const_element_iterator
Mesh::not_level_elements_begin (const unsigned int level) const
{
  Predicates::NotLevel<const_elem_iterator_imp> p(level);
  return const_element_iterator(_elements.begin(), _elements.end(), p);
}



// const pid begin() accessor
Mesh::const_element_iterator
Mesh::pid_elements_begin (const unsigned int proc_id) const
{
  Predicates::PID<const_elem_iterator_imp> p(proc_id);
  return const_element_iterator(_elements.begin(), _elements.end(), p);
}


// const type begin() accessor
Mesh::const_element_iterator
Mesh::type_elements_begin (const ElemType type) const
{
  Predicates::Type<const_elem_iterator_imp> p(type);
  return const_element_iterator(_elements.begin(), _elements.end(), p);
}



// const active type begin() accessor

Mesh::const_element_iterator
Mesh::active_type_elements_begin (const ElemType type) const
{
  Predicates::ActiveType<const_elem_iterator_imp> p(type);
  return const_element_iterator(_elements.begin(), _elements.end(), p);
}




// const active pid elements begin() accessor
Mesh::const_element_iterator
Mesh::active_pid_elements_begin (const unsigned int proc_id) const
{
  Predicates::ActivePID<const_elem_iterator_imp> p(proc_id);
  return const_element_iterator(_elements.begin(), _elements.end(), p);
}









// default end() accessor
Mesh::element_iterator
Mesh::elements_end ()
{
  Predicates::NotNull<elem_iterator_imp> p;
  return element_iterator(_elements.end(), _elements.end(), p);
}



// active end() accessor
Mesh::element_iterator
Mesh::active_elements_end ()
{
  Predicates::Active<elem_iterator_imp> p;
  return element_iterator(_elements.end(), _elements.end(), p);
}



// active end() accessor
Mesh::element_iterator
Mesh::not_active_elements_end ()
{
  Predicates::NotActive<elem_iterator_imp> p;
  return element_iterator(_elements.end(), _elements.end(), p);
}



// subactive end() accessor
Mesh::element_iterator
Mesh::subactive_elements_end ()
{
  Predicates::SubActive<elem_iterator_imp> p;
  return element_iterator(_elements.end(), _elements.end(), p);
}



// subactive end() accessor
Mesh::element_iterator
Mesh::not_subactive_elements_end ()
{
  Predicates::NotSubActive<elem_iterator_imp> p;
  return element_iterator(_elements.end(), _elements.end(), p);
}



// local end() accessor
Mesh::element_iterator
Mesh::local_elements_end ()
{
  Predicates::Local<elem_iterator_imp> p;
  return element_iterator(_elements.end(), _elements.end(), p);
}


// active local end() accessor
Mesh::element_iterator
Mesh::active_local_elements_end ()
{
  Predicates::ActiveLocal<elem_iterator_imp> p;
  return element_iterator(_elements.end(), _elements.end(), p);
}



// level end() accessor
Mesh::element_iterator
Mesh::level_elements_end (const unsigned int level)
{
  Predicates::Level<elem_iterator_imp> p(level);
  return element_iterator(_elements.end(), _elements.end(), p);
}



// not level end() accessor
Mesh::element_iterator
Mesh::not_level_elements_end (const unsigned int level)
{
  Predicates::NotLevel<elem_iterator_imp> p(level);
  return element_iterator(_elements.end(), _elements.end(), p);
}




// pid end() accessor
Mesh::element_iterator
Mesh::pid_elements_end (const unsigned int proc_id)
{
  Predicates::PID<elem_iterator_imp> p(proc_id);
  return element_iterator(_elements.end(), _elements.end(), p);
}



// type end() accessor
Mesh::element_iterator
Mesh::type_elements_end (const ElemType type)
{
  Predicates::Type<elem_iterator_imp> p(type);
  return element_iterator(_elements.end(), _elements.end(), p);
}



// active type end() accessor
Mesh::element_iterator
Mesh::active_type_elements_end (const ElemType type)
{
  Predicates::ActiveType<elem_iterator_imp> p(type);
  return element_iterator(_elements.end(), _elements.end(), p);
}


// active PID end() accessor
Mesh::element_iterator
Mesh::active_pid_elements_end (const unsigned int proc_id)
{
  Predicates::ActivePID<elem_iterator_imp> p(proc_id);
  return element_iterator(_elements.end(), _elements.end(), p);
}












// default const end() accessor
Mesh::const_element_iterator
Mesh::elements_end () const
{
  Predicates::NotNull<const_elem_iterator_imp> p;
  return const_element_iterator(_elements.end(), _elements.end(), p);
}




// active const end() accessor
Mesh::const_element_iterator
Mesh::active_elements_end () const
{
  Predicates::Active<const_elem_iterator_imp> p;
  return const_element_iterator(_elements.end(), _elements.end(), p);
}



// not active const end() accessor
Mesh::const_element_iterator
Mesh::not_active_elements_end () const
{
  Predicates::NotActive<const_elem_iterator_imp> p;
  return const_element_iterator(_elements.end(), _elements.end(), p);
}



// subactive const end() accessor
Mesh::const_element_iterator
Mesh::subactive_elements_end () const
{
  Predicates::SubActive<const_elem_iterator_imp> p;
  return const_element_iterator(_elements.end(), _elements.end(), p);
}



// not subactive const end() accessor
Mesh::const_element_iterator
Mesh::not_subactive_elements_end () const
{
  Predicates::NotSubActive<const_elem_iterator_imp> p;
  return const_element_iterator(_elements.end(), _elements.end(), p);
}




// local const end() accessor
Mesh::const_element_iterator
Mesh::local_elements_end () const
{
  Predicates::Local<const_elem_iterator_imp> p;
  return const_element_iterator(_elements.end(), _elements.end(), p);
}


// local active const end() accessor
Mesh::const_element_iterator
Mesh::active_local_elements_end () const
{
  Predicates::ActiveLocal<const_elem_iterator_imp> p;
  return const_element_iterator(_elements.end(), _elements.end(), p);
}




// level const end() accessor
Mesh::const_element_iterator
Mesh::level_elements_end (const unsigned int level) const
{
  Predicates::Level<const_elem_iterator_imp> p(level);
  return const_element_iterator(_elements.end(), _elements.end(), p);
}




// not level const end() accessor
Mesh::const_element_iterator
Mesh::not_level_elements_end (const unsigned int level) const
{
  Predicates::NotLevel<const_elem_iterator_imp> p(level);
  return const_element_iterator(_elements.end(), _elements.end(), p);
}





// pid const end() accessor
Mesh::const_element_iterator
Mesh::pid_elements_end (const unsigned int proc_id) const
{
  Predicates::PID<const_elem_iterator_imp> p(proc_id);
  return const_element_iterator(_elements.end(), _elements.end(), p);
}



// type const end() accessor
Mesh::const_element_iterator
Mesh::type_elements_end (const ElemType type) const
{
  Predicates::Type<const_elem_iterator_imp> p(type);
  return const_element_iterator(_elements.end(), _elements.end(), p);
}



// active type const end() accessor
Mesh::const_element_iterator
Mesh::active_type_elements_end (const ElemType type) const
{
  Predicates::ActiveType<const_elem_iterator_imp> p(type);
  return const_element_iterator(_elements.end(), _elements.end(), p);
}


// active PID end() accessor
Mesh::const_element_iterator
Mesh::active_pid_elements_end (const unsigned int proc_id) const
{
  Predicates::ActivePID<const_elem_iterator_imp> p(proc_id);
  return const_element_iterator(_elements.end(), _elements.end(), p);
}









// default nodes begin() accessor
Mesh::node_iterator
Mesh::nodes_begin ()
{
  Predicates::NotNull<node_iterator_imp> p;
  return node_iterator(_nodes.begin(), _nodes.end(), p);
}




// active nodes begin() accessor
Mesh::node_iterator
Mesh::active_nodes_begin ()
{
  Predicates::Active<node_iterator_imp> p;
  return node_iterator(_nodes.begin(), _nodes.end(), p);
}



// local nodes begin() accessor
Mesh::node_iterator
Mesh::local_nodes_begin ()
{
  Predicates::Local<node_iterator_imp> p;
  return node_iterator(_nodes.begin(), _nodes.end(), p);
}



// pid nodes begin() accessor
Mesh::node_iterator
Mesh::pid_nodes_begin (const unsigned int proc_id)
{
  Predicates::PID<node_iterator_imp> p(proc_id);
  return node_iterator(_nodes.begin(), _nodes.end(), p);
}



// default const nodes begin() accessor
Mesh::const_node_iterator
Mesh::nodes_begin () const
{
  Predicates::NotNull<const_node_iterator_imp> p;
  return const_node_iterator(_nodes.begin(), _nodes.end(), p);
}



// active const nodes begin() accessor
Mesh::const_node_iterator
Mesh::active_nodes_begin () const
{
  Predicates::Active<const_node_iterator_imp> p;
  return const_node_iterator(_nodes.begin(), _nodes.end(), p);
}



// local const nodes begin() accessor
Mesh::const_node_iterator
Mesh::local_nodes_begin () const
{
  Predicates::Local<const_node_iterator_imp> p;
  return const_node_iterator(_nodes.begin(), _nodes.end(), p);
}


// pid const nodes begin() accessor
Mesh::const_node_iterator
Mesh::pid_nodes_begin (const unsigned int proc_id) const
{
  Predicates::PID<const_node_iterator_imp> p(proc_id);
  return const_node_iterator(_nodes.begin(), _nodes.end(), p);
}





// default nodes end() accessor
Mesh::node_iterator
Mesh::nodes_end ()
{
  Predicates::NotNull<node_iterator_imp> p;
  return node_iterator(_nodes.end(), _nodes.end(), p);
}




// active nodes end() accessor
Mesh::node_iterator
Mesh::active_nodes_end ()
{
  Predicates::Active<node_iterator_imp> p;
  return node_iterator(_nodes.end(), _nodes.end(), p);
}



// local nodes end() accessor
Mesh::node_iterator
Mesh::local_nodes_end ()
{
  Predicates::Local<node_iterator_imp> p;
  return node_iterator(_nodes.end(), _nodes.end(), p);
}


// pid nodes end() accessor
Mesh::node_iterator
Mesh::pid_nodes_end (const unsigned int proc_id)
{
  Predicates::PID<node_iterator_imp> p(proc_id);
  return node_iterator(_nodes.end(), _nodes.end(), p);
}



// default const nodes end() accessor
Mesh::const_node_iterator
Mesh::nodes_end () const
{
  Predicates::NotNull<const_node_iterator_imp> p;
  return const_node_iterator(_nodes.end(), _nodes.end(), p);
}


// const active nodes end() accessor
Mesh::const_node_iterator
Mesh::active_nodes_end () const
{
  Predicates::Active<const_node_iterator_imp> p;
  return const_node_iterator(_nodes.end(), _nodes.end(), p);
}


// local const nodes end() accessor
Mesh::const_node_iterator
Mesh::local_nodes_end () const
{
  Predicates::Local<const_node_iterator_imp> p;
  return const_node_iterator(_nodes.end(), _nodes.end(), p);
}


// pid const nodes end() accessor
Mesh::const_node_iterator
Mesh::pid_nodes_end (const unsigned int proc_id) const
{
  Predicates::PID<const_node_iterator_imp> p(proc_id);
  return const_node_iterator(_nodes.end(), _nodes.end(), p);
}
