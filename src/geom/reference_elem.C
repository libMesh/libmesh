// The libMesh Finite Element Library.
// Copyright (C) 2002-2013 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/node.h"
#include "libmesh/elem.h"
#include "libmesh/reference_elem.h"
#include "libmesh/libmesh_singleton.h"
#include "libmesh/threads.h"

// C++ includes
#include <map>
#include <sstream>



//-----------------------------------------------
// anonymous namespace for implementation details
namespace
{
  using namespace libMesh;

  namespace ElemDataStrings
  {
#include "reference_elem.data"
  }

  typedef Threads::spin_mutex InitMutex;

  // Mutex for thread safety.
  InitMutex init_mtx;

  // map from ElemType to reference element file system object name
  typedef std::map<ElemType, const char*> FileMapType;
  FileMapType ref_elem_file;
  Elem* ref_elem_map[INVALID_ELEM];



  class SingletonCache : public libMesh::Singleton
  {
  public:
    ~SingletonCache()
    {
      for (unsigned int e=0; e<elem_list.size(); e++)
	if (elem_list[e])
	  {
	    delete elem_list[e];
	    elem_list[e] = NULL;
	  }

      elem_list.clear();

      for (unsigned int n=0; n<node_list.size(); n++)
	if (node_list[n])
	  {
	    delete node_list[n];
	    node_list[n] = NULL;
	  }

      node_list.clear();
    }

    std::vector<Node*> node_list;
    std::vector<Elem*> elem_list;
  };

  // singleton object, dynamically created and then
  // removed at program exit
  SingletonCache *singleton_cache = NULL;



  Elem* read_ref_elem (const ElemType Type,
		       std::istream &in)
  {
    libmesh_assert (singleton_cache != NULL);

    static const unsigned int comm_len = 1024;
    char comm[comm_len];

    std::string foo;
    unsigned int n_elem, n_nodes, elem_type, nn;
    double x, y, z;

    in >> foo;
    in >> n_elem;  /**/ in.getline (comm, comm_len); libmesh_assert_equal_to (n_elem, 1);
    in >> n_nodes; /**/ in.getline (comm, comm_len);
    in >> foo;     /**/ in.getline (comm, comm_len);
    in >> foo;     /**/ in.getline (comm, comm_len);
    in >> foo;     /**/ in.getline (comm, comm_len);
    in >> foo;     /**/ in.getline (comm, comm_len);
    in >> n_elem;  /**/ in.getline (comm, comm_len); libmesh_assert_equal_to (n_elem, 1);

    in >> elem_type;

    libmesh_assert_less (elem_type, INVALID_ELEM);
    libmesh_assert_equal_to (elem_type, Type);
    libmesh_assert_equal_to (n_nodes, Elem::type_to_n_nodes_map[elem_type]);

    // Construct the elem
    Elem *elem = Elem::build(static_cast<ElemType>(elem_type)).release();

    // We are expecing an identity map, so assert it!
    for (unsigned int n=0; n<n_nodes; n++)
      {
	in >> nn;
	libmesh_assert_equal_to (n,nn);
      }

    for (unsigned int n=0; n<n_nodes; n++)
      {
	in >> x >> y >> z;

	Node *node = new Node(x,y,z,n);
	singleton_cache->node_list.push_back(node);

	elem->set_node(n) = node;
      }


    // it is entirely possible we ran out of file or encountered
    // another error.  If so, cleanly abort.
    if (!in)
      {
	delete elem;
	elem = NULL;
	libMesh::err << "ERROR while creating element singleton!\n";
	libmesh_error();
      }

    else
      singleton_cache->elem_list.push_back (elem);

    ref_elem_map[Type] = elem;

    return elem;
  }



  void init_ref_elem_table()
  {
    // ouside mutex - if this pointer is set, we can trust it.
    if (singleton_cache != NULL) return;

    // playing with fire here - lock before touching shared
    // data structures
    InitMutex::scoped_lock lock(init_mtx);

    // inside mutex - pointer may have changed while waiting
    // for the lock to acquire, check it again.
    if (singleton_cache != NULL) return;

    // OK, if we get here we have the lock and we are not
    // initialized.  populate singleton.
    singleton_cache = new SingletonCache;

    // initialize the reference file table
    {
      ref_elem_file.clear();

      // // 1D elements
      ref_elem_file[EDGE2]    = ElemDataStrings::one_edge;
      ref_elem_file[EDGE3]    = ElemDataStrings::one_edge3;
      ref_elem_file[EDGE4]    = ElemDataStrings::one_edge4;

      // 2D elements
      ref_elem_file[TRI3]     = ElemDataStrings::one_tri;
      ref_elem_file[TRI6]     = ElemDataStrings::one_tri6;

      ref_elem_file[QUAD4]    = ElemDataStrings::one_quad;
      ref_elem_file[QUAD8]    = ElemDataStrings::one_quad8;
      ref_elem_file[QUAD9]    = ElemDataStrings::one_quad9;

      // 3D elements
      ref_elem_file[HEX8]     = ElemDataStrings::one_hex;
      ref_elem_file[HEX20]    = ElemDataStrings::one_hex20;
      ref_elem_file[HEX27]    = ElemDataStrings::one_hex27;

      ref_elem_file[TET4]     = ElemDataStrings::one_tet;
      ref_elem_file[TET10]    = ElemDataStrings::one_tet10;

      ref_elem_file[PRISM6]   = ElemDataStrings::one_prism;
      ref_elem_file[PRISM15]  = ElemDataStrings::one_prism15;
      ref_elem_file[PRISM18]  = ElemDataStrings::one_prism18;

      ref_elem_file[PYRAMID5] = ElemDataStrings::one_pyramid;
      ref_elem_file[PYRAMID14] = ElemDataStrings::one_pyramid14;
    }

    // Read'em
    for (FileMapType::const_iterator it=ref_elem_file.begin();
	 it != ref_elem_file.end(); ++it)
      {
	std::istringstream stream(it->second);

	read_ref_elem(it->first,
		      stream);
      }
  }


  // no reason to do this at startup -
  // data structures will get initialized *if*
  // ReferenceElem::get() is ever called.
  // // Class to setup singleton data
  // class ReferenceElemSetup : public Singleton::Setup
  // {
  //   void setup ()
  //   {
  //     init_ref_elem_table();
  //   }
  // } reference_elem_setup;

} // anonymous namespace



//----------------------------------------------------------------------------
// external API Implementation
namespace libMesh
{
  namespace ReferenceElem
  {
    const Elem & get (const ElemType Type)
    {
      libmesh_assert_less (Type, INVALID_ELEM);

      init_ref_elem_table();

      libmesh_assert (ref_elem_map[Type] != NULL);

      return *ref_elem_map[Type];
    }
  } // namespace ReferenceElem
} // namespace libMesh
