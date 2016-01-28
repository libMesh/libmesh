// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_REMOTE_ELEM_H
#define LIBMESH_REMOTE_ELEM_H

// Local includes
#include "libmesh/elem.h"
#include "libmesh/libmesh_singleton.h"

// C++ includes
#include <cstddef>
#include <limits>

namespace libMesh
{

/**
 * In parallel meshes where a ghost element has neighbors which do
 * not exist on the local processor, the ghost element's neighbors
 * are set to point to the singleton RemoteElement instead.
 * Library code can then distinguish between such elements and
 * boundary elements (with NULL neighbors).
 *
 * \author Roy H. Stogner
 * \date 2007
 */
class RemoteElem : public Elem,
                   public Singleton
{
public:

  /**
   * A unique \p id to distinguish remote element links
   */
  static const dof_id_type remote_elem_id = static_cast<dof_id_type>(-2);

  /**
   * Constructor. Private to force use of the \p create() member.
   */
private:
  RemoteElem () : Elem(0,
                       0,
                       libmesh_nullptr,
                       _elemlinks_data,
                       libmesh_nullptr)
  { this->set_id(remote_elem_id); }

public:
  /**
   * Destructor.
   */
  virtual ~RemoteElem();

  /**
   * Return a reference to the global \p RemoteElem
   * sigleton object.
   */
  static const Elem & create ();

  virtual Point master_point (const unsigned int /*i*/) const libmesh_override
  { libmesh_not_implemented(); return Point(); }

  virtual Node * & set_node (const unsigned int i) libmesh_override
  { libmesh_not_implemented(); return Elem::set_node(i); }

  /**
   * Don't hide Elem::key() defined in the base class.
   */
  using Elem::key;

  virtual dof_id_type key (const unsigned int) const libmesh_override
  { libmesh_not_implemented(); return 0; }

  virtual bool is_remote () const libmesh_override
  { return true; }

  virtual void connectivity(const unsigned int,
                            const IOPackage,
                            std::vector<dof_id_type> &) const libmesh_override
  { libmesh_not_implemented(); }

  virtual ElemType type () const libmesh_override
  { return REMOTEELEM; }

  virtual unsigned int dim () const libmesh_override
  { libmesh_not_implemented(); return 0; }

  virtual unsigned int n_nodes () const libmesh_override
  { libmesh_not_implemented(); return 0; }

  virtual unsigned int n_sides () const libmesh_override
  { libmesh_not_implemented(); return 0; }

  virtual unsigned int n_vertices () const libmesh_override
  { libmesh_not_implemented(); return 0; }

  virtual unsigned int n_edges () const libmesh_override
  { libmesh_not_implemented(); return 0; }

  virtual unsigned int n_faces () const libmesh_override
  { libmesh_not_implemented(); return 0; }

  virtual unsigned int n_children () const libmesh_override
  { libmesh_not_implemented(); return 0; }

  virtual bool is_vertex(const unsigned int) const libmesh_override
  { libmesh_not_implemented(); return false; }

  virtual bool is_edge(const unsigned int) const libmesh_override
  { libmesh_not_implemented(); return false; }

  virtual bool is_face(const unsigned int) const libmesh_override
  { libmesh_not_implemented(); return false; }

  virtual bool is_node_on_side(const unsigned int,
                               const unsigned int) const libmesh_override
  { libmesh_not_implemented(); return false; }

  virtual bool is_child_on_side(const unsigned int,
                                const unsigned int) const libmesh_override
  { libmesh_not_implemented(); return false; }

  virtual bool is_edge_on_side(const unsigned int,
                               const unsigned int) const libmesh_override
  { libmesh_not_implemented(); return false; }

  virtual bool is_node_on_edge(const unsigned int,
                               const unsigned int) const libmesh_override
  { libmesh_not_implemented(); return false; }

  virtual unsigned int n_sub_elem () const libmesh_override
  { libmesh_not_implemented(); return 0; }

  virtual UniquePtr<Elem> side (const unsigned int) const libmesh_override
  { libmesh_not_implemented(); return UniquePtr<Elem>(); }

  virtual UniquePtr<Elem> build_side (const unsigned int,
                                      bool) const libmesh_override
  { libmesh_not_implemented(); return UniquePtr<Elem>(); }

  virtual UniquePtr<Elem> build_edge (const unsigned int) const libmesh_override
  { libmesh_not_implemented(); return UniquePtr<Elem>(); }

  virtual Order default_order () const libmesh_override
  { libmesh_not_implemented(); return FIRST; }

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

  virtual bool infinite () const libmesh_override
  { libmesh_not_implemented(); return false; }

#endif

#ifdef LIBMESH_ENABLE_AMR

  /**
   * Matrix that transforms the parents nodes into the children's
   * nodes
   */
  virtual float embedding_matrix (const unsigned int,
                                  const unsigned int,
                                  const unsigned int) const libmesh_override
  { libmesh_not_implemented(); return 0.; }

  LIBMESH_ENABLE_TOPOLOGY_CACHES;

#endif // LIBMESH_ENABLE_AMR


protected:

  /**
   * Data for link to (NULL!) parent
   */
  Elem * _elemlinks_data[1];
};

// Singleton RemoteElem
extern const RemoteElem * remote_elem;

} // namespace libMesh

#endif // LIBMESH_REMOTE_ELEM_H
