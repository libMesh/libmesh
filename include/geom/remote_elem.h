// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include <limits>

namespace libMesh
{

/**
 * In parallel meshes where a ghost element has neighbors which do
 * not exist on the local processor, the ghost element's neighbors
 * are set to point to the singleton RemoteElement instead.
 * Library code can then distinguish between such elements and
 * boundary elements (with nullptr neighbors).
 *
 * \author Roy H. Stogner
 * \date 2007
 * \brief Used by ParallelMesh to represent an Elem owned by another processor.
 */
template <typename RealType = Real>
class RemoteElemTempl : public ElemTempl<RealType>,
                        public Singleton
{
public:
  typedef RemoteElemTempl<RealType> RemoteElem;
  typedef ElemTempl<RealType> Elem;
  typedef PointTempl<RealType> Point;
  typedef NodeTempl<RealType> Node;

  /**
   * A unique \p id to distinguish remote element links
   */
  static const dof_id_type remote_elem_id = static_cast<dof_id_type>(-2);

  /**
   * Constructor. Private to force use of the \p create() member.
   */
private:
  RemoteElemTempl () : Elem(0,
                       0,
                       nullptr,
                       _elemlinks_data,
                       nullptr)
  { this->set_id(remote_elem_id); }

  static RemoteElem * remote_elem;

public:

  RemoteElemTempl (RemoteElem &&) = delete;
  RemoteElemTempl (const RemoteElem &) = delete;
  RemoteElem & operator= (const RemoteElem &) = delete;
  RemoteElem & operator= (RemoteElem &&) = delete;

  /**
   * Sets remote_elem to nullptr.
   */
  virtual ~RemoteElemTempl();

  /**
   * Return a reference to the global \p RemoteElem
   * singleton object.
   */
  static const Elem & create ();

  /**
   * Return a reference to the global \p RemoteElem
   * singleton object.
   */
  static const RemoteElem * get_instance ()
    {
      if (!remote_elem)
        create();

      return remote_elem;
    }

  virtual Point master_point (const unsigned int /*i*/) const override
  { libmesh_not_implemented(); return Point(); }

  virtual Node * & set_node (const unsigned int i) override
  { libmesh_not_implemented(); return Elem::set_node(i); }

  /**
   * Don't hide Elem::key() defined in the base class.
   */
  using Elem::key;

  virtual dof_id_type key (const unsigned int) const override
  { libmesh_not_implemented(); return 0; }

  virtual unsigned int which_node_am_i(unsigned int /*side*/,
                                       unsigned int /*side_node*/) const override
  { libmesh_not_implemented(); return 0; }

  virtual bool is_remote () const override
  { return true; }

  virtual void connectivity(const unsigned int,
                            const IOPackage,
                            std::vector<dof_id_type> &) const override
  { libmesh_not_implemented(); }

  virtual ElemType type () const override
  { return REMOTEELEM; }

  virtual unsigned short dim () const override
  { libmesh_not_implemented(); return 0; }

  virtual unsigned int n_nodes () const override
  { libmesh_not_implemented(); return 0; }

  virtual unsigned int n_sides () const override
  { libmesh_not_implemented(); return 0; }

  virtual unsigned int n_vertices () const override
  { libmesh_not_implemented(); return 0; }

  virtual unsigned int n_edges () const override
  { libmesh_not_implemented(); return 0; }

  virtual unsigned int n_faces () const override
  { libmesh_not_implemented(); return 0; }

  virtual unsigned int n_children () const override
  { libmesh_not_implemented(); return 0; }

  virtual bool is_vertex(const unsigned int) const override
  { libmesh_not_implemented(); return false; }

  virtual bool is_edge(const unsigned int) const override
  { libmesh_not_implemented(); return false; }

  virtual bool is_face(const unsigned int) const override
  { libmesh_not_implemented(); return false; }

  virtual bool is_node_on_side(const unsigned int,
                               const unsigned int) const override
  { libmesh_not_implemented(); return false; }

  virtual std::vector<unsigned int> nodes_on_side(const unsigned int) const override
  {
    libmesh_not_implemented();
    return {0};
  }

  virtual bool is_child_on_side(const unsigned int,
                                const unsigned int) const override
  { libmesh_not_implemented(); return false; }

  virtual bool is_edge_on_side(const unsigned int,
                               const unsigned int) const override
  { libmesh_not_implemented(); return false; }

  virtual bool is_node_on_edge(const unsigned int,
                               const unsigned int) const override
  { libmesh_not_implemented(); return false; }

  virtual unsigned int n_sub_elem () const override
  { libmesh_not_implemented(); return 0; }

  virtual std::unique_ptr<Elem> side_ptr (const unsigned int) override
  { libmesh_not_implemented(); return std::unique_ptr<Elem>(); }

  virtual void side_ptr (std::unique_ptr<Elem> &,
                         const unsigned int) override
  { libmesh_not_implemented(); }

  virtual std::unique_ptr<Elem> build_side_ptr (const unsigned int,
                                                bool) override
  { libmesh_not_implemented(); return std::unique_ptr<Elem>(); }

  virtual void build_side_ptr (std::unique_ptr<Elem> &,
                               const unsigned int) override
  { libmesh_not_implemented(); }

  virtual std::unique_ptr<Elem> build_edge_ptr (const unsigned int) override
  { libmesh_not_implemented(); return std::unique_ptr<Elem>(); }

  virtual Order default_order () const override
  { libmesh_not_implemented(); return static_cast<Order>(1); }

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

  virtual bool infinite () const override
  { libmesh_not_implemented(); return false; }

#endif

#ifdef LIBMESH_ENABLE_AMR

  /**
   * Matrix that transforms the parents nodes into the children's
   * nodes.
   */
  virtual float embedding_matrix (const unsigned int,
                                  const unsigned int,
                                  const unsigned int) const override
  { libmesh_not_implemented(); return 0.; }

  LIBMESH_ENABLE_TOPOLOGY_CACHES;

#endif // LIBMESH_ENABLE_AMR


protected:

  /**
   * Data for link to (nullptr!) parent.
   */
  Elem * _elemlinks_data[1];
};

template <typename RealType>
RemoteElemTempl<RealType> *
RemoteElemTempl<RealType>::remote_elem = nullptr;


typedef RemoteElemTempl<Real> RemoteElem;

} // namespace libMesh

#endif // LIBMESH_REMOTE_ELEM_H
