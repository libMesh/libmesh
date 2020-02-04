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



#ifndef LIBMESH_PATCH_H
#define LIBMESH_PATCH_H

// Local includes
#include "libmesh/id_types.h"

// C++ includes
#include <vector>
#include <string>
#include <set>

namespace libMesh
{

// Forward Declarations
template <typename> class ElemTempl;
typedef ElemTempl<Real> Elem;

/**
 * This class implements useful utility functions for a patch of
 * elements
 *
 * \author Varis Carey
 * \author Benjamin S. Kirk
 * \date 2004
 * \author Roy H. Stogner
 * \date 2007
 */
class Patch : public std::set<const Elem *>
{
public:

  /**
   * Constructor. Requires the processor ID to be interpreted as "local".
   */
  Patch(const processor_id_type my_procid = static_cast<processor_id_type>(-1)) :
    _my_procid(my_procid)
  {}

  /**
   * Destructor.
   */
  ~Patch() {}

  /**
   * This function finds all elements which touch the current patch at
   * a face, and adds them to the patch.
   */
  void add_face_neighbors();

  /**
   * This function finds all elements on the current processor which
   * touch the current patch at a face, and adds them to the patch.
   */
  void add_local_face_neighbors();

  /**
   * This function finds all elements which touch the current patch at
   * a face and which touch one of our processor's elements at any
   * point, and it adds them to the patch.
   */
  void add_semilocal_face_neighbors();

  /**
   * This function finds all elements which touch the current patch at
   * any point, and adds them to the patch.
   */
  void add_point_neighbors();

  /**
   * This function finds all elements on the current processor which
   * touch the current patch at any point, and adds them to the patch.
   */
  void add_local_point_neighbors();

  /**
   * This function finds all elements which touch the current patch at
   * any point and which touch one of our processor's elements at any
   * point, and it adds them to the patch.
   */
  void add_semilocal_point_neighbors();

  /**
   * Pointer to Member Function typedef
   */
  typedef void (Patch::*PMF)();

  /**
   * Erases any elements in the current patch, then builds a new patch
   * containing element \p elem by repeated addition of neighbors on
   * the current processor.  This procedure is repeated until the
   * number of elements meets or exceeds \p target_patch_size, or
   * until the patch has no more local neighbors.
   */
  void build_around_element(const Elem * elem,
                            const unsigned int target_patch_size = 10,
                            PMF patchtype = &Patch::add_local_face_neighbors);

protected:

  /**
   * This function finds all elements which
   * touch the current patch at a face
   */
  void find_face_neighbors(std::set<const Elem *> & neighbor_set);

  /**
   * This function finds all elements which
   * touch the current patch at any point
   */
  void find_point_neighbors(std::set<const Elem *> & neighbor_set);

  const processor_id_type _my_procid;
};

} // namespace libMesh


#endif // LIBMESH_PATCH_H
