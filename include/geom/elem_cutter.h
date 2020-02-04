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



#ifndef LIBMESH_ELEM_CUTTER_H
#define LIBMESH_ELEM_CUTTER_H

#include "libmesh/libmesh_config.h"

#if defined(LIBMESH_HAVE_TRIANGLE) && defined(LIBMESH_HAVE_TETGEN)

// libMesh includes
#include "libmesh/libmesh_common.h"
#include "libmesh/point.h"

// TIMPI includes
#include "timpi/communicator.h"

// C++ includes
#include <vector>
#include <memory>

namespace libMesh
{

// Forward declarations
template <typename> class ElemTempl;
typedef ElemTempl<Real> Elem;
class ReplicatedMesh;
class TriangleInterface;
class TetGenMeshInterface;

// This is for backwards compatibility, but if your code relies on
// forward declarations in our headers then fix it.
class SerialMesh;

/**
 * This class implements cutting a single element into a collection
 * of subelements.  This class depends on libmesh's Triangle and Tetgen
 * interfaces, the former of which is only defined if libmesh is configured
 * with --disable-strict-lgpl.
 *
 * \author Benjamin S. Kirk
 * \date 2013
 * \brief Subdivides a single element using a mesh generator.
 */
class ElemCutter
{
public:

  /**
   * Constructor.  Initializes pointer data
   * without requiring a full \p ReplicatedMesh in this header file.
   */
  ElemCutter();

  /**
   * Destructor.
   */
  ~ElemCutter();

  /**
   * \returns \p true if the element is completely inside the
   * interface defined implicitly by the vertex values of the signed
   * \p vertex_distance_func.
   */
  bool is_inside (const Elem & elem,
                  const std::vector<Real> & vertex_distance_func) const;

  /**
   * \returns \p true if the element is completely outside the
   * interface defined implicitly by the vertex values of the signed
   * \p vertex_distance_func.
   */
  bool is_outside (const Elem & elem,
                   const std::vector<Real> & vertex_distance_func) const;

  /**
   * \returns \p true if the element is cut by the interface defined
   * implicitly by the vertex values of the signed
   * \p vertex_distance_func.
   */
  bool is_cut (const Elem & elem,
               const std::vector<Real> & vertex_distance_func) const;

  /**
   * This function implements cutting an element by a signed distance
   * function. The input array \p vertex_distance_func contains the
   * vertex values of a signed distance function, from which the cutting
   * interface is inferred from the 0 level set.  If all vertex values
   * are positive, the element is outside the cutting surface and is not cut.
   * Likewise if all vertex values are negative, the element is inside the
   * cutting surface and is not cut.
   */
  void operator()(const Elem & elem_in,
                  const std::vector<Real> & vertex_distance_func);

  /**
   * \returns A list of general element pieces considered inside the
   * cutting surface.  These are subelements whose geometric union
   * defines the spatial domain of the inside portion of the cut element.
   */
  const std::vector<Elem const *> & inside_elements () const
  { return _inside_elem; }

  /**
   * \returns A list of general element pieces considered outside the
   * cutting surface.  These are subelements whose geometric union
   * defines the spatial domain of the outside portion of the cut element.
   */
  const std::vector<Elem const *> & outside_elements() const
  { return _outside_elem; }

protected:

  /**
   * Finds the points where the cutting surface
   * intersects the element edges.
   */
  void find_intersection_points(const Elem & elem,
                                const std::vector<Real> & vertex_distance_func);

  /**
   * cutting algorithm in 1D.
   */
  void cut_1D(const Elem & elem,
              const std::vector<Real> & vertex_distance_func);

  /**
   * cutting algorithm in 2D.
   */
  void cut_2D(const Elem & elem,
              const std::vector<Real> & vertex_distance_func);

  /**
   * cutting algorithm in 3D.
   */
  void cut_3D(const Elem & elem,
              const std::vector<Real> & vertex_distance_func);

  std::vector<Elem const *> _inside_elem;
  std::vector<Elem const *> _outside_elem;

  Parallel::Communicator _comm_self; // defaults to MPI_COMM_SELF

  std::unique_ptr<ReplicatedMesh> _inside_mesh_2D;
  std::unique_ptr<TriangleInterface> _triangle_inside;

  std::unique_ptr<ReplicatedMesh> _outside_mesh_2D;
  std::unique_ptr<TriangleInterface> _triangle_outside;

  std::unique_ptr<ReplicatedMesh> _inside_mesh_3D;
  std::unique_ptr<TetGenMeshInterface> _tetgen_inside;

  std::unique_ptr<ReplicatedMesh> _outside_mesh_3D;
  std::unique_ptr<TetGenMeshInterface> _tetgen_outside;

  std::vector<Point> _intersection_pts;
};


} // namespace libMesh

#endif // LIBMESH_HAVE_TRIANGLE && LIBMESH_HAVE_TETGEN
#endif // LIBMESH_ELEM_CUTTER_H
