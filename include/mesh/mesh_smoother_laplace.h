// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_MESH_SMOOTHER_LAPLACE_H
#define LIBMESH_MESH_SMOOTHER_LAPLACE_H

// C++ Includes
#include <vector>

// Local Includes
#include "libmesh/libmesh.h"
#include "libmesh/mesh_smoother.h"

namespace libMesh
{

/**
 * This class defines the data structures necessary for Laplace
 * smoothing.
 *
 * \note This is a simple averaging smoother, which does \e not
 * guarantee that points will be smoothed to valid locations, e.g.
 * locations inside the boundary!  This aspect could use work.
 *
 * \author John W. Peterson
 * \date 2002-2007
 */
class LaplaceMeshSmoother : public MeshSmoother
{
public:
  /**
   * Constructor.  Sets the constant mesh reference
   * in the protected data section of the class.
   * @param n_iterations The number of smoothing iterations to be performed.
   */
  explicit LaplaceMeshSmoother(UnstructuredMesh &mesh,
                               const unsigned int n_iterations=1);

#ifdef LIBMESH_ENABLE_DEPRECATED
  /**
   * Constructor. Sets the constant mesh reference
   * in the protected data section of the class.
   *
   * \deprecated This constructor has been deprecated in favor of the
   * (UnstructuredMesh, const unsigned int) constructor. By specifying the
   * number of smoothing iterations in the constructor, there is no need to pass
   * this parameter to the smooth method. As such, the smooth method's signature
   * will match that of the parent class and the sister VariationalMeshSmoother
   * class.
   */
  explicit
  LaplaceMeshSmoother(UnstructuredMesh & mesh);
#endif

  /**
   * Destructor.
   */
  virtual ~LaplaceMeshSmoother() = default;

  /**
   * Redefinition of the smooth function from the
   * base class.
   */
  virtual void smooth() override;

#ifdef LIBMESH_ENABLE_DEPRECATED
  /**
   * The actual smoothing function, gets called whenever
   * the user specifies an actual number of smoothing
   * iterations.
   *
   * \deprecated The number of iterations should be set in the class constructor.
   * The parameterless smooth() override should be used to smooth.
   */
  void smooth(unsigned int n_iterations);
#endif

  /**
   * Initialization for the Laplace smoothing routine
   * is basically identical to building an "L-graph"
   * which is expensive.  It's provided separately from
   * the constructor since you may or may not want
   * to build the L-graph on construction.
   */
  void init();

  /**
   * Mainly for debugging, this function will print
   * out the connectivity graph which has been created.
   */
  void print_graph(std::ostream & out_stream = libMesh::out) const;

private:
  /**
   * This function allgather's the (local) graph after
   * it is computed on each processor by the init() function.
   */
  void allgather_graph();

  /**
   * True if the L-graph has been created, false otherwise.
   */
  bool _initialized;

  /**
   * Data structure for holding the L-graph
   */
  std::vector<std::vector<dof_id_type>> _graph;

  /**
  * Number of smoothing iterations to perform
  */
  unsigned int _n_iterations;
};


} // namespace libMesh

#endif // LIBMESH_MESH_SMOOTHER_LAPLACE_H
