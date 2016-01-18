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



#ifndef LIBMESH_MESH_SMOOTHER_VSMOOTHER_H
#define LIBMESH_MESH_SMOOTHER_VSMOOTHER_H

#include "libmesh/libmesh_config.h"
#ifdef LIBMESH_ENABLE_VSMOOTHER

// Local Includes -----------------------------------
#include "libmesh/libmesh_common.h"
#include "libmesh/mesh_smoother.h"

// C++ Includes   -----------------------------------
#include <cstddef>
#include <vector>
#include <map>
#include <fstream>

namespace libMesh
{

// Forward declarations
class UnstructuredMesh;

/**
 * This is an implementation of Larisa Branets' smoothing algorithms.
 * The initial implementation was done by her, the adaptation to
 * libmesh was completed by Derek Gaston.  The code was heavily
 * refactored into something more closely resembling C++ by John
 * Peterson in 2014.
 *
 * Here are the relevant publications:
 * 1) L. Branets, G. Carey, "Extension of a mesh quality metric for
 * elements with a curved boundary edge or surface",
 * Journal of Computing and Information Science in Engineering, vol. 5(4), pp.302-308, 2005.
 *
 * 2) L. Branets, G. Carey, "A local cell quality metric and variational grid
 * smoothing algorithm", Engineering with Computers, vol. 21, pp.19-28, 2005.
 *
 * 3) L. Branets, "A variational grid optimization algorithm based on a local
 * cell quality metric", Ph.D. thesis, The University of Texas at Austin, 2005.
 *
 * \author Derek R. Gaston
 * \date 2006
 */
class VariationalMeshSmoother : public MeshSmoother
{
public:

  /**
   * Simple constructor to use for smoothing purposes
   */
  VariationalMeshSmoother(UnstructuredMesh & mesh,
                          double theta=0.5,
                          unsigned miniter=2,
                          unsigned maxiter=5,
                          unsigned miniterBC=5);

  /**
   * Slightly more complicated constructor for mesh redistribution based on adapt_data
   */
  VariationalMeshSmoother(UnstructuredMesh & mesh,
                          std::vector<float> * adapt_data,
                          double theta=0.5,
                          unsigned miniter=2,
                          unsigned maxiter=5,
                          unsigned miniterBC=5,
                          double percent_to_move=1);

  /**
   * Even more complicated constructor for mesh redistribution based on adapt_data with an
   * area of interest
   */
  VariationalMeshSmoother(UnstructuredMesh & mesh,
                          const UnstructuredMesh * area_of_interest,
                          std::vector<float> * adapt_data,
                          double theta=0.5,
                          unsigned miniter=2,
                          unsigned maxiter=5,
                          unsigned miniterBC=5,
                          double percent_to_move=1);

  enum MetricType
    {
      UNIFORM = 1,
      VOLUMETRIC = 2,
      DIRECTIONAL = 3
    };

  enum AdaptType
    {
      CELL = -1,
      NONE = 0,
      NODE = 1
    };

  /**
   * Destructor.
   */
  virtual ~VariationalMeshSmoother() {}

  /**
   * Redefinition of the smooth function from the
   * base class.  All this does is call the smooth
   * function in this class which takes an int, using
   * a default value of 1.
   */
  virtual void smooth() libmesh_override { _distance = this->smooth(1); }

  /**
   * The actual smoothing function, gets called whenever
   * the user specifies an actual number of smoothing
   * iterations.
   */
  double smooth(unsigned int n_iterations);

  /**
   * @return max distance a node moved during the last smooth.
   */
  double distance_moved() const { return _distance; }

  /**
   * Allow user to control whether the metric is generated from the initial mesh.
   */
  void set_generate_data(bool b) { _generate_data = b; }

  /**
   * Allow user to control the smoothing metric used.
   */
  void set_metric(MetricType t) { _metric = t; }

private:

  /**
   * Max distance of the last set of movement.
   */
  double _distance;

  /**
   * Dampening factor
   */
  const double _percent_to_move;

  /**
   * Records a relative "distance moved"
   */
  double _dist_norm;

  /**
   * Map for hanging_nodes
   */
  std::map<dof_id_type, std::vector<dof_id_type> > _hanging_nodes;

  /**
   * Vector for holding adaptive data
   */
  std::vector<float> * _adapt_data;

  /**
   * Smoother control variables
   */
  const unsigned _dim;
  const unsigned _miniter;
  const unsigned _maxiter;
  const unsigned _miniterBC;
  MetricType _metric;
  AdaptType _adaptive_func;
  const double _theta;
  bool _generate_data;

  /**
   * The number of nodes in the Mesh at the time of smoothing.
   * Not set until smooth() is actually called to mimic the
   * original code's behavior.
   */
  dof_id_type _n_nodes;

  /**
   * The number of active elements in the Mesh at the time of smoothing.
   * Not set until smooth() is actually called to mimic the
   * original code's behavior.
   */
  dof_id_type _n_cells;

  /**
   * The number of hanging node edges in the Mesh at the time of smoothing.
   * Not set until smooth() is actually called to mimic the
   * original code's behavior.
   */
  dof_id_type _n_hanging_edges;

  /**
   * All output (including debugging) is sent to the _logfile.
   */
  std::ofstream _logfile;

  /**
   * Area of Interest Mesh
   */
  const UnstructuredMesh * _area_of_interest;

  void adjust_adapt_data();
  float adapt_minimum() const;

  /**
   * 2D array type for interfacing with C APIs.
   */
  template <typename T>
  struct Array2D
  {
    Array2D(unsigned nx, unsigned ny) :
      _data(nx, std::vector<T>(ny)) {}

    // Accessors
    std::vector<T> & operator[](unsigned i) {return _data[i];}
    const std::vector<T> & operator[](unsigned i) const {return _data[i];}

  private:
    std::vector<std::vector<T> > _data;
  };



  /**
   * 3D array type for interfacing with C APIs.
   */
  template <typename T>
  struct Array3D
  {
    Array3D(unsigned nx, unsigned ny, unsigned nz)
    {
      _data.resize(nx, Array2D<T>(ny,nz));
    }

    // Accessors
    Array2D<T> & operator[](unsigned i) {return _data[i];}
    const Array2D<T> & operator[](unsigned i) const {return _data[i];}

  private:
    std::vector<Array2D<T> > _data;
  };


  int writegr(const Array2D<double> & R);

  int readgr(Array2D<double> & R,
             std::vector<int> & mask,
             Array2D<int> & cells,
             std::vector<int> & mcells,
             std::vector<int> & edges,
             std::vector<int> & hnodes);

  int readmetr(std::string name,
               Array3D<double> & H);

  int read_adp(std::vector<double> & afun);

  double jac3(double x1, double y1, double z1,
              double x2, double y2, double z2,
              double x3, double y3, double z3);

  double jac2(double x1, double y1,
              double x2, double y2);

  int basisA(Array2D<double> & Q,
             int nvert,
             const std::vector<double> & K,
             const Array2D<double> & H,
             int me);

  void adp_renew(const Array2D<double> & R,
                 const Array2D<int> & cells,
                 std::vector<double> & afun,
                 int adp);

  void full_smooth(Array2D<double> & R,
                   const std::vector<int> & mask,
                   const Array2D<int> & cells,
                   const std::vector<int> & mcells,
                   const std::vector<int> & edges,
                   const std::vector<int> & hnodes,
                   double w,
                   const std::vector<int> & iter,
                   int me,
                   const Array3D<double> & H,
                   int adp,
                   int gr);

  double maxE(Array2D<double> & R,
              const Array2D<int> & cells,
              const std::vector<int> & mcells,
              int me,
              const Array3D<double> & H,
              double v,
              double epsilon,
              double w,
              std::vector<double> & Gamma,
              double & qmin);

  double minq(const Array2D<double> & R,
              const Array2D<int> & cells,
              const std::vector<int> & mcells,
              int me,
              const Array3D<double> & H,
              double & vol,
              double & Vmin);

  double minJ(Array2D<double> & R,
              const std::vector<int> & mask,
              const Array2D<int> & cells,
              const std::vector<int> & mcells,
              double epsilon,
              double w,
              int me,
              const Array3D<double> & H,
              double vol,
              const std::vector<int> & edges,
              const std::vector<int> & hnodes,
              int msglev,
              double & Vmin,
              double & emax,
              double & qmin,
              int adp,
              const std::vector<double> & afun);

  double minJ_BC(Array2D<double> & R,
                 const std::vector<int> & mask,
                 const Array2D<int> & cells,
                 const std::vector<int> & mcells,
                 double epsilon,
                 double w,
                 int me,
                 const Array3D<double> & H,
                 double vol,
                 int msglev,
                 double & Vmin,
                 double & emax,
                 double & qmin,
                 int adp,
                 const std::vector<double> & afun,
                 int NCN);

  double localP(Array3D<double> & W,
                Array2D<double> & F,
                Array2D<double> & R,
                const std::vector<int> & cell_in,
                const std::vector<int> & mask,
                double epsilon,
                double w,
                int nvert,
                const Array2D<double> & H,
                int me,
                double vol,
                int f,
                double & Vmin,
                double & qmin,
                int adp,
                const std::vector<double> & afun,
                std::vector<double> & Gloc);

  double avertex(const std::vector<double> & afun,
                 std::vector<double> & G,
                 const Array2D<double> & R,
                 const std::vector<int> & cell_in,
                 int nvert,
                 int adp);

  double vertex(Array3D<double> & W,
                Array2D<double> & F,
                const Array2D<double> & R,
                const std::vector<int> & cell_in,
                double epsilon,
                double w,
                int nvert,
                const std::vector<double> & K,
                const Array2D<double> & H,
                int me,
                double vol,
                int f,
                double & Vmin,
                int adp,
                const std::vector<double> & g,
                double sigma);

  void metr_data_gen(std::string grid,
                     std::string metr,
                     int me);

  int solver(int n,
             const std::vector<int> & ia,
             const std::vector<int> & ja,
             const std::vector<double> & a,
             std::vector<double> & x,
             const std::vector<double> & b,
             double eps,
             int maxite,
             int msglev);

  int pcg_ic0(int n,
              const std::vector<int> & ia,
              const std::vector<int> & ja,
              const std::vector<double> & a,
              const std::vector<double> & u,
              std::vector<double> & x,
              const std::vector<double> & b,
              std::vector<double> & r,
              std::vector<double> & p,
              std::vector<double> & z,
              double eps,
              int maxite,
              int msglev);

  int pcg_par_check(int n,
                    const std::vector<int> & ia,
                    const std::vector<int> & ja,
                    const std::vector<double> & a,
                    double eps,
                    int maxite,
                    int msglev);

  void gener(char grid[], int n);
};

} // namespace libMesh

#endif // LIBMESH_ENABLE_VSMOOTHER

#endif // LIBMESH_MESH_SMOOTHER_VSMOOTHER_H
