// The libMesh Finite Element Library.
// Copyright (C) 2002-2014 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/mesh_smoother.h"
#include "libmesh/unstructured_mesh.h"

// C++ Includes   -----------------------------------
#include <cstddef>
#include <vector>
#include <map>

namespace libMesh
{

/**
 * This is an implementation of Larisa Branets' smoothing
 * algorithms.  The initial implementation was done by her,
 * the adaptation to libmesh was completed by Derek Gaston.
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
  VariationalMeshSmoother(UnstructuredMesh& mesh,
                          double theta=0.5,
                          unsigned miniter=2,
                          unsigned maxiter=5,
                          unsigned miniterBC=5);

  /**
   * Slightly more complicated constructor for mesh redistribution based on adapt_data
   */
  VariationalMeshSmoother(UnstructuredMesh& mesh,
                          std::vector<float>* adapt_data,
                          double theta=0.5,
                          unsigned miniter=2,
                          unsigned maxiter=5,
                          unsigned miniterBC=5,
                          double percent_to_move=1);

  /**
   * Even more complicated constructor for mesh redistribution based on adapt_data with an
   * area of interest
   */
  VariationalMeshSmoother(UnstructuredMesh& mesh,
                          const UnstructuredMesh* area_of_interest,
                          std::vector<float>* adapt_data,
                          double theta=0.5,
                          unsigned miniter=2,
                          unsigned maxiter=5,
                          unsigned miniterBC=5,
                          double percent_to_move=1);

  enum metric_type
    {
      uniform=1,
      volumetric=2,
      directional=3
    };

  enum adapt_type
    {
      cell=-1,
      none=0,
      node=1
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
  virtual void smooth() { _distance = this->smooth(1); }

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
  const unsigned _dim,_miniter,_maxiter,_miniterBC;
  const metric_type _metric;
  const adapt_type _adaptive_func;
  const double _theta;
  const bool _generate_data;

  /**
   * Area of Interest Mesh
   */
  const UnstructuredMesh * _area_of_interest;

  void adjust_adapt_data();
  float adapt_minimum() const;

  /**
   * Imported stuff
   */

  /**
   * allocate memory
   */
  double* alloc_d_n1(int m1)              { return((double*)malloc(m1*sizeof(double))); }
  int* alloc_i_n1(int m1)                 { return((int*)malloc(m1*sizeof(int))); }
  int** alloc_i_n1_n2(int m1, int)    { return((int**)malloc(m1*sizeof(int*))); }
  double** alloc_d_n1_n2(int m1, int) { return((double**)malloc(m1*sizeof(double*))); }
  double*** alloc_d_n1_n2_n3(int m1, int, int) { return((double***)malloc(m1*sizeof(double**))); }

  int writegr(int n, int N, double** R, int* mask, int ncells, int** cells,
              int* mcells, int nedges, int* edges, int* hnodes, const char grid[],
              int me, const char grid_old[], FILE *sout);

  int readgr(int n, int N, double** R, int* mask, int ncells,
             int** cells, int* mcells, int nedges, int* edges, int* hnodes, FILE *sout);

  int readmetr(char *name, double*** H, int ncells, int n, FILE *sout);

  int read_adp(double* afun, char *adap, FILE *sout);

  double jac3(double x1, double y1, double z1, double x2, double y2,
              double z2, double x3, double y3, double z3);

  double jac2(double x1, double y1, double x2, double y2);

  int basisA(int n, double** Q, int nvert, double* K, double** H, int me);

  void adp_renew(int n, int N, double** R, int ncells, int** cells,
                 double* afun, int adp, FILE *sout);

  void full_smooth(int n, int N, double** R, int* mask, int ncells, int** cells, int* mcells,
                   int nedges, int* edges, int* hnodes, double w, int* iter, int me,
                   double*** H, int adp, char *adap, int gr, FILE *sout);

  double maxE(int n, int N, double** R, int ncells, int** cells, int* mcells,
              int me, double*** H, double v, double epsilon, double w, double* Gamma,
              double *qmin, FILE *sout);

  double minq(int n, int N, double** R, int ncells, int** cells, int* mcells,
              int me, double*** H, double *vol, double *Vmin, FILE *sout);

  double minJ(int n, int N, double** R, int* mask, int ncells, int** cells, int* mcells,
              double epsilon, double w, int me, double*** H, double vol, int nedges,
              int* edges, int* hnodes, int msglev, double *Vmin, double *emax, double *qmin,
              int adp, double* afun, FILE *sout);

  double minJ_BC(int N, double** R, int* mask, int ncells, int** cells, int* mcells,
                 double epsilon, double w, int me, double*** H, double vol, int msglev,
                 double *Vmin, double *emax, double *qmin, int adp, double* afun, int NCN, FILE *sout);

  double localP(int n, double*** W, double** F, double** R, int* cell, int* mask, double epsilon,
                double w, int nvert, double** H, int me, double vol, int f, double *Vmin,
                double *qmin, int adp, double* afun, double* Gloc, FILE *sout);

  double avertex(int n, double* afun, double* G, double** R, int* cell, int nvert, int adp, FILE *sout);

  double vertex(int n, double*** W, double** F, double** R, int* cell,
                double epsilon, double w, int nvert, double* K,
                double** H, int me, double vol, int f, double *Vmin, int adp,
                double* G, double sigma, FILE *sout);

  void metr_data_gen(char grid[], char metr[], int n, int me, FILE *sout);

  int solver(int n, int* ia, int* ja, double* a, double* x, double* b, double eps,
             int maxite, int msglev, FILE *sout);

  int pcg_ic0(int n, int* ia, int* ja, double* a, double* u, double* x, double* b, double* r,
              double* p, double* z, double eps, int maxite, int msglev, FILE *sout);

  int pcg_par_check(int n, int* ia, int* ja, double* a, double eps, int maxite, int msglev, FILE *sout);

  void gener(char grid[], int n, FILE *sout);

  void local_sweep(int n, int N, double** R, int* mask, int ncells, int** cells, int* mcells,
                   int nedges, int* edges, int* hnodes, double w, int* iter, int me,
                   double*** H, int adp, int OPT, FILE *sout);

  double minJ_l(int n, int N, double** R, int* mask, int ncells, int** cells, int* mcells,
                double epsilon, double w, int me, double*** H, double vol, int nedges,
                int* edges, int* hnodes, int msglev, double *Vmin, double *emax, double *qmin,
                int adp, double* afun, FILE *sout);
};

} // namespace libMesh

#endif // LIBMESH_ENABLE_VSMOOTHER

#endif // LIBMESH_MESH_SMOOTHER_VSMOOTHER_H
