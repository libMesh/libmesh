// $Id: mesh_smoother_vsmoother.h,v 1.1 2007-04-05 18:49:00 friedmud Exp $

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



#ifndef __mesh_smoother_vsmoother_h__
#define __mesh_smoother_vsmoother_h__

// C++ Includes   -----------------------------------
#include <vector>
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

//jj
typedef double     * LPDOUBLE;
typedef LPDOUBLE   * LPLPDOUBLE;
typedef LPLPDOUBLE   * LPLPLPDOUBLE;
typedef void     * LPVOID;
typedef LPVOID   * LPLPVOID;
typedef int  * LPINT;
typedef LPINT  * LPLPINT;

// Local Includes -----------------------------------
#include "mesh_smoother.h"
#include "mesh.h"

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
 * This is an implementation of Larisa Branets smoothing
 * algorithms.  The initial implementation was done by her,
 * the adaptation to libmesh was completed by Derek Gaston.
 *
 * \author Derek R. Gaston
 * \date 2006
 * \version $Revision: 1.1 $
 */


// ------------------------------------------------------------
// VariationalMeshSmoother class definition
class VariationalMeshSmoother : public MeshSmoother
{
public:

  /**
   * Simple constructor to use for smoothing purposes
   */
  VariationalMeshSmoother(Mesh& mesh, const double& theta=0.5, const uint& miniter=2,
			  const uint& maxiter=5, const uint& miniterBC=5)
    :MeshSmoother(mesh),
     _dim(mesh.mesh_dimension()),
     _theta(theta),
     _miniter(miniter),
     _maxiter(maxiter),
     _miniterBC(miniterBC),
     
     _metric(uniform),
     _generate_data(false),
     _adaptive_func(none),

     _area_of_interest(NULL),
     _adapt_data(NULL),
     _percent_to_move(1)
  {}

  /**
   * Slightly more complicated constructor for mesh redistribution based on adapt_data
   */
  VariationalMeshSmoother(Mesh& mesh, std::vector<float>* adapt_data, const double& theta=0.5,
			  const uint& miniter=2, const uint& maxiter=5, const uint& miniterBC=5,
			  const double& percent_to_move=1)
    :MeshSmoother(mesh),
     _dim(mesh.mesh_dimension()),
     _adapt_data(adapt_data),
     _theta(theta),
     _miniter(miniter),
     _maxiter(maxiter),
     _miniterBC(miniterBC),
     _percent_to_move(percent_to_move),
     
     _metric(uniform),
     _generate_data(false),
     _adaptive_func(cell),

     _area_of_interest(NULL)
  {}
  
  /**
   * Even more complicated constructor for mesh redistribution based on adapt_data with an
   * area of interest
   */
  VariationalMeshSmoother(Mesh& mesh, const Mesh* area_of_interest, std::vector<float>* adapt_data,
			  const double& theta=0.5, const uint& miniter=2, const uint& maxiter=5,
			  const uint& miniterBC=5, const double& percent_to_move=1)
    :MeshSmoother(mesh),
     _dim(mesh.mesh_dimension()),
     _area_of_interest(area_of_interest),
     _adapt_data(adapt_data),
     _theta(theta),
     _miniter(miniter),
     _maxiter(maxiter),
     _miniterBC(miniterBC),
     _percent_to_move(percent_to_move),

     _metric(uniform),
     _generate_data(false),
     _adaptive_func(cell)
  {}

/* Old constructors... will eventually be removed.
  VariationalMeshSmoother(Mesh& mesh, uint metric, bool generate_data, int adaptive_func, double theta,
                          uint miniter, uint maxiter, uint miniterBC, std::vector<float> &adapt_data, const double percent_to_move)
  :MeshSmoother(mesh),
  _dim(mesh.mesh_dimension()),
  _metric(metric),
  _generate_data(generate_data),
  _adaptive_func(adaptive_func),
  _theta(theta),
  _miniter(miniter),
  _maxiter(maxiter),
  _miniterBC(miniterBC),
  _adapt_data(adapt_data),
  _area_of_interest(NULL),
  _percent_to_move(percent_to_move),
  _dist_norm(0){}
  
  VariationalMeshSmoother(Mesh& mesh, uint metric, bool generate_data, int adaptive_func, double theta,
                          uint miniter, uint maxiter, uint miniterBC, std::vector<float> &adapt_data, Mesh& area_of_interest)
  :MeshSmoother(mesh),
  _dim(mesh.mesh_dimension()),
  _metric(metric),
  _generate_data(generate_data),
  _adaptive_func(adaptive_func),
  _theta(theta),
  _miniter(miniter),
  _maxiter(maxiter),
  _miniterBC(miniterBC),
  _adapt_data(adapt_data),
  _area_of_interest(&area_of_interest){}
    
  VariationalMeshSmoother(Mesh& mesh,
                          uint metric,
                          bool generate_data,
                          int adaptive_func,
                          double theta,
                          uint miniter,
                          uint maxiter,
                          uint miniterBC)
                          :MeshSmoother(mesh),_dim(mesh.mesh_dimension()),_metric(metric),_generate_data(generate_data),
                          _adaptive_func(adaptive_func),_theta(theta),_miniter(miniter),_maxiter(maxiter),
                          _miniterBC(miniterBC),_adapt_data(*(new std::vector<float>())),
                          _area_of_interest(NULL){}
*/
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
   * @brief Member function <code>distanceMoved</code>
   *
   * @return a <code>double</code> max distance a node moved during the last smooth.
   */
  double distanceMoved() const {return _distance;}

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
  std::map<unsigned int, std::vector<unsigned int> > _hanging_nodes;
  
  /**
   * Vector for holding adaptive data
   */
  std::vector<float> * _adapt_data;

  /**
   * Smoother control variables
   */
  const uint& _dim,_miniter,_maxiter,_miniterBC;
  const metric_type& _metric;
  const adapt_type& _adaptive_func;
  const double& _theta;
  const bool& _generate_data;

  /**
   * Area of Interest Mesh
   */
  const Mesh * _area_of_interest;
  
  void adjust_adapt_data();
  float adapt_minimum() const;
  
  /**
   *Imported stuff: ..........
   */

  /*-- allocate memory --*/
  LPDOUBLE alloc_d_n1(int m1)              { return((LPDOUBLE)malloc(m1*sizeof(double))); }
  LPINT alloc_i_n1(int m1)                 { return((LPINT)malloc(m1*sizeof(int))); }
  LPLPINT alloc_i_n1_n2(int m1, int m2)    { return((LPLPINT)malloc(m1*sizeof(LPINT))); }
  LPLPDOUBLE alloc_d_n1_n2(int m1, int m2) { return((LPLPDOUBLE)malloc(m1*sizeof(LPDOUBLE))); }
  LPLPLPDOUBLE alloc_d_n1_n2_n3(int m1, int m2, int m3) { return((LPLPLPDOUBLE)malloc(m1*sizeof(LPLPDOUBLE))); }

  int writegr(int n, int N, LPLPDOUBLE R, LPINT mask, int ncells, LPLPINT cells,
              LPINT mcells, int nedges, LPINT edges, LPINT hnodes, char grid[],
              int me, char grid_old[], FILE *sout);

  int readgr(int n, int N, LPLPDOUBLE R, LPINT mask, int ncells,
             LPLPINT cells, LPINT mcells, int nedges, LPINT edges, LPINT hnodes, FILE *sout);

  int readmetr(char *name, LPLPLPDOUBLE H, int ncells, int n, FILE *sout);

  int read_adp(LPDOUBLE afun, char *adap, FILE *sout);

  double jac3(double x1,double y1,double z1,double x2,double y2,
              double z2,double x3,double y3,double z3);

  double jac2(double x1,double y1,double x2,double y2);

  int basisA(int n, LPLPDOUBLE Q, int nvert, LPDOUBLE K, LPLPDOUBLE H, int me);

  void adp_renew(int n, int N, LPLPDOUBLE R, int ncells, LPLPINT cells,
                 LPDOUBLE afun, int adp, FILE *sout);

  void full_smooth(int n, int N, LPLPDOUBLE R, LPINT mask, int ncells, LPLPINT cells, LPINT mcells, 
                   int nedges, LPINT edges, LPINT hnodes, double w, LPINT iter, int me, 
                   LPLPLPDOUBLE H, int adp, char *adap, int gr, FILE *sout);

  double maxE(int n, int N, LPLPDOUBLE R, int ncells, LPLPINT cells, LPINT mcells, 
              int me, LPLPLPDOUBLE H, double v, double epsilon, double w, LPDOUBLE Gamma, 
              double *qmin, FILE *sout);

  double minq(int n, int N, LPLPDOUBLE R, int ncells, LPLPINT cells, LPINT mcells, 
              int me, LPLPLPDOUBLE H, double *vol, double *Vmin, FILE *sout);

  double minJ(int n, int N, LPLPDOUBLE R, LPINT mask, int ncells, LPLPINT cells, LPINT mcells,
              double epsilon, double w, int me, LPLPLPDOUBLE H, double vol, int nedges, 
              LPINT edges, LPINT hnodes, int msglev, double *Vmin, double *emax, double *qmin, 
              int adp, LPDOUBLE afun, FILE *sout);

  double minJ_BC(int N, LPLPDOUBLE R, LPINT mask, int ncells, LPLPINT cells, LPINT mcells,
                double epsilon, double w, int me, LPLPLPDOUBLE H, double vol, int msglev, 
                double *Vmin, double *emax, double *qmin, int adp, LPDOUBLE afun, int NCN, FILE *sout);

  double localP(int n, LPLPLPDOUBLE W, LPLPDOUBLE F, LPLPDOUBLE R, LPINT cell, LPINT mask, double epsilon, 
                double w, int nvert, LPLPDOUBLE H, int me, double vol, int f, double *Vmin, 
                double *qmin, int adp, LPDOUBLE afun, LPDOUBLE Gloc, FILE *sout);

  double avertex(int n, LPDOUBLE afun, LPDOUBLE G, LPLPDOUBLE R, LPINT cell, int nvert, int adp, FILE *sout);

  double vertex(int n, LPLPLPDOUBLE W, LPLPDOUBLE F, LPLPDOUBLE R, LPINT cell, 
                double epsilon, double w, int nvert, LPDOUBLE K, 
                LPLPDOUBLE H, int me, double vol, int f, double *Vmin, int adp, 
                LPDOUBLE G, double sigma, FILE *sout);

  void metr_data_gen(char grid[], char metr[], int n, int me, FILE *sout);

  int solver(int n, LPINT ia, LPINT ja, LPDOUBLE a, LPDOUBLE x, LPDOUBLE b, double eps,
             int maxite, int msglev, FILE *sout);

  int pcg_ic0(int n, LPINT ia, LPINT ja, LPDOUBLE a, LPDOUBLE u, LPDOUBLE x, LPDOUBLE b, LPDOUBLE r, 
              LPDOUBLE p, LPDOUBLE z, double eps, int maxite, int msglev, FILE *sout);

  int pcg_par_check(int n, LPINT ia, LPINT ja, LPDOUBLE a, double eps, int maxite, int msglev, FILE *sout);

  void gener(char grid[], int n, FILE *sout);
  
  void local_sweep(int n, int N, LPLPDOUBLE R, LPINT mask, int ncells, LPLPINT cells, LPINT mcells, 
                  int nedges, LPINT edges, LPINT hnodes, double w, LPINT iter, int me, 
                      LPLPLPDOUBLE H, int adp, int OPT, FILE *sout);
  
  double minJ_l(int n, int N, LPLPDOUBLE R, LPINT mask, int ncells, LPLPINT cells, LPINT mcells,
              double epsilon, double w, int me, LPLPLPDOUBLE H, double vol, int nedges, 
              LPINT edges, LPINT hnodes, int msglev, double *Vmin, double *emax, double *qmin, 
              int adp, LPDOUBLE afun, FILE *sout);
};


#endif
