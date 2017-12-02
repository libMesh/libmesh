#ifndef BIHARMONIC_H
#define BIHARMONIC_H

#include "libmesh/equation_systems.h"
#include "libmesh/replicated_mesh.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/mesh_refinement.h"

// Bring in bits from the libMesh namespace.
// Just the bits we're using, since this is a header.
using libMesh::EquationSystems;
using libMesh::ExodusII_IO;
using libMesh::Point;
using libMesh::Real;
using libMesh::ReplicatedMesh;

#ifdef LIBMESH_ENABLE_AMR
using libMesh::MeshRefinement;
#endif

/**
 * The Biharmonic class encapsulates most of the data structures
 * necessary to calculate the biharmonic residual and Jacobian,
 * auxiliary quantities, to take a timestep, and to output the state --
 * biharmonic solution and vectors of auxiliary quantities.
 *
 * The main reason for this design is to have a data structure that
 * has all of the necessary data in one place, where all of the
 * calculation subroutines can access these data. Currently these data
 * are split up among several interdependent objects with no clear
 * hierarchy between them: mesh, equation system, equation system
 * bundle, residual/Jacobian calculator.
 *
 * Since no object contains all others and the data are distributed
 * among many objects, the natural control and data flow resides outside
 * of these objects and is typically implemented in main().  We,
 * however, would like to split the calculation into natural chunks --
 * subroutines -- while retaining these subroutines access to the common
 * necessary data -- biharmonic parameters, mesh and time interval
 * sizes, etc. Thus, class Biharmonic.  Finally, making Biharmonic
 * inherit from EquationSystems makes it possible to include it in the
 * most common callbacks that do not pass back a user context, but only
 * an EquationSystems object.
 */
class Biharmonic : public EquationSystems
{
public:
  // Initial state enumeration
  enum InitialStateEnum {STRIP = 0,
                         ROD   = 1,
                         BALL  = 2};

  // Free energy enumeration
  enum FreeEnergyEnum {DOUBLE_WELL         = 1,
                       DOUBLE_OBSTACLE     = 2,
                       LOG_DOUBLE_WELL     = 3,
                       LOG_DOUBLE_OBSTACLE = 4};

  /**
   * Constructor retrieves command-line options, setting  defaults, if necessary.
   * It then builds the mesh using these options, then the equations systems around it,
   * and, finally, sets up the output.
   * We recommend that this be used through the factory Create function, which allocates
   * the mesh. In that case don't forget to call Destroy at the end, to free the mesh up.
   */
  Biharmonic(ReplicatedMesh & mesh);


  /**
   * Destructor
   */
  ~Biharmonic()
  {
    // delete _meshRefinement;
  }


  // Misc. getters
  bool verbose()         { return _verbose; }
  Real dt0()             { return _dt0; }
  Real dt()              { return _dt; }


  // Public interface functions
  void viewParameters();
  void init();
  void step(const Real & dt = -1.0);
  void output(int timestep, const Real & t, Real & o_t, bool force = false);
  void run();

private:
  unsigned int  _dim, _N;
  Real _kappa, _theta, _theta_c;
  Real _tol;
  bool _growth, _degenerate, _cahn_hillard, _netforce;
  FreeEnergyEnum  _energy;
  int _log_truncation;
  bool _verbose;
  InitialStateEnum  _initialState;
  Point _initialCenter;
  Real _initialWidth;
  Real _dt0, _dt, _t0, _T, _t1;
  Real _cnWeight;
  //
  std::string  _ofile_base, _ofile;
  std::unique_ptr<ExodusII_IO> _exio;
  Real    _o_dt;
  int     _o_count;
  //
  friend class JR;
  class JR;       // forward
  ReplicatedMesh & _mesh;
  JR * _jr;
};

#endif // BIHARMONIC_H
