#include "libmesh/fem_system.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_function.h"
#include "libmesh/fem_physics.h"
#include "libmesh/system.h"

// Using namespace libMesh
using namespace libMesh;

// FEMSystem, TimeSolver and  NewtonSolver will handle most tasks,
// but we must specify element residuals
class SigmaPhysics : public FEMPhysics
{
public:
  // Constructor
  SigmaPhysics() : FEMPhysics() {}

  // Special functions
  SigmaPhysics (SigmaPhysics &&) = default;
  SigmaPhysics (const SigmaPhysics &) = default;
  SigmaPhysics & operator= (const SigmaPhysics &) = default;
  SigmaPhysics & operator= (SigmaPhysics &&) = default;
  virtual ~SigmaPhysics () = default;

  Real & k() { return _k; }
  bool & analytic_jacobians() { return _analytic_jacobians; }

  // System initialization
  virtual void init_data (System & sys);

  virtual std::unique_ptr<DifferentiablePhysics> clone_physics() override;

protected:

  // Context initialization
  virtual void init_context (DiffContext & context) override;

  // Element residual and jacobian calculations
  // Time dependent parts
  virtual bool element_time_derivative (bool request_jacobian,
                                        DiffContext & context) override;

  // The parameters to solve for
  Real _k;

  // Variables to hold the computed QoIs
  Number computed_QoI[1];

  // The FE type to use
  std::string _fe_family;
  unsigned int _fe_order;

  // Index for T variable
  unsigned int T_var;

  // Calculate Jacobians analytically or not?
  bool _analytic_jacobians;
};
