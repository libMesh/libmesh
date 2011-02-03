#include "enum_fe_family.h"
#include "fem_system.h"
#include "system.h"

#include "parameter_vector.h"

#include "trilinos_epetra_matrix.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// FEMSystem, TimeSolver and  NewtonSolver will handle most tasks,
// but we must specify element residuals
class LaplaceSystem : public FEMSystem
{
public:
  // Constructor
  LaplaceSystem(EquationSystems& es,
               const std::string& name,
               const unsigned int number)
  : FEMSystem(es, name, number),
    _fe_family("LAGRANGE"), _fe_order(1),
    _analytic_jacobians(true) { qoi.resize(1); }

  std::string & fe_family() { return _fe_family;  }
  unsigned int & fe_order() { return _fe_order;  }
  bool & analytic_jacobians() { return _analytic_jacobians; }
    
  // Postprocessing function which we are going to override for this application

  virtual void postprocess(void);

  Number &get_QoI_value(std::string type, unsigned int QoI_index)
    {
      if(type == "exact")
	{
	  return exact_QoI[QoI_index];
	}
      else
	{
	  return computed_QoI[QoI_index];
	}
    }

  Number &get_parameter_value(unsigned int parameter_index)
    {
      return parameters[parameter_index];
    }

  ParameterVector &get_parameter_vector()
    {
      parameter_vector.resize(parameters.size());
      for(unsigned int i = 0; i != parameters.size(); ++i)
	{
	  parameter_vector[i] = &parameters[i];
	}
      
      return parameter_vector;
    }
 
  protected:
  // System initialization
  virtual void init_data ();

  // Context initialization
  virtual void init_context (DiffContext &context);

  // Element residual and jacobian calculations
  // Time dependent parts
  virtual bool element_time_derivative (bool request_jacobian,
					DiffContext &context);

  // Constraint parts
  virtual bool side_constraint (bool request_jacobian,
				DiffContext &context);

  // Overloading the qoi function on elements

  virtual void element_qoi_derivative
    (DiffContext &context,
     const QoISet & qois);  

  virtual void element_qoi
    (DiffContext &context,
     const QoISet & qois);  
 
  Number exact_solution (const Point&);

  // Parameters associated with the system
  std::vector<Number> parameters;
  
  // The ParameterVector object that will contain pointers to
  // the system parameters
  ParameterVector parameter_vector;

  // Variables to hold the computed QoIs
  
  Number computed_QoI[2];

  // Variables to read in the exact QoIs from l-shaped.in
 
  Number exact_QoI[2];
  
  // The FE type to use
  std::string _fe_family;
  unsigned int _fe_order;

  // Calculate Jacobians analytically or not?
  bool _analytic_jacobians;
};
