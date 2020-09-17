// General libMesh includes
#include "libmesh/libmesh_common.h"
#include "libmesh/elem.h"
#include "libmesh/fe_base.h"
#include "libmesh/fem_context.h"
#include "libmesh/point.h"
#include "libmesh/quadrature.h"

// Local includes
#include "heatsystem.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Define the postprocess function to compute QoI 0, the weighted flux

void HeatSystem::element_qoi (DiffContext & context, const QoISet & /* qois */)
{
  FEMContext & c = cast_ref<FEMContext &>(context);

  FEBase * elem_fe = nullptr;
  c.get_element_fe( 0, elem_fe );

  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real> & JxW = elem_fe->get_JxW();

  // The number of local degrees of freedom in each variable

  unsigned int n_qpoints = c.get_element_qrule().n_points();

  // A reference to the system context is built with
  const System & sys = c.get_system();

  Number dQoI_0 = 0.;
  Number dQoI_1 = 0.;

  // Smooth QoI weight
  Number q = 0.0;

  // Loop over quadrature points
  for (unsigned int qp = 0; qp != n_qpoints; qp++)
  {
    Number u = c.interior_value(0, qp);

    // The final time QoI
    if(!std::isnan(dynamic_cast<const HeatSystem &>(sys).time_solver->get_final_time()))
    {
      if(std::abs(sys.time - dynamic_cast<const HeatSystem &>(sys).time_solver->get_final_time()) < TOLERANCE)
      {
        dQoI_0 += JxW[qp] * ( u * 1.0 );
      }
    }

    // Smooth QoI contribution
    q = 1.0;

    dQoI_1 += (JxW[qp] * ( u * q ));
  }
  // End loop over quadrature points

  // Update the computed value of the global functional R, by adding the contribution from this element
  c.get_qois()[0] = c.get_qois()[0] + dQoI_0;
  c.get_qois()[1] = c.get_qois()[1] + dQoI_1;

}
