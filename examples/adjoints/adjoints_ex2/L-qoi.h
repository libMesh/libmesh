
#ifndef L_QOI_H
#define L_QOI_H

#include "libmesh/libmesh_common.h"
#include "libmesh/elem.h"
#include "libmesh/fe_base.h"
#include "libmesh/fem_context.h"
#include "libmesh/point.h"
#include "libmesh/quadrature.h"
#include "libmesh/diff_qoi.h"
#include "libmesh/auto_ptr.h" // libmesh_make_unique

// Bring in everything from the libMesh namespace
using namespace libMesh;

class LaplaceQoI : public DifferentiableQoI
{
public:
  LaplaceQoI(){}
  virtual ~LaplaceQoI(){}

  virtual void init_qoi(std::vector<Number> & sys_qoi);

  // Context initialization
  virtual void init_context (DiffContext & context);

  virtual void postprocess() {}

  virtual void element_qoi_derivative(DiffContext & context, const QoISet & qois);

  virtual void element_qoi (DiffContext & context, const QoISet & qois);

  virtual std::unique_ptr<DifferentiableQoI> clone()
  {
    return libmesh_make_unique<LaplaceQoI>(*this);
  }

};
#endif // L_QOI_H
